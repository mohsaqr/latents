# Inference, tidy tables and S3 methods of the growth mixture model
# (mixture_regression() with `random =`).

#' Covariance of the growth mixture estimates
#'
#' The observed information is the numerical Jacobian of the analytic total
#' score; `"robust"` is the person-level sandwich, the only valid choice
#' under sampling weights.
#'
#' @param fit A `latents_growth_mixture` fit.
#' @param vcov_type `"observed"`, `"robust"` or `"opg"`.
#' @return A list with `theta`, `vcov`, `vcov_type`, `information` and
#'   `score_norm` (the largest absolute total score at the estimate).
#' @noRd
.growth_inference <- function(fit, vcov_type = c("observed", "robust", "opg")) {
  vcov_type <- match.arg(vcov_type)
  spec <- fit$spec
  stats <- fit$stats
  params <- fit$params
  theta <- .growth_pack(spec, params)
  # One weight per independent unit (person, or cluster); unweighted, one.
  weights <- spec$sampling_weights %||% 1
  scores <- .growth_scores(spec, stats, params)
  total <- function(v) {
    colSums(.growth_scores(spec, stats, .growth_unpack(spec, v, params)) * weights)
  }
  # Weighted OPG, or central differences of the weighted total score with
  # h = 1e-5 * (1 + |theta|) (the service's "central"/"plus1" variant).
  information <- if (identical(vcov_type, "opg")) crossprod(scores * sqrt(weights)) else
    .inference_score_information(total, theta, 1e-5, scheme = "central",
                                 h_rule = "plus1")
  dimnames(information) <- list(names(theta), names(theta))
  inverse <- .mixture_invert(information)
  vcov <- if (identical(vcov_type, "robust")) {
    .inference_sandwich(inverse, crossprod(scores * weights))
  } else inverse
  dimnames(vcov) <- list(names(theta), names(theta))
  list(theta = theta, vcov = vcov, vcov_type = vcov_type, information = information,
       score_norm = max(abs(colSums(scores * weights))))
}

#' The inference of a growth fit at a requested covariance type
#' @noRd
.growth_resolve_inference <- function(x, vcov_type = NULL) {
  stored <- x$inference
  wanted <- vcov_type %||% stored$vcov_type %||% if (.latents_is_weighted(x))
    "robust" else "observed"
  wanted <- .latents_weighted_vcov(.latents_is_weighted(x), wanted, is.null(vcov_type))
  if (!is.null(stored) && identical(stored$vcov_type, wanted)) return(stored)
  .growth_inference(x, wanted)
}

#' Class names of a growth fit
#' @noRd
.growth_class_names <- function(x) paste0("class_", seq_len(x$spec$n_classes))

#' Delta-method table of a transformation of the estimates
#'
#' Intervals are formed on `scale` (`"identity"`, `"log"` for variances and
#' standard deviations, `"atanh"` for correlations, `"logit"` for shares) and
#' carried back, so they respect the parameter's range.
#' @noRd
.growth_delta_table <- function(inference, transform, scale, level) {
  link <- switch(scale, identity = identity, log = log, atanh = atanh,
                 logit = stats::qlogis)
  inverse <- switch(scale, identity = identity, log = exp, atanh = tanh,
                    logit = stats::plogis)
  natural <- .mixture_delta(inference$theta, inference$vcov, transform)
  linked <- .mixture_delta(inference$theta, inference$vcov,
                           function(v) link(transform(v)))
  z <- .inference_critical(level, "one_minus")
  data.frame(estimate = natural$estimate, std_error = natural$std_error,
             conf_low = inverse(linked$estimate - z * linked$std_error),
             conf_high = inverse(linked$estimate + z * linked$std_error))
}

#' Regression coefficients, one row per class and term
#' @noRd
.growth_coefficient_table <- function(x, level, inference) {
  theta <- inference$theta
  se <- sqrt(pmax(diag(inference$vcov), 0))
  picked <- startsWith(names(theta), "coefficient.")
  parts <- strsplit(sub("^coefficient\\.", "", names(theta)[picked]), ".", fixed = TRUE)
  class <- vapply(parts, `[`, character(1), 1L)
  term <- vapply(parts, function(p) paste(p[-1L], collapse = "."), character(1))
  data.frame(class = class, term = term,
             .inference_wald(unname(theta[picked]), unname(se[picked]), level),
             stringsAsFactors = FALSE)
}

#' Random-effect variances, standard deviations and correlations, per class
#' @noRd
.growth_random_table <- function(x, level, inference) {
  spec <- x$spec
  terms <- colnames(spec$random_design)
  n_random <- length(terms)
  classes <- .growth_class_names(x)
  template <- x$params
  covariance_of <- function(v, k) {
    .growth_covariances(spec, .growth_unpack(spec, v, template))[[k]]
  }
  pairs <- which(upper.tri(diag(n_random)), arr.ind = TRUE)
  # A covariance shared by every class is one set of parameters: reported once.
  shared <- identical(spec$random_covariance, "equal")
  if (shared) classes <- "all"
  do.call(rbind, lapply(seq_along(classes), function(k) {
    variances <- lapply(seq_len(n_random), function(j) {
      cbind(data.frame(class = classes[k], parameter = "variance", term = terms[j],
                       stringsAsFactors = FALSE),
            .growth_delta_table(inference, function(v) covariance_of(v, k)[j, j],
                                "log", level))
    })
    deviations <- lapply(seq_len(n_random), function(j) {
      cbind(data.frame(class = classes[k], parameter = "sd", term = terms[j],
                       stringsAsFactors = FALSE),
            .growth_delta_table(inference,
                                function(v) sqrt(covariance_of(v, k)[j, j]), "log", level))
    })
    correlations <- if (isTRUE(spec$random_diagonal)) list() else
      lapply(seq_len(nrow(pairs)), function(p) {
        a <- pairs[p, 1L]
        b <- pairs[p, 2L]
        cbind(data.frame(class = classes[k], parameter = "correlation",
                         term = paste(terms[a], terms[b], sep = ":"),
                         stringsAsFactors = FALSE),
              .growth_delta_table(inference, function(v) {
                covariance <- covariance_of(v, k)
                covariance[a, b] / sqrt(covariance[a, a] * covariance[b, b])
              }, "atanh", level))
      })
    rows <- do.call(rbind, c(variances, deviations, correlations))
    # On a degenerate covariance the likelihood is flat: no Wald error or
    # interval is meaningful there, so none is reported.
    rows$boundary <- .growth_class_is_degenerate(x, k)
    rows[rows$boundary, c("std_error", "conf_low", "conf_high")] <- NA_real_
    rows
  }))
}

#' Whether class k's random-effect covariance is degenerate
#' @noRd
.growth_class_is_degenerate <- function(x, k) {
  degenerate <- x$degenerate_classes %||% character()
  identical(degenerate, "all classes") || sprintf("class %d", k) %in% degenerate
}

#' Class sizes, residual standard deviations and classification quality
#' @noRd
.growth_class_table <- function(x, level, inference) {
  spec <- x$spec
  template <- x$params
  classes <- .growth_class_names(x)
  posterior <- x$expectation$posterior
  weights <- spec$sampling_weights %||% rep(1, nrow(posterior))
  modal <- max.col(posterior, ties.method = "first")
  shares <- do.call(rbind, lapply(seq_along(classes), function(k) {
    .growth_delta_table(inference, function(v) {
      prior <- .growth_class_prior(spec, .growth_unpack(spec, v, template))
      sum(weights * prior[, k]) / sum(weights)
    }, "logit", level)
  }))
  sigma <- do.call(rbind, lapply(seq_along(classes), function(k) {
    .growth_delta_table(inference, function(v) {
      sqrt(.growth_unpack(spec, v, template)$sigma2[k])
    }, "log", level)
  }))
  average <- vapply(seq_along(classes), function(k) {
    if (!any(modal == k)) NA_real_ else mean(posterior[modal == k, k])
  }, numeric(1))
  data.frame(class = classes, share = shares$estimate, share_std_error = shares$std_error,
             share_conf_low = shares$conf_low, share_conf_high = shares$conf_high,
             expected_persons = colSums(posterior * weights),
             assigned_persons = tabulate(modal, length(classes)),
             average_posterior = average, residual_sd = sigma$estimate,
             residual_sd_std_error = sigma$std_error, stringsAsFactors = FALSE)
}

#' Each person's prior class probabilities
#'
#' The membership logits of their covariates; in a multilevel growth mixture,
#' averaged over the group classes of their cluster with the cluster's
#' group-class probabilities.
#'
#' @param spec The growth specification.
#' @param params A parameter list.
#' @return A persons-by-classes matrix whose rows sum to one.
#' @noRd
.growth_class_prior <- function(spec, params) {
  structure <- .growth_structure(spec)
  if (identical(structure$nesting, "observation")) {
    return(exp(.mixture_log_softmax(spec$w, params$gamma)))
  }
  group_prior <- exp(.mixture_log_softmax(structure$v, params$delta))
  Reduce(`+`, lapply(seq_len(structure$n_group_classes), function(h) {
    exp(.mixture_log_softmax(.mixture_two_level_design(structure, h),
                             params$class_logits)) *
      group_prior[structure$group_index, h]
  }))
}

#' Membership logits, one row per class (against class 1) and term
#' @noRd
.growth_membership_table <- function(x, level, inference) {
  theta <- inference$theta
  se <- sqrt(pmax(diag(inference$vcov), 0))
  # Group-class logits of a multilevel growth mixture (against group class
  # 1) come first, labelled by their group class.
  picked <- startsWith(names(theta), "membership.") |
    startsWith(names(theta), "group_membership.")
  if (!any(picked)) {
    return(data.frame(class = character(), term = character(), estimate = numeric(),
                      std_error = numeric(), statistic = numeric(), p_value = numeric(),
                      conf_low = numeric(), conf_high = numeric()))
  }
  parts <- strsplit(sub("^(group_)?membership\\.", "", names(theta)[picked]), ".",
                    fixed = TRUE)
  data.frame(class = vapply(parts, `[`, character(1), 1L),
             term = vapply(parts, function(p) paste(p[-1L], collapse = "."), character(1)),
             .inference_wald(unname(theta[picked]), unname(se[picked]), level),
             stringsAsFactors = FALSE)
}

#' Relative entropy of the person posteriors
#' @noRd
.growth_entropy <- function(posterior) {
  if (ncol(posterior) < 2L) return(NA_real_)
  logs <- ifelse(posterior > 0, log(posterior), 0)
  1 - sum(-posterior * logs) / (nrow(posterior) * log(ncol(posterior)))
}

#' The random effects in words: "intercept", "intercept + time", "time"
#' @noRd
.growth_random_description <- function(spec) {
  if (is.null(spec$random)) return("none")
  paste(sub("^\\(Intercept\\)$", "intercept", colnames(spec$random_design)),
        collapse = " + ")
}

#' Classification entropy of the person posteriors (for ICL)
#' @noRd
.growth_entropy_sum <- function(x) {
  posterior <- x$expectation$posterior
  -sum(posterior[posterior > 0] * log(posterior[posterior > 0]))
}

#' One-row summary of the fit
#' @noRd
.growth_fit_table <- function(x) {
  spec <- x$spec
  n_persons <- length(x$stats$n)
  k <- x$n_parameters
  log_likelihood <- x$log_likelihood
  if (!is.null(spec$cluster)) return(.growth_cluster_fit_table(x))
  data.frame(
    n_classes = spec$n_classes, n_persons = n_persons, n_observations = spec$n,
    random = .growth_random_description(spec), random_covariance = spec$random_covariance,
    random_diagonal = isTRUE(spec$random_diagonal), residual_variance = spec$variance,
    n_parameters = k, log_likelihood = log_likelihood,
    aic = -2 * log_likelihood + 2 * k,
    bic = -2 * log_likelihood + log(n_persons) * k,
    bic_rows = -2 * log_likelihood + log(spec$n) * k,
    sabic = -2 * log_likelihood + log((n_persons + 2) / 24) * k,
    icl = -2 * log_likelihood + log(n_persons) * k + 2 * .growth_entropy_sum(x),
    entropy = .growth_entropy(x$expectation$posterior),
    smallest_share = min(colSums(x$expectation$posterior *
                                   (spec$sampling_weights %||% 1)) /
                           sum(spec$sampling_weights %||% rep(1, n_persons))),
    random_boundary = isTRUE(x$random_boundary),
    converged = isTRUE(x$converged), n_best_replicated = x$n_best_replicated,
    weights = x$weights %||% NA_character_, stringsAsFactors = FALSE)
}

#' Fit of a multilevel growth mixture
#'
#' The clusters are the independent units, so `bic` counts clusters, as for
#' the two-level regression mixture and as Lukociene, Varriale and Vermunt
#' (2010) recommend for choosing the number of group classes; `bic_persons`
#' and `bic_rows` count persons and observations instead.
#' @noRd
.growth_cluster_fit_table <- function(x) {
  spec <- x$spec
  n_persons <- length(x$stats$n)
  n_clusters <- spec$n_clusters
  k <- x$n_parameters
  log_likelihood <- x$log_likelihood
  rho <- x$expectation$structure$rho
  group_entropy <- -sum(rho[rho > 0] * log(rho[rho > 0]))
  bic <- -2 * log_likelihood + log(n_clusters) * k
  data.frame(
    n_classes = spec$n_classes, n_group_classes = spec$n_cluster_classes,
    n_clusters = n_clusters, n_persons = n_persons, n_observations = spec$n,
    random = .growth_random_description(spec), random_covariance = spec$random_covariance,
    random_diagonal = isTRUE(spec$random_diagonal), residual_variance = spec$variance,
    n_parameters = k, log_likelihood = log_likelihood,
    aic = -2 * log_likelihood + 2 * k, bic = bic,
    bic_persons = -2 * log_likelihood + log(n_persons) * k,
    bic_rows = -2 * log_likelihood + log(spec$n) * k,
    sabic = -2 * log_likelihood + log((n_clusters + 2) / 24) * k,
    icl = bic + 2 * (.growth_entropy_sum(x) + group_entropy),
    entropy = .growth_entropy(x$expectation$posterior),
    group_entropy = .growth_entropy(rho),
    smallest_share = min(colMeans(x$expectation$posterior)),
    random_boundary = isTRUE(x$random_boundary),
    converged = isTRUE(x$converged), n_best_replicated = x$n_best_replicated,
    stringsAsFactors = FALSE)
}

#' Group classes of a multilevel growth mixture: each group class's share of
#' the clusters and its trajectory-class probabilities
#' @noRd
.growth_group_class_table <- function(x, level, inference) {
  spec <- x$spec
  structure <- .growth_structure(spec)
  if (identical(structure$nesting, "observation")) {
    return(data.frame(group_class = character(), class = character(),
                      probability = numeric(), group_share = numeric()))
  }
  rho <- x$expectation$structure$rho
  n_group_classes <- structure$n_group_classes
  group_names <- paste0("group_class_", seq_len(n_group_classes))
  classes <- .growth_class_names(x)
  # Model-implied class probabilities within each group class, averaged over
  # the persons of the clusters it is responsible for.
  probabilities <- function(params, cluster_weights) {
    vapply(seq_len(n_group_classes), function(h) {
      prior <- exp(.mixture_log_softmax(.mixture_two_level_design(structure, h),
                                        params$class_logits))
      weight <- cluster_weights[structure$group_index, h]
      colSums(prior * weight) / sum(weight)
    }, numeric(spec$n_classes)) |> matrix(spec$n_classes)
  }
  estimate <- probabilities(x$params, rho)
  standard_error <- if (ncol(structure$w) == 0L && !is.null(inference)) {
    .mixture_delta(inference$theta, inference$vcov, function(v) {
      as.vector(probabilities(.growth_unpack(spec, v, x$params), rho))
    })$std_error
  } else rep(NA_real_, length(estimate))
  data.frame(
    group_class = rep(group_names, each = spec$n_classes),
    class = rep(classes, n_group_classes),
    probability = as.vector(estimate), std_error = standard_error,
    group_share = rep(colMeans(rho), each = spec$n_classes),
    assigned_clusters = rep(tabulate(max.col(rho, ties.method = "first"),
                                     n_group_classes), each = spec$n_classes),
    stringsAsFactors = FALSE)
}

#' Each cluster's group-class posteriors and modal group class
#' @noRd
.growth_cluster_table <- function(x) {
  spec <- x$spec
  if (is.null(spec$cluster)) {
    return(data.frame(cluster = character(), group_class = character()))
  }
  rho <- x$expectation$structure$rho
  modal <- max.col(rho, ties.method = "first")
  out <- data.frame(cluster = spec$cluster_levels,
                    n_persons = tabulate(spec$cluster_index, spec$n_clusters),
                    group_class = paste0("group_class_", modal),
                    probability = rho[cbind(seq_along(modal), modal)],
                    stringsAsFactors = FALSE)
  names(out)[1L] <- spec$cluster
  posterior <- as.data.frame(rho)
  names(posterior) <- paste0("posterior_group_class_", seq_len(ncol(rho)))
  cbind(out, posterior)
}

#' Each person's posterior class probabilities and modal class
#' @noRd
.growth_assignment_table <- function(x) {
  view <- .trajectory_view(x)
  posterior <- view$posterior
  out <- data.frame(id = view$spec$group_levels, class = view$classes[view$modal],
                    probability = apply(posterior, 1L, max), stringsAsFactors = FALSE)
  names(out)[1L] <- view$spec$id
  if (!is.null(view$spec$cluster)) {
    out <- cbind(out[1L], stats::setNames(
      data.frame(view$spec$cluster_levels[view$spec$cluster_index]), view$spec$cluster),
      out[-1L])
  }
  posterior_frame <- as.data.frame(posterior)
  names(posterior_frame) <- paste0("posterior_", view$classes)
  cbind(out, posterior_frame)
}

#' Each person's predicted random effects (posterior means, modal class)
#' @noRd
.growth_random_effects_table <- function(x) {
  posterior <- x$expectation$posterior
  modal <- max.col(posterior, ties.method = "first")
  terms <- colnames(x$spec$random_design)
  effects <- t(vapply(seq_along(modal), function(i) {
    x$expectation$means[[modal[i]]][i, ]
  }, numeric(length(terms)))) |> matrix(length(modal), length(terms))
  out <- data.frame(id = x$spec$group_levels, class = .growth_class_names(x)[modal],
                    effects, stringsAsFactors = FALSE)
  names(out) <- c(x$spec$id, "class", terms)
  out
}

#' The variable trajectories are drawn along
#'
#' The caller's choice, or the first numeric variable of `random` (the time
#' scale of a random slope), or else the first numeric predictor.
#' @noRd
.growth_time_variable <- function(x, time = NULL) {
  data <- x$spec$model_data
  if (!is.null(time)) {
    if (!is.character(time) || length(time) != 1L || !time %in% names(data) ||
        !is.numeric(data[[time]])) {
      stop(errorCondition("`time` must name a numeric variable of the model.",
                          class = "latents_bad_argument", call = NULL))
    }
    return(time)
  }
  candidates <- c(all.vars(x$spec$random), all.vars(stats::delete.response(x$spec$terms)))
  candidates <- setdiff(unique(candidates), x$spec$id)
  numeric_ones <- candidates[vapply(candidates, function(v) is.numeric(data[[v]]),
                                    logical(1))]
  if (length(numeric_ones) == 0L) {
    stop(errorCondition("The model has no numeric predictor to draw trajectories along.",
                        class = "latents_bad_argument", call = NULL))
  }
  numeric_ones[1L]
}

#' Regression design of new rows, in the fitted column order
#' @noRd
.growth_design <- function(x, newdata) {
  spec <- x$spec
  model_terms <- stats::delete.response(spec$terms)
  frame <- stats::model.frame(model_terms, newdata, xlev = spec$xlevels,
                              na.action = stats::na.fail)
  full <- stats::model.matrix(model_terms, frame, contrasts.arg = spec$contrasts)
  cbind(full[, colnames(spec$x), drop = FALSE], full[, colnames(spec$z), drop = FALSE])
}

#' Mean trajectory of every class over a time grid, with Wald bands
#'
#' Other predictors are held at their mean (numeric) or first level.
#' @noRd
.growth_trajectory_table <- function(x, level, inference, time = NULL, n_points = 60L) {
  time <- .growth_time_variable(x, time)
  data <- x$spec$model_data
  grid <- seq(min(data[[time]]), max(data[[time]]), length.out = n_points)
  predictors <- setdiff(all.vars(stats::delete.response(x$spec$terms)), x$spec$id)
  rows <- .mixture_typical_rows(data[predictors], time, grid)
  band <- .trajectory_band(x, rows, inference, level)
  out <- data.frame(class = band$class, time = rep(grid, x$spec$n_classes),
                    band[c("estimate", "std_error", "conf_low", "conf_high")],
                    stringsAsFactors = FALSE)
  names(out)[names(out) == "time"] <- time
  out
}

#' Every observation with its class trajectory and its person's own curve
#'
#' `class_mean` is the modal class's fixed trajectory; `predicted` adds the
#' person's predicted random effects (the conditional, empirical Bayes fit).
#' @noRd
.growth_individual_table <- function(x) {
  view <- .trajectory_view(x)
  spec <- view$spec
  modal <- view$modal
  person <- spec$group_index
  design <- cbind(spec$x, spec$z)
  linear <- vapply(seq_len(spec$n), function(r) {
    sum(design[r, ] * .trajectory_coefficients(x, modal[person[r]]))
  }, numeric(1)) + spec$offset
  # With random effects, a person's own curve adds their predicted effects
  # under their most likely class; without, it is the class trajectory.
  own <- if (view$growth) {
    effects <- t(vapply(seq_len(spec$n), function(r) {
      x$expectation$means[[modal[person[r]]]][person[r], ]
    }, numeric(ncol(spec$random_design)))) |> matrix(spec$n)
    linear + rowSums(spec$random_design * effects)
  } else linear
  out <- data.frame(id = spec$model_data[[spec$id]], class = view$classes[modal[person]],
                    observed = .trajectory_observed(view),
                    class_mean = view$response_mean(linear, modal[person]),
                    predicted = view$response_mean(own, modal[person]),
                    stringsAsFactors = FALSE)
  names(out)[1L] <- spec$id
  time <- .growth_time_variable(x)
  out[[time]] <- spec$model_data[[time]]
  out[c(spec$id, time, "class", "observed", "class_mean", "predicted")]
}

.growth_tables <- function(x = NULL) {
  c("coefficients", "classes", "random", "membership", "fit", "assignments",
    "random_effects", "trajectories", "individual", "starts",
    if (!is.null(x$spec$cluster)) c("group_classes", "clusters"))
}

#' @noRd
.growth_table <- function(x, what, level = 0.95, vcov_type = NULL, time = NULL,
                          data = NULL, truth = NULL) {
  needs_inference <- what %in% c("coefficients", "classes", "random", "membership",
                                 "trajectories", "group_classes")
  inference <- if (needs_inference) .growth_resolve_inference(x, vcov_type)
  table <- switch(what,
    coefficients = .growth_coefficient_table(x, level, inference),
    classes = .growth_class_table(x, level, inference),
    random = .growth_random_table(x, level, inference),
    membership = .growth_membership_table(x, level, inference),
    fit = .growth_fit_table(x),
    assignments = .growth_assignment_table(x),
    random_effects = .growth_random_effects_table(x),
    trajectories = .growth_trajectory_table(x, level, inference, time),
    individual = .growth_individual_table(x),
    starts = x$starts,
    group_classes = .growth_group_class_table(x, level, inference),
    clusters = .growth_cluster_table(x),
    # The persons' posteriors in the place the regression table reads them.
    recovery = .mixture_recovery_table(
      list(spec = x$spec, expectation = list(group_tau = x$expectation$posterior)),
      data, truth, "class"))
  .growth_present(table, what, level, x)
}

#' Attach the printed layout of a growth table
#'
#' The data are untouched; this only chooses the title, columns, labels and
#' number formats `print()` uses.
#' @noRd
.growth_present <- function(table, what, level, x) {
  column <- .latents_column
  interval <- function(low, high, from) {
    column(sprintf("%s%% CI", format(100 * level)), c(low, high), "ci", from)
  }
  percent <- format(100 * level)
  id <- x$spec$id
  time <- if (what %in% c("trajectories", "individual")) .growth_time_variable(x) else NULL
  # The recovery table's second column is the caller's truth column.
  truth_label <- if (identical(what, "recovery")) names(table)[2L] else NULL
  plan <- switch(what,
    coefficients = list(
      title = sprintf("Trajectory coefficients (%s%% CI)", percent),
      display = list(column("Class", "class", "label"), column("Term", "term", "label"),
                     column("Estimate", "estimate"),
                     interval("conf_low", "conf_high", "estimate"),
                     column("p", "p_value", "p"))),
    classes = list(
      title = "Classes",
      display = list(column("Class", "class", "label"), column("Share", "share"),
                     interval("share_conf_low", "share_conf_high", "share"),
                     column("Persons", "assigned_persons", "integer"),
                     column("Avg. posterior", "average_posterior"),
                     column("Residual SD", "residual_sd"))),
    random = list(
      title = sprintf("Random effects (%s%% CI)", percent),
      display = list(column("Class", "class", "label"),
                     column("Parameter", "parameter", "label"),
                     column("Term", "term", "label"), column("Estimate", "estimate"),
                     interval("conf_low", "conf_high", "estimate"))),
    membership = list(
      title = sprintf("Class membership (log odds against Class 1, %s%% CI)", percent),
      display = list(column("Class", "class", "label"), column("Term", "term", "label"),
                     column("Estimate", "estimate"),
                     interval("conf_low", "conf_high", "estimate"),
                     column("p", "p_value", "p"))),
    fit = list(
      title = "Model fit",
      display = list(column("Classes", "n_classes", "integer"),
                     if (!is.null(x$spec$cluster))
                       column("Group classes", "n_group_classes", "integer"),
                     if (!is.null(x$spec$cluster))
                       column("Clusters", "n_clusters", "integer"),
                     column("Persons", "n_persons", "integer"),
                     column("Observations", "n_observations", "integer"),
                     column("Random effects", "random", "text"),
                     column("Covariance", "random_covariance", "text"),
                     column("Parameters", "n_parameters", "integer"),
                     column("Log likelihood", "log_likelihood"),
                     column("AIC", "aic"), column("BIC", "bic"), column("SABIC", "sabic"),
                     column("ICL", "icl"), column("Entropy", "entropy"),
                     column("Smallest class", "smallest_share", "percent"),
                     column("Converged", "converged", "logical"),
                     column("Boundary", "random_boundary", "logical"))),
    assignments = list(
      title = "Class assignments",
      display = list(column(id, id, "text"),
                     if (!is.null(x$spec$cluster))
                       column(x$spec$cluster, x$spec$cluster, "text"),
                     column("Class", "class", "label"),
                     column("Probability", "probability"))),
    random_effects = list(
      title = "Predicted random effects (most likely class)",
      display = c(list(column(id, id, "text"), column("Class", "class", "label")),
                  lapply(colnames(x$spec$random_design), function(term) {
                    column(.latents_label(term), term)
                  }))),
    trajectories = list(
      title = sprintf("Class trajectories (%s%% CI)", percent),
      display = list(column("Class", "class", "label"), column(time, time),
                     column("Estimate", "estimate"),
                     interval("conf_low", "conf_high", "estimate"))),
    individual = list(
      title = "Observations with fitted trajectories",
      display = list(column(id, id, "text"), column(time, time),
                     column("Class", "class", "label"), column("Observed", "observed"),
                     column("Class mean", "class_mean"), column("Predicted", "predicted"))),
    starts = list(title = "Starts", display = NULL),
    group_classes = list(
      title = "Group classes: trajectory classes within each",
      display = list(column("Group class", "group_class", "label"),
                     column("Class", "class", "label"),
                     column("Probability", "probability", "share"),
                     column("SE", "std_error", "share"),
                     column("Share of clusters", "group_share", "share"),
                     column("Clusters", "assigned_clusters", "integer"))),
    clusters = list(
      title = "Cluster group classes",
      display = list(column(x$spec$cluster %||% "cluster", x$spec$cluster %||% "cluster",
                            "text"),
                     column("Persons", "n_persons", "integer"),
                     column("Group class", "group_class", "label"),
                     column("Probability", "probability"))),
    recovery = list(
      title = "Recovery of a known classification",
      display = list(column("Assigned", "assigned", "label"),
                     column(truth_label, truth_label, "text"),
                     column("Persons", "n", "integer"),
                     column("Share of assigned", "share", "share"))))
  # Columns a single-level fit does not have are NULL in the plan.
  display <- if (is.null(plan$display)) NULL else Filter(Negate(is.null), plan$display)
  table <- .latents_table(table, level, plan$title, display)
  if (identical(what, "fit")) attr(table, "card") <- TRUE
  if (identical(what, "random")) {
    attr(table, "marks") <- ifelse(table$boundary, "  boundary: no interval", "")
  }
  table
}

#' Tidy results of a growth mixture model
#'
#' Every table of a `mixture_regression()` fit with `random =` is a base
#' `data.frame` with one row per observation unit of that table.
#'
#' @param x A `latents_growth_mixture` fit.
#' @param what Which table: `"coefficients"` (one row per class and term,
#'   with Wald statistics; shared terms have class `"common"`), `"classes"`
#'   (share with interval, expected and assigned persons, average posterior,
#'   residual standard deviation), `"random"` (per class: random-effect
#'   variances, standard deviations and correlations with intervals formed on
#'   the log and Fisher-z scales), `"membership"` (membership logits against
#'   class 1), `"fit"` (one row: log likelihood, information criteria with
#'   persons as the sample size, entropy, convergence), `"assignments"` (one
#'   row per person: modal class and posteriors), `"random_effects"` (one row
#'   per person: predicted random effects under the modal class),
#'   `"trajectories"` (one row per class and time point: the mean trajectory
#'   with a confidence band), `"individual"` (one row per observation:
#'   observed value, class trajectory and the person's own predicted curve),
#'   `"starts"` (one row per start), `"recovery"` (needs `data` and `truth`:
#'   the modal classes cross-tabulated against a known classification that is
#'   constant within persons, with counts and shares of each assigned class),
#'   or `"all"` (a named list of every table but `"recovery"`). A multilevel
#'   growth mixture (`cluster`) adds `"group_classes"` (one row per group
#'   class and class: the class probabilities within the group class, their
#'   standard errors without membership covariates, the group class's share
#'   of the clusters and its assigned clusters) and `"clusters"` (one row per
#'   cluster: persons, modal group class and posteriors); its `"membership"`
#'   table holds the group-class logits (against group class 1) and the class
#'   logits with one intercept per group class, `"assignments"` gains each
#'   person's cluster, and `"fit"` counts clusters in `bic` (persons in
#'   `bic_persons`).
#' @param level Confidence level of the intervals.
#' @param vcov_type `NULL` (the stored type), `"observed"`, `"robust"` or
#'   `"opg"`. A weighted fit allows `"robust"` only.
#' @param time For `"trajectories"` and `"individual"`, the numeric variable
#'   to draw along; by default the first numeric variable of `random`.
#' @param data,truth For `"recovery"`: the data frame the model was fitted to
#'   and the name of its column holding the known classes.
#' @param ... Unused.
#' @return A base `data.frame` (or, for `"all"`, a named list of them).
#' @examples
#' \donttest{
#' fit <- mixture_regression(score ~ wave, growth_scores, n_classes = 3,
#'                           id = "student", class_level = "group",
#'                           random = "wave", random_covariance = "equal",
#'                           n_starts = 3, seed = 1)
#' get_results(fit, "classes")
#' get_results(fit, "random")
#' get_results(fit, "recovery", data = growth_scores, truth = "trajectory")
#' }
#' @export
get_results.latents_growth_mixture <- function(x, what = "coefficients", level = 0.95,
                                               vcov_type = NULL, time = NULL,
                                               data = NULL, truth = NULL, ...) {
  what <- match.arg(what, c(.growth_tables(x), "recovery", "all"))
  stopifnot("`level` must be a single number in (0, 1)" =
              is.numeric(level) && length(level) == 1L && is.finite(level) &&
              level > 0 && level < 1)
  if (identical(what, "all")) {
    return(stats::setNames(lapply(.growth_tables(x), function(w) {
      .growth_table(x, w, level, vcov_type, time)
    }), .growth_tables(x)))
  }
  .growth_table(x, what, level, vcov_type, time, data, truth)
}

#' @rdname get_results.latents_growth_mixture
#' @param row.names,optional Unused; for the generic.
#' @export
as.data.frame.latents_growth_mixture <- function(x, row.names = NULL, optional = FALSE,
                                                 what = "coefficients", ...) {
  get_results(x, what = what, ...)
}

#' Print a growth mixture model
#' @param x A `latents_growth_mixture` fit.
#' @param digits Significant digits.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.latents_growth_mixture <- function(x, digits = 4L, ...) {
  fit <- .growth_fit_table(x)
  random <- paste("random", .growth_random_description(x$spec))
  covariance <- switch(fit$random_covariance, varying = "class-specific",
                       equal = "shared across classes",
                       proportional = "proportional across classes")
  if (is.null(x$spec$cluster)) {
    cat("Growth mixture model\n")
    cat(sprintf("  %d classes, %d persons, %d observations\n", fit$n_classes,
                fit$n_persons, fit$n_observations))
  } else {
    cat("Multilevel growth mixture model\n")
    cat(sprintf("  %d classes and %d group classes; %d clusters, %d persons, %d observations\n",
                fit$n_classes, fit$n_group_classes, fit$n_clusters, fit$n_persons,
                fit$n_observations))
  }
  cat(sprintf("  %s%s, covariance %s; residual variance %s\n",
              toupper(substr(random, 1L, 1L)), substring(random, 2L), covariance,
              if (identical(fit$residual_variance, "equal")) "shared" else "by class"))
  .latents_print_weights(x)
  cat(sprintf("  Log likelihood %.2f, BIC %.2f, entropy %.2f, %s\n", fit$log_likelihood,
              fit$bic, fit$entropy, if (fit$converged) "converged" else "NOT converged"))
  cat("\n")
  print(get_results(x, "classes"))
  if (!is.null(x$spec$cluster)) {
    cat("\n")
    print(get_results(x, "group_classes"))
  }
  cat("\n")
  print(get_results(x, "coefficients"))
  if (isTRUE(x$random_boundary)) {
    cat("\nNote: a random-effect covariance is degenerate (see the warning at fit).\n")
  }
  cat("\nMore: get_results(x, \"random\"), plot(x), summary(x)\n")
  invisible(x)
}

#' Summarize a growth mixture model
#' @param object A `latents_growth_mixture` fit.
#' @param level Confidence level.
#' @param vcov_type As in [get_results.latents_growth_mixture()].
#' @param ... Unused.
#' @return An object of class `summary_latents_growth_mixture`: a list of the
#'   `fit`, `classes`, `coefficients`, `random` and `membership` tables.
#' @export
summary.latents_growth_mixture <- function(object, level = 0.95, vcov_type = NULL, ...) {
  tables <- c("fit", "classes", if (!is.null(object$spec$cluster)) "group_classes",
              "coefficients", "random", "membership")
  structure(stats::setNames(lapply(tables, function(w) {
    .growth_table(object, w, level, vcov_type)
  }), tables), class = "summary_latents_growth_mixture")
}

#' @rdname summary.latents_growth_mixture
#' @param x A `summary_latents_growth_mixture` object.
#' @param digits Significant digits.
#' @export
print.summary_latents_growth_mixture <- function(x, digits = 4L, ...) {
  invisible(lapply(x, function(table) {
    print(table)
    cat("\n")
  }))
  invisible(x)
}

#' @rdname get_results.latents_growth_mixture
#' @param object A `latents_growth_mixture` fit.
#' @export
coef.latents_growth_mixture <- function(object, ...) .growth_pack(object$spec, object$params)

#' @rdname get_results.latents_growth_mixture
#' @param type As `vcov_type`.
#' @export
vcov.latents_growth_mixture <- function(object, type = NULL, ...) {
  .growth_resolve_inference(object, type)$vcov
}

#' @rdname get_results.latents_growth_mixture
#' @param parm Unused.
#' @export
confint.latents_growth_mixture <- function(object, parm, level = 0.95, ...) {
  inference <- .growth_resolve_inference(object)
  table <- .inference_wald(unname(inference$theta),
                         sqrt(pmax(diag(inference$vcov), 0)), level)
  data.frame(parameter = names(inference$theta), table[c("estimate", "conf_low", "conf_high")],
             stringsAsFactors = FALSE)
}

#' @rdname get_results.latents_growth_mixture
#' @export
logLik.latents_growth_mixture <- function(object, ...) {
  structure(object$log_likelihood, df = object$n_parameters,
            nobs = length(object$stats$n), class = "logLik")
}

#' @rdname get_results.latents_growth_mixture
#' @export
nobs.latents_growth_mixture <- function(object, ...) length(object$stats$n)
