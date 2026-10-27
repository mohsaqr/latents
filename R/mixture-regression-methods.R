# Tables, accessors and S3 methods for mixture_regression() fits.

#' The table names a mixture_regression fit offers
#' @noRd
.mixture_tables <- function() {
  c("coefficients", "classes", "membership", "group_classes", "fit",
    "assignments", "groups", "fitted", "starts", "classification",
    "recovery")
}

#' Inference at the requested covariance type, reusing the stored one
#' @noRd
.mixture_resolve_inference <- function(x, vcov_type) {
  if (is.null(vcov_type)) {
    if (!is.null(x$inference)) return(x$inference)
    vcov_type <- if (.latents_is_weighted(x)) "robust" else "observed"
  }
  vcov_type <- .latents_weighted_vcov(.latents_is_weighted(x),
    match.arg(vcov_type, c("observed", "robust", "opg")), FALSE)
  if (!is.null(x$inference) && identical(x$inference$vcov_type, vcov_type)) {
    return(x$inference)
  }
  .mixture_inference(x, vcov_type)
}

#' Mean response of each class on the response scale
#'
#' For the ordinal family, the expected category score (categories scored
#' 1..C in order).
#' @noRd
.mixture_class_means <- function(spec, params) {
  eta <- .mixture_linear_predictors(spec, params)
  switch(spec$family,
    gaussian = eta,
    binomial = spec$trials * stats::plogis(eta),
    poisson = ,
    negative_binomial = exp(eta),
    ordinal = .ordinal_expected_scores(eta, params$thresholds))
}

#' Coefficient table
#' @noRd
.mixture_coefficient_table <- function(x, level, inference) {
  spec <- x$spec
  theta <- inference$theta
  std_error <- sqrt(pmax(diag(inference$vcov), 0))
  names(std_error) <- names(theta)
  regression <- grepl("^coefficient\\.", names(theta))
  pieces <- strsplit(sub("^coefficient\\.", "", names(theta)[regression]),
                     ".", fixed = TRUE)
  class_label <- vapply(pieces, `[`, character(1), 1L)
  term <- vapply(pieces, function(piece) paste(piece[-1L], collapse = "."),
                 character(1))
  table <- cbind(
    data.frame(class = class_label, term = term),
    .inference_wald(theta[regression], std_error[regression], level))
  ordinal <- identical(spec$family, "ordinal")
  if (ordinal) {
    # Each class's thresholds lead its rows, as an intercept would.
    table <- rbind(.ordinal_threshold_table(spec, inference, level), table)
    rank <- match(table$class,
                  c(paste0("class_", seq_len(spec$n_classes)), "common"))
    table <- table[order(rank, seq_len(nrow(table))), , drop = FALSE]
    rownames(table) <- NULL
  }
  thresholds <- startsWith(table$term, "threshold:") & ordinal
  slopes <- table$term != "(Intercept)" & !thresholds
  table$p_adjusted <- NA_real_
  table$p_adjusted[slopes] <- stats::p.adjust(table$p_value[slopes],
                                              method = "BH")
  if (!identical(spec$family, "gaussian")) {
    table$exp_estimate <- exp(table$estimate)
    table$exp_estimate[thresholds] <- NA_real_
  }
  table
}

#' Class table: sizes, residual spread and classification quality
#' @noRd
.mixture_class_table <- function(x, level, inference) {
  spec <- x$spec
  tau <- x$expectation$tau
  unit_tau <- if (identical(spec$nesting, "group")) x$expectation$group_tau else tau
  unit_weights <- if (identical(spec$nesting, "group")) spec$sampling_weights else
    spec$row_weights
  weighted_tau <- if (is.null(unit_weights)) unit_tau else unit_tau * unit_weights
  # A unit of weight zero is not in the sample (ML-1): it is assigned to no
  # class (review 2026-10-26, round 2, N7).
  if (!is.null(unit_weights)) {
    kept <- unit_weights > 0
    unit_tau <- unit_tau[kept, , drop = FALSE]
    unit_weights <- unit_weights[kept]
  }
  modal <- max.col(unit_tau, ties.method = "first")
  class_names <- paste0("class_", seq_len(spec$n_classes))
  table <- data.frame(
    class = class_names,
    share = colSums(weighted_tau) / sum(weighted_tau),
    count = colSums(weighted_tau),
    n_assigned = tabulate(modal, spec$n_classes),
    mean_posterior = vapply(seq_len(spec$n_classes), function(k) {
      assigned <- modal == k
      if (any(assigned)) stats::weighted.mean(unit_tau[assigned, k],
        (unit_weights %||% rep(1, nrow(unit_tau)))[assigned]) else NA_real_
    }, numeric(1)),
    row.names = NULL)
  if (identical(spec$family, "gaussian")) {
    sigma_names <- if (identical(spec$variance, "equal")) {
      rep("log_sigma.all", spec$n_classes)
    } else paste0("log_sigma.", class_names)
    theta <- inference$theta
    delta <- .mixture_delta(theta, inference$vcov, function(t) {
      exp(t[sigma_names])
    })
    table$sigma <- delta$estimate
    table$sigma_std_error <- delta$std_error
  }
  if (identical(spec$family, "negative_binomial")) {
    dispersion_names <- if (identical(spec$variance, "equal")) {
      rep("log_dispersion.all", spec$n_classes)
    } else paste0("log_dispersion.", class_names)
    delta <- .mixture_delta(inference$theta, inference$vcov, function(t) {
      exp(t[dispersion_names])
    })
    table$dispersion <- delta$estimate
    table$dispersion_std_error <- delta$std_error
  }
  if (!identical(spec$nesting, "two-level") && spec$intercept_only) {
    delta <- .mixture_delta(inference$theta, inference$vcov, function(t) {
      params <- .mixture_unpack(spec, t, x$params)
      drop(exp(.mixture_log_softmax(spec$w[1L, , drop = FALSE],
                                   params$gamma)))
    })
    table$prior <- delta$estimate
    table$prior_std_error <- delta$std_error
  }
  table
}

#' Membership logit table
#' @noRd
.mixture_membership_table <- function(x, level, inference) {
  theta <- inference$theta
  std_error <- sqrt(pmax(diag(inference$vcov), 0))
  names(std_error) <- names(theta)
  mixing <- grepl("^(group_)?membership\\.", names(theta))
  if (!any(mixing)) {
    return(data.frame(model = character(), class = character(),
                      term = character(), estimate = numeric(),
                      std_error = numeric(), statistic = numeric(),
                      p_value = numeric(), conf_low = numeric(),
                      conf_high = numeric(), odds_ratio = numeric()))
  }
  labels <- names(theta)[mixing]
  model <- ifelse(startsWith(labels, "group_membership."), "group_class",
                  "class")
  rest <- sub("^(group_)?membership\\.", "", labels)
  class_label <- sub("\\..*$", "", rest)
  term <- sub("^[^.]*\\.", "", rest)
  table <- cbind(data.frame(model = model, class = class_label, term = term,
                            reference = ifelse(model == "class", "class_1",
                                               "group_class_1")),
                 .inference_wald(theta[mixing], std_error[mixing], level))
  table$odds_ratio <- exp(table$estimate)
  table
}

#' Group-class table for the two-level model
#' @noRd
.mixture_group_class_table <- function(x, level, inference) {
  spec <- x$spec
  if (!identical(spec$nesting, "two-level")) {
    return(data.frame(group_class = character(), class = character(),
                      probability = numeric(), group_share = numeric()))
  }
  rho <- x$expectation$rho
  weighted_rho <- if (is.null(spec$sampling_weights)) rho else rho * spec$sampling_weights
  group_names <- paste0("group_class_", seq_len(spec$n_group_classes))
  class_names <- paste0("class_", seq_len(spec$n_classes))
  # Model-implied class probabilities within each group class, averaged over
  # the rows the group class is responsible for.
  probability <- vapply(seq_len(spec$n_group_classes), function(h) {
    prior <- exp(x$expectation$log_prior_by_group_class[[h]])
    weight <- weighted_rho[spec$group_index, h]
    colSums(prior * weight) / sum(weight)
  }, numeric(spec$n_classes))
  probability <- matrix(probability, spec$n_classes)
  table <- data.frame(
    group_class = rep(group_names, each = spec$n_classes),
    class = rep(class_names, spec$n_group_classes),
    probability = as.vector(probability),
    group_share = rep(colSums(weighted_rho) / sum(weighted_rho), each = spec$n_classes),
    n_groups_assigned = rep(tabulate(max.col(rho, ties.method = "first"),
                                     spec$n_group_classes),
                            each = spec$n_classes))
  if (ncol(spec$w) == 0L) {
    delta <- .mixture_delta(inference$theta, inference$vcov, function(t) {
      params <- .mixture_unpack(spec, t, x$params)
      as.vector(t(exp(.mixture_log_softmax(
        diag(1, spec$n_group_classes), params$class_logits))))
    })
    table$std_error <- delta$std_error
  }
  table
}

#' Relative entropy of a posterior matrix (1 = perfectly separated)
#' @noRd
.mixture_relative_entropy <- function(posterior, weights = NULL) {
  .multilpa_relative_entropy(posterior, weights)
}

#' Complete classification entropy of a regression mixture
#' @param x A regression mixture fit.
#' @return Sampling-weighted joint entropy of every latent class.
#' @noRd
.mixture_complete_entropy <- function(x) {
  stopifnot("`x` must be a fitted regression mixture" =
              inherits(x, "latents_mixture_regression"))
  spec <- x$spec
  if (identical(spec$nesting, "two-level")) {
    within <- .latents_row_conditional_entropy(
      x$expectation$tau_by_group_class, x$expectation$rho, spec$group_index)
    return(.multilpa_entropy_sum(x$expectation$rho, spec$sampling_weights) +
             sum(within * (spec$row_weights %||% 1)))
  }
  posterior <- if (identical(spec$nesting, "group")) x$expectation$group_tau else
    x$expectation$tau
  weights <- if (identical(spec$nesting, "group")) spec$sampling_weights else
    spec$row_weights
  .multilpa_entropy_sum(posterior, weights)
}

#' Fit-statistics table
#' @noRd
.mixture_fit_table <- function(x) {
  spec <- x$spec
  unit_posterior <- switch(spec$nesting,
    observation = x$expectation$tau,
    group = x$expectation$group_tau,
    "two-level" = x$expectation$tau)
  sample <- .mixture_weighted_sample(spec)
  n_units <- sample$units
  k <- x$n_parameters
  log_lik <- x$log_likelihood
  # Classification likelihood and reporting use the same sampling weights
  # as the fitted likelihood; two-level ICL includes both latent levels.
  entropy_sum <- .mixture_complete_entropy(x)
  bic <- -2 * log_lik + k * log(n_units)
  unit_weights <- if (identical(spec$nesting, "group")) spec$sampling_weights else
    spec$row_weights
  weighted_posterior <- unit_posterior * (unit_weights %||% 1)
  data.frame(
    family = spec$family, nesting = spec$nesting,
    n_classes = spec$n_classes, n_group_classes = spec$n_group_classes,
    n_observations = spec$n,
    n_groups = if (is.null(spec$id)) NA_integer_ else spec$n_groups,
    log_likelihood = log_lik, n_parameters = k,
    aic = -2 * log_lik + 2 * k,
    bic = bic,
    bic_rows = -2 * log_lik + k * log(sample$rows),
    sabic = -2 * log_lik + k * log((n_units + 2) / 24),
    icl = bic + 2 * entropy_sum,
    entropy = .mixture_relative_entropy(unit_posterior, unit_weights),
    group_entropy = if (identical(spec$nesting, "two-level"))
      .mixture_relative_entropy(x$expectation$rho, spec$sampling_weights) else NA_real_,
    smallest_share = min(colSums(weighted_posterior) / sum(weighted_posterior)),
    converged = x$converged, iterations = x$iterations,
    n_starts = nrow(x$starts), n_best_replicated = x$n_best_replicated,
    vcov_type = x$inference$vcov_type %||% NA_character_)
}

#' Per-row assignment table
#' @noRd
.mixture_assignment_table <- function(x) {
  spec <- x$spec
  tau <- x$expectation$tau
  colnames(tau) <- paste0("probability_class_", seq_len(spec$n_classes))
  table <- data.frame(row = spec$kept_rows)
  if (!is.null(spec$id)) table[[spec$id]] <- spec$group_levels[spec$group_index]
  table$class <- max.col(tau, ties.method = "first")
  table$posterior <- tau[cbind(seq_len(spec$n), table$class)]
  if (identical(spec$nesting, "two-level")) {
    rho <- x$expectation$rho
    group_class <- max.col(rho, ties.method = "first")
    table$group_class <- group_class[spec$group_index]
    table$group_posterior <- rho[cbind(spec$group_index,
                                       table$group_class)]
  }
  cbind(table, as.data.frame(tau))
}

#' Per-group table
#' @noRd
.mixture_group_table <- function(x) {
  spec <- x$spec
  if (is.null(spec$id) || identical(spec$nesting, "observation")) {
    return(data.frame(group = character(), n_rows = integer()))
  }
  posterior <- if (identical(spec$nesting, "group")) x$expectation$group_tau else
    x$expectation$rho
  prefix <- if (identical(spec$nesting, "group")) "class" else "group_class"
  colnames(posterior) <- paste0("probability_", prefix, "_",
                                seq_len(ncol(posterior)))
  table <- data.frame(group = spec$group_levels,
                      n_rows = tabulate(spec$group_index, spec$n_groups))
  names(table)[1L] <- spec$id
  table[[prefix]] <- max.col(posterior, ties.method = "first")
  table$posterior <- posterior[cbind(seq_len(nrow(posterior)),
                                     table[[prefix]])]
  cbind(table, as.data.frame(posterior))
}

#' Fitted values table
#' @noRd
.mixture_fitted_table <- function(x) {
  spec <- x$spec
  means <- .mixture_class_means(spec, x$params)
  tau <- x$expectation$tau
  modal <- max.col(tau, ties.method = "first")
  observed <- spec$y
  fitted <- rowSums(tau * means)
  data.frame(row = spec$kept_rows, observed = observed,
             fitted = fitted, fitted_modal = means[cbind(seq_len(spec$n), modal)],
             residual = observed - fitted, class = modal)
}

#' Average posterior probabilities by modal class
#'
#' The classification-quality matrix of latent class reporting: for the units
#' assigned to each class, the mean posterior probability of every class. A
#' diagonal near one means the classes are well separated.
#' @noRd
.mixture_classification_table <- function(x) {
  spec <- x$spec
  posterior <- if (identical(spec$nesting, "group")) x$expectation$group_tau else
    x$expectation$tau
  unit_weights <- if (identical(spec$nesting, "group")) spec$sampling_weights else
    spec$row_weights
  # Units of weight zero are not in the sample (ML-1, round 2 N7).
  if (!is.null(unit_weights)) {
    kept <- unit_weights > 0
    posterior <- posterior[kept, , drop = FALSE]
    unit_weights <- unit_weights[kept]
  }
  assigned <- max.col(posterior, ties.method = "first")
  k <- spec$n_classes
  sizes <- tabulate(assigned, k)
  occupied <- which(sizes > 0L)
  means <- matrix(NA_real_, k, k)
  weights <- unit_weights %||% rep(1, nrow(posterior))
  means[occupied, ] <- rowsum(posterior * weights, assigned) /
    as.vector(rowsum(weights, assigned))
  data.frame(assigned = rep(paste0("class_", seq_len(k)), times = k),
             class = rep(paste0("class_", seq_len(k)), each = k),
             mean_posterior = as.vector(means),
             n_assigned = rep(tabulate(assigned, k), times = k))
}

#' Cross-tabulate assignments against a known classification
#'
#' @param x A fit.
#' @param data The data frame the fit was made from, carrying `truth`.
#' @param truth Name of the column holding the known classification.
#' @param by `"class"` or `"group_class"`.
#' @return A long data frame: `assigned`, the truth column, `n`, `share` of
#'   the assigned class.
#' @noRd
.mixture_recovery_table <- function(x, data, truth, by) {
  spec <- x$spec
  if (is.null(data) || is.null(truth)) {
    stop(errorCondition(paste(
      "`what = \"recovery\"` needs `data`, the data frame the model was",
      "fitted to, and `truth`, the column holding the known classes."),
      class = "latents_bad_argument", call = NULL))
  }
  stopifnot(
    "`data` must be a data frame" = is.data.frame(data),
    "`truth` must be a single column name of `data`" =
      is.character(truth) && length(truth) == 1L && truth %in% names(data))
  if (max(spec$kept_rows) > nrow(data)) {
    stop(errorCondition(paste(
      "`data` has fewer rows than the fit used; pass the data frame the",
      "model was fitted to."), class = "latents_bad_data", call = NULL))
  }
  values <- data[[truth]][spec$kept_rows]
  group_level <- identical(by, "group_class") ||
    identical(spec$nesting, "group")
  if (identical(by, "group_class") && !identical(spec$nesting, "two-level")) {
    stop(errorCondition("`by = \"group_class\"` needs a two-level fit.",
                        class = "latents_bad_argument", call = NULL))
  }
  if (group_level) {
    first <- match(seq_len(spec$n_groups), spec$group_index)
    reference <- values[first][spec$group_index]
    # A missing known class is a value of its own: it must be missing for the
    # whole group, and is tabulated (as NA) rather than dropped.
    same <- (is.na(values) & is.na(reference)) |
      (!is.na(values) & !is.na(reference) & values == reference)
    if (!all(same)) {
      stop(errorCondition(sprintf(paste(
        "`%s` varies within groups, so it cannot be compared with a",
        "group-level classification."), truth),
        class = "latents_bad_data", call = NULL))
    }
    values <- values[first]
    posterior <- if (identical(by, "group_class")) x$expectation$rho else
      x$expectation$group_tau
  } else {
    posterior <- x$expectation$tau
  }
  assigned <- max.col(posterior, ties.method = "first")
  prefix <- if (identical(by, "group_class")) "group_class_" else "class_"
  counts <- as.data.frame(table(
    assigned = factor(paste0(prefix, assigned),
                      levels = paste0(prefix, seq_len(ncol(posterior)))),
    truth = values, useNA = "ifany"), responseName = "n", stringsAsFactors = FALSE)
  names(counts)[names(counts) == "truth"] <- truth
  totals <- tapply(counts$n, counts$assigned, sum)
  counts$share <- counts$n / pmax(totals[counts$assigned], 1)
  counts
}

#' Build one table by name
#' @noRd
.mixture_table <- function(x, what, level = 0.95, vcov_type = NULL,
                          data = NULL, truth = NULL, by = "class") {
  needs_inference <- what %in% c("coefficients", "classes", "membership",
                                 "group_classes")
  inference <- if (needs_inference) .mixture_resolve_inference(x, vcov_type)
  switch(what,
    coefficients = .mixture_coefficient_table(x, level, inference),
    classes = .mixture_class_table(x, level, inference),
    membership = .mixture_membership_table(x, level, inference),
    group_classes = .mixture_group_class_table(x, level, inference),
    fit = .mixture_fit_table(x),
    assignments = .mixture_assignment_table(x),
    groups = .mixture_group_table(x),
    fitted = .mixture_fitted_table(x),
    starts = x$starts,
    classification = .mixture_classification_table(x),
    recovery = .mixture_recovery_table(x, data, truth, by)) |>
    .mixture_present(what, level, x)
}

#' Attach the printed layout of a mixture-regression table
#'
#' The data are untouched; this chooses the title, columns, labels and
#' number formats `print()` uses. Columns a fit does not have (`sigma` for
#' a binomial fit, the two-level columns) are left out.
#' @noRd
.mixture_present <- function(table, what, level, x) {
  column <- .latents_column
  percent <- format(100 * level)
  interval <- function(low, high, from) {
    column(sprintf("%s%% CI", percent), c(low, high), "ci", from)
  }
  units <- if (identical(x$spec$nesting, "observation")) "Rows" else "Groups"
  plan <- switch(what,
    coefficients = list(
      title = sprintf("Regression coefficients (%s%% CI)", percent),
      display = list(column("Class", "class", "label"), column("Term", "term", "label"),
                     column("Estimate", "estimate"),
                     interval("conf_low", "conf_high", "estimate"),
                     column("p", "p_value", "p"),
                     column(if (x$spec$family %in% c("binomial", "ordinal")) "Odds ratio" else
                       "Rate ratio", "exp_estimate"))),
    classes = list(
      title = "Classes",
      display = list(column("Class", "class", "label"), column("Share", "share", "share"),
                     column("Expected", "count"),
                     column(sprintf("%s assigned", units), "n_assigned", "integer"),
                     column("Avg. posterior", "mean_posterior", "share"),
                     column("Residual SD", "sigma"),
                     column("Dispersion", "dispersion"))),
    membership = list(
      title = sprintf("Class membership (log odds, %s%% CI)", percent),
      display = list(column("Model", "model", "label"), column("Class", "class", "label"),
                     column("Term", "term", "label"), column("Estimate", "estimate"),
                     interval("conf_low", "conf_high", "estimate"),
                     column("p", "p_value", "p"), column("Odds ratio", "odds_ratio"))),
    group_classes = list(
      title = "Group classes",
      display = list(column("Group class", "group_class", "label"),
                     column("Class", "class", "label"),
                     column("Probability", "probability", "share"),
                     column("Group share", "group_share", "share"),
                     column("Groups assigned", "n_groups_assigned", "integer"))),
    fit = list(
      title = "Model fit",
      display = list(column("Family", "family", "text"), column("Nesting", "nesting", "text"),
                     column("Classes", "n_classes", "integer"),
                     column("Group classes", "n_group_classes", "integer"),
                     column("Observations", "n_observations", "integer"),
                     column("Groups", "n_groups", "integer"),
                     column("Parameters", "n_parameters", "integer"),
                     column("Log likelihood", "log_likelihood"), column("AIC", "aic"),
                     column("BIC", "bic"), column("SABIC", "sabic"), column("ICL", "icl"),
                     column("Entropy", "entropy", "share"),
                     column("Smallest class", "smallest_share", "percent"),
                     column("Converged", "converged", "logical"))),
    # The second column is the caller's truth column.
    recovery = list(
      title = "Recovery of a known classification",
      display = list(column("Assigned", "assigned", "label"),
                     column(names(table)[2L], names(table)[2L], "text"),
                     column(units, "n", "integer"),
                     column("Share of assigned", "share", "share"))),
    list(title = switch(what, assignments = "Class assignments",
                        groups = "Groups", fitted = "Fitted values", starts = "Starts",
                        classification = "Classification (average posteriors)", NULL),
         display = NULL))
  display <- Filter(function(spec) all(spec$column %in% names(table)), plan$display)
  if (!identical(x$spec$nesting, "two-level")) {
    display <- Filter(function(spec) !identical(spec$column, "model"), display)
  }
  table <- .latents_table(table, level, plan$title,
                          if (length(display) > 0L) display)
  if (identical(what, "fit")) attr(table, "card") <- TRUE
  table
}

#' @noRd
#' @export
get_results.latents_mixture_regression <- function(x, what = "coefficients", level = 0.95,
                                       vcov_type = NULL, data = NULL,
                                       truth = NULL,
                                       by = c("class", "group_class"), ...) {
  what <- match.arg(what, c(.mixture_tables(), "all"))
  by <- match.arg(by)
  stopifnot("`level` must be a single number in (0, 1)" =
              is.numeric(level) && length(level) == 1L && level > 0 &&
              level < 1)
  if (identical(what, "all")) {
    tables <- setdiff(.mixture_tables(), "recovery")
    return(stats::setNames(lapply(tables, function(name) {
      .mixture_table(x, name, level, vcov_type)
    }), tables))
  }
  .mixture_table(x, what, level, vcov_type, data, truth, by)
}

#' @noRd
#' @param row.names,optional Unused; part of the generic.
#' @export
as.data.frame.latents_mixture_regression <- function(x, row.names = NULL, optional = FALSE,
                                         what = "coefficients", level = 0.95,
                                         vcov_type = NULL, ...) {
  get_results.latents_mixture_regression(x, what = what, level = level,
                                         vcov_type = vcov_type)
}

#' @noRd
#' @export
print.latents_mixture_regression <- function(x, digits = 4L, ...) {
  spec <- x$spec
  fit <- .mixture_fit_table(x)
  nesting <- switch(spec$nesting,
    observation = "one class per row",
    group = sprintf("one class per `%s`", spec$id),
    "two-level" = sprintf("two-level: %d group classes of `%s`",
                          spec$n_group_classes, spec$id))
  cat(sprintf("Mixture of %s regressions: %d classes (%s)\n", spec$family,
              spec$n_classes, nesting))
  cat(sprintf("%s ~ %s\n", spec$response_name,
              deparse1(spec$formula[[3L]])))
  cat(sprintf(paste0("%d rows%s | log likelihood %.", digits, "f | ",
                     "BIC %.2f | entropy %.3f\n"),
              spec$n, if (is.null(spec$id)) "" else
                sprintf(" in %d groups", spec$n_groups),
              fit$log_likelihood, fit$bic, fit$entropy))
  .latents_print_weights(x)
  cat(sprintf(paste("Converged: %s | iterations: %d | best likelihood reached",
                    "by %d of %d completed starts (%d run)\n\n"),
              fit$converged, fit$iterations, fit$n_best_replicated,
              sum(x$starts$stage == "completed"), fit$n_starts))
  if (!is.null(x$inference)) {
    print(.mixture_table(x, "classes"))
    cat("\n")
    print(.mixture_table(x, "coefficients"))
  } else if (identical(spec$family, "ordinal")) {
    print(rbind(x$params$thresholds, x$params$beta), digits = digits)
  } else {
    print(x$params$beta, digits = digits)
  }
  cat("\nEvery table: get_results(x, what = ), e.g. \"classes\", \"membership\", \"fit\", \"assignments\".\n")
  invisible(x)
}

#' @noRd
#' @export
summary.latents_mixture_regression <- function(object, level = 0.95, vcov_type = NULL, ...) {
  tables <- c("fit", "coefficients", "classes", "membership", "group_classes")
  structure(stats::setNames(lapply(tables, function(name) {
    .mixture_table(object, name, level, vcov_type)
  }), tables), class = "summary_latents_mixture_regression")
}

#' @noRd
#' @export
print.summary_latents_mixture_regression <- function(x, digits = 4L, ...) {
  # Each table carries its own title and layout.
  invisible(lapply(x, function(table) {
    if (nrow(table) == 0L) return(NULL)
    print(table)
    cat("\n")
  }))
  invisible(x)
}

#' @noRd
#' @param object A fit from [mixture_regression()].
#' @export
coef.latents_mixture_regression <- function(object, ...) {
  .mixture_pack(object$spec, object$params)
}

#' @noRd
#' @param type Covariance type for `vcov()`; as `vcov_type`.
#' @export
vcov.latents_mixture_regression <- function(object, type = NULL, ...) {
  .mixture_resolve_inference(object, type)$vcov
}

#' @noRd
#' @param parm Parameter names or positions; all by default.
#' @export
confint.latents_mixture_regression <- function(object, parm, level = 0.95, ...) {
  inference <- .mixture_resolve_inference(object, NULL)
  theta <- inference$theta
  std_error <- sqrt(pmax(diag(inference$vcov), 0))
  z <- .inference_critical(level, "one_minus")
  bounds <- cbind(theta - z * std_error, theta + z * std_error)
  percent <- paste(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2),
                          trim = TRUE, scientific = FALSE, digits = 3), "%")
  dimnames(bounds) <- list(names(theta), percent)
  if (!missing(parm)) bounds <- bounds[.mixture_confint_rows(parm, names(theta)), , drop = FALSE]
  bounds
}

#' Rows of a `confint()` result picked by `parm`
#' @param parm Parameter names, or positions among `names`.
#' @param names The names of the estimates, in order.
#' @return Integer positions; raises `latents_bad_argument` for a name or
#'   position that is not an estimate.
#' @noRd
.mixture_confint_rows <- function(parm, names) {
  rows <- if (is.character(parm)) match(parm, names) else if (is.numeric(parm) &&
    all(is.finite(parm)) && all(abs(parm - round(parm)) < sqrt(.Machine$double.eps)))
    as.integer(round(parm)) else NA_integer_
  if (length(rows) == 0L || anyNA(rows) || any(rows < 1L | rows > length(names))) {
    stop(errorCondition(paste(
      "`parm` must name estimates of the fit (see names(coef(x))) or give",
      "their positions."), class = "latents_bad_argument", call = NULL))
  }
  rows
}

#' @noRd
#' @export
logLik.latents_mixture_regression <- function(object, ...) {
  structure(object$log_likelihood, df = object$n_parameters,
            nobs = .mixture_weighted_sample(object$spec)$units, class = "logLik")
}

#' The units and rows that carry sampling weight
#'
#' A unit of weight zero contributes nothing to the pseudo likelihood, so it is
#' in no sample size, criterion or entropy (review 2026-10-26, ML-1), as
#' `nobs.glm` counts only non-zero weights. The data counts
#' (`n_observations`, `n_groups`) still count it.
#'
#' @param spec A regression mixture specification.
#' @return A list with `units` and `rows` (counts) and `unit_kept` (one
#'   logical per unit of the likelihood: a row in the observation nesting, a
#'   group otherwise).
#' @noRd
.mixture_weighted_sample <- function(spec) {
  n_units <- if (identical(spec$nesting, "observation")) spec$n else spec$n_groups
  row_kept <- if (is.null(spec$row_weights)) rep(TRUE, spec$n) else spec$row_weights > 0
  # The observation nesting's units are its rows, with or without `id`: the
  # weights are then given per `id` group but carried by each row (review
  # 2026-10-26, round 2, N3).
  unit_kept <- if (identical(spec$nesting, "observation")) row_kept else
    if (is.null(spec$sampling_weights)) rep(TRUE, n_units) else spec$sampling_weights > 0
  list(units = sum(unit_kept), rows = sum(row_kept), unit_kept = unit_kept)
}

#' @noRd
#' @export
nobs.latents_mixture_regression <- function(object, ...) {
  .mixture_weighted_sample(object$spec)$units
}
