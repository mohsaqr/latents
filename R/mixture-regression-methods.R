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
  modal <- max.col(unit_tau, ties.method = "first")
  class_names <- paste0("class_", seq_len(spec$n_classes))
  table <- data.frame(
    class = class_names,
    share = colSums(weighted_tau) / sum(weighted_tau),
    count = colSums(weighted_tau),
    n_assigned = tabulate(modal, spec$n_classes),
    mean_posterior = vapply(seq_len(spec$n_classes), function(k) {
      assigned <- modal == k
      if (any(assigned)) mean(unit_tau[assigned, k]) else NA_real_
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
.mixture_relative_entropy <- function(posterior) {
  if (ncol(posterior) < 2L) return(NA_real_)
  positive <- posterior > 0
  entropy <- -sum(posterior[positive] * log(posterior[positive]))
  1 - entropy / (nrow(posterior) * log(ncol(posterior)))
}

#' Fit-statistics table
#' @noRd
.mixture_fit_table <- function(x) {
  spec <- x$spec
  unit_posterior <- switch(spec$nesting,
    observation = x$expectation$tau,
    group = x$expectation$group_tau,
    "two-level" = x$expectation$tau)
  n_units <- if (identical(spec$nesting, "observation")) spec$n else spec$n_groups
  k <- x$n_parameters
  log_lik <- x$log_likelihood
  entropy_sum <- -sum(unit_posterior[unit_posterior > 0] *
                        log(unit_posterior[unit_posterior > 0]))
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
    bic_rows = -2 * log_lik + k * log(spec$n),
    sabic = -2 * log_lik + k * log((n_units + 2) / 24),
    icl = bic + 2 * entropy_sum,
    entropy = .mixture_relative_entropy(unit_posterior),
    group_entropy = if (identical(spec$nesting, "two-level"))
      .mixture_relative_entropy(x$expectation$rho) else NA_real_,
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
  assigned <- max.col(posterior, ties.method = "first")
  k <- spec$n_classes
  means <- rowsum(posterior, factor(assigned, levels = seq_len(k))) /
    pmax(tabulate(assigned, k), 1)
  means[tabulate(assigned, k) == 0L, ] <- NA_real_
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
    if (any(values != values[first][spec$group_index])) {
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
    truth = values), responseName = "n", stringsAsFactors = FALSE)
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

#' Tables of a mixture-of-regressions fit
#'
#' @param x A fit from [mixture_regression()].
#' @param what The table to return; see *Tables*. `"all"` returns every
#'   table in a named list.
#' @param level Confidence level for the intervals.
#' @param vcov_type `NULL` for the covariance stored with the fit, or
#'   `"observed"`, `"robust"` or `"opg"` to recompute it.
#' @param data,truth For `what = "recovery"`: the data frame the model was
#'   fitted to and the name of its column holding a known classification.
#' @param by For `what = "recovery"`: `"class"`, or `"group_class"` for the
#'   group classes of a two-level fit.
#' @param ... Unused.
#'
#' @section Tables:
#' \describe{
#'   \item{`coefficients`}{One row per class and regression term: `class`
#'     (`"common"` for a shared term), `term`, `estimate`, `std_error`,
#'     `statistic` (Wald z), `p_value`, `conf_low`, `conf_high`,
#'     `p_adjusted` (Benjamini-Hochberg across the non-intercept rows;
#'     `NA` for intercepts) and, for the binomial and Poisson families,
#'     `exp_estimate` (odds or rate ratio). For the ordinal family each
#'     class's thresholds lead its rows as terms `"threshold:<lower>|<upper>"`
#'     (delta-method errors; no `p_adjusted` or `exp_estimate`), and
#'     `exp_estimate` of a slope is the cumulative odds ratio of a higher
#'     category.}
#'   \item{`classes`}{One row per class: `share` (mean posterior), `count`
#'     (summed posterior), `n_assigned` (modal assignment), `mean_posterior`
#'     (average posterior of the units assigned to it); `sigma` and
#'     `sigma_std_error` for the Gaussian family; `prior` and
#'     `prior_std_error` when membership has no covariates. Units are groups
#'     under `class_level = "group"`.}
#'   \item{`membership`}{One row per multinomial-logit coefficient: `model`
#'     (`"class"` or `"group_class"`), `class`, `term`, `reference`, the Wald
#'     columns and `odds_ratio`. In the two-level model the class intercepts
#'     are per group class, `(Intercept):group_class_h`.}
#'   \item{`group_classes`}{Two-level model: one row per group class and
#'     class, with the model-implied class `probability`, the `group_share`,
#'     `n_groups_assigned`, and `std_error` of the probability when
#'     membership has no covariates.}
#'   \item{`fit`}{One row: log likelihood, parameter count, `aic`, `bic`
#'     (penalized by the number of independent units: rows, or groups when
#'     `id` defines them), `bic_rows` (always by rows, as flexmix reports it),
#'     `sabic`, `icl`, relative `entropy` (and `group_entropy`), the smallest
#'     class share, convergence and start replication.}
#'   \item{`assignments`}{One row per data row used: `row` (its position in
#'     the supplied data), the `id` column, modal `class`, its `posterior`,
#'     and one `probability_class_k` column per class; two-level fits add
#'     `group_class` and `group_posterior`.}
#'   \item{`groups`}{One row per group (`class_level = "group"` or the
#'     two-level model): its size, modal class or group class, and
#'     posterior probabilities.}
#'   \item{`fitted`}{One row per data row: `observed`, the posterior-weighted
#'     `fitted` mean, `fitted_modal` (the modal class's mean), `residual`
#'     and `class`.}
#'   \item{`starts`}{One row per start: its log likelihood, convergence,
#'     iterations, `stage` (`"screened"`: stopped after the 50-iteration
#'     screening because better starts existed; `"completed"`: run to
#'     convergence or `max_iter`), whether it degenerated, and which one was
#'     selected. `n_best_replicated` counts completed starts only.}
#'   \item{`classification`}{One row per assigned (modal) class and class:
#'     the `mean_posterior` probability of `class` among the units assigned to
#'     `assigned`, and `n_assigned`. Its diagonal is each class's average
#'     posterior probability; values near one mean well-separated classes.}
#'   \item{`recovery`}{Needs `data` and `truth`: the modal classes
#'     cross-tabulated against a known classification, one row per assigned
#'     class and true value, with the count `n` and its `share` of the
#'     assigned class. Group-level under `class_level = "group"` or
#'     `by = "group_class"`, where `truth` must be constant within groups.}
#' }
#' @return A base `data.frame` as described under *Tables*, or a named list
#'   of them for `what = "all"`.
#' @examples
#' fit <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2,
#'                           n_starts = 3, seed = 1)
#' get_results(fit, "coefficients")
#' get_results(fit, "fit")
#' get_results(fit, "coefficients", vcov_type = "robust")
#' get_results(fit, "recovery", data = study_hours, truth = "strategy")
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

#' @rdname get_results.latents_mixture_regression
#' @param row.names,optional Unused; part of the generic.
#' @export
as.data.frame.latents_mixture_regression <- function(x, row.names = NULL, optional = FALSE,
                                         what = "coefficients", level = 0.95,
                                         vcov_type = NULL, ...) {
  get_results.latents_mixture_regression(x, what = what, level = level,
                                         vcov_type = vcov_type)
}

#' Print a mixture-of-regressions fit
#' @param x A fit from [mixture_regression()].
#' @param digits Significant digits.
#' @param ... Unused.
#' @return `x`, invisibly.
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

#' Summarize a mixture-of-regressions fit
#' @param object A fit from [mixture_regression()].
#' @param level Confidence level.
#' @param vcov_type As in [get_results.latents_mixture_regression()].
#' @param ... Unused.
#' @return An object of class `summary_latents_mixture_regression`: a named list of the
#'   `fit`, `coefficients`, `classes`, `membership` and `group_classes`
#'   tables, printed by its print method.
#' @export
summary.latents_mixture_regression <- function(object, level = 0.95, vcov_type = NULL, ...) {
  tables <- c("fit", "coefficients", "classes", "membership", "group_classes")
  structure(stats::setNames(lapply(tables, function(name) {
    .mixture_table(object, name, level, vcov_type)
  }), tables), class = "summary_latents_mixture_regression")
}

#' @rdname summary.latents_mixture_regression
#' @param x A `summary_latents_mixture_regression` object.
#' @param digits Significant digits.
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

#' @rdname get_results.latents_mixture_regression
#' @param object A fit from [mixture_regression()].
#' @export
coef.latents_mixture_regression <- function(object, ...) {
  .mixture_pack(object$spec, object$params)
}

#' @rdname get_results.latents_mixture_regression
#' @param type Covariance type for `vcov()`; as `vcov_type`.
#' @export
vcov.latents_mixture_regression <- function(object, type = NULL, ...) {
  .mixture_resolve_inference(object, type)$vcov
}

#' @rdname get_results.latents_mixture_regression
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
  if (!missing(parm)) bounds <- bounds[parm, , drop = FALSE]
  bounds
}

#' @rdname get_results.latents_mixture_regression
#' @export
logLik.latents_mixture_regression <- function(object, ...) {
  units <- if (identical(object$spec$nesting, "observation")) object$spec$n else
    object$spec$n_groups
  structure(object$log_likelihood, df = object$n_parameters, nobs = units,
            class = "logLik")
}

#' @rdname get_results.latents_mixture_regression
#' @export
nobs.latents_mixture_regression <- function(object, ...) {
  if (identical(object$spec$nesting, "observation")) object$spec$n else
    object$spec$n_groups
}
