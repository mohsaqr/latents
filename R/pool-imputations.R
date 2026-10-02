# Fitting a model to multiply imputed data and pooling it with Rubin's rules.
#
# Missing indicators are integrated out by `missing = "fiml"`; missing
# covariates cannot be, because the model conditions on them and has no
# distribution to integrate over. Multiple imputation is the standard answer:
# impute the covariates (with the indicators in the imputation model), fit each
# completed data set, and pool. Pooling a mixture has one step a regression
# does not: the imputations label their profiles and group classes
# independently, so they are aligned to the first before anything is averaged.

#' Fit a model to multiply imputed data and pool it with Rubin's rules
#'
#' Fits [multilpa()] to every completed data set, aligns each fit's profile and
#' group-class labels to the first one, and pools the estimates with Rubin's
#' rules. This is the route for missing *covariates*, which
#' `missing = "fiml"` cannot integrate out: the model conditions on
#' covariates, so it has no distribution for them to integrate over. Impute
#' them first (for example with the mice package, including the indicators and
#' the group identifier in the imputation model), then pass the completed data
#' here.
#'
#' Pooling is on the reported (natural) scale: for each parameter the pooled
#' estimate is the mean over imputations, and its variance is the mean
#' within-imputation variance plus `(1 + 1/m)` times the between-imputation
#' variance. Degrees of freedom follow Rubin (1987) with a large complete-data
#' sample, which is what maximum-likelihood inference assumes.
#' Pooled intervals are symmetric t intervals on this scale; unlike the
#' transformed Wald intervals of [parameter_inference()], they can extend
#' outside the parameter space for probabilities or variances.
#'
#' @section Choosing the imputation model:
#' Rubin's rules are only as good as the imputations. The covariate and the
#' indicators are related through the latent classes, so an imputation model
#' that is linear in the indicators (mice's defaults) misses that structure:
#' in `validation/pool-imputations-coverage.R` (a membership slope with 30%
#' of its covariate missing at random given the indicators, 200 replications,
#' ten imputations) 95% intervals covered 0.925 under predictive mean
#' matching and 0.935 under Bayesian linear regression, against 0.995 when the
#' covariate was approximately drawn from its conditional distribution using
#' known generating parameters. That oracle arm omits uncertainty in the
#' imputation parameters and is not a validation of proper multiple
#' imputation or nominal interval coverage. The study is one design, not a
#' general guarantee for any method. Use an imputation model that accounts for
#' the latent structure and clustering, includes the indicators and relevant
#' nonlinearities, and propagates parameter uncertainty. Imputing separately
#' within assigned classes alone ignores classification uncertainty.
#' Complete-case coverage was 0.955 in that design, but selecting on the
#' indicators changes their within-class distributions; the usual logistic
#' case-control slope argument does not guarantee unbiased slopes when the
#' outcome is latent and its measurement model is estimated.
#'
#' @param imputed The completed data sets: a list of at least two data frames
#'   with the same columns and rows in the same order, or a `mids` object from
#'   the mice package (which is then completed here).
#' @param vars,id,n_profiles As in [multilpa()].
#' @param ... Further arguments for [multilpa()], such as `n_group_classes`,
#'   `profile_covariates`, `group_covariates`, `categorical`, `n_starts` and
#'   `seed`. The same arguments fit every imputation.
#' @param level Confidence level of the pooled intervals.
#' @param vcov_type,adjust As in [parameter_inference()], applied to every
#'   imputation's fit.
#' @return An object of class `latents_pooled`. [as.data.frame()] and
#'   `get_results(x, "estimates")` give one row per parameter, with the columns
#'   of [parameter_inference()] (`level`, `outcome`, `term`, `parameter`,
#'   `estimate`, `standard_error`, `statistic`, `p_value`, `p_adjusted`,
#'   `conf_low`, `conf_high`) plus `df` (Rubin's degrees of freedom),
#'   `within` and `between` (the two variance components), `riv` (the
#'   relative increase in variance due to missingness) and `fmi` (the fraction
#'   of missing information). `statistic` and `p_value` are `NA` where
#'   [parameter_inference()] reports no test (variances and probabilities).
#'   `get_results(x, "imputations")` gives each imputation's aligned estimates,
#'   one row per imputation and parameter, and `get_results(x, "fits")` one row
#'   per imputation with its log likelihood, convergence and the label
#'   permutation that aligned it.
#' @section Conditions:
#'   `latents_bad_data` when `imputed` is not a list of at least two data frames
#'   with identical columns and row counts; `latents_pooling_failed` when an
#'   imputation's fit or inference fails (its message names the imputation and
#'   carries the original reason), since Rubin's rules need every imputation;
#'   `latents_missing_package` for a `mids` object without mice installed.
#'   Warnings from the individual fits are passed on as they are raised.
#' @references
#' Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys*.
#' Wiley.
#'
#' Meng, X.-L. (1994). Multiple-imputation inferences with uncongenial
#' sources of input. *Statistical Science*, 9(4), 538--558.
#'
#' van Buuren, S. (2018). *Flexible Imputation of Missing Data* (2nd ed.),
#' section 5.2. Chapman & Hall/CRC.
#' @examples
#' set.seed(1)
#' # Two stand-in completed data sets. With real missing covariates, pass the
#' # `mids` object mice::mice() returns instead.
#' completed <- lapply(1:2, function(i) {
#'   within(subset(course_engagement, student <= 40),
#'          previous_grade <- previous_grade + rnorm(length(previous_grade), sd = 0.1))
#' })
#' pooled <- pool_imputations(completed, c("browse", "lectures", "forum_read"),
#'                            "student", n_profiles = 2, n_group_classes = 1,
#'                            profile_covariates = "previous_grade",
#'                            n_starts = 2, seed = 1)
#' pooled
#' get_results(pooled, "fits")
#' @export
pool_imputations <- function(imputed, vars, id, n_profiles, ..., level = 0.95,
                             vcov_type = c("observed", "robust", "opg"),
                             adjust = .multilpa_p_adjust_methods) {
  stopifnot(
    "`level` must be a single number in (0, 1)" =
      is.numeric(level) && length(level) == 1L && is.finite(level) &&
      level > 0 && level < 1)
  vcov_defaulted <- missing(vcov_type)
  vcov_type <- match.arg(vcov_type)
  adjust <- match.arg(adjust)
  completed <- .multilpa_completed_data(imputed)
  call <- match.call()
  fits <- lapply(seq_along(completed), function(index) {
    .multilpa_pool_step(index, "fitting", {
      multilpa(completed[[index]], vars, id, n_profiles, ...)
    })
  })
  reference <- fits[[1L]]
  vcov_type <- .latents_weighted_vcov(.latents_is_weighted(reference), vcov_type,
                                      vcov_defaulted)
  alignments <- lapply(fits, .multilpa_pool_align, reference = reference)
  aligned <- lapply(alignments, `[[`, "fit")
  tables <- lapply(seq_along(aligned), function(index) {
    .multilpa_pool_step(index, "inference", {
      parameter_inference(aligned[[index]], vcov_type = vcov_type,
                          adjust = "none")
    })
  })
  keys <- lapply(tables, .multilpa_pool_keys)
  if (!all(vapply(keys, identical, logical(1), keys[[1L]]))) {
    stop(errorCondition(paste(
      "The imputations' fits do not report the same parameters, so they cannot",
      "be pooled. This happens when a categorical indicator has a category in",
      "one completed data set that another lacks."),
      class = "latents_pooling_failed", call = NULL))
  }
  pooled <- .multilpa_rubin(tables, level, adjust)
  structure(list(
    estimates = pooled,
    imputations = .multilpa_pool_long(tables),
    fits = .multilpa_pool_fit_frame(fits, alignments),
    covariance = .multilpa_rubin_covariance(tables, aligned, vcov_type),
    m = length(completed), level = level, vcov_type = vcov_type,
    adjust = adjust, model = class(reference)[1L],
    n_profiles = reference$n_profiles,
    n_group_classes = reference$n_group_classes, call = call),
    class = "latents_pooled")
}

#' The completed data sets, as a checked list of data frames
#' @param imputed A list of data frames or a `mids` object.
#' @return A list of at least two data frames.
#' @noRd
.multilpa_completed_data <- function(imputed) {
  bad <- function(message) {
    stop(errorCondition(message, class = "latents_bad_data", call = NULL))
  }
  if (inherits(imputed, "mids")) {
    if (!requireNamespace("mice", quietly = TRUE)) {
      stop(errorCondition(paste(
        "A `mids` object needs the mice package to be completed.",
        "Install it with install.packages(\"mice\")."),
        class = "latents_missing_package", call = NULL))
    }
    imputed <- mice::complete(imputed, action = "all")
  }
  if (!is.list(imputed) || is.data.frame(imputed) ||
      !all(vapply(imputed, is.data.frame, logical(1)))) {
    bad("`imputed` must be a list of completed data frames, or a `mids` object.")
  }
  if (length(imputed) < 2L) {
    bad("`imputed` must hold at least two completed data sets; one is not pooling.")
  }
  first <- imputed[[1L]]
  same_shape <- vapply(imputed, function(frame) {
    identical(names(frame), names(first)) && nrow(frame) == nrow(first)
  }, logical(1))
  if (!all(same_shape)) {
    bad(sprintf(paste(
      "Every completed data set must have the same columns and rows; set%s %s",
      "differ%s from the first."),
      if (sum(!same_shape) > 1L) "s" else "",
      paste(which(!same_shape), collapse = ", "),
      if (sum(!same_shape) > 1L) "" else "s"))
  }
  unname(imputed)
}

#' Run one imputation's step, naming the imputation if it fails
#'
#' Rubin's rules need every imputation, so a failure is not dropped: it is
#' raised again as `latents_pooling_failed`, with the imputation's number and
#' the original condition's class and message.
#' @param index The imputation's number.
#' @param stage `"fitting"` or `"inference"`.
#' @param expr The step.
#' @return The step's value.
#' @noRd
.multilpa_pool_step <- function(index, stage, expr) {
  tryCatch(expr, error = function(error) {
    stop(errorCondition(sprintf(
      "Imputation %d failed at %s (%s): %s", index, stage,
      class(error)[1L], conditionMessage(error)),
      class = "latents_pooling_failed", call = NULL))
  })
}

#' Align one imputation's labels to the reference fit's
#' @param fit A fit from one imputation.
#' @param reference The first imputation's fit.
#' @return A list with the aligned `fit` and the `profile_order` and
#'   `group_order` that aligned it.
#' @noRd
.multilpa_pool_align <- function(fit, reference) {
  if (inherits(fit, "multilpa_covariates")) {
    profile_order <- .multilpa_match_order(
      .multilpa_cov_profile_signature(reference, reference),
      .multilpa_cov_profile_signature(fit, reference))
    fit <- .multilpa_cov_permute(fit, profile_order, seq_len(fit$n_group_classes))
    group_order <- .multilpa_match_order(.multilpa_cov_group_signature(reference),
                                         .multilpa_cov_group_signature(fit))
    fit <- .multilpa_cov_permute(fit, seq_len(fit$n_profiles), group_order)
    return(list(fit = fit, profile_order = profile_order, group_order = group_order))
  }
  profile_order <- .multilpa_match_order(
    .multilpa_profile_signature(reference, reference),
    .multilpa_profile_signature(fit, reference))
  fit <- .multilpa_permute_profiles(fit, profile_order)
  group_order <- if (reference$n_group_classes == 1L) 1L else {
    signature <- function(x) cbind(x$profile_probabilities, x$group_probabilities)
    .multilpa_match_order(signature(reference), signature(fit))
  }
  fit <- .multilpa_permute_group_classes(fit, group_order)
  list(fit = fit, profile_order = profile_order, group_order = group_order)
}

#' What a covariate fit's profile looks like, for matching
#'
#' As `.multilpa_profile_signature()`, with the profile's overall share taken
#' from the posteriors, since a covariate fit has no constant mixing
#' probabilities.
#' @param x,reference `multilpa_covariates` fits.
#' @return A profiles-by-features matrix.
#' @noRd
.multilpa_cov_profile_signature <- function(x, reference = x) {
  blocks <- list()
  if (length(.multilpa_continuous_names(x)) > 0L) {
    spread <- sqrt(colMeans(reference$variances))
    spread[!is.finite(spread) | spread <= 0] <- 1
    blocks$means <- sweep(unname(x$means), 2L, spread, "/")
  }
  responses <- x$response_probabilities %||% list()
  if (length(responses) > 0L) {
    blocks$response <- do.call(cbind, lapply(responses, unname))
  }
  blocks <- c(blocks, .multilpa_extra_signature(x, reference))
  blocks$mixing <- matrix(colMeans(x$subject_posteriors), ncol = 1L)
  do.call(cbind, blocks)
}

#' What a covariate fit's group class looks like, for matching
#'
#' The profile mixture a group class implies, averaged over the observed
#' covariates, next to the class's share of the groups. Read after the profiles
#' are aligned, since the mixture is only comparable once they agree.
#' @param x A `multilpa_covariates` fit.
#' @return A group-classes-by-features matrix.
#' @noRd
.multilpa_cov_group_signature <- function(x) {
  mixtures <- t(vapply(x$profile_design, function(design) {
    colMeans(.multilpa_softmax(design, x$profile_coefficients))
  }, numeric(x$n_profiles)))
  cbind(matrix(mixtures, x$n_group_classes), colMeans(x$group_posteriors))
}

#' Relabel a covariate fit's profiles and group classes
#'
#' New profile `j` is old profile `profile_order[j]`, and likewise for group
#' classes. Measurement blocks and posteriors are reordered. The membership
#' logits are reparameterized, not merely reordered: they are contrasts against
#' the last profile (or group class), and relabelling may move a different one
#' into that reference position. The full logits (the reference's being zero)
#' are reordered and the new reference's are subtracted. The profile design
#' needs no change: group class `h`'s design is the same matrix whichever
#' fitted class carries that label.
#'
#' @param x A `multilpa_covariates` fit.
#' @param profile_order,group_order Integer permutations.
#' @return The relabelled fit, whose likelihood is unchanged.
#' @noRd
.multilpa_cov_permute <- function(x, profile_order, group_order) {
  stopifnot(inherits(x, "multilpa_covariates"),
            identical(sort(as.integer(profile_order)), seq_len(x$n_profiles)),
            identical(sort(as.integer(group_order)), seq_len(x$n_group_classes)))
  x <- .multilpa_permute_extra(x, profile_order)
  n_profiles <- x$n_profiles
  n_groups <- x$n_group_classes
  profiles <- paste0("profile_", seq_len(n_profiles))
  reorder_rows <- function(block) {
    if (is.null(block)) return(NULL)
    names_kept <- colnames(block)
    reordered <- block[profile_order, , drop = FALSE]
    dimnames(reordered) <- list(if (!is.null(rownames(block))) profiles, names_kept)
    reordered
  }
  x$means <- reorder_rows(x$means)
  x$variances <- reorder_rows(x$variances)
  x$standard_deviations <- reorder_rows(x$standard_deviations)
  if (!is.null(x$covariances)) {
    x$covariances <- x$covariances[, , profile_order, drop = FALSE]
  }
  if (!is.null(x$response_probabilities)) {
    x$response_probabilities <- lapply(x$response_probabilities, reorder_rows)
  }
  ## Baseline-category logits: complete them with the reference's zeros,
  ## reorder, and take contrasts against the new last category.
  recontrast <- function(coefficients, order) {
    full <- cbind(coefficients, 0)[, order, drop = FALSE]
    result <- full[, -ncol(full), drop = FALSE] - full[, ncol(full)]
    dimnames(result) <- dimnames(coefficients)
    result
  }
  beta <- recontrast(x$profile_coefficients, profile_order)
  ## Rows: one intercept per group class, then the slopes, which come once
  ## (shared) or in one block per group class.
  n_slopes <- nrow(beta) - n_groups
  slope_rows <- if (identical(x$profile_slopes %||% "shared", "group_class")) {
    width <- n_slopes / n_groups
    unlist(lapply(group_order, function(h) n_groups + (h - 1L) * width + seq_len(width)))
  } else n_groups + seq_len(n_slopes)
  kept_names <- dimnames(beta)
  beta <- beta[c(group_order, slope_rows), , drop = FALSE]
  dimnames(beta) <- kept_names
  x$profile_coefficients <- beta
  x$group_coefficients <- recontrast(x$group_coefficients, group_order)
  ## Posteriors and assignments follow the labels.
  x$subject_posteriors <- x$subject_posteriors[, profile_order, drop = FALSE]
  x$group_posteriors <- x$group_posteriors[, group_order, drop = FALSE]
  x$subject_profiles <- match(x$subject_profiles, profile_order)
  x$group_classes <- match(x$group_classes, group_order)
  x$effective_profile_counts <- x$effective_profile_counts[profile_order]
  x$effective_group_counts <- x$effective_group_counts[group_order]
  if (!is.null(x$group_priors)) {
    x$group_priors <- x$group_priors[, group_order, drop = FALSE]
  }
  if (!is.null(x$profile_priors)) {
    x$profile_priors <- lapply(x$profile_priors[group_order], function(prior) {
      prior[, profile_order, drop = FALSE]
    })
  }
  x
}

#' The identifying columns of an inference table, as one key per row
#' @param table A [parameter_inference()] table.
#' @return A character vector.
#' @noRd
.multilpa_pool_keys <- function(table) {
  paste(table$level, table$outcome, table$term, table$parameter, sep = "\r")
}

#' Pool inference tables with Rubin's rules
#'
#' Rubin (1987): the pooled estimate is the mean `Q`, the total variance
#' `T = U + (1 + 1/m) B` with `U` the mean within-imputation variance and `B`
#' the between-imputation variance, `r=(1 + 1/m) * B / U`, degrees of freedom
#' `(m - 1)(1 + 1/r)^2` (infinite when `B` is zero), and the fraction of
#' missing information `(r + 2 / (df + 3)) / (r + 1)`.
#'
#' @param tables Aligned [parameter_inference()] tables, one per imputation.
#' @param level Interval level.
#' @param adjust A [stats::p.adjust()] method for `p_adjusted`.
#' @return A data frame, one row per parameter.
#' @noRd
.multilpa_rubin <- function(tables, level, adjust = "none") {
  m <- length(tables)
  estimates <- vapply(tables, `[[`, numeric(nrow(tables[[1L]])), "estimate")
  variances <- vapply(tables, function(table) table$standard_error^2,
                      numeric(nrow(tables[[1L]])))
  estimates <- matrix(estimates, ncol = m)
  variances <- matrix(variances, ncol = m)
  pooled <- rowMeans(estimates)
  within <- rowMeans(variances)
  between <- apply(estimates, 1L, stats::var)
  total <- within + (1 + 1 / m) * between
  riv <- (1 + 1 / m) * between / within
  df <- ifelse(between > 0, (m - 1) * (1 + 1 / riv)^2, Inf)
  fmi <- (riv + 2 / (df + 3)) / (riv + 1)
  standard_error <- sqrt(total)
  ## A test is reported only where the imputations' tables report one.
  tested <- !is.na(tables[[1L]]$statistic)
  statistic <- ifelse(tested, pooled / standard_error, NA_real_)
  p_value <- ifelse(tested, 2 * stats::pt(-abs(statistic), df), NA_real_)
  critical <- stats::qt((1 + level) / 2, df)
  data.frame(
    tables[[1L]][c("level", "outcome", "term", "parameter")],
    estimate = pooled, standard_error = standard_error,
    statistic = statistic, p_value = p_value,
    p_adjusted = stats::p.adjust(p_value, method = adjust),
    conf_low = pooled - critical * standard_error,
    conf_high = pooled + critical * standard_error,
    df = df, within = within, between = between, riv = riv, fmi = fmi,
    row.names = NULL, stringsAsFactors = FALSE)
}

#' Pool the imputations' covariance matrices
#'
#' A covariate-free inference table carries its covariance; a covariate fit's
#' comes from [vcov()] with the same `vcov_type`. Either way its diagonal must
#' reproduce the table's squared standard errors, which is what ties its rows
#' to the table's rows.
#' @param tables Aligned [parameter_inference()] tables.
#' @param fits The aligned fits the tables came from.
#' @param vcov_type As in [parameter_inference()].
#' @return The total covariance `U + (1 + 1/m) B`, named by parameter.
#' @noRd
.multilpa_rubin_covariance <- function(tables, fits, vcov_type) {
  m <- length(tables)
  covariances <- lapply(seq_len(m), function(index) {
    attr(tables[[index]], "covariance") %||%
      vcov(fits[[index]], vcov_type = vcov_type)
  })
  invisible(lapply(seq_len(m), function(index) {
    stopifnot("each covariance must match its inference table row for row" =
                isTRUE(all.equal(sqrt(diag(covariances[[index]])),
                                 tables[[index]]$standard_error,
                                 check.attributes = FALSE, tolerance = 1e-8)))
  }))
  within <- Reduce(`+`, covariances) / m
  estimates <- do.call(rbind, lapply(tables, `[[`, "estimate"))
  total <- within + (1 + 1 / m) * stats::cov(estimates)
  dimnames(total) <- dimnames(covariances[[1L]])
  total
}

#' Every imputation's aligned estimates, one row per imputation and parameter
#' @param tables Aligned [parameter_inference()] tables.
#' @return A data frame.
#' @noRd
.multilpa_pool_long <- function(tables) {
  do.call(rbind, lapply(seq_along(tables), function(index) {
    data.frame(imputation = index,
               tables[[index]][c("level", "outcome", "term", "parameter",
                                 "estimate", "standard_error")],
               row.names = NULL, stringsAsFactors = FALSE)
  }))
}

#' One row per imputation: its fit and the relabelling that aligned it
#' @param fits The unaligned fits.
#' @param alignments The alignments from `.multilpa_pool_align()`.
#' @return A data frame.
#' @noRd
.multilpa_pool_fit_frame <- function(fits, alignments) {
  data.frame(
    imputation = seq_along(fits),
    log_likelihood = vapply(fits, `[[`, numeric(1), "log_likelihood"),
    n_parameters = vapply(fits, function(fit) as.numeric(fit$n_parameters),
                          numeric(1)),
    converged = vapply(fits, function(fit) isTRUE(fit$converged), logical(1)),
    profile_order = vapply(alignments, function(alignment) {
      paste(alignment$profile_order, collapse = ",")
    }, character(1)),
    group_order = vapply(alignments, function(alignment) {
      paste(alignment$group_order, collapse = ",")
    }, character(1)),
    stringsAsFactors = FALSE)
}

#' Tables of a pooled multiply imputed fit
#'
#' @param x A `latents_pooled` object from [pool_imputations()].
#' @param what `"estimates"` (the pooled table), `"imputations"` (every
#'   imputation's aligned estimates) or `"fits"` (one row per imputation).
#' @param ... Unused.
#' @return A base `data.frame`; see [pool_imputations()] for the columns.
#' @examples
#' set.seed(1)
#' # Two stand-in completed data sets. With real missing covariates, pass the
#' # `mids` object mice::mice() returns instead.
#' completed <- lapply(1:2, function(i) {
#'   within(subset(course_engagement, student <= 40),
#'          previous_grade <- previous_grade + rnorm(length(previous_grade), sd = 0.1))
#' })
#' pooled <- pool_imputations(completed, c("browse", "lectures", "forum_read"),
#'                            "student", n_profiles = 2, n_group_classes = 1,
#'                            profile_covariates = "previous_grade",
#'                            n_starts = 2, seed = 1)
#' get_results(pooled, "imputations")
#' @export
get_results.latents_pooled <- function(x, what = c("estimates", "imputations", "fits"),
                                       ...) {
  what <- match.arg(what)
  x[[what]]
}

#' @rdname get_results.latents_pooled
#' @param row.names,optional Unused; part of the generic.
#' @export
as.data.frame.latents_pooled <- function(x, row.names = NULL, optional = FALSE,
                                         what = "estimates", ...) {
  get_results.latents_pooled(x, what = what)
}

#' Print a pooled multiply imputed fit
#' @param x A `latents_pooled` object.
#' @param digits Significant digits.
#' @param rows How many rows of the pooled table to show.
#' @param ... Unused.
#' @return `x`, invisibly. Called for the side effect of printing the model,
#'   the number of imputations and the pooled table.
#' @export
print.latents_pooled <- function(x, digits = 4L, rows = 20L, ...) {
  cat(sprintf(paste(
    "Pooled over %d imputations (Rubin's rules): %s, %d profiles,",
    "%d group classes\n"), x$m, x$model, x$n_profiles, x$n_group_classes))
  cat(sprintf("Largest fraction of missing information: %.3f\n",
              max(x$estimates$fmi)))
  shown <- utils::head(x$estimates[c("level", "outcome", "term", "parameter",
                                     "estimate", "standard_error", "conf_low",
                                     "conf_high", "fmi")], rows)
  print(shown, digits = digits, row.names = FALSE)
  if (nrow(x$estimates) > rows) {
    cat(sprintf("... %d more rows: as.data.frame(x)\n",
                nrow(x$estimates) - rows))
  }
  invisible(x)
}

#' Summarize a pooled multiply imputed fit
#' @param object A `latents_pooled` object.
#' @param ... Unused.
#' @return A one-row base `data.frame`: the number of imputations, the model,
#'   whether every imputation converged, the range of their log likelihoods,
#'   how many needed relabelling, and the largest and median fraction of
#'   missing information.
#' @export
summary.latents_pooled <- function(object, ...) {
  fits <- object$fits
  identity_profiles <- paste(seq_len(object$n_profiles), collapse = ",")
  identity_groups <- paste(seq_len(object$n_group_classes), collapse = ",")
  data.frame(
    m = object$m, model = object$model, n_profiles = object$n_profiles,
    n_group_classes = object$n_group_classes,
    all_converged = all(fits$converged),
    log_likelihood_min = min(fits$log_likelihood),
    log_likelihood_max = max(fits$log_likelihood),
    relabelled = sum(fits$profile_order != identity_profiles |
                       fits$group_order != identity_groups),
    fmi_max = max(object$estimates$fmi),
    fmi_median = stats::median(object$estimates$fmi),
    stringsAsFactors = FALSE)
}

#' Pooled covariance of a multiply imputed fit
#' @param object A `latents_pooled` object.
#' @param ... Unused.
#' @return The total covariance matrix `U + (1 + 1/m) B` of the natural-scale
#'   parameters, named as [vcov()] names them for a single fit, or `NULL` when
#'   the imputations' fits carried no covariance.
#' @export
vcov.latents_pooled <- function(object, ...) {
  object$covariance
}

#' Plot a pooled multiply imputed fit
#'
#' One row per parameter: each imputation's aligned estimate as an open circle
#' and the pooled estimate as a filled circle with its interval, so the
#' between-imputation spread is visible next to the total uncertainty.
#' @param x A `latents_pooled` object.
#' @param level Which parameter level to show: `"profile"` (the profile
#'   membership coefficients or probabilities), `"group"` or `"measurement"`.
#' @param main Plot title; `NULL` for the default.
#' @param ... Unused.
#' @return `x`, invisibly. Called for the side effect of drawing.
#' @export
plot.latents_pooled <- function(x, level = c("profile", "group", "measurement"),
                                main = NULL, ...) {
  level <- match.arg(level)
  estimates <- x$estimates[x$estimates$level %in% level, , drop = FALSE]
  imputations <- x$imputations[x$imputations$level %in% level, , drop = FALSE]
  if (nrow(estimates) == 0L) {
    stop(errorCondition(sprintf("This fit has no `%s` parameters to plot.", level),
                        class = "latents_bad_argument", call = NULL))
  }
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE, after = FALSE)
  label <- function(table) paste(table$outcome, table$term, table$parameter)
  labels <- unique(label(estimates))
  position <- rev(seq_along(labels))
  pooled_y <- position[match(label(estimates), labels)]
  imputed_y <- position[match(label(imputations), labels)]
  colours <- .multilpa_palette(2L)
  ## Sized from the label lengths: measuring text needs an open plot.
  graphics::par(mar = c(4, min(0.45 * max(nchar(labels)) + 1, 20), 3, 1))
  range_x <- range(c(estimates$conf_low, estimates$conf_high,
                     imputations$estimate), na.rm = TRUE)
  graphics::plot(imputations$estimate, imputed_y, xlim = range_x,
                 ylim = c(0.5, length(labels) + 0.5), yaxt = "n", pch = 1,
                 col = colours[2L],
                 xlab = sprintf("Estimate with pooled %g%% interval", 100 * x$level),
                 ylab = "", main = main %||% sprintf("Pooled over %d imputations", x$m))
  graphics::segments(estimates$conf_low, pooled_y, estimates$conf_high, pooled_y,
                     col = colours[1L], lwd = 2)
  graphics::points(estimates$estimate, pooled_y, pch = 16, col = colours[1L],
                   cex = 1.3)
  graphics::axis(2, at = position, labels = labels, las = 1, cex.axis = 0.8)
  graphics::legend("topright", legend = c("pooled", "one imputation"),
                   pch = c(16, 1), col = colours, bty = "n")
  invisible(x)
}
