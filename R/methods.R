#' Print a fitted multilevel latent profile model
#' @param x An `multilpa` model.
#' @param ... Reserved for compatibility with `print()`.
#' @return The input model, invisibly.
#' @examples
#' # Printing a fitted model displays its convergence and information criteria.
#' @export
print.multilpa <- function(x, ...) {
  stopifnot(inherits(x, "multilpa"))
  cat(sprintf("Two-level latent profile analysis: %d profiles, %d group classes\n",
              x$n_profiles, x$n_group_classes))
  covariance_model <- x$covariance_model %||% "diagonal"
  cat(sprintf("%d individuals in %d groups; %s %s residual covariance\n",
              x$n_observations, x$n_groups, x$variance_model, covariance_model))
  cat(sprintf("Log likelihood: %.6f | AIC: %.3f | BIC (groups): %.3f\n",
              x$log_likelihood, x$aic, x$bic))
  cat(sprintf("Converged: %s | iterations: %d | best start: %d/%d\n",
              x$converged, x$iterations, x$best_start, nrow(x$starts)))
  if (!is.null(x$n_informative) && x$n_informative < x$n_observations) {
    cat(sprintf("%d row(s) carry no observed indicator; the individual-level BIC uses %d.\n",
                x$n_observations - x$n_informative, x$n_informative))
  }
  if (length(x$fixed %||% character()) > 0L) {
    cat(sprintf("Held fixed, not estimated here: %s%s\n",
                paste(x$fixed, collapse = ", "),
                if (isTRUE(x$staged)) " (first stage)" else ""))
    cat(sprintf("Parameters estimated here: %d; with the held measurement: %d\n",
                x$n_parameters,
                x$n_parameters_with_measurement %||% x$n_parameters))
  }
  if (x$boundary) cat("A variance is at its specified lower bound.\n")
  if (x$small_classes) cat("An effective class membership is below one.\n")
  invisible(x)
}

#' Summarize a fitted multilevel latent profile model
#'
#' The summary collects the estimates, the effective class counts, the
#' information criteria and the restart diagnostics of a fit. Every field is
#' built explicitly, and a field the fit does not carry takes its documented
#' default rather than being silently absent, so the summary of a diagonal fit
#' and the summary of a staged full-covariance fit have exactly the same names.
#' Read the numbers with [as.data.frame()], which returns them as tidy tables;
#' the printed form is a human-facing report.
#'
#' @param object An `multilpa` model.
#' @param ... Reserved for compatibility with `summary()`.
#' @return A `summary_multilpa` object: a named list, never carrying an `NA`
#'   name, with the fit's dimensions (`n_observations`, `n_informative`,
#'   `n_groups`, `n_profiles`, `n_group_classes`), its specification
#'   (`variance_model`, `covariance_model`, `missing`, `continuous`, `fixed`,
#'   `staged`), its estimates (`means`, `variances`, `standard_deviations`,
#'   `covariances`, `response_probabilities`, `profile_probabilities`,
#'   `group_probabilities`), the effective class counts at both levels, the
#'   likelihood, parameter counts and information criteria, and the restart
#'   diagnostics. `covariances` is `NULL` under the diagonal parameterization
#'   and `response_probabilities` is `NULL` when no indicator is categorical;
#'   both keep their names in either case. Use
#'   [as.data.frame.summary_multilpa()] rather than reading the fields.
#' @seealso [as.data.frame.summary_multilpa()] for the tidy tables.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                 n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' summary(fit)
#' @export
summary.multilpa <- function(object, ...) {
  stopifnot("`object` must be an `multilpa` fit" = inherits(object, "multilpa"))
  # Built field by field rather than by subsetting the fit with a vector of
  # names: subsetting a list by a name it does not carry yields a `NULL`
  # element keyed `NA`, so an absent field became an unnamed hole in the
  # summary instead of a documented default.
  result <- list(
    call = object$call,
    n_observations = object$n_observations,
    n_informative = object$n_informative %||% object$n_observations,
    n_groups = object$n_groups,
    n_profiles = object$n_profiles,
    n_group_classes = object$n_group_classes,
    variance_model = object$variance_model,
    covariance_model = object$covariance_model %||% "diagonal",
    missing = object$missing %||% "listwise",
    continuous = .multilpa_continuous_names(object),
    fixed = object$fixed %||% character(),
    staged = isTRUE(object$staged),
    means = object$means,
    variances = object$variances,
    standard_deviations = object$standard_deviations,
    covariances = object$covariances,
    response_probabilities = object$response_probabilities,
    profile_probabilities = object$profile_probabilities,
    group_probabilities = object$group_probabilities,
    effective_profile_counts = colSums(object$subject_posteriors),
    effective_group_counts = colSums(object$group_posteriors),
    log_likelihood = object$log_likelihood,
    n_parameters = object$n_parameters,
    n_parameters_with_measurement =
      object$n_parameters_with_measurement %||% object$n_parameters,
    aic = object$aic,
    bic = object$bic,
    bic_individual = object$bic_individual,
    converged = object$converged,
    iterations = object$iterations,
    boundary = object$boundary,
    min_variance = object$min_variance,
    small_classes = object$small_classes,
    starts = object$starts,
    best_start = object$best_start,
    n_best_replicated = object$n_best_replicated,
    replication_tolerance = object$replication_tolerance)
  stopifnot("the summary must not carry an unnamed field" =
              !anyNA(names(result)) && !any(names(result) == ""))
  class(result) <- "summary_multilpa"
  result
}

#' Tidy a multilevel LPA summary
#'
#' The printed summary is a report for a reader; this returns the same content
#' as tidy tables, so nothing has to be read out of the summary object by hand.
#'
#' @param x A `summary_multilpa` object from [summary.multilpa()].
#' @param row.names Passed to `data.frame()`; `NULL` gives default row names.
#' @param optional Ignored, present for generic compatibility.
#' @param what Which table to return. `"profiles"` gives the Gaussian
#'   measurement model, `"responses"` the categorical measurement model,
#'   `"profile_probabilities"` the profile prevalence within each group class,
#'   `"counts"` the effective class memberships at both levels,
#'   `"covariances"` the within-profile residual covariance matrices,
#'   `"fit"` the one-row fit summary, and `"starts"` the restart diagnostics.
#' @param ... Ignored.
#' @return A base `data.frame` whose columns depend on `what`.
#'   `"profiles"` has one row per profile and continuous indicator, with columns
#'   `profile`, `indicator`, `mean`, `variance` and `standard_deviation`, and
#'   zero rows when every indicator is categorical. `"responses"` has one row
#'   per profile, categorical indicator and category, with columns `profile`,
#'   `indicator`, `category`, `probability` and `threshold`, and zero rows when
#'   no indicator is categorical. `"profile_probabilities"` has one row per
#'   group class and profile, with columns `group_class`, `profile`,
#'   `probability` and `group_class_probability`. `"counts"` has one row per
#'   class at each level, with columns `level` (`"individuals"` or `"groups"`),
#'   `class`, `effective_count` (the summed posteriors) and
#'   `effective_proportion`. `"covariances"` has one row per profile and
#'   ordered pair of continuous indicators, with columns `profile`, `indicator`,
#'   `indicator_2` and `covariance`; under the diagonal parameterization the
#'   off-diagonal covariances are zero by assumption rather than missing.
#'   `"fit"` has exactly one row, carrying the fit's dimensions, likelihood,
#'   parameter counts, information criteria and convergence diagnostics.
#'   `"starts"` has one row per EM start.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                 n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' as.data.frame(summary(fit))
#' as.data.frame(summary(fit), what = "counts")
#' as.data.frame(summary(fit), what = "fit")
#' @export
as.data.frame.summary_multilpa <- function(x, row.names = NULL,
                                           optional = FALSE,
                                           what = c("profiles", "responses",
                                                    "profile_probabilities",
                                                    "counts", "covariances",
                                                    "fit", "starts"), ...) {
  stopifnot("`x` must be a `summary_multilpa` object" =
              inherits(x, "summary_multilpa"))
  what <- match.arg(what)
  result <- switch(what,
    profiles = .multilpa_profile_frame(x),
    responses = .multilpa_response_frame(x),
    profile_probabilities = .multilpa_probability_frame(x),
    counts = .multilpa_count_frame(x),
    covariances = .multilpa_covariance_frame(x),
    fit = .multilpa_fit_frame(x),
    starts = x$starts)
  row.names(result) <- row.names
  result
}

#' Effective class memberships at both levels as a tidy table
#' @param x A `summary_multilpa` object.
#' @return One row per class at each level.
#' @noRd
.multilpa_count_frame <- function(x) {
  counts <- list(individuals = x$effective_profile_counts,
                 groups = x$effective_group_counts)
  rows <- lapply(names(counts), function(level) {
    count <- unname(counts[[level]])
    data.frame(level = level, class = seq_along(count),
               effective_count = count,
               effective_proportion = count / sum(count))
  })
  do.call(rbind, rows)
}

#' Within-profile residual covariances as a tidy table
#'
#' The diagonal parameterization carries no covariance array, but it does state
#' a covariance matrix: the variances on the diagonal and exact zeros off it.
#' Reporting that explicitly means one shape answers both parameterizations.
#'
#' @param x A `summary_multilpa` object.
#' @return One row per profile and ordered pair of continuous indicators.
#' @noRd
.multilpa_covariance_frame <- function(x) {
  indicators <- .multilpa_continuous_names(x)
  n_indicators <- length(indicators)
  n_profiles <- x$n_profiles
  if (n_indicators == 0L) {
    return(data.frame(profile = integer(), indicator = character(),
                      indicator_2 = character(), covariance = numeric()))
  }
  blocks <- lapply(seq_len(n_profiles), function(profile) {
    block <- if (is.null(x$covariances)) {
      diag(x$variances[profile, ], n_indicators, n_indicators)
    } else {
      matrix(x$covariances[, , profile], n_indicators, n_indicators)
    }
    data.frame(profile = profile,
               indicator = rep(indicators, times = n_indicators),
               indicator_2 = rep(indicators, each = n_indicators),
               covariance = as.vector(block))
  })
  do.call(rbind, blocks)
}

#' The one-row fit summary as a tidy table
#' @param x A `summary_multilpa` object.
#' @return One row, describing the whole fit.
#' @noRd
.multilpa_fit_frame <- function(x) {
  data.frame(n_observations = x$n_observations,
             n_informative = x$n_informative,
             n_groups = x$n_groups,
             n_profiles = x$n_profiles,
             n_group_classes = x$n_group_classes,
             n_parameters = x$n_parameters,
             n_parameters_with_measurement = x$n_parameters_with_measurement,
             log_likelihood = x$log_likelihood,
             aic = x$aic, bic = x$bic, bic_individual = x$bic_individual,
             converged = x$converged, iterations = x$iterations,
             boundary = x$boundary, small_classes = x$small_classes,
             best_start = x$best_start,
             n_best_replicated = x$n_best_replicated,
             row.names = NULL)
}

#' Print a multilevel LPA summary
#'
#' A human-facing report of the fit. The same content is available as tidy
#' tables from [as.data.frame.summary_multilpa()].
#'
#' @param x A `summary_multilpa` object.
#' @param digits Number of printed significant digits.
#' @param ... Additional arguments passed to matrix printing.
#' @return The summary, invisibly.
#' @seealso [as.data.frame.summary_multilpa()] for the same content as data.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                 n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' print(summary(fit), digits = 3)
#' @export
print.summary_multilpa <- function(x, digits = 4L, ...) {
  stopifnot("`x` must be a `summary_multilpa` object" =
              inherits(x, "summary_multilpa"),
            "`digits` must be a single number between 1 and 22" =
              is.numeric(digits) && length(digits) == 1L &&
              is.finite(digits) && digits >= 1 && digits <= 22)
  cat(sprintf("Multilevel LPA: %d profiles and %d group classes\n",
              x$n_profiles, x$n_group_classes))
  cat(sprintf("Individuals: %d; groups: %d; parameters: %d; converged: %s\n",
              x$n_observations, x$n_groups, x$n_parameters, x$converged))
  if (length(x$fixed %||% character()) > 0L) {
    cat(sprintf("Held fixed%s: %s; parameters with the held measurement: %d\n",
                if (isTRUE(x$staged)) " from the first stage" else "",
                paste(x$fixed, collapse = ", "),
                x$n_parameters_with_measurement %||% x$n_parameters))
  }
  # An all-categorical fit has a zero-column Gaussian block, which printed as
  # an empty matrix under a heading promising means.
  if (length(.multilpa_continuous_names(x)) > 0L) {
    cat("\nProfile means:\n")
    print(x$means, digits = digits, ...)
    cat("\nProfile standard deviations:\n")
    print(x$standard_deviations, digits = digits, ...)
    if (identical(x$covariance_model, "full")) {
      cat("\nProfile residual covariance matrices:\n")
      print(x$covariances, digits = digits, ...)
    }
  }
  if (length(x$response_probabilities) > 0L) {
    cat("\nCategorical response probabilities:\n")
    print(x$response_probabilities, digits = digits, ...)
  }
  cat("\nProfile probabilities within each group class:\n")
  print(x$profile_probabilities, digits = digits, ...)
  cat("\nGroup-class probabilities:\n")
  print(x$group_probabilities, digits = digits, ...)
  cat("\nEffective individual memberships:\n")
  print(x$effective_profile_counts, digits = digits, ...)
  cat("\nEffective group memberships:\n")
  print(x$effective_group_counts, digits = digits, ...)
  cat(sprintf("\nLog likelihood: %.6f; AIC: %.3f\nBIC (groups): %.3f; BIC (individuals): %.3f\n",
              x$log_likelihood, x$aic, x$bic, x$bic_individual))
  cat(sprintf("Best likelihood replicated in %d/%d starts (absolute tolerance %.3g).\n",
              x$n_best_replicated, nrow(x$starts), x$replication_tolerance))
  if (!x$converged) cat("WARNING: best start did not converge.\n")
  if (x$boundary) cat("WARNING: at least one variance is at min_variance.\n")
  if (x$small_classes) cat("WARNING: an effective class membership is below one.\n")
  cat("\nStart diagnostics:\n")
  print(x$starts, digits = digits, row.names = FALSE)
  invisible(x)
}

#' Extract the multilevel model log likelihood
#' @param object An `multilpa` model.
#' @param ... Reserved for compatibility with `logLik()`.
#' @return A `logLik` object with parameter count `df` and the number of observed
#'   groups as `nobs`. Thus `stats::BIC()` uses group-count BIC. For the
#'   individual-count alternative, and every other criterion, call
#'   [information_criteria()], which reports both conventions side by side.
#' @examples
#' # After fitting: logLik(fit); AIC(fit); BIC(fit)
#' @export
#' @importFrom stats logLik
logLik.multilpa <- function(object, ...) {
  stopifnot(inherits(object, "multilpa"))
  structure(object$log_likelihood, df = object$n_parameters,
            nobs = object$n_groups, class = "logLik")
}

#' Extract the number of independent groups
#' @param object An `multilpa` model.
#' @param ... Reserved for compatibility with `nobs()`.
#' @return Number of observed groups, which are the independent units of the
#'   two-level likelihood. For the individual count alongside every other
#'   sample-size-dependent quantity, call [information_criteria()].
#' @examples
#' # After fitting: nobs(fit)
#' @export
#' @importFrom stats nobs
nobs.multilpa <- function(object, ...) {
  stopifnot(inherits(object, "multilpa"))
  object$n_groups
}
