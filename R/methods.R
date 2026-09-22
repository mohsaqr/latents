#' Print a fitted multilevel latent profile model
#' @param x An `multilpa` model.
#' @param rows How many rows of the printed table to show before truncating.
#' @param ... Reserved for compatibility with `print()`.
#' @return The input model, invisibly. Called for the side effect of printing
#'   the class counts, the sample sizes and covariance specification, the log
#'   likelihood with AIC and group-level BIC, the convergence and restart
#'   diagnostics, and any blocks the fit held fixed.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                 n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' print(fit)
#' @export
print.multilpa <- function(x, rows = 20L, ...) {
  stopifnot(inherits(x, "multilpa"))
  ## A single-level fit has one observation per unit and one group class, so
  ## naming either would describe machinery rather than the model asked for.
  if (isTRUE(x$single_level)) {
    cat(sprintf("Latent profile analysis: %d profile%s\n",
                x$n_profiles, if (x$n_profiles == 1L) "" else "s"))
  } else {
    cat(sprintf("Two-level latent profile analysis: %d profile%s, %d group class%s\n",
                x$n_profiles, if (x$n_profiles == 1L) "" else "s",
                x$n_group_classes, if (x$n_group_classes == 1L) "" else "es"))
  }
  covariance_model <- x$covariance_model %||% "diagonal"
  cat(sprintf("%d %s; %s %s residual covariance%s\n",
              x$n_observations,
              if (isTRUE(x$single_level)) "observations" else
                sprintf("individuals in %d groups", x$n_groups),
              x$variance_model, covariance_model,
              if (is.null(x$covariance_structure)) "" else
                sprintf(" (%s)", x$covariance_structure)))
  if (!identical(x$centering %||% "none", "none")) {
    cat(sprintf("Indicators %s-centred, so the profiles are profiles of change\n",
                x$centering))
  }
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
  .multilpa_print_primary(x, rows = rows)
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
#'   [get_data()] rather than reading the fields.
#' @seealso [get_data()] for the tidy tables.
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
  # Every table the fit can produce, built once here, so `get_data()` on the
  # summary serves the same tables the fit would and `print()` can show them
  # all without recomputing anything.
  result$tables <- get_data(object, "all")
  class(result) <- "summary_multilpa"
  result
}

#' Coerce a multilevel LPA summary to its primary table
#'
#' Plain coercion, as the base generic means it: one object, one data frame.
#' A summary carries every table the object it describes can produce, and
#' [get_data()] names them.
#'
#' @param x An object of class `summary_multilpa`.
#' @param row.names Passed to `data.frame()`; `NULL` gives default row names.
#' @param optional Ignored, present for generic compatibility.
#' @param ... Must be empty. An argument here raises `multilpa_bad_argument`
#'   naming it, rather than being dropped.
#' @return A base `data.frame`: one row per profile and continuous indicator.
#' @seealso [get_data()] for every other table this summary holds.
#' @examples
#' fit <- multilpa(
#'   course_engagement,
#'   vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
#'   id = "student", n_profiles = 2, n_group_classes = 2, n_starts = 4,
#'   seed = 1
#' )
#' as.data.frame(summary(fit))
#' get_data(summary(fit), what = "counts")
#' @export
as.data.frame.summary_multilpa <- function(x, row.names = NULL, optional = FALSE, ...) {
  stopifnot("`x` must be an object of class `summary_multilpa`" = inherits(x, "summary_multilpa"))
  .multilpa_coerce(x, row.names, list(...))
}

#' Effective class memberships at both levels as a tidy table
#'
#' Read from the summary's stored counts where there are any and summed from
#' the posteriors otherwise, so the one definition serves a fit and its summary
#' alike. A family with no discrete group classes contributes one level rather
#' than an empty second one.
#'
#' @param x A fitted model of this package, or its summary.
#' @return One row per class at each level the object has.
#' @noRd
.multilpa_count_frame <- function(x) {
  individuals <- x$effective_profile_counts %||% colSums(x$subject_posteriors)
  groups <- x$effective_group_counts
  if (is.null(groups) && !is.null(x$group_posteriors)) {
    groups <- colSums(x$group_posteriors)
  }
  counts <- Filter(Negate(is.null),
                   list(individuals = individuals, groups = groups))
  rows <- lapply(names(counts), function(level) {
    count <- unname(counts[[level]])
    data.frame(level = level, class = seq_along(count),
               effective_count = count,
               effective_proportion = count / sum(count))
  })
  result <- do.call(rbind, rows)
  row.names(result) <- NULL
  result
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
  vars <- .multilpa_continuous_names(x)
  n_indicators <- length(vars)
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
               indicator = rep(vars, times = n_indicators),
               indicator_2 = rep(vars, each = n_indicators),
               covariance = as.vector(block))
  })
  do.call(rbind, blocks)
}

#' The one-row description of a fit as a tidy table
#'
#' Family-shaped rather than one column set for every family: a random
#' intercept records its integration diagnostics and a covariate model its
#' predictor counts, and padding those onto a two-level fit would make the
#' common case mostly `NA`. Comparing fits across families is what the
#' information criteria are for, and that table does have one column set.
#'
#' @param x A fitted model of this package.
#' @return One row, describing the whole fit.
#' @noRd
.multilpa_fit_frame <- function(x) {
  if (inherits(x, "multilpa_covariates")) {
    return(.multilpa_covariate_fit_frame(x))
  }
  # A field an older fit or another family does not carry is absent, not
  # short: `NULL` and a length-two vector would both corrupt the row.
  data.frame(n_observations = .multilpa_one(x$n_observations, NA_integer_),
             n_informative = .multilpa_one(x$n_informative %||% x$n_observations,
                                           NA_integer_),
             n_groups = .multilpa_one(x$n_groups, NA_integer_),
             n_profiles = .multilpa_one(x$n_profiles, NA_integer_),
             n_group_classes = .multilpa_one(x$n_group_classes, NA_integer_),
             centering = .multilpa_one(x$centering %||% "none", NA_character_),
             covariance_structure = .multilpa_one(x$covariance_structure,
                                                  NA_character_),
             n_parameters = .multilpa_one(x$n_parameters, NA_integer_),
             n_parameters_with_measurement =
               .multilpa_one(x$n_parameters_with_measurement %||%
                               x$n_parameters, NA_integer_),
             log_likelihood = .multilpa_one(x$log_likelihood, NA_real_),
             aic = .multilpa_one(x$aic, NA_real_),
             bic_groups = .multilpa_one(x$bic, NA_real_),
             bic_individual = .multilpa_one(x$bic_individual, NA_real_),
             converged = .multilpa_one(x$converged, NA),
             iterations = .multilpa_one(x$iterations, NA_integer_),
             boundary = .multilpa_one(x$boundary, NA),
             small_classes = .multilpa_one(x$small_classes, NA),
             best_start = .multilpa_one(x$best_start, NA_integer_),
             n_best_replicated = .multilpa_one(x$n_best_replicated, NA_integer_),
             row.names = NULL)
}

#' One scalar, or a stated absence
#'
#' A field an older fit or another family does not carry is `NULL`, and one
#' recorded per class is longer than one. Either would corrupt a one-row frame,
#' so both become the absent value of the right type.
#'
#' @param value The field.
#' @param absent What to report when it is not a single value.
#' @return A length-one vector.
#' @noRd
.multilpa_one <- function(value, absent) {
  if (is.null(value) || length(value) != 1L) absent else value
}

#' The one-row description of a covariate fit
#' @param x A fitted `multilpa_covariates` model.
#' @return One row.
#' @noRd
.multilpa_covariate_fit_frame <- function(x) {
  data.frame(
    n_observations = x$n_observations,
    n_groups = x$n_groups,
    n_profiles = x$n_profiles,
    n_group_classes = x$n_group_classes,
    n_profile_covariates = length(x$profile_covariates),
    n_group_covariates = length(x$group_covariates),
    variance_model = x$variance_model,
    covariance_model = x$covariance_model %||% "diagonal",
    n_parameters = x$n_parameters,
    log_likelihood = x$log_likelihood,
    aic = x$aic,
    bic_groups = x$bic,
    bic_individual = x$bic_individual,
    converged = x$converged,
    boundary = x$boundary,
    extreme_logits = x$extreme_logits,
    n_starts = nrow(x$starts),
    row.names = NULL, stringsAsFactors = FALSE)
}

#' Print a multilevel LPA summary
#'
#' A human-facing report of the fit. The same content is available as tidy
#' tables from [get_data()].
#'
#' @param x A `summary_multilpa` object.
#' @param digits Number of printed significant digits.
#' @param rows How many rows of each table to print. A longer table is shown
#'   to that depth, with its remaining row count and the `get_data()` call that
#'   returns it whole.
#' @param ... Additional arguments passed to matrix printing.
#' @return The summary, invisibly. Called for the side effect of printing the
#'   estimates block by block, the effective class memberships at both levels,
#'   the likelihood and information criteria, any convergence or boundary
#'   warnings, and the restart diagnostics.
#' @seealso [get_data()] for the same content as data.
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
print.summary_multilpa <- function(x, digits = 4L, rows = 10L, ...) {
  stopifnot("`x` must be a `summary_multilpa` object" =
              inherits(x, "summary_multilpa"))
  .multilpa_check_print_arguments(digits, rows)
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
  cat(sprintf("Log likelihood: %.6f; AIC: %.3f\nBIC (groups): %.3f; BIC (individuals): %.3f\n",
              x$log_likelihood, x$aic, x$bic, x$bic_individual))
  cat(sprintf("Best likelihood replicated in %d/%d starts (absolute tolerance %.3g).\n",
              x$n_best_replicated, nrow(x$starts), x$replication_tolerance))
  if (!x$converged) cat("WARNING: best start did not converge.\n")
  if (x$boundary) cat("WARNING: at least one variance is at min_variance.\n")
  if (x$small_classes) cat("WARNING: an effective class membership is below one.\n")
  .multilpa_print_tables(x$tables, rows = rows, digits = digits)
  .multilpa_print_table_footer(x$tables)
  invisible(x)
}

#' Extract the multilevel model log likelihood
#' @param object An `multilpa` model.
#' @param ... Reserved for compatibility with `logLik()`.
#' @return A `logLik` object with parameter count `df` and the number of observed
#'   groups as `nobs`. Thus `stats::BIC()` uses group-count BIC. For the
#'   individual-count alternative, and every other criterion, call
#'   `get_data(x, "information_criteria")`, which reports both conventions
#'   side by side.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                 n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' logLik(fit)
#' AIC(fit)
#' BIC(fit)
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
#' @return A single integer: the number of observed groups, which are the
#'   independent units of the two-level likelihood. For the individual count
#'   alongside every other sample-size-dependent quantity, call
#'   `get_data(x, "information_criteria")`.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                 n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' nobs(fit)
#' @export
#' @importFrom stats nobs
nobs.multilpa <- function(object, ...) {
  stopifnot(inherits(object, "multilpa"))
  object$n_groups
}
