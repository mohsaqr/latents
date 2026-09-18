#' Tidy a fitted multilevel latent profile model
#'
#' Returns any part of a fitted model as a tidy `data.frame`, so that results
#' can be printed, joined, and written out without reaching into the fitted
#' object. Every table is one row per observation of the thing it describes.
#'
#' @param x A fitted `multilpa` model.
#' @param row.names Passed to `data.frame()`; `NULL` gives default row names.
#' @param optional Ignored, present for generic compatibility.
#' @param what Which table to return. `"profiles"` gives the Gaussian
#'   measurement model, `"responses"` the categorical measurement model,
#'   `"profile_probabilities"` the profile prevalence within each group class,
#'   `"posteriors"` the individual posteriors, `"group_posteriors"` the group
#'   posteriors, `"starts"` the restart diagnostics,
#'   `"information_criteria"` the information criteria, `"classification"` the
#'   classification quality, and `"entropy"` the entropy summary.
#' @param ... Passed to the underlying accessor; `classification` accepts
#'   `level` and `detail`.
#' @return A base `data.frame` whose columns depend on `what`.
#'   `"profiles"` has one row per profile and continuous indicator, with columns
#'   `profile`, `indicator`, `mean`, `variance`, and `standard_deviation`; it has
#'   zero rows when every indicator is categorical.
#'   `"responses"` has one row per profile, categorical indicator and category,
#'   with columns `profile`, `indicator`, `category`, `probability`, and
#'   `threshold`. The threshold is `qlogis(P(y <= category))`, the cumulative
#'   logit that mixture software reports, and is `NA_real_` for the final
#'   category of each indicator, where the cumulative probability is one.
#'   `"profile_probabilities"` has one row per group class and profile, with
#'   columns `group_class`, `profile`, `probability`, and
#'   `group_class_probability`. `"posteriors"` has one row per input individual,
#'   with columns `row`, `group`, `profile` (the modal profile), and one
#'   `posterior_profile_*` column per profile. `"group_posteriors"` has one row
#'   per observed group, with columns `group`, `group_size`, `group_class`, the
#'   group's `log_likelihood`, and one `posterior_group_class_*` column per
#'   group class. `"starts"` has one row per EM start. The remaining values
#'   return exactly [information_criteria()], [classification_table()], and
#'   [entropy_table()].
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                   n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' as.data.frame(fit)
#' as.data.frame(fit, what = "profile_probabilities")
#' @export
as.data.frame.multilpa <- function(x, row.names = NULL, optional = FALSE,
                                 what = c("profiles", "responses",
                                          "profile_probabilities",
                                          "posteriors", "group_posteriors",
                                          "starts", "information_criteria",
                                          "classification", "entropy"), ...) {
  stopifnot("`x` must be an `multilpa` fit" = inherits(x, "multilpa"))
  what <- match.arg(what)
  result <- switch(what,
    profiles = .multilpa_profile_frame(x),
    responses = .multilpa_response_frame(x),
    profile_probabilities = .multilpa_probability_frame(x),
    posteriors = .multilpa_posterior_frame(x),
    group_posteriors = .multilpa_group_posterior_frame(x),
    starts = x$starts,
    information_criteria = information_criteria(x),
    classification = classification_table(x, ...),
    entropy = entropy_table(x))
  row.names(result) <- row.names
  result
}

#' Names of the Gaussian indicators in a fit
#'
#' Fits made before categorical indicators existed carry no `continuous` field,
#' so their indicators are all Gaussian.
#'
#' @param x A fitted `multilpa` model.
#' @return Character vector of continuous indicator names.
#' @noRd
.multilpa_continuous_names <- function(x) {
  x$continuous %||% x$indicators
}

#' Measurement parameters as a tidy table
#' @param x A fitted `multilpa` model.
#' @return One row per profile and continuous indicator.
#' @noRd
.multilpa_profile_frame <- function(x) {
  indicators <- .multilpa_continuous_names(x)
  n_profiles <- x$n_profiles
  data.frame(
    profile = rep(seq_len(n_profiles), each = length(indicators)),
    indicator = rep(indicators, times = n_profiles),
    mean = as.vector(t(x$means)),
    variance = as.vector(t(x$variances)),
    # Derived rather than read from the fit, because not every result class
    # stores a standard-deviation matrix.
    standard_deviation = sqrt(as.vector(t(x$variances))))
}

#' Categorical response probabilities as a tidy table
#' @param x A fitted `multilpa` model.
#' @return One row per profile, categorical indicator and category.
#' @noRd
.multilpa_response_frame <- function(x) {
  stopifnot("`x` must be an `multilpa` fit" = inherits(x, "multilpa"))
  blocks <- x$response_probabilities
  if (is.null(blocks) || length(blocks) == 0L) {
    return(data.frame(profile = integer(), indicator = character(),
                      category = character(), probability = numeric(),
                      threshold = numeric()))
  }
  rows <- lapply(names(blocks), function(indicator) {
    block <- blocks[[indicator]]
    thresholds <- .multilpa_categorical_thresholds(block)
    # The final category has cumulative probability one, so its threshold is
    # infinite and is reported as missing rather than as a number.
    thresholds <- cbind(thresholds, NA_real_)
    data.frame(
      profile = rep(seq_len(nrow(block)), times = ncol(block)),
      indicator = indicator,
      category = rep(colnames(block), each = nrow(block)),
      probability = as.vector(block),
      threshold = as.vector(thresholds))
  })
  do.call(rbind, rows)
}

#' Profile prevalence within group classes as a tidy table
#' @param x A fitted `multilpa` model.
#' @return One row per group class and profile.
#' @noRd
.multilpa_probability_frame <- function(x) {
  stopifnot("`x` must be an `multilpa` fit" = inherits(x, "multilpa"))
  n_types <- x$n_group_classes
  n_profiles <- x$n_profiles
  data.frame(
    group_class = rep(seq_len(n_types), each = n_profiles),
    profile = rep(seq_len(n_profiles), times = n_types),
    probability = as.vector(t(x$profile_probabilities)),
    group_class_probability = rep(unname(x$group_probabilities), each = n_profiles))
}

#' Individual posteriors as a tidy table
#' @param x A fitted `multilpa` model.
#' @return One row per input individual, in the original row order.
#' @noRd
.multilpa_posterior_frame <- function(x) {
  posteriors <- as.data.frame(unname(x$subject_posteriors))
  names(posteriors) <- sprintf("posterior_profile_%d", seq_len(x$n_profiles))
  cbind(data.frame(row = seq_len(x$n_observations),
                   group = x$group_values[x$group_index],
                   profile = max.col(x$subject_posteriors,
                                     ties.method = "first")),
        posteriors)
}

#' Group posteriors as a tidy table
#' @param x A fitted `multilpa` model.
#' @return One row per observed group, in first-occurrence order.
#' @noRd
.multilpa_group_posterior_frame <- function(x) {
  posteriors <- as.data.frame(unname(x$group_posteriors))
  names(posteriors) <- sprintf("posterior_group_class_%d", seq_len(x$n_group_classes))
  # Sizes and modal classes are derived rather than read from the fit, because
  # not every result class stores them as fields.
  cbind(data.frame(group = x$group_values,
                   group_size = tabulate(x$group_index,
                                         nbins = length(x$group_values)),
                   group_class = max.col(x$group_posteriors,
                                         ties.method = "first"),
                   log_likelihood = unname(x$group_log_likelihood)),
        posteriors)
}

#' Tidy a class-enumeration grid
#'
#' @param x An `multilpa_enumeration` result from [enumerate_classes()].
#' @param row.names Passed to `data.frame()`; `NULL` gives default row names.
#' @param optional Ignored, present for generic compatibility.
#' @param ... Ignored.
#' @return A base `data.frame` with one row per candidate model in the grid,
#'   carrying its class counts, log likelihood, parameter count, every
#'   information criterion under both sample-size conventions, both entropies,
#'   and the convergence, boundary, replication, warning, and error diagnostics.
#'   Failed candidates are retained with `NA` estimates and their error text.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' candidates <- enumerate_classes(example_data, c("score_a", "score_b"), "school",
#'                                profiles = 1:2, group_classes = 1, n_starts = 2,
#'                                seed = 1)
#' as.data.frame(candidates)
#' @export
as.data.frame.multilpa_enumeration <- function(x, row.names = NULL,
                                             optional = FALSE, ...) {
  stopifnot("`x` must be an `multilpa_enumeration` result" =
              inherits(x, "multilpa_enumeration"))
  result <- x$table
  row.names(result) <- row.names
  result
}


#' Identify mixing-probability coefficient names
#' @param parameter Character vector of natural-scale parameter names.
#' @return Logical vector marking probability parameters.
#' @noRd
.multilpa_is_probability <- function(parameter) {
  stopifnot("`parameter` must be character" = is.character(parameter))
  grepl("^(profile|group)_probability\\[", parameter)
}

#' Print a class-enumeration grid
#' @param x An `multilpa_enumeration` result.
#' @param ... Passed to the underlying `data.frame` printing.
#' @return The input, invisibly.
#' @examples
#' # After enumerating: print(candidates)
#' @export
print.multilpa_enumeration <- function(x, ...) {
  stopifnot("`x` must be an `multilpa_enumeration` result" =
              inherits(x, "multilpa_enumeration"))
  table <- x$table
  cat(sprintf("Class enumeration: %d candidate models\n", nrow(table)))
  columns <- intersect(c("n_profiles", "n_group_classes", "log_likelihood",
                         "n_parameters", "bic_individual", "sabic_individual",
                         "profile_entropy", "converged"), names(table))
  print(table[, columns, drop = FALSE], row.names = FALSE, ...)
  failures <- sum(!table$converged)
  if (failures > 0L) {
    cat(sprintf("%d candidate(s) did not converge; see as.data.frame(x).\n", failures))
  }
  cat("No candidate is selected automatically. Compare one convention consistently.\n")
  invisible(x)
}


#' Build starting values for a multilevel latent profile fit
#'
#' Returns the starting-value list [multilpa()] accepts, taken from a fitted
#' model or from any object that already carries the parameter blocks, such as a
#' retained reference solution. This exists so that callers never assemble the
#' list by hand, and so that starting values are validated where they are built
#' rather than deep inside the fitting loop.
#'
#' @param object A fitted `multilpa` model, or a list carrying `means`,
#'   `variances` or `covariances`, `profile_probabilities`, and
#'   `group_probabilities`. Any other elements are ignored.
#' @param covariance `"auto"` keeps `covariances` when the object carries them,
#'   `"drop"` always returns the diagonal parameterization, and `"keep"`
#'   requires `covariances` and fails when they are absent.
#' @return A list with elements `means`, `variances`, `profile_probabilities`,
#'   and `group_probabilities`, plus `covariances` when the full-covariance
#'   parameterization is returned. Dimension names are dropped, matching what
#'   [multilpa()] expects of `start`.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                   n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' refit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                     n_profiles = 2, n_group_classes = 1, n_starts = 1,
#'                     start = starting_values(fit))
#' logLik(refit)
#' @export
starting_values <- function(object, covariance = c("auto", "drop", "keep")) {
  stopifnot("`object` must be a list or an `multilpa` fit" = is.list(object))
  covariance <- match.arg(covariance)
  required <- c("means", "profile_probabilities", "group_probabilities")
  missing_fields <- setdiff(required, names(object))
  if (length(missing_fields) > 0L) {
    stop(errorCondition(sprintf("`object` is missing starting values for %s.",
                                paste(missing_fields, collapse = ", ")),
                        class = "multilpa_bad_start", call = NULL))
  }
  has_covariances <- !is.null(object$covariances)
  if (identical(covariance, "keep") && !has_covariances) {
    stop(errorCondition("`object` carries no covariances to keep.",
                        class = "multilpa_bad_start", call = NULL))
  }
  use_covariances <- has_covariances && !identical(covariance, "drop")
  variances <- object$variances
  if (is.null(variances)) {
    if (!has_covariances) {
      stop(errorCondition("`object` is missing starting values for variances.",
                          class = "multilpa_bad_start", call = NULL))
    }
    dimension <- dim(object$covariances)[1L]
    variances <- t(vapply(seq_len(dim(object$covariances)[3L]), function(profile) {
      diag(matrix(object$covariances[, , profile], dimension, dimension))
    }, numeric(dimension)))
  }
  start <- list(means = unname(as.matrix(object$means)),
                variances = unname(as.matrix(variances)),
                profile_probabilities = unname(as.matrix(object$profile_probabilities)),
                group_probabilities = unname(as.vector(object$group_probabilities)))
  if (use_covariances) {
    start$covariances <- unname(object$covariances)
    start$variances <- NULL
  }
  start
}

#' Tidy a one-step membership-covariate fit
#'
#' @param x An `multilpa_covariates` model from [fit_covariates()].
#' @param row.names Passed to `data.frame()`; `NULL` gives default row names.
#' @param optional Ignored, present for generic compatibility.
#' @param what Which table to return. `"profiles"` gives the measurement model,
#'   `"coefficients"` the membership regressions at both levels, `"posteriors"`
#'   the individual posteriors, and `"group_posteriors"` the group posteriors.
#' @param ... Ignored.
#' @return A base `data.frame`. `"profiles"` has one row per profile and
#'   indicator, with columns `profile`, `indicator`, `mean`, `variance`, and
#'   `standard_deviation`. `"coefficients"` has one row per estimated
#'   membership coefficient, with columns `level` (`"profile"` or `"group"`),
#'   `outcome` (the class the coefficient predicts), `term`, and `estimate`;
#'   these models carry no standard errors, so none is reported. `"posteriors"`
#'   and `"group_posteriors"` match the corresponding [as.data.frame.multilpa()]
#'   tables.
#' @examples
#' # After fitting: as.data.frame(with_predictors, what = "coefficients")
#' @export
as.data.frame.multilpa_covariates <- function(x, row.names = NULL, optional = FALSE,
                                            what = c("profiles", "coefficients",
                                                     "posteriors",
                                                     "group_posteriors"), ...) {
  stopifnot("`x` must be an `multilpa_covariates` fit" =
              inherits(x, "multilpa_covariates"))
  what <- match.arg(what)
  result <- switch(what,
    profiles = .multilpa_profile_frame(x),
    coefficients = .multilpa_coefficient_frame(x),
    posteriors = .multilpa_posterior_frame(x),
    group_posteriors = .multilpa_group_posterior_frame(x))
  row.names(result) <- row.names
  result
}

#' Membership regression coefficients as a tidy table
#' @return One row per estimated coefficient at either level.
#' @noRd
.multilpa_coefficient_frame <- function(x) {
  blocks <- list(profile = x$profile_coefficients, group = x$group_coefficients)
  rows <- lapply(names(blocks), function(level) {
    block <- blocks[[level]]
    if (is.null(block) || length(block) == 0L) return(NULL)
    data.frame(level = level,
               outcome = rep(colnames(block), each = nrow(block)),
               term = rep(rownames(block), times = ncol(block)),
               estimate = as.vector(block))
  })
  result <- do.call(rbind, rows)
  if (is.null(result)) {
    return(data.frame(level = character(), outcome = character(),
                      term = character(), estimate = numeric()))
  }
  result
}

#' Tidy a continuous group random-intercept fit
#'
#' @param x An `multilpa_random_intercept` model from
#'   [fit_random_intercept()].
#' @param row.names Passed to `data.frame()`; `NULL` gives default row names.
#' @param optional Ignored, present for generic compatibility.
#' @param what Which table to return. `"profiles"` gives the measurement model,
#'   `"posteriors"` the individual posteriors, and `"random_intercepts"` the
#'   posterior group effects.
#' @param ... Ignored.
#' @return A base `data.frame`. `"profiles"` has one row per profile and
#'   indicator, with columns `profile`, `indicator`, `mean`, `variance`, and
#'   `standard_deviation`. `"random_intercepts"` has one row per observed group,
#'   with columns `group`, `group_size`, `mean` and `sd`, the posterior mean and
#'   standard deviation of that group's scalar intercept. `"posteriors"` has one
#'   row per individual, with columns `row`, `group`, `profile`, and one
#'   `posterior_profile_*` column per profile.
#' @examples
#' # After fitting: as.data.frame(random_intercept, what = "random_intercepts")
#' @export
as.data.frame.multilpa_random_intercept <- function(x, row.names = NULL,
                                                  optional = FALSE,
                                                  what = c("profiles",
                                                           "posteriors",
                                                           "random_intercepts"),
                                                  ...) {
  stopifnot("`x` must be an `multilpa_random_intercept` fit" =
              inherits(x, "multilpa_random_intercept"))
  what <- match.arg(what)
  result <- switch(what,
    profiles = .multilpa_profile_frame(x),
    posteriors = .multilpa_posterior_frame(x),
    random_intercepts = data.frame(
      group = x$group_values,
      group_size = unname(x$group_sizes),
      mean = unname(x$random_intercept_mean),
      sd = unname(x$random_intercept_sd)))
  row.names(result) <- row.names
  result
}
