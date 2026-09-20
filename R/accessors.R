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
#'   `"stages"` the estimation stages and what each held fixed,
#'   `"information_criteria"` the information criteria, `"classification"` the
#'   classification quality, and `"entropy"` the entropy summary.
#' @param scale For `what = "profiles"`, `"raw"` reports the estimates in input
#'   units and `"standardized"` divides each indicator's deviation from its
#'   grand mean by that indicator's observed standard deviation, which is what
#'   `plot(x, scale = "standardized")` draws. Supplying it with any other
#'   `what` is an error rather than a silent no-op.
#' @param format For the posterior tables, `"long"` (the default) gives one row
#'   per unit and class, and `"wide"` gives one row per unit with one column per
#'   class. The wide column set grows with the number of classes, so no caller
#'   can address it generically; it stays available because joining posteriors
#'   back onto the fitting data wants one row per unit. Supplying it with any
#'   other `what` is an error rather than a silent no-op.
#' @param ... Passed to the underlying accessor. `classification` accepts
#'   `level`, and `information_criteria` accepts `definitions`. The
#'   per-assigned-class breakdown that `classification_table(detail = TRUE)`
#'   once returned is now its own verb, [average_posteriors()]. `profiles` and `responses` accept `data`, the data
#'   frame the model was fitted to, which adds a standard error beside every
#'   estimate so the measurement model can be read complete from one call
#'   rather than joined to [parameter_inference()] by hand.
#' @return A base `data.frame` whose columns depend on `what`.
#'   `"profiles"` has one row per profile and continuous indicator, with columns
#'   `profile`, `indicator`, `mean`, `variance`, and `standard_deviation`; it has
#'   zero rows when every indicator is categorical. Given `data`, it also carries
#'   `mean_standard_error` and `variance_standard_error`, and `"responses"`
#'   carries `probability_standard_error`. Supplying `data` runs
#'   [parameter_inference()], so it raises whatever that would raise for a fit
#'   whose standard errors are unavailable, rather than returning empty columns.
#'   With `scale = "standardized"` every estimate in that table, and every
#'   standard error beside it, is divided by the indicator's observed standard
#'   deviation, so the whole row is on one scale; it raises
#'   `multilpa_no_indicator_data` for a fit that kept no indicators and
#'   `multilpa_bad_scale` for an indicator with no spread.
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
#'   group class. `"starts"` has one row per EM start, with columns `start`,
#'   `log_likelihood`, `converged`, `iterations`, `error` and `boundary`.
#'   `"stages"` has one row per estimation stage, with columns `stage`,
#'   `group_classes`, `fixed`, `log_likelihood`, `parameters`,
#'   `parameters_with_measurement` and `converged`; an ordinary fit has a single
#'   `"joint"` row and a [fit_staged()] result has two. `fixed` names the blocks
#'   that stage held fixed, comma separated, and is `NA_character_` when the
#'   stage estimated every block, which is the case for an ordinary fit and for
#'   the measurement stage of a staged one; it is never a sentinel such as
#'   `"none"`. The remaining values
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
#' as.data.frame(fit, scale = "standardized")
#' as.data.frame(fit, what = "profile_probabilities")
#' @export
as.data.frame.multilpa <- function(x, row.names = NULL, optional = FALSE,
                                 what = c("profiles", "responses",
                                          "profile_probabilities",
                                          "posteriors", "group_posteriors",
                                          "starts", "stages",
                                          "information_criteria",
                                          "classification", "entropy"),
                                 scale = c("raw", "standardized"),
                                 format = c("long", "wide"), ...) {
  stopifnot("`x` must be an `multilpa` fit" = inherits(x, "multilpa"))
  what <- match.arg(what)
  scale <- match.arg(scale)
  format <- match.arg(format)
  stopifnot("`scale` applies only to `what = \"profiles\"`" =
              identical(scale, "raw") || identical(what, "profiles"),
            "`format` applies only to the posterior tables" =
              identical(format, "long") ||
              what %in% c("posteriors", "group_posteriors"))
  result <- switch(what,
    profiles = .multilpa_profile_frame(x, scale = scale, ...),
    responses = .multilpa_response_frame(x, ...),
    profile_probabilities = .multilpa_probability_frame(x),
    posteriors = .multilpa_posterior_frame(x, format),
    group_posteriors = .multilpa_group_posterior_frame(x, format),
    starts = x$starts,
    stages = .multilpa_stage_frame(x),
    information_criteria = information_criteria(x, ...),
    classification = classification_table(x, ...),
    entropy = entropy_table(x))
  row.names(result) <- row.names
  result
}

#' The stages of a fit as a tidy table
#'
#' An ordinary fit has one stage, estimated jointly; a staged fit has two, the
#' second holding the first's measurement solution fixed. Both are described by
#' the same columns, so a caller need not branch on which kind of fit it holds.
#'
#' @param x A fitted `multilpa` model.
#' @return One row per estimation stage.
#' @noRd
.multilpa_stage_frame <- function(x) {
  describe <- function(stage, fit, held) {
    data.frame(stage = stage, group_classes = fit$n_group_classes,
               fixed = if (length(held) == 0L) NA_character_ else
                 paste(held, collapse = ", "),
               log_likelihood = fit$log_likelihood,
               parameters = fit$n_parameters,
               parameters_with_measurement =
                 fit$n_parameters_with_measurement %||% fit$n_parameters,
               converged = fit$converged, row.names = NULL)
  }
  if (is.null(x$stage_one)) {
    return(describe("joint", x, x$fixed %||% character()))
  }
  rbind(describe("measurement", x$stage_one, character()),
        describe("membership", x, x$fixed %||% character()))
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
  x$continuous %||% x$vars
}

#' Measurement parameters as a tidy table
#' @param x A fitted `multilpa` model.
#' @return One row per profile and continuous indicator.
#' @noRd
.multilpa_profile_frame <- function(x, data = NULL,
                                    scale = c("raw", "standardized")) {
  scale <- match.arg(scale)
  vars <- .multilpa_continuous_names(x)
  n_profiles <- x$n_profiles
  # `spread` of one and `centre` of zero leave the raw scale untouched, so the
  # two scales share one code path and cannot drift apart.
  basis <- if (identical(scale, "standardized") && length(vars) > 0L) {
    .multilpa_standardization(x, vars)
  } else {
    list(centre = rep(0, length(vars)),
         spread = rep(1, length(vars)))
  }
  means <- sweep(sweep(x$means, 2L, basis$centre, "-"), 2L, basis$spread, "/")
  variances <- sweep(x$variances, 2L, basis$spread^2, "/")
  frame <- data.frame(
    profile = rep(seq_len(n_profiles), each = length(vars)),
    indicator = rep(vars, times = n_profiles),
    mean = as.vector(t(means)),
    variance = as.vector(t(variances)),
    # Derived rather than read from the fit, because not every result class
    # stores a standard-deviation matrix.
    standard_deviation = sqrt(as.vector(t(variances))))
  if (is.null(data)) return(frame)
  errors <- .multilpa_measurement_errors(x, data)
  row_spread <- rep(basis$spread, times = n_profiles)
  # Standardizing divides each estimate by a constant, so its standard error
  # divides by the same constant; a table must not mix the two scales.
  frame$mean_standard_error <- .multilpa_match_error(errors, "mean",
    frame$profile, frame$indicator) / row_spread
  frame$variance_standard_error <- .multilpa_match_error(errors, "variance",
    frame$profile, frame$indicator) / row_spread^2
  frame
}

#' The observed centre and spread each indicator is standardized by
#'
#' The same grand mean and observed standard deviation that
#' `plot(scale = "standardized")` divides by, so the plotted shape and the
#' table of numbers behind it cannot disagree. The spread is the observed
#' standard deviation of the indicator, not the within-profile residual one,
#' so the standardized means are comparable across indicators without being
#' effect sizes.
#'
#' @param x A fitted model of this package.
#' @param vars The continuous indicator names of that fit.
#' @return A list of the numeric vectors `centre` and `spread`, one element per
#'   continuous indicator.
#' @noRd
.multilpa_standardization <- function(x, vars) {
  observed <- x$indicator_data
  if (is.null(observed) || ncol(observed) != length(vars)) {
    stop(errorCondition(
      "This fit did not retain indicator data, so it cannot be standardized.",
      class = "multilpa_no_indicator_data", call = NULL))
  }
  spread <- vapply(seq_along(vars), function(index) {
    stats::sd(observed[, index], na.rm = TRUE)
  }, numeric(1))
  if (any(!is.finite(spread)) || any(spread <= 0)) {
    stop(errorCondition(
      "An indicator has zero or undefined standard deviation.",
      class = "multilpa_bad_scale", call = NULL))
  }
  list(centre = unname(colMeans(observed, na.rm = TRUE)), spread = spread)
}

#' Standard errors for the measurement parameters of a fit
#'
#' [parameter_inference()] reports one row per free parameter, which is the
#' right shape for reading coefficients and the wrong one for reading a
#' measurement model. This fetches them so the measurement tables can carry
#' their own errors in the shape they already have.
#'
#' @param x A fitted model of this package.
#' @param data The data frame the model was fitted to.
#' @return The inference table, restricted to measurement parameters.
#' @noRd
.multilpa_measurement_errors <- function(x, data) {
  inference <- parameter_inference(x, data)
  inference[inference$level == "measurement", , drop = FALSE]
}

#' Look up one measurement standard error per requested cell
#' @param errors Measurement rows of an inference table.
#' @param parameter Which kind of parameter to take.
#' @param profile Integer profile indices wanted, one per output row.
#' @param term Term labels wanted, one per output row.
#' @return A numeric vector, `NA_real_` where the fit reports no such parameter.
#' @noRd
.multilpa_match_error <- function(errors, parameter, profile, term) {
  wanted <- errors[errors$parameter == parameter, , drop = FALSE]
  if (nrow(wanted) == 0L) return(rep(NA_real_, length(profile)))
  key <- paste(sub("^profile_", "", wanted$outcome), wanted$term, sep = "\r")
  wanted$standard_error[match(paste(profile, term, sep = "\r"), key)]
}

#' Categorical response probabilities as a tidy table
#' @param x A fitted `multilpa` model.
#' @return One row per profile, categorical indicator and category.
#' @noRd
.multilpa_response_frame <- function(x, data = NULL) {
  stopifnot("`x` must be a fitted model of this package, or its summary" =
              .multilpa_any_fit(x) || inherits(x, "summary_multilpa"))
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
  frame <- do.call(rbind, rows)
  if (is.null(data)) return(frame)
  errors <- .multilpa_measurement_errors(x, data)
  frame$probability_standard_error <- .multilpa_match_error(errors, "response",
    frame$profile, paste(frame$indicator, frame$category, sep = ":"))
  frame
}

#' Profile prevalence within group classes as a tidy table
#' @param x A fitted `multilpa` model.
#' @return One row per group class and profile.
#' @noRd
.multilpa_probability_frame <- function(x) {
  stopifnot("`x` must be an `multilpa` fit, or its summary" =
              inherits(x, "multilpa") || inherits(x, "summary_multilpa"))
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
.multilpa_posterior_frame <- function(x, format = "long") {
  stopifnot("`format` must be \"long\" or \"wide\"" =
              format %in% c("long", "wide"))
  modal <- max.col(x$subject_posteriors, ties.method = "first")
  if (identical(format, "wide")) {
    posteriors <- as.data.frame(unname(x$subject_posteriors))
    names(posteriors) <- sprintf("posterior_profile_%d", seq_len(x$n_profiles))
    return(cbind(data.frame(row = seq_len(x$n_observations),
                            group = x$group_values[x$group_index],
                            profile = modal),
                 posteriors))
  }
  # One row per individual and profile. The wide form's column set grows with the
  # number of profiles, so no caller can address it generically; it stays
  # available because joining posteriors back onto the fitting data wants one row
  # per individual, but it is no longer what the accessor hands back by default.
  data.frame(
    row = rep(seq_len(x$n_observations), times = x$n_profiles),
    group = rep(x$group_values[x$group_index], times = x$n_profiles),
    profile = rep(seq_len(x$n_profiles), each = x$n_observations),
    posterior = as.vector(x$subject_posteriors),
    modal = rep(seq_len(x$n_profiles), each = x$n_observations) ==
      rep(modal, times = x$n_profiles))
}

#' Group posteriors as a tidy table
#' @param x A fitted `multilpa` model.
#' @return One row per observed group, in first-occurrence order.
#' @noRd
.multilpa_group_posterior_frame <- function(x, format = "long") {
  stopifnot("`format` must be \"long\" or \"wide\"" =
              format %in% c("long", "wide"))
  modal <- max.col(x$group_posteriors, ties.method = "first")
  sizes <- tabulate(x$group_index, nbins = length(x$group_values))
  if (identical(format, "long")) {
    n_groups <- length(x$group_values)
    return(data.frame(
      group = rep(x$group_values, times = x$n_group_classes),
      group_size = rep(sizes, times = x$n_group_classes),
      log_likelihood = rep(unname(x$group_log_likelihood),
                           times = x$n_group_classes),
      group_class = rep(seq_len(x$n_group_classes), each = n_groups),
      posterior = as.vector(x$group_posteriors),
      modal = rep(seq_len(x$n_group_classes), each = n_groups) ==
        rep(modal, times = x$n_group_classes)))
  }
  posteriors <- as.data.frame(unname(x$group_posteriors))
  names(posteriors) <- sprintf("posterior_group_class_%d", seq_len(x$n_group_classes))
  # Sizes and modal classes are derived rather than read from the fit, because
  # not every result class stores them as fields.
  cbind(data.frame(group = x$group_values,
                   group_size = sizes,
                   group_class = modal,
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
#'                                n_profiles = 1:2, n_group_classes = 1, n_starts = 2,
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
#' @param x A fitted `multilpa` model, or a list carrying
#'   `profile_probabilities`, `group_probabilities`, and whichever measurement
#'   blocks the model uses: `means` with `variances` or `covariances` for
#'   continuous indicators, `response_probabilities` for categorical ones. An
#'   all-categorical model needs no Gaussian block. Any other elements are
#'   ignored.
#' @param what `"all"` returns the measurement and mixing blocks.
#'   `"measurement"` omits `profile_probabilities` and `group_probabilities`,
#'   which are sized for the model that produced them and must not be carried
#'   into a staged fit that changes the number of group classes. Pass the
#'   result as `start` alongside `fixed` in [multilpa()].
#' @param covariance `"auto"` keeps `covariances` when the object carries them,
#'   `"drop"` always returns the diagonal parameterization, and `"keep"`
#'   requires `covariances` and fails when they are absent.
#' @return An object of class `multilpa_start`: the list [multilpa()] accepts as
#'   `start`, with elements `means`, `variances`, `profile_probabilities`, and
#'   `group_probabilities`, plus `covariances` when the full-covariance
#'   parameterization is returned and `response_probabilities` when the object
#'   carries a categorical measurement model. `what = "measurement"` drops
#'   `profile_probabilities` and `group_probabilities`. Dimension names are
#'   dropped, matching what [multilpa()] expects of `start`, so the values are
#'   indexed by position; read them with [as.data.frame()], which labels the
#'   positions, rather than out of the list. The list payload is unchanged by
#'   the class, so a `multilpa_start` can be passed straight back as `start`.
#' @seealso [as.data.frame.multilpa_start()] for the values as tidy tables.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                   n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' starting_values(fit)
#' as.data.frame(starting_values(fit))
#' refit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                     n_profiles = 2, n_group_classes = 1, n_starts = 1,
#'                     start = starting_values(fit))
#' logLik(refit)
#' @export
starting_values <- function(x, covariance = c("auto", "drop", "keep"),
                            what = c("all", "measurement")) {
  stopifnot("`object` must be a list or an `multilpa` fit" = is.list(x))
  covariance <- match.arg(covariance)
  what <- match.arg(what)
  has_responses <- !is.null(x$response_probabilities)
  # A model with only categorical indicators carries no Gaussian block.
  required <- c("means", "profile_probabilities", "group_probabilities")
  if (has_responses && is.null(x$means)) {
    required <- setdiff(required, "means")
  }
  missing_fields <- setdiff(required, names(x))
  if (length(missing_fields) > 0L) {
    stop(errorCondition(sprintf("`object` is missing starting values for %s.",
                                paste(missing_fields, collapse = ", ")),
                        class = "multilpa_bad_start", call = NULL))
  }
  has_covariances <- !is.null(x$covariances)
  if (identical(covariance, "keep") && !has_covariances) {
    stop(errorCondition("`object` carries no covariances to keep.",
                        class = "multilpa_bad_start", call = NULL))
  }
  use_covariances <- has_covariances && !identical(covariance, "drop")
  n_profiles <- ncol(as.matrix(x$profile_probabilities))
  means <- x$means %||% matrix(0, n_profiles, 0L)
  variances <- x$variances
  if (is.null(variances)) {
    if (has_responses && ncol(as.matrix(means)) == 0L) {
      variances <- matrix(0, n_profiles, 0L)
    } else if (!has_covariances) {
      stop(errorCondition("`object` is missing starting values for variances.",
                          class = "multilpa_bad_start", call = NULL))
    } else {
    dimension <- dim(x$covariances)[1L]
    variances <- t(matrix(vapply(seq_len(dim(x$covariances)[3L]), function(profile) {
      diag(matrix(x$covariances[, , profile], dimension, dimension))
    }, numeric(dimension)), nrow = dimension, ncol = n_profiles))
    }
  }
  start <- list(means = unname(as.matrix(means)),
                variances = unname(as.matrix(variances)),
                profile_probabilities = unname(as.matrix(x$profile_probabilities)),
                group_probabilities = unname(as.vector(x$group_probabilities)))
  if (use_covariances) {
    start$covariances <- unname(x$covariances)
    start$variances <- NULL
  }
  if (has_responses) {
    start$response_probabilities <- unname(lapply(x$response_probabilities,
                                                  \(block) unname(as.matrix(block))))
  }
  # The mixing blocks are sized for the model that produced them, so a staged
  # fit that changes the number of group classes must not carry them across.
  if (identical(what, "measurement")) {
    start[c("profile_probabilities", "group_probabilities")] <- NULL
  }
  # Classed so that the blocks can be printed and tabulated; the payload is
  # the plain list `multilpa(start = )` validates, and a class attribute does
  # not change how that validation reads it.
  class(start) <- "multilpa_start"
  start
}

#' Print a set of starting values
#'
#' @param x A `multilpa_start` object from [starting_values()].
#' @param ... Reserved for compatibility with `print()`.
#' @return The input, invisibly.
#' @seealso [as.data.frame.multilpa_start()] for the values as tidy tables.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                   n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' print(starting_values(fit))
#' @export
print.multilpa_start <- function(x, ...) {
  stopifnot("`x` must be a `multilpa_start` object" =
              inherits(x, "multilpa_start"))
  n_profiles <- nrow(as.matrix(x$means))
  n_indicators <- ncol(as.matrix(x$means))
  cat(sprintf("Starting values: %d profile(s), %d continuous indicator(s)\n",
              n_profiles, n_indicators))
  cat(sprintf("Residual covariance: %s\n",
              if (is.null(x$covariances)) "diagonal" else "full"))
  cat(sprintf("Categorical indicators: %d\n",
              length(x$response_probabilities)))
  cat(sprintf("Mixing blocks: %s\n",
              if (is.null(x$profile_probabilities)) "none (measurement only)"
              else sprintf("%d group class(es)",
                           length(x$group_probabilities))))
  cat("Pass this to multilpa(start = ) as it is.\n")
  invisible(x)
}

#' Tidy a set of starting values
#'
#' [starting_values()] drops dimension names, because that is what [multilpa()]
#' expects of `start`. This labels the positions again, so the values can be
#' read and compared without indexing the list.
#'
#' @param x A `multilpa_start` object from [starting_values()].
#' @param row.names Passed to `data.frame()`; `NULL` gives default row names.
#' @param optional Ignored, present for generic compatibility.
#' @param what Which table to return. `"profiles"` gives the Gaussian
#'   measurement block, `"responses"` the categorical measurement block,
#'   `"profile_probabilities"` the mixing block, and `"covariances"` the
#'   within-profile residual covariance matrices.
#' @param ... Ignored.
#' @return A base `data.frame`. `"profiles"` has one row per profile and
#'   continuous indicator, with columns `profile`, `indicator` (the indicator's
#'   position, since the names are dropped), `mean`, `variance` and
#'   `standard_deviation`, and zero rows when the block is empty. `"responses"`
#'   has one row per profile, categorical indicator and category, with columns
#'   `profile`, `indicator`, `category` and `probability`, and zero rows when
#'   the object carries no categorical block. `"profile_probabilities"` has one
#'   row per group class and profile, with columns `group_class`, `profile`,
#'   `probability` and `group_class_probability`, and zero rows for a
#'   measurement-only object, which carries no mixing block. `"covariances"`
#'   has one row per profile and ordered pair of continuous indicators, with
#'   columns `profile`, `indicator`, `indicator_2` and `covariance`; under the
#'   diagonal parameterization the off-diagonal covariances are zero by
#'   assumption rather than missing.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                   n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' as.data.frame(starting_values(fit))
#' as.data.frame(starting_values(fit), what = "profile_probabilities")
#' @export
as.data.frame.multilpa_start <- function(x, row.names = NULL, optional = FALSE,
                                         what = c("profiles", "responses",
                                                  "profile_probabilities",
                                                  "covariances"), ...) {
  stopifnot("`x` must be a `multilpa_start` object" =
              inherits(x, "multilpa_start"))
  what <- match.arg(what)
  result <- switch(what,
    profiles = .multilpa_start_profile_frame(x),
    responses = .multilpa_start_response_frame(x),
    profile_probabilities = .multilpa_start_probability_frame(x),
    covariances = .multilpa_start_covariance_frame(x))
  row.names(result) <- row.names
  result
}

#' The Gaussian block of a starting-value set as a tidy table
#' @param x A `multilpa_start` object.
#' @return One row per profile and continuous indicator.
#' @noRd
.multilpa_start_profile_frame <- function(x) {
  means <- as.matrix(x$means)
  variances <- .multilpa_start_variances(x)
  if (ncol(means) == 0L) {
    return(data.frame(profile = integer(), indicator = integer(),
                      mean = numeric(), variance = numeric(),
                      standard_deviation = numeric()))
  }
  data.frame(profile = rep(seq_len(nrow(means)), each = ncol(means)),
             indicator = rep(seq_len(ncol(means)), times = nrow(means)),
             mean = as.vector(t(means)),
             variance = as.vector(t(variances)),
             standard_deviation = sqrt(as.vector(t(variances))))
}

#' The variances a starting-value set implies
#'
#' The full-covariance parameterization carries no `variances` block, but its
#' covariance array states them on its diagonal.
#'
#' @param x A `multilpa_start` object.
#' @return A profiles by indicators matrix.
#' @noRd
.multilpa_start_variances <- function(x) {
  if (!is.null(x$variances)) return(as.matrix(x$variances))
  n_indicators <- dim(x$covariances)[1L]
  t(vapply(seq_len(dim(x$covariances)[3L]), function(profile) {
    diag(matrix(x$covariances[, , profile], n_indicators, n_indicators))
  }, numeric(n_indicators)))
}

#' The categorical block of a starting-value set as a tidy table
#' @param x A `multilpa_start` object.
#' @return One row per profile, categorical indicator and category.
#' @noRd
.multilpa_start_response_frame <- function(x) {
  blocks <- x$response_probabilities
  if (length(blocks) == 0L) {
    return(data.frame(profile = integer(), indicator = integer(),
                      category = integer(), probability = numeric()))
  }
  rows <- lapply(seq_along(blocks), function(indicator) {
    block <- as.matrix(blocks[[indicator]])
    data.frame(profile = rep(seq_len(nrow(block)), times = ncol(block)),
               indicator = indicator,
               category = rep(seq_len(ncol(block)), each = nrow(block)),
               probability = as.vector(block))
  })
  do.call(rbind, rows)
}

#' The mixing block of a starting-value set as a tidy table
#' @param x A `multilpa_start` object.
#' @return One row per group class and profile.
#' @noRd
.multilpa_start_probability_frame <- function(x) {
  # `what = "measurement"` deliberately drops the mixing blocks, so an empty
  # table is the honest answer rather than an error.
  if (is.null(x$profile_probabilities)) {
    return(data.frame(group_class = integer(), profile = integer(),
                      probability = numeric(),
                      group_class_probability = numeric()))
  }
  probabilities <- as.matrix(x$profile_probabilities)
  data.frame(group_class = rep(seq_len(nrow(probabilities)),
                               each = ncol(probabilities)),
             profile = rep(seq_len(ncol(probabilities)),
                           times = nrow(probabilities)),
             probability = as.vector(t(probabilities)),
             group_class_probability = rep(unname(x$group_probabilities),
                                           each = ncol(probabilities)))
}

#' The residual covariances a starting-value set implies
#' @param x A `multilpa_start` object.
#' @return One row per profile and ordered pair of continuous indicators.
#' @noRd
.multilpa_start_covariance_frame <- function(x) {
  variances <- .multilpa_start_variances(x)
  n_indicators <- ncol(variances)
  if (n_indicators == 0L) {
    return(data.frame(profile = integer(), indicator = integer(),
                      indicator_2 = integer(), covariance = numeric()))
  }
  blocks <- lapply(seq_len(nrow(variances)), function(profile) {
    block <- if (is.null(x$covariances)) {
      diag(variances[profile, ], n_indicators, n_indicators)
    } else {
      matrix(x$covariances[, , profile], n_indicators, n_indicators)
    }
    data.frame(profile = profile,
               indicator = rep(seq_len(n_indicators), times = n_indicators),
               indicator_2 = rep(seq_len(n_indicators), each = n_indicators),
               covariance = as.vector(block))
  })
  do.call(rbind, blocks)
}

#' Tidy a one-step membership-covariate fit
#'
#' @param x An `multilpa_covariates` model from [fit_covariates()].
#' @param row.names Passed to `data.frame()`; `NULL` gives default row names.
#' @param optional Ignored, present for generic compatibility.
#' @param what Which table to return. `"profiles"` gives the measurement model,
#'   `"coefficients"` the membership regressions at both levels, `"posteriors"`
#'   the individual posteriors, and `"group_posteriors"` the group posteriors.
#' @param scale For `what = "profiles"`, `"raw"` reports the estimates in input
#'   units and `"standardized"` divides each indicator's deviation from its
#'   grand mean by that indicator's observed standard deviation, matching
#'   `plot(x, scale = "standardized")`.
#' @param format For the posterior tables, `"long"` (the default) gives one row
#'   per unit and class, and `"wide"` gives one row per unit with one column per
#'   class. The wide column set grows with the number of classes, so no caller
#'   can address it generically; it stays available because joining posteriors
#'   back onto the fitting data wants one row per unit. Supplying it with any
#'   other `what` is an error rather than a silent no-op.
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
                                                     "group_posteriors",
                                                     "starts"),
                                            scale = c("raw", "standardized"),
                                            format = c("long", "wide"), ...) {
  stopifnot("`x` must be an `multilpa_covariates` fit" =
              inherits(x, "multilpa_covariates"))
  what <- match.arg(what)
  scale <- match.arg(scale)
  format <- match.arg(format)
  stopifnot("`scale` applies only to `what = \"profiles\"`" =
              identical(scale, "raw") || identical(what, "profiles"),
            "`format` applies only to the posterior tables" =
              identical(format, "long") ||
              what %in% c("posteriors", "group_posteriors"))
  result <- switch(what,
    profiles = .multilpa_profile_frame(x, scale = scale),
    coefficients = .multilpa_coefficient_frame(x),
    posteriors = .multilpa_posterior_frame(x, format),
    group_posteriors = .multilpa_group_posterior_frame(x, format),
    starts = x$starts)
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
               # These are multinomial logits, not probabilities. Naming the
               # parameter is what stops the table being read as the latter, and
               # it matches the decomposition parameter_inference() reports.
               parameter = "logit",
               estimate = as.vector(block))
  })
  result <- do.call(rbind, rows)
  if (is.null(result)) {
    return(data.frame(level = character(), outcome = character(),
                      term = character(), parameter = character(),
                      estimate = numeric()))
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
#' @param scale For `what = "profiles"`, `"raw"` reports the estimates in input
#'   units and `"standardized"` divides each indicator's deviation from its
#'   grand mean by that indicator's observed standard deviation.
#' @param format For the posterior tables, `"long"` (the default) gives one row
#'   per unit and class, and `"wide"` gives one row per unit with one column per
#'   class. The wide column set grows with the number of classes, so no caller
#'   can address it generically; it stays available because joining posteriors
#'   back onto the fitting data wants one row per unit. Supplying it with any
#'   other `what` is an error rather than a silent no-op.
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
                                                           "random_intercepts",
                                                           "starts"),
                                                  scale = c("raw",
                                                            "standardized"),
                                                  format = c("long", "wide"),
                                                  ...) {
  stopifnot("`x` must be an `multilpa_random_intercept` fit" =
              inherits(x, "multilpa_random_intercept"))
  what <- match.arg(what)
  scale <- match.arg(scale)
  format <- match.arg(format)
  stopifnot("`scale` applies only to `what = \"profiles\"`" =
              identical(scale, "raw") || identical(what, "profiles"),
            "`format` applies only to the posterior tables" =
              identical(format, "long") || identical(what, "posteriors"))
  result <- switch(what,
    profiles = .multilpa_profile_frame(x, scale = scale),
    posteriors = .multilpa_posterior_frame(x, format),
    random_intercepts = data.frame(
      group = x$group_values,
      group_size = unname(x$group_sizes),
      mean = unname(x$random_intercept_mean),
      sd = unname(x$random_intercept_sd)),
    starts = x$starts)
  row.names(result) <- row.names
  result
}
