#' Coerce a fitted multilevel latent profile model to its primary table
#'
#' Plain coercion, as the base generic means it: one object, one data frame.
#' The other tables are named rather than positional, so they belong to
#' [get_results()], which takes `what` and refuses a name this object has not.
#'
#' @param x A fitted multilevel latent profile model.
#' @param row.names Passed to `data.frame()`; `NULL` gives default row names.
#' @param optional Ignored, present for generic compatibility.
#' @param ... Must be empty. An argument here raises `multilpa_bad_argument`
#'   naming it, rather than being dropped, because `what =` used to live on
#'   this generic and silently returning the primary table instead of the one
#'   that was asked for is the one outcome worth refusing.
#' @return A base `data.frame`: the Gaussian measurement model, one row per
#'   profile and continuous indicator.
#' @seealso [get_results()] for every other table this object holds.
#' @examples
#' fit <- multilpa(
#'   course_engagement,
#'   vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
#'   id = "student", n_profiles = 2, n_group_classes = 2, n_starts = 4,
#'   seed = 1
#' )
#' as.data.frame(fit)
#' @export
as.data.frame.multilpa <- function(x, row.names = NULL, optional = FALSE, ...) {
  stopifnot("`x` must be an object of class `multilpa`" =
              inherits(x, "multilpa"))
  .multilpa_coerce(x, row.names, list(...))
}

#' Every observation with the class it was assigned to
#'
#' The join the caller would otherwise write by hand. Comparing an assignment
#' with anything else -- a known label, an outcome, a covariate -- means putting
#' them in the same row, and doing that with two separate objects silently
#' assumes they share an order. Here the alignment is the verb's, not the
#' reader's, and every column of `data` that the fit can recognise is checked
#' against what the fit holds for that row.
#'
#' Documented on `?get_results`, which is where a caller reaches this table from,
#' including the recovery cross-tabulation `truth` returns and the rule that
#' pairs each truth column with the level it describes.
#'
#' @param x A fitted model of this package.
#' @param data A data frame with one row per observation of the fit, in the
#'   order it was fitted in, or `NULL` for the columns the fit carries.
#' @param truth Column names of `data` holding known labels, or `NULL`.
#' @return A base `data.frame`: one row per observation, or, with `truth`, one
#'   row per class and truth value.
#' @noRd
.multilpa_assignments <- function(x, data = NULL, truth = NULL) {
  stopifnot("`x` must be a fitted model of this package" = .multilpa_any_fit(x))
  frame <- .multilpa_resolve_data(x, data)
  .multilpa_check_row_count(x, frame)
  .multilpa_check_alignment(x, frame)
  posteriors <- as.data.frame(unname(x$subject_posteriors))
  names(posteriors) <- sprintf("posterior_profile_%d", seq_len(x$n_profiles))
  # How much of this row's membership the modal profile did not account for.
  # Modal assignment throws it away, so the table that reports the assignment
  # reports what reporting it cost.
  uncertainty <- 1 - apply(x$subject_posteriors, 1L, max)
  labels <- data.frame(profile = x$subject_profiles)
  if (!is.null(x$group_classes)) {
    labels$group_class <- x$group_classes[x$group_index]
  }
  added <- cbind(labels, uncertainty = uncertainty, posteriors)
  clashes <- intersect(names(frame), names(added))
  if (length(clashes) > 0L) {
    stop(errorCondition(sprintf(
      "`data` already has a column named %s, which the assignments would overwrite.",
      paste(sprintf("`%s`", clashes), collapse = ", ")),
      class = "multilpa_bad_data", call = NULL))
  }
  result <- cbind(frame, added)
  row.names(result) <- NULL
  if (is.null(truth)) return(result)
  .multilpa_recovery_frame(result, truth, x)
}

#' Cross-tabulate the model's labels against known ones
#'
#' The recovery check, as a table rather than as a ritual the reader assembles
#' from an assignment frame and `xtabs()`.
#'
#' @param joined The per-observation frame, already aligned and joined.
#' @param truth Column names of `joined` holding known labels.
#' @param x The fit, for its group index and whether it has group classes.
#' @return One row per class and truth value.
#' @noRd
.multilpa_recovery_frame <- function(joined, truth, x) {
  stopifnot(
    "`truth` must be a character vector of column names" =
      is.character(truth) && length(truth) > 0L && !anyNA(truth),
    "`truth` must not be duplicated" = !anyDuplicated(truth))
  absent <- setdiff(truth, names(joined))
  if (length(absent) > 0L) {
    stop(errorCondition(sprintf(
      "`truth` names %s, which `data` has not.",
      paste(sprintf("`%s`", absent), collapse = ", ")),
      class = "multilpa_bad_data", call = NULL))
  }
  assignment_columns <- intersect(c("profile", "group_class"), truth)
  if (length(assignment_columns) > 0L) {
    stop(errorCondition(sprintf(
      "`truth` names %s, which is the model's own label, not a known one.",
      paste(sprintf("`%s`", assignment_columns), collapse = ", ")),
      class = "multilpa_bad_data", call = NULL))
  }
  rows <- lapply(truth, .multilpa_recovery_rows, joined = joined, x = x)
  result <- do.call(rbind, rows)
  row.names(result) <- NULL
  result
}

#' One truth column cross-tabulated against the level it describes
#' @param name The truth column's name.
#' @param joined The per-observation frame.
#' @param x The fit.
#' @return One row per class and truth value.
#' @noRd
.multilpa_recovery_rows <- function(name, joined, x) {
  values <- joined[[name]]
  # Which level a truth column describes is read off the data, not guessed: a
  # column taking one value within every group is a property of the group, and
  # one that varies inside any group cannot be.
  #
  # That reading needs groups with something inside them. Where every group
  # holds one observation -- an `id = NULL` fit, or a two-level fit with no
  # repeated units -- "constant within group" is true of every column by
  # construction, and the test would send every truth column to the group level
  # and cross-tabulate it against a group class that is the same for all of
  # them. There is no group-level question to ask of singleton groups, so the
  # answer is the profile.
  grouped <- !is.null(x$group_classes) && max(x$group_sizes) > 1L &&
    .multilpa_constant_within(values, x$group_index)
  assignment <- if (grouped) "group_class" else "profile"
  # A group-level truth labels groups, not the observations inside them.
  # Counting every row would give larger groups more weight in recovery.
  units <- if (grouped) !duplicated(x$group_index) else
    rep(TRUE, length(values))
  values <- values[units]
  assigned <- joined[[assignment]][units]
  # `useNA = "ifany"` so a missing label is a visible value rather than a row
  # that quietly leaves the table.
  tabulation <- table(class = assigned, value = values, useNA = "ifany")
  # Factors may retain unused levels. They have no units and hence no
  # defined within-truth proportion; only observed truth values are rows.
  tabulation <- tabulation[, colSums(tabulation) > 0L, drop = FALSE]
  counts <- as.data.frame(tabulation, responseName = "n",
                          stringsAsFactors = FALSE)
  data.frame(
    assignment = assignment,
    class = as.integer(counts$class),
    truth = name,
    value = counts$value,
    n = counts$n,
    # Within a truth value, so the column reads as the share of the units
    # carrying that label which the model put in each class.
    proportion = as.vector(sweep(tabulation, 2L, colSums(tabulation), "/")),
    row.names = NULL, stringsAsFactors = FALSE)
}

#' Does a column take one value within every group?
#' @param values The column.
#' @param group_index One group index per observation.
#' @return A single logical.
#' @noRd
.multilpa_constant_within <- function(values, group_index) {
  distinct <- tapply(as.character(values), group_index,
                     function(within) length(unique(within)))
  all(distinct == 1L)
}

#' The columns of the original data the model was fitted to
#'
#' Rebuilt from what the fit already carries rather than stored a second time:
#' the continuous indicators are kept in `indicator_data` as the model saw
#' them, with centring offsets and input storage types recorded. Categorical
#' indicators are codes in `categorical_data` with one original typed value
#' per code. The identifier is `group_values[group_index]` and the occasion is
#' `time_values`. Continuous values restored after centring agree to floating-
#' point precision; other columns retain their input types. An inference verb
#' therefore no longer has to ask for data the fit already holds.
#'
#' Columns a model never saw -- an outcome, a covariate -- are not here, which
#' is why [three_step()] and [r3step()] still take `data`.
#'
#' @param x A fitted model of this package.
#' @return A base `data.frame`, one row per observation in input order, with the
#'   identifier, the occasion where there is one, and the indicators under their
#'   original names and in their original order.
#' @noRd
.multilpa_model_frame <- function(x) {
  if (is.null(x$indicator_data) && is.null(x$categorical_data)) {
    stop(errorCondition(
      "This fit carries no indicator data, so the data it was fitted to cannot be rebuilt.",
      class = "multilpa_incomplete_fit", call = NULL))
  }
  columns <- list()
  columns[[x$id]] <- x$group_values[x$group_index]
  if (!is.null(x$time) && !is.null(x$time_values)) columns[[x$time]] <- x$time_values
  raw <- .multilpa_uncentered_indicators(x)
  continuous <- if (is.null(raw)) list() else
    stats::setNames(lapply(colnames(raw), function(name) {
      value <- raw[, name]
      # Combining integer and double indicators in a matrix promotes the
      # integers to doubles. Put the supplied type back in the data table;
      # rounding first removes centring's floating-point residue.
      if (identical(x$continuous_types[[name]], "integer")) {
        as.integer(round(value))
      } else value
    }),
                    colnames(raw))
  categorical <- if (is.null(x$categorical_data)) list() else
    stats::setNames(lapply(colnames(x$categorical_data), function(name) {
      original_values <- x$categorical_values[[name]]
      if (!is.null(original_values)) {
        return(original_values[x$categorical_data[, name]])
      }
      # Fits saved before typed category values were retained still expose
      # their labels, even though their original storage class is unknown.
      levels_for <- x$categorical_levels[[name]]
      if (is.null(levels_for)) x$categorical_data[, name] else
        levels_for[x$categorical_data[, name]]
    }), colnames(x$categorical_data))
  indicators <- c(continuous, categorical)
  # Back in the order the caller named them, not the order the model stored them.
  ordered <- x$vars[x$vars %in% names(indicators)]
  frame <- c(columns, indicators[ordered])
  as.data.frame(frame, stringsAsFactors = FALSE, check.names = FALSE)
}

#' The data a verb should use, given what the caller supplied
#'
#' `NULL` means the fit's own columns. Anything else is used as given, so the
#' verbs that need a column the model never saw keep working unchanged.
#'
#' @param x A fitted model of this package.
#' @param data The caller's `data`, possibly `NULL`.
#' @return A base `data.frame`.
#' @noRd
.multilpa_resolve_data <- function(x, data) {
  if (is.null(data)) return(.multilpa_model_frame(x))
  stopifnot("`data` must be a data frame" = is.data.frame(data))
  data
}

#' Do two vectors hold the same values in the same positions
#'
#' Numbers are compared to tolerance rather than bit for bit, so a frame that
#' has been written out and read back still aligns; everything else is compared
#' as text, which is what makes an identifier stored as a factor and supplied as
#' a character vector the same identifier.
#'
#' @param supplied The caller's column.
#' @param stored What the fit holds for those rows.
#' @return `TRUE` when the two agree position by position, `FALSE` otherwise.
#' @noRd
.multilpa_same_values <- function(supplied, stored) {
  if (length(supplied) != length(stored)) return(FALSE)
  absent <- unname(is.na(stored))
  if (!identical(unname(is.na(supplied)), absent)) return(FALSE)
  observed <- !absent
  if (is.numeric(stored) && is.numeric(supplied)) {
    return(isTRUE(all.equal(as.numeric(supplied[observed]),
                            as.numeric(stored[observed]))))
  }
  identical(as.character(supplied[observed]), as.character(stored[observed]))
}

#' Refuse data with a different number of rows from the fit
#'
#' Every verb that reads `data` alongside a fit's posteriors pairs them row by
#' row, so a different row count is the first way data can fail to reproduce
#' the fit, and it is reported with the same class as the others.
#'
#' @param x A fitted model of this package.
#' @param data The frame supplied with it.
#' @return `NULL`, invisibly; raises `multilpa_bad_inference_data` otherwise.
#' @noRd
.multilpa_check_row_count <- function(x, data) {
  if (!identical(as.integer(NROW(data)), as.integer(x$n_observations))) {
    stop(errorCondition(sprintf(
      "`data` must have one row per observation of the fit: %d rows supplied, %d expected.",
      NROW(data), x$n_observations),
      class = "multilpa_bad_inference_data", call = NULL))
  }
  invisible(NULL)
}

#' Check that a frame really is in the fit's own row order
#'
#' Row count alone does not establish alignment: a reversed, resorted or
#' otherwise permuted frame has exactly the right number of rows and pairs every
#' assignment with the wrong observation. The fit already holds, per row, the
#' group index, the indicators and the occasion, so each of those columns the
#' caller supplies can be compared against the fit rather than trusted.
#'
#' A frame with no shared column, or with only repeated identifiers that cannot
#' distinguish rows within a group, does not establish its order. That warns
#' rather than passing silently.
#'
#' @param x A fitted model of this package.
#' @param data The frame to check, already known to have the right row count.
#' @return The names of the columns that were checked, invisibly.
#' @noRd
.multilpa_check_alignment <- function(x, data) {
  checked <- character()
  disagreeing <- character()
  record <- function(name, agrees) {
    checked <<- c(checked, name)
    if (!isTRUE(agrees)) disagreeing <<- c(disagreeing, name)
  }
  if (!is.null(x$id) && x$id %in% names(data) &&
      !is.null(x$group_values) && !is.null(x$group_index)) {
    supplied_index <- match(data[[x$id]], x$group_values)
    record(x$id, identical(as.integer(supplied_index), as.integer(x$group_index)))
  }
  if (!is.null(x$time) && x$time %in% names(data) && !is.null(x$time_values)) {
    record(x$time, .multilpa_same_values(data[[x$time]], x$time_values))
  }
  raw <- .multilpa_uncentered_indicators(x)
  continuous <- intersect(colnames(raw), names(data))
  lapply(continuous, function(name)
    record(name, .multilpa_same_values(data[[name]], raw[, name])))
  categorical <- intersect(colnames(x$categorical_data), names(data))
  lapply(categorical, function(name) {
    levels_for <- x$categorical_levels[[name]]
    codes <- if (is.null(levels_for)) data[[name]] else
      match(as.character(data[[name]]), levels_for)
    record(name, .multilpa_same_values(codes, x$categorical_data[, name]))
  })
  if (length(disagreeing) > 0L) {
    stop(errorCondition(sprintf(
      "`data` is not in the order the model was fitted in: %s disagree%s with the fit row by row.",
      paste(sprintf("`%s`", disagreeing), collapse = ", "),
      if (length(disagreeing) == 1L) "s" else ""),
      class = "multilpa_bad_inference_data", call = NULL))
  }
  # A fit that carries no row-level record of its own has nothing to compare
  # against. Otherwise the supplied columns establish order if their observed
  # combinations uniquely identify rows, or if they contain the complete
  # fitted input (indistinguishable duplicate input rows then share a fitted
  # posterior). A repeated group identifier alone cannot distinguish a swap
  # of two different observations within that group.
  checkable <- c(if (!is.null(x$group_values) && !is.null(x$group_index)) x$id,
                 if (!is.null(x$time_values)) x$time,
                 colnames(x$indicator_data), colnames(x$categorical_data))
  measurements <- c(colnames(x$indicator_data), colnames(x$categorical_data))
  unique_rows <- length(checked) > 0L && !anyDuplicated(data[checked])
  complete_inputs <- length(measurements) > 0L &&
    all(checkable %in% checked)
  if (length(checkable) > 0L && !unique_rows && !complete_inputs) {
    warning(warningCondition(sprintf(
      "`data` does not establish its row order uniquely against the fit; include more fitted columns (for example %s) to have it verified.",
      paste(sprintf("`%s`", checkable), collapse = ", ")),
      class = "multilpa_unverified_alignment"))
  }
  invisible(checked)
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

#' Coerce a class-enumeration grid to its primary table
#'
#' Plain coercion, as the base generic means it: one object, one data frame.
#' The other tables are named rather than positional, so they belong to
#' [get_results()], which takes `what` and refuses a name this object has not.
#'
#' @param x An object of class `multilpa_enumeration`.
#' @param row.names Passed to `data.frame()`; `NULL` gives default row names.
#' @param optional Ignored, present for generic compatibility.
#' @param ... Must be empty. An argument here raises `multilpa_bad_argument`
#'   naming it, rather than being dropped, because `what =` used to live on
#'   this generic and silently returning the primary table instead of the one
#'   that was asked for is the one outcome worth refusing.
#' @return A base `data.frame`: one row per candidate model in the grid.
#' @seealso [get_results()] for every other table this object holds.
#' @examples
#' candidates <- enumerate_classes(
#'   course_engagement,
#'   vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
#'   id = "student", n_profiles = 1:2, n_group_classes = 1, n_starts = 2,
#'   seed = 1
#' )
#' as.data.frame(candidates)
#' @export
as.data.frame.multilpa_enumeration <- function(x, row.names = NULL, optional = FALSE, ...) {
  stopifnot("`x` must be an object of class `multilpa_enumeration`" = inherits(x, "multilpa_enumeration"))
  .multilpa_coerce(x, row.names, list(...))
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
#' @return The input, invisibly. Called for the side effect of printing the
#'   candidate grid: one line per candidate with its class counts, covariance
#'   model, log likelihood, parameter count, AIC, BIC under both sample-size
#'   conventions, ICL counted over individuals, entropy at both levels, and the
#'   diagnostics needed to trust a candidate (convergence, a bound reached, and
#'   how many starts reached the best likelihood), then a count of the
#'   candidates that did not converge. [summary()] and [as.data.frame()] give
#'   every criterion.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' candidates <- enumerate_classes(
#'   example_data, c("score_a", "score_b"), "school", n_profiles = 1:2,
#'   n_group_classes = 1, n_starts = 2, seed = 1
#' )
#' print(candidates)
#' @export
print.multilpa_enumeration <- function(x, ...) {
  stopifnot("`x` must be an `multilpa_enumeration` result" =
              inherits(x, "multilpa_enumeration"))
  table <- x$table
  cat(sprintf("Class enumeration: %d candidate models\n", nrow(table)))
  columns <- intersect(c("n_profiles", "n_group_classes", "model",
                         "log_likelihood", "n_parameters", "aic", "bic_groups",
                         "bic_individual", "icl_individual", "profile_entropy",
                         "group_entropy", "converged", "boundary",
                         "n_best_replicated"),
                       names(table))
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
#'   `profile_probabilities` and `group_probabilities`. The Gaussian blocks are
#'   unnamed and indexed by position, matching what [multilpa()] expects of
#'   `start`; read them with [as.data.frame()], which labels the positions,
#'   rather than out of the list. `response_probabilities` is the exception: the
#'   list keeps each indicator's name and each block keeps its category names,
#'   because those labels are what lets [multilpa()] attach a block to the item
#'   it was estimated for rather than to whichever item happens to sit in that
#'   position. A fit that names its categorical indicators in another order, or
#'   that has other categories, therefore refuses the start with
#'   `multilpa_bad_start` instead of reporting a likelihood for a measurement
#'   model it never held. The list payload is unchanged by the class, so a
#'   `multilpa_start` can be passed straight back as `start`.
#' @seealso [get_results()] for the values as tidy tables.
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
    # The item and category labels are the only thing that tells a later fit
    # which block belongs to which indicator. Stripping them left the blocks to
    # be matched by position, so a fit that names its categorical indicators in
    # another order attached a whole item's distribution to the wrong item and
    # reported the result as an ordinary likelihood. `multilpa(start = )` aligns
    # on these labels and refuses ones it cannot place; the profile row names
    # carry no such meaning and are still dropped.
    start$response_probabilities <- lapply(x$response_probabilities,
                                           function(block) {
      block <- as.matrix(block)
      dimnames(block) <- list(NULL, colnames(block))
      block
    })
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
#' @seealso [get_results()] for the values as tidy tables.
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

#' Coerce a starting-value set to its primary table
#'
#' Plain coercion, as the base generic means it: one object, one data frame.
#' The other tables are named rather than positional, so they belong to
#' [get_results()], which takes `what` and refuses a name this object has not.
#'
#' @param x An object of class `multilpa_start`.
#' @param row.names Passed to `data.frame()`; `NULL` gives default row names.
#' @param optional Ignored, present for generic compatibility.
#' @param ... Must be empty. An argument here raises `multilpa_bad_argument`
#'   naming it, rather than being dropped, because `what =` used to live on
#'   this generic and silently returning the primary table instead of the one
#'   that was asked for is the one outcome worth refusing.
#' @return A base `data.frame`: one row per profile and continuous indicator.
#' @seealso [get_results()] for every other table this object holds.
#' @examples
#' fit <- multilpa(
#'   course_engagement,
#'   vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
#'   id = "student", n_profiles = 2, n_group_classes = 2, n_starts = 4,
#'   seed = 1
#' )
#' as.data.frame(starting_values(fit))
#' @export
as.data.frame.multilpa_start <- function(x, row.names = NULL, optional = FALSE, ...) {
  stopifnot("`x` must be an object of class `multilpa_start`" = inherits(x, "multilpa_start"))
  .multilpa_coerce(x, row.names, list(...))
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
  ## A start now carries the indicator and category labels when it came from
  ## `starting_values()`, and those are what a reader wants to see -- the fitted
  ## object's own response table reports them. An unlabelled start, which is
  ## still a supported input, has only positions to report, so fall back to them
  ## rather than inventing names.
  block_names <- names(blocks)
  rows <- lapply(seq_along(blocks), function(indicator) {
    block <- as.matrix(blocks[[indicator]])
    categories <- colnames(block)
    data.frame(profile = rep(seq_len(nrow(block)), times = ncol(block)),
               indicator = if (is.null(block_names)) indicator else
                 block_names[indicator],
               category = if (is.null(categories))
                 rep(seq_len(ncol(block)), each = nrow(block)) else
                 rep(categories, each = nrow(block)),
               probability = as.vector(block),
               stringsAsFactors = FALSE)
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

#' Coerce a fitted covariate model to its primary table
#'
#' Plain coercion, as the base generic means it: one object, one data frame.
#' The other tables are named rather than positional, so they belong to
#' [get_results()], which takes `what` and refuses a name this object has not.
#'
#' @param x A fitted covariate model.
#' @param row.names Passed to `data.frame()`; `NULL` gives default row names.
#' @param optional Ignored, present for generic compatibility.
#' @param ... Must be empty. An argument here raises `multilpa_bad_argument`
#'   naming it, rather than being dropped, because `what =` used to live on
#'   this generic and silently returning the primary table instead of the one
#'   that was asked for is the one outcome worth refusing.
#' @return A base `data.frame`: the Gaussian measurement model, one row per
#'   profile and continuous indicator.
#' @seealso [get_results()] for every other table this object holds.
#' @examples
#' fit <- multilpa(
#'   course_engagement,
#'   vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
#'   id = "student", n_profiles = 2, n_group_classes = 2,
#'   profile_covariates = "sequence", n_starts = 4, seed = 1
#' )
#' as.data.frame(fit)
#' @export
as.data.frame.multilpa_covariates <- function(x, row.names = NULL, optional = FALSE, ...) {
  stopifnot("`x` must be an object of class `multilpa_covariates`" =
              inherits(x, "multilpa_covariates"))
  .multilpa_coerce(x, row.names, list(...))
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
