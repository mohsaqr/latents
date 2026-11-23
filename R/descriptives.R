#' Describe each variable of a data frame or a fitted model, with its variation
#' within and between groups
#'
#' @description
#' `descriptives()` describes each variable of a data frame, or each indicator
#' of a fitted model: the number of observed and missing values, the mean,
#' standard deviation, minimum and maximum of a numeric variable, and the number
#' of distinct values. When the rows are nested in groups it adds the number of
#' groups, the standard deviations within and between groups, and the intraclass
#' correlation, the share of a variable's variance that lies between groups.
#'
#' The table can be computed for the whole data set, or once per level of an
#' observed grouping column (`by`). For a fitted model it can be computed
#' once per assigned profile or group class. Use it to inspect variability,
#' ranges and missingness before fitting, and to describe the observations
#' assigned to each class afterwards.
#'
#' @details
#' Summaries by fitted class use assigned membership and do not incorporate
#' classification uncertainty. They describe the observations assigned to each
#' class and should be distinguished from the model's estimated class-specific
#' parameters. Rows excluded from classification are retained in a separate
#' stratum with a missing class label, following the fit's `drop_empty` setting.
#' These are unweighted descriptions of the recorded rows, including rows with
#' sampling weight zero. Describing a weighted fit warns explicitly; the ICC
#' is the unweighted one-way ANOVA quantity, not a survey-weighted estimate.
#'
#' An intraclass correlation (ICC) near zero indicates little between-group
#' variation in that variable's mean. It does not rule out group classes that
#' differ in other variables, in profile prevalence, or in other features of
#' the distributions.
#'
#' @param x What to describe:
#'   * a data frame, with `vars`, `id` and `by` naming its columns;
#'   * a profile or class fit from [lpa()], [lca()], [multilpa()] or
#'     [multilca()], with or without membership covariates;
#'   * an [lta()] fit with its default homogeneous first-order transitions,
#'     invariant measurement and continuous or categorical indicators.
#'
#'   A fit is described on its own indicators and group identifier, so
#'   `vars` and `id` are not given with it. A single-level fit has no
#'   groups, so its table has no `n_groups`, `sd_within`, `sd_between` or
#'   `icc`. Anything else raises
#'   `latents_not_a_fit`, including the other `multilpa()` families
#'   (additive, dispersion and cross-level), an `lta()` fit with ordinal or
#'   count indicators, covariates, occasion-specific or second-order
#'   transitions, occasion-specific measurement, a stayer class or weights,
#'   and a fit from [mixture_regression()]; describe the data frame it was
#'   fitted to instead.
#' @param vars For a data frame: character vector of column names to
#'   describe; repeated names are described once. `NULL` (default) takes
#'   every numeric column other than `id` and `by`. Columns must be vectors:
#'   matrix-valued columns raise `latents_bad_data`.
#' @param id For a data frame: optional name of the column identifying the
#'   groups the rows are nested in (for example the person in repeated
#'   measures). Supplying it adds `n_groups`, `sd_within`, `sd_between` and
#'   `icc`.
#' @param by Optional stratifier; each variable is then described once per
#'   level of it. `NULL` (default) describes all rows together.
#'   * On a data frame: the name of one of its columns.
#'   * On a fitted model, `"profile"`: the profile (class) each observation
#'     is assigned to.
#'   * On a fitted model, `"group_class"`: the group class the observation's
#'     group is assigned to; a single-level fit has one group class, so one
#'     stratum.
#'
#'   Rows whose stratifier is `NA` form their own stratum, labelled `NA`
#'   and reported last, so every input row is counted in exactly one
#'   stratum. On a fitted model any other value raises
#'   `latents_bad_argument`.
#' @param ... Not used. For a data frame, any further argument raises a
#'   `latents_unused_argument` warning; for a fit, it raises
#'   `latents_bad_argument`.
#' @return A base `data.frame`, one row per variable (or per variable and `by`
#'   level), with the columns
#'   \describe{
#'     \item{`variable`}{character: the column described.}
#'     \item{`n`}{integer: observations with a value.}
#'     \item{`n_missing`}{integer: observations without one.}
#'     \item{`mean`,`sd`,`min`,`max`}{numeric: the mean, standard deviation,
#'       minimum and maximum of the observed values; `NA_real_` for a
#'       variable that is not numeric, such as a categorical indicator.}
#'     \item{`n_distinct`}{integer: distinct observed values.}
#'     \item{`n_groups`}{integer, present only when `id` is known: the groups
#'       with an observed value of the variable.}
#'     \item{`sd_within`}{numeric, present only when `id` is known: the
#'       pooled within-group standard deviation, how much the values vary
#'       around their own group's mean (the square root of the within-group
#'       mean square). `NA` when no group holds more than one value, and for
#'       a variable that is not numeric.}
#'     \item{`sd_between`}{numeric, present only when `id` is known: the
#'       standard deviation of the group means, each group counted once, how
#'       much the groups differ from each other. `NA` with fewer than two
#'       groups, and for a variable that is not numeric.}
#'     \item{`icc`}{numeric, present only when `id` is known: the one-way
#'       random-effects intraclass correlation,
#'       `(MSB - MSW) / (MSB + (n0 - 1) * MSW)`, from the between-group and
#'       within-group mean squares of a one-way analysis of variance, with
#'       `n0` the average group size corrected for unequal sizes. It can be
#'       negative when the between-group mean square is smaller than the
#'       within-group mean square; this does not represent a negative
#'       population variance. It is `NA` for a variable that is not numeric,
#'       with fewer than two groups, with no group holding more than one
#'       value, or constant over all rows. Rows with a missing group
#'       identifier contribute to neither statistic.}
#'   }
#'   With `by`, a first column named after the stratifier holds each stratum's
#'   value (`NA` for the stratum of missing values). For every variable,
#'   `n + n_missing` summed over the strata equals the number of input rows.
#'   When describing a fit by `"profile"` or `"group_class"`, the model's
#'   assignment is used as the stratifier without replacing an indicator of
#'   the same name.
#'
#'   Errors: `latents_not_a_fit` for an `x` that is neither a data frame
#'   nor a supported fit; `latents_bad_data` for a `vars` name the data
#'   frame does not have, a matrix-valued column, or a `by` named like one
#'   of the summary's own columns; `latents_nothing_to_describe` when there
#'   are no variables to describe; `latents_bad_argument` for any other
#'   invalid argument. Describing a fit with sampling weights warns with
#'   `latents_unweighted_descriptives`.
#' @seealso [multilpa()] for model estimation and [diagnostics()] for
#'   diagnostics of a fitted model.
#' @examples
#' descriptives(
#'   course_engagement,
#'   vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
#'   id = "student"
#' )
#'
#' # Missing stratification values form a separate stratum.
#' descriptives(data.frame(cohort = c("A", "A", "B", NA), score = 1:4),
#'              vars = "score", by = "cohort")
#'
#' fit <- multilpa(
#'   course_engagement,
#'   vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
#'   id = "student", n_profiles = 2, n_group_classes = 2, n_starts = 4,
#'   seed = 1
#' )
#' descriptives(fit)
#' descriptives(fit, by = "profile")
#' @export
descriptives <- function(x, ...) {
  UseMethod("descriptives")
}

#' @noRd
#' @export
descriptives.default <- function(x, ...) {
  # Any other object reached UseMethod()'s unclassed "no applicable method";
  # a fit of a family without indicators to describe is told what to pass.
  .latents_require_fit(x, substitute(x), expected = paste(
    "a data frame, or", .latents_profile_fit_text),
    hint = "Describe the data frame it was fitted to: descriptives(data, vars = ).")
}

#' @rdname descriptives
#' @export
descriptives.data.frame <- function(x, vars = NULL, id = NULL, by = NULL, ...) {
  absent <- if (is.character(vars)) setdiff(vars, names(x)) else character()
  if (length(absent) > 0L) {
    stop(errorCondition(sprintf("`vars` names %s, which `x` does not have.%s",
                                .latents_name_list(absent),
                                .latents_did_you_mean(absent, names(x))),
                        class = "latents_bad_data", call = NULL))
  }
  .latents_require("`vars` must be column names of `x`" =
              is.null(vars) || (is.character(vars) && all(vars %in% names(x))),
            "`id` must be a single column name of `x`" =
              is.null(id) || (is.character(id) && length(id) == 1L && id %in% names(x)),
            "`by` must be a single column name of `x`" =
              is.null(by) || (is.character(by) && length(by) == 1L && by %in% names(x)))
  # Released versions dropped anything else without a word; it still describes
  # as before, but says what it did not use.
  .latents_warn_unused_arguments(list(...), "descriptives()",
                                 "It takes `vars`, `id` and `by`.")
  if (is.null(vars)) {
    numeric_columns <- names(x)[vapply(x, is.numeric, logical(1))]
    vars <- setdiff(numeric_columns, c(id, by))
  }
  vars <- unique(vars)
  if (length(vars) == 0L) {
    stop(errorCondition("There are no variables to describe.",
                        class = "latents_nothing_to_describe", call = NULL))
  }
  .multilpa_check_descriptive_columns(x, unique(c(vars, id, by)))
  if (is.null(by)) return(.multilpa_describe(x, vars, id))
  .multilpa_describe_by(x, vars, id, x[[by]], by)
}

#' @noRd
#' @export
descriptives.multilpa <-function(x, by = NULL, ...) {
  if (!(is.null(by) || (is.character(by) && length(by) == 1L &&
                        by %in% c("profile", "group_class")))) {
    # A column name here is the usual slip: the fit's indicators are already
    # its own, and only the grouping is chosen.
    stop(errorCondition(sprintf(paste(
      "For a fitted model, `by` must be \"profile\" or \"group_class\", not %s.",
      "The fit describes its own indicators; to describe chosen columns, pass the",
      "data frame: descriptives(data, vars = ...)."),
      if (is.character(by)) sprintf("\"%s\"", paste(by, collapse = "\", \"")) else
        .latents_describe(by)),
      class = "latents_bad_argument", call = NULL))
  }
  # `vars` and `id` are the fit's own, so there is nothing further to pass on;
  # a supplied one would otherwise collide with them or be dropped unread.
  .multilpa_reject_extra_arguments(
    list(...), "descriptives()",
    "A fit supplies its own `vars` and `id`; describe a data frame to choose them.")
  frame <- .multilpa_model_frame(x)
  if (.latents_is_weighted(x)) {
    warning(warningCondition(paste(
      "descriptives() summarizes the recorded rows without sampling weights,",
      "including zero-weight rows; its ICC is also unweighted. Describe the",
      "data frame directly when that unweighted summary is intended."),
      class = "latents_unweighted_descriptives", call = NULL))
  }
  # A single-level fit numbers its rows as one-row groups of its own making.
  # Its grouping summaries would describe that numbering, not the data.
  id <- if (isTRUE(x$single_level)) NULL else x$id
  if (is.null(by)) {
    return(descriptives.data.frame(frame, vars = x$vars, id = id))
  }
  # The assignment is metadata, not a measurement: it is carried beside the
  # frame rather than written into it, so an indicator the caller named
  # `profile` or `group_class` keeps its own values.
  .multilpa_describe_by(frame, x$vars, id,
                        .multilpa_assignment_values(x, by), by)
}

#' @noRd
#' @export
descriptives.multilpa_transitions <-descriptives.multilpa

#' @noRd
#' @export
descriptives.multilpa_covariates <-descriptives.multilpa

#' The class a fit assigned to each of its observations
#'
#' Read from the fit rather than written into the frame it describes, so that
#' the assignment can never displace a column of the caller's own.
#'
#' @param x A fitted model of this package.
#' @param by `"profile"` or `"group_class"`.
#' @return A vector with one element per observation of the fit.
#' @noRd
.multilpa_assignment_values <- function(x, by) {
  values <- if (identical(by, "profile")) x$subject_profiles else
    x$group_classes[x$group_index]
  if (is.null(values) || length(values) != x$n_observations) {
    stop(errorCondition(sprintf(
      "This fit carries no `%s` assignment for each of its observations.", by),
      class = "latents_incomplete_fit", call = NULL))
  }
  included <- if (identical(by, "profile")) .multilpa_informative_rows(x) else
    .multilpa_informative_groups(x)[x$group_index]
  values[!included] <- NA_integer_
  values
}

#' Describe each variable once per stratum
#'
#' The stratifier arrives as a vector beside the frame, not as a column of it,
#' so describing by a class the model assigned cannot overwrite an indicator of
#' the same name.
#'
#' Rows are keyed to their stratum with `match()`, which returns `NA` for a
#' missing stratifier. Selecting with `which()` on that key therefore drops the
#' missing rows from every observed stratum instead of turning each selector
#' into a vector of `NA`s, which is what `x[[by]] == level` did: `NA` selectors
#' fabricate rows that are missing on every variable and lose the observed
#' values of the rows that carry them.
#'
#' @param data A data frame.
#' @param vars Column names to describe.
#' @param id Optional grouping column name.
#' @param by_values One stratifier value per row of `data`.
#' @param by_label The name the stratum column takes in the result.
#' @return A base `data.frame`, one row per stratum and variable, the stratum of
#'   rows with a missing stratifier last.
#' @noRd
.multilpa_describe_by <- function(data, vars, id, by_values, by_label) {
  summary_columns <- c("variable", "n", "n_missing", "mean", "sd", "min", "max",
                       "n_distinct", "n_groups", "sd_within", "sd_between",
                       "icc")
  stopifnot(
    "`vars` must be column names of the frame being described" =
      all(vars %in% names(data)),
    "`by` must have one value for each row described" =
      length(by_values) == nrow(data))
  if (by_label %in% summary_columns) {
    stop(errorCondition(sprintf(
      "`by` names the column `%s`, which the summary itself uses, so the two could not be told apart in the result.",
      by_label), class = "latents_bad_data", call = NULL))
  }
  # Strata in a fixed order -- a factor's own level order, otherwise sorted --
  # so the table reads class 1, 2, ... whatever order the rows came in.
  present <- by_values[!is.na(by_values)]
  observed_levels <- if (is.factor(by_values)) {
    factor(intersect(levels(by_values), as.character(present)),
           levels = levels(by_values))
  } else sort(unique(present))
  keys <- match(by_values, observed_levels)
  missing_stratum <- anyNA(keys)
  strata <- c(lapply(seq_along(observed_levels), \(index) which(keys == index)),
              if (missing_stratum) list(which(is.na(keys))))
  # A length-one vector of the stratifier's own type, so a factor keeps its
  # levels and a numeric stratifier is not coerced to character by its label.
  labels <- c(observed_levels, if (missing_stratum) by_values[NA_integer_])
  if (length(strata) == 0L) {
    empty <- cbind(stats::setNames(data.frame(by_values[NA_integer_]), by_label),
                   .multilpa_describe(data, vars, id))
    return(empty[0L, , drop = FALSE])
  }
  described <- lapply(seq_along(strata), function(index) {
    cbind(stats::setNames(data.frame(labels[index], stringsAsFactors = FALSE),
                          by_label),
          .multilpa_describe(data[strata[[index]], , drop = FALSE], vars, id))
  })
  result <- do.call(rbind, described)
  row.names(result) <- NULL
  accounted <- vapply(split(result$n + result$n_missing, result$variable),
                      sum, integer(1))
  stopifnot(
    "Every described row must fall in exactly one stratum" =
      identical(sum(lengths(strata)), nrow(data)),
    "Each variable must be counted once for every described row" =
      all(accounted == nrow(data)))
  result
}

#' Describe each variable of one frame
#' @param data A data frame.
#' @param vars Column names to describe.
#' @param id Optional grouping column name.
#' @return A base `data.frame`, one row per variable.
#' @noRd
.multilpa_describe <- function(data, vars, id) {
  .multilpa_check_descriptive_columns(data, unique(c(vars, id)))
  described <- lapply(vars, function(name) {
    values <- data[[name]]
    observed <- values[!is.na(values)]
    row <- data.frame(
      variable = name,
      n = length(observed),
      n_missing = sum(is.na(values)),
      mean = if (is.numeric(values)) mean(observed) else NA_real_,
      sd = if (is.numeric(values) && length(observed) > 1L) stats::sd(observed) else NA_real_,
      min = if (is.numeric(values) && length(observed) > 0L) min(observed) else NA_real_,
      max = if (is.numeric(values) && length(observed) > 0L) max(observed) else NA_real_,
      n_distinct = length(unique(observed)),
      stringsAsFactors = FALSE)
    if (is.null(id)) return(row)
    grouped <- !is.na(values) & !is.na(data[[id]])
    groups <- data[[id]][grouped]
    row$n_groups <- length(unique(groups))
    spread <- if (is.numeric(values)) .multilpa_level_sds(values[grouped], groups) else
      c(within = NA_real_, between = NA_real_)
    row$sd_within <- spread[["within"]]
    row$sd_between <- spread[["between"]]
    row$icc <- if (is.numeric(values)) .multilpa_icc(values[grouped], groups) else NA_real_
    row
  })
  result <- do.call(rbind, described)
  row.names(result) <- NULL
  result
}

#' Require one descriptive value per data-frame row
#' @param data A data frame.
#' @param columns The variables, identifier and stratifier to be used.
#' @return `NULL`, invisibly; raises `latents_bad_data` for matrix/array columns.
#' @noRd
.multilpa_check_descriptive_columns <- function(data, columns) {
  stopifnot(is.data.frame(data), is.character(columns),
            all(columns %in% names(data)))
  invalid <- vapply(data[columns], function(value) {
    !is.null(dim(value)) || length(value) != nrow(data)
  }, logical(1))
  if (any(invalid)) {
    stop(errorCondition(sprintf(
      "Descriptive columns must be vectors with one value per row: %s.",
      paste(sprintf("`%s`", columns[invalid]), collapse = ", ")),
      class = "latents_bad_data", call = NULL))
  }
  invisible(NULL)
}

#' One-way random-effects intraclass correlation
#'
#' The share of variance lying between groups, from the one-way ANOVA mean
#' squares. `n0` is the size-corrected average group size, which equals the
#' common size when the groups are balanced and is smaller when they are not, so
#' a few large groups cannot inflate the estimate.
#'
#' The estimate is returned as computed. It can fall slightly below zero when
#' the between-group mean square is smaller than the within-group one, which is
#' the estimator saying there is no group structure; truncating it at zero would
#' hide that rather than report it.
#'
#' @param values Numeric observations, already complete.
#' @param groups The group each observation belongs to.
#' @return A single number, or `NA_real_` where the design cannot support one.
#' @noRd
#' Within-group and between-group standard deviations
#'
#' The within-group SD is the square root of the pooled within-group
#' variance, the within-group mean square of a one-way analysis of variance:
#' each value's deviation from its own group's mean, summed over groups and
#' divided by n - k. The between-group SD is the standard deviation of the
#' group means, each group counted once whatever its size.
#' @param values Numeric values, none missing.
#' @param groups The group of each value, none missing.
#' @return A named numeric vector, `within` and `between`; `NA` where
#'   undefined (no group with two values for `within`, fewer than two groups
#'   for `between`).
#' @noRd
.multilpa_level_sds <- function(values, groups) {
  stopifnot("one group per value" = length(values) == length(groups),
            "no missing values or groups" = !anyNA(values) && !anyNA(groups))
  index <- match(groups, unique(groups))
  k <- max(c(0L, index))
  n <- length(values)
  if (k == 0L) return(c(within = NA_real_, between = NA_real_))
  group_means <- vapply(split(values, index), mean, numeric(1))
  within <- if (n > k) sqrt(sum((values - group_means[index])^2) / (n - k)) else NA_real_
  between <- if (k >= 2L) stats::sd(group_means) else NA_real_
  c(within = within, between = between)
}

.multilpa_icc <- function(values, groups) {
  index <- match(groups, unique(groups))
  sizes <- tabulate(index)
  k <- length(sizes)
  n <- length(values)
  if (k < 2L || n <= k) return(NA_real_)
  grand_mean <- mean(values)
  # Integer indices preserve native identifiers that have identical printed
  # labels (for example 0.3 and 0.1 + 0.2).
  group_means <- vapply(split(values, index), mean, numeric(1))
  between <- sum(sizes * (group_means - grand_mean)^2) / (k - 1L)
  within <- sum((values - group_means[index])^2) / (n - k)
  if (!is.finite(between) || !is.finite(within)) return(NA_real_)
  n0 <- (n - sum(sizes^2) / n) / (k - 1L)
  denominator <- between + (n0 - 1) * within
  # Zero within-group variance is not a degenerate case: every observation
  # equals its group's own value, so all of the variance is between groups and
  # the ratio is exactly one. Only a variable that is constant everywhere has
  # no ratio to report.
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  (between - within) / denominator
}
