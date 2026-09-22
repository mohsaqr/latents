#' Describe the variables a model uses, or will use
#'
#' One row per variable, with the summaries you would look at before fitting and
#' the one that decides whether a two-level model has anything to find: the
#' intraclass correlation, the share of each variable's variance that lies
#' *between* groups rather than within them.
#'
#' An ICC near zero suggests little between-group mean variation in that
#' variable. It does not rule out group classes that differ in other variables,
#' profile prevalence, or other features of their distributions.
#'
#' @param x A data frame, or a fitted model of this package. A fit carries the
#'   columns it was built from, so `vars` and `id` are read from it.
#' @param vars Character vector of column names to describe. Repeated names
#'   are described once. For a data frame, defaults to every numeric column
#'   that is not `id`.
#' @param id Optional. The column identifying the groups rows are nested in.
#'   Supplying it adds `n_groups` and `icc`; omitting it leaves them out.
#' @param by Optional. `"profile"` or `"group_class"` on a fitted model, or a
#'   column name on a data frame: describe each variable once per level of it,
#'   rather than once overall. Rows whose stratifier is `NA` are neither dropped
#'   nor folded into another level: they form their own stratum, labelled `NA`
#'   in the `by` column and reported last, so that the strata account for every
#'   input row exactly once. The stratifier is never compared with `==`, so a
#'   missing value cannot select rows it does not belong to.
#' @param ... Passed between methods.
#' @return A base `data.frame`, one row per variable, or one row per variable
#'   and `by` level, with the columns
#'   \describe{
#'     \item{`variable`}{character: the column described.}
#'     \item{`n`}{integer: observations with a value.}
#'     \item{`n_missing`}{integer: observations without one.}
#'     \item{`mean`,`sd`,`min`,`max`}{numeric, and `NA_real_` for a variable
#'       that is not numeric.}
#'     \item{`n_distinct`}{integer: distinct observed values, which is the
#'       useful summary where a mean is not.}
#'     \item{`n_groups`,`icc`}{present only when `id` is known. `icc` is the
#'       one-way random-effects estimate, and may be slightly negative when the
#'       between-group mean square falls below the within-group one; that is the
#'       estimator reporting no group structure, not an error. Rows with a
#'       missing group ID do not contribute to either statistic.}
#'   }
#'   A `by` column comes first, named after what it splits on and holding that
#'   level's own value, `NA` for the stratum of rows whose stratifier is
#'   missing. For every variable, `n + n_missing` summed over the strata equals
#'   the number of input rows; that invariant is asserted before the table is
#'   returned. Describing by a class the model assigned keeps that assignment
#'   out of the described frame, so an indicator of the caller's own that is
#'   named `profile` or `group_class` is summarized as itself, not replaced by
#'   the assignment. A `by` whose name is one of the summary's own columns
#'   cannot be told apart from them and raises `multilpa_bad_data`.
#' @references Bliese, P. D. (2000). Within-group agreement, non-independence,
#'   and reliability. In K. J. Klein & S. W. J. Kozlowski (Eds.), *Multilevel
#'   theory, research, and methods in organizations*.
#' @seealso [multilpa()] to fit, [diagnostics()] for the after-the-fit
#'   counterpart.
#' @examples
#' descriptives(
#'   course_engagement,
#'   vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
#'   id = "student"
#' )
#'
#' # A missing stratifier is a stratum of its own, not a silent loss and not
#' # missingness invented in a variable that has none.
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

#' @rdname descriptives
#' @export
descriptives.data.frame <- function(x, vars = NULL, id = NULL, by = NULL, ...) {
  stopifnot("`vars` must be column names of `x`" =
              is.null(vars) || (is.character(vars) && all(vars %in% names(x))),
            "`id` must be a single column name of `x`" =
              is.null(id) || (is.character(id) && length(id) == 1L && id %in% names(x)),
            "`by` must be a single column name of `x`" =
              is.null(by) || (is.character(by) && length(by) == 1L && by %in% names(x)))
  if (is.null(vars)) {
    numeric_columns <- names(x)[vapply(x, is.numeric, logical(1))]
    vars <- setdiff(numeric_columns, c(id, by))
  }
  vars <- unique(vars)
  if (length(vars) == 0L) {
    stop(errorCondition("There are no variables to describe.",
                        class = "multilpa_nothing_to_describe", call = NULL))
  }
  if (is.null(by)) return(.multilpa_describe(x, vars, id))
  .multilpa_describe_by(x, vars, id, x[[by]], by)
}

#' @rdname descriptives
#' @export
descriptives.multilpa <- function(x, by = NULL, ...) {
  stopifnot("`by` must be \"profile\" or \"group_class\"" =
              is.null(by) || (is.character(by) && length(by) == 1L &&
                                by %in% c("profile", "group_class")))
  frame <- .multilpa_model_frame(x)
  if (is.null(by)) {
    return(descriptives.data.frame(frame, vars = x$vars, id = x$id, ...))
  }
  # The assignment is metadata, not a measurement: it is carried beside the
  # frame rather than written into it, so an indicator the caller named
  # `profile` or `group_class` keeps its own values.
  .multilpa_describe_by(frame, x$vars, x$id,
                        .multilpa_assignment_values(x, by), by)
}

#' @rdname descriptives
#' @export
descriptives.multilpa_transitions <- descriptives.multilpa

#' @rdname descriptives
#' @export
descriptives.multilpa_covariates <- descriptives.multilpa

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
      class = "multilpa_incomplete_fit", call = NULL))
  }
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
                       "n_distinct", "n_groups", "icc")
  stopifnot(
    "`vars` must be column names of the frame being described" =
      all(vars %in% names(data)),
    "`by` must have one value for each row described" =
      length(by_values) == nrow(data))
  if (by_label %in% summary_columns) {
    stop(errorCondition(sprintf(
      "`by` names the column `%s`, which the summary itself uses, so the two could not be told apart in the result.",
      by_label), class = "multilpa_bad_data", call = NULL))
  }
  observed_levels <- unique(by_values[!is.na(by_values)])
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
    row$icc <- if (is.numeric(values)) .multilpa_icc(values[grouped], groups) else NA_real_
    row
  })
  result <- do.call(rbind, described)
  row.names(result) <- NULL
  result
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
.multilpa_icc <- function(values, groups) {
  sizes <- tabulate(match(groups, unique(groups)))
  k <- length(sizes)
  n <- length(values)
  if (k < 2L || n <= k) return(NA_real_)
  grand_mean <- mean(values)
  group_means <- vapply(split(values, groups), mean, numeric(1))
  between <- sum(sizes * (group_means[match(unique(groups), names(group_means))] -
                            grand_mean)^2) / (k - 1L)
  within <- sum((values - group_means[match(groups, names(group_means))])^2) / (n - k)
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
