#' Describe the variables a model uses, or will use
#'
#' One row per variable, with the summaries you would look at before fitting and
#' the one that decides whether a two-level model has anything to find: the
#' intraclass correlation, the share of each variable's variance that lies
#' *between* groups rather than within them.
#'
#' An ICC near zero says the groups do not differ on that variable, so a model
#' whose whole purpose is to tell groups apart will be fitting noise. It is the
#' cheapest check available and the one most often skipped.
#'
#' @param x A data frame, or a fitted model of this package. A fit carries the
#'   columns it was built from, so `vars` and `id` are read from it.
#' @param vars Character vector of column names to describe. For a data frame,
#'   defaults to every numeric column that is not `id`.
#' @param id Optional. The column identifying the groups rows are nested in.
#'   Supplying it adds `n_groups` and `icc`; omitting it leaves them out.
#' @param by Optional. `"profile"` or `"group_class"` on a fitted model, or a
#'   column name on a data frame: describe each variable once per level of it,
#'   rather than once overall.
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
#'       estimator reporting no group structure, not an error.}
#'   }
#'   A `by` column comes first, named after what it splits on.
#' @references Bliese, P. D. (2000). Within-group agreement, non-independence,
#'   and reliability. In K. J. Klein & S. W. J. Kozlowski (Eds.), *Multilevel
#'   theory, research, and methods in organizations*.
#' @seealso [multilpa()] to fit, [diagnostics()] for the after-the-fit
#'   counterpart.
#' @examples
#' descriptives(school_engagement,
#'              vars = c("homework_hours", "participation", "interest"),
#'              id = "school")
#'
#' fit <- multilpa(school_engagement,
#'                 c("homework_hours", "participation", "interest"),
#'                 id = "school", n_profiles = 2, n_group_classes = 2,
#'                 n_starts = 4, seed = 1)
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
  if (length(vars) == 0L) {
    stop(errorCondition("There are no variables to describe.",
                        class = "multilpa_nothing_to_describe", call = NULL))
  }
  if (is.null(by)) return(.multilpa_describe(x, vars, id))
  levels_of <- unique(x[[by]])
  described <- lapply(levels_of, function(level) {
    rows <- x[x[[by]] == level, , drop = FALSE]
    cbind(stats::setNames(data.frame(level, stringsAsFactors = FALSE), by),
          .multilpa_describe(rows, vars, id))
  })
  result <- do.call(rbind, described)
  row.names(result) <- NULL
  result
}

#' @rdname descriptives
#' @export
descriptives.multilpa <- function(x, by = NULL, ...) {
  stopifnot("`by` must be \"profile\" or \"group_class\"" =
              is.null(by) || (is.character(by) && length(by) == 1L &&
                                by %in% c("profile", "group_class")))
  frame <- .multilpa_model_frame(x)
  if (!is.null(by)) {
    frame[[by]] <- if (identical(by, "profile")) x$subject_profiles else
      x$group_classes[x$group_index]
  }
  descriptives.data.frame(frame, vars = x$vars, id = x$id, by = by, ...)
}

#' @rdname descriptives
#' @export
descriptives.multilpa_transitions <- descriptives.multilpa

#' @rdname descriptives
#' @export
descriptives.multilpa_covariates <- descriptives.multilpa

#' @rdname descriptives
#' @export
descriptives.multilpa_random_intercept <- function(x, by = NULL, ...) {
  stopifnot("`by` must be \"profile\"" =
              is.null(by) || identical(by, "profile"))
  frame <- .multilpa_model_frame(x)
  if (!is.null(by)) frame[[by]] <- x$subject_profiles
  descriptives.data.frame(frame, vars = x$vars, id = x$id, by = by, ...)
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
    numeric_values <- if (is.numeric(values)) values else NA_real_
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
    groups <- data[[id]][!is.na(values)]
    row$n_groups <- length(unique(groups))
    row$icc <- if (is.numeric(values)) .multilpa_icc(observed, groups) else NA_real_
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
