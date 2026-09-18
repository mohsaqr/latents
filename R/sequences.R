#' The time column of a fit, or a classed error
#'
#' @param object A fitted `multilpa` model.
#' @return The stored time values, in input row order.
#' @noRd
.multilpa_require_time <- function(object) {
  stopifnot("`object` must be a fitted model of this package" =
              .multilpa_any_fit(object))
  if (is.null(object$time_values)) {
    stop(errorCondition(
      paste("This fit carries no ordering. Refit with `time =` naming the column",
            "that gives each observation's position within its group."),
      class = "multilpa_no_time", call = NULL))
  }
  object$time_values
}

#' Assignments in sequence order
#'
#' Returns which profile the model gave each observation, laid out in the order
#' the observations occur within their group, together with the group's own
#' latent class. The model itself ignores the ordering, so this is the view that
#' shows whether a group's assignments drift, alternate or hold steady --
#' something the profile prevalences alone cannot say.
#'
#' @param object A fitted `multilpa` model, fitted with `time =`.
#' @param format `"long"` (the default) gives one row per observation.
#'   `"wide"` gives one row per group and one column per position, which is the
#'   shape sequence-plotting functions expect.
#' @return For `format = "long"`, a base `data.frame` with one row per
#'   observation and the columns `group`, `group_class`, `time` and `profile`,
#'   ordered by group and then by time. For `format = "wide"`, a base
#'   `data.frame` with one row per group and one column per distinct position,
#'   holding the profile as a factor and `NA` where the group has no observation
#'   at that position; row names are the group identifiers and the columns are
#'   in increasing time order.
#' @details Profile and group-class labels are arbitrary and are returned as
#'   integers; two fits must have their labels aligned before their sequences
#'   are compared.
#' @seealso [sequence_summary()] for one row per group class, and
#'   `plot(object, what = "sequences")` to draw them.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10), wave = rep(seq_len(10), times = 12),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                 n_profiles = 2, n_group_classes = 2, n_starts = 2,
#'                 seed = 1, time = "wave")
#' head(sequences(fit))
#' head(sequences(fit, format = "wide"))
#' @export
sequences <- function(object, format = c("long", "wide")) {
  time_values <- .multilpa_require_time(object)
  format <- match.arg(format)
  long <- data.frame(
    group = object$group_values[object$group_index],
    group_class = object$group_classes[object$group_index],
    time = time_values,
    profile = object$subject_profiles,
    row.names = NULL, stringsAsFactors = FALSE
  )
  long <- long[order(long$group, long$time), , drop = FALSE]
  row.names(long) <- NULL
  if (identical(format, "long")) return(long)

  positions <- sort(unique(long$time))
  groups <- unique(long$group)
  cells <- matrix(NA_integer_, length(groups), length(positions),
                  dimnames = list(as.character(groups), as.character(positions)))
  cells[cbind(match(long$group, groups), match(long$time, positions))] <- long$profile
  wide <- as.data.frame(lapply(seq_along(positions), function(column) {
    factor(cells[, column], levels = seq_len(object$n_profiles))
  }), col.names = as.character(positions), optional = TRUE)
  names(wide) <- as.character(positions)
  row.names(wide) <- as.character(groups)
  wide
}

#' Sequence lengths by group class
#'
#' Summarizes how much data each latent group class actually contributes. A
#' class whose groups are systematically shorter is a different finding from a
#' class whose groups simply differ in profile mix, and the prevalence table
#' cannot show it.
#'
#' @param object A fitted `multilpa` model, fitted with `time =`.
#' @return A base `data.frame` with one row per group class and the columns
#'   `group_class`, `groups` (how many groups are assigned to it),
#'   `observations` (their total number of observations), `mean_length`,
#'   `median_length`, `shortest` and `complete` (how many of its groups are
#'   observed at every position seen anywhere in the data).
#' @seealso [sequences()].
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10), wave = rep(seq_len(10), times = 12),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                 n_profiles = 2, n_group_classes = 2, n_starts = 2,
#'                 seed = 1, time = "wave")
#' sequence_summary(fit)
#' @export
sequence_summary <- function(object) {
  long <- sequences(object, format = "long")
  positions <- length(unique(long$time))
  lengths_by_group <- tapply(long$time, long$group, length)
  class_by_group <- tapply(long$group_class, long$group, function(v) v[[1L]])
  classes <- seq_len(object$n_group_classes)
  summarize <- function(class) {
    taken <- lengths_by_group[class_by_group == class]
    data.frame(
      group_class = class,
      groups = length(taken),
      observations = sum(taken),
      mean_length = if (length(taken)) mean(taken) else NA_real_,
      median_length = if (length(taken)) stats::median(taken) else NA_real_,
      shortest = if (length(taken)) min(taken) else NA_integer_,
      complete = sum(taken == positions),
      row.names = NULL
    )
  }
  do.call(rbind, lapply(classes, summarize))
}
