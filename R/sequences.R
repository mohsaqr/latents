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
      class = "latents_no_time", call = NULL))
  }
  object$time_values
}

#' Name of the column the ordering came from
#'
#' Used to build self-describing occasion column names. A fit always stores the
#' name beside the values, but a hand-assembled fixture may carry only the
#' values, so a neutral fallback keeps the names syntactic either way.
#'
#' @param object A fitted model carrying an ordering.
#' @return A single character string.
#' @noRd
.multilpa_time_name <- function(object) {
  if (is.null(object$time) || !is.character(object$time) ||
      length(object$time) != 1L || is.na(object$time)) "time" else object$time
}

#' Assigned profiles as a groups-by-occasions integer matrix
#'
#' The rectangular layout the drawing code needs. A matrix is the right
#' structure inside a function body and the wrong one to return to a caller,
#' so this stays internal and the sequence table is a data frame instead.
#'
#' @param object A fitted model carrying an ordering.
#' @param long The long form, if the caller already has it; recomputed
#'   otherwise, so that an outside caller needs only the fit.
#' @return An integer matrix with one row per group and one column per position
#'   seen anywhere in the data, holding the assigned profile or `NA`. Row names
#'   are the group identifiers made unique, in increasing group order; column
#'   names are `<time>_<position>`.
#' @noRd
.multilpa_sequence_matrix <- function(object,
                                      long = .multilpa_sequences(object,
                                                       format = "long")) {
  positions <- sort(unique(long$time))
  groups <- unique(long$group)
  group_labels <- make.unique(as.character(object$group_values))[
    match(groups, object$group_values)]
  cells <- matrix(
    NA_integer_, length(groups), length(positions),
    dimnames = list(group_labels, .multilpa_occasion_names(object, positions)))
  cells[cbind(match(long$group, groups), match(long$time, positions))] <-
    long$profile
  cells
}

#' Syntactic, self-describing names for the occasion columns
#'
#' @param object A fitted model carrying an ordering.
#' @param positions The distinct positions, in increasing order.
#' @return A character vector of syntactic, unique column names.
#' @noRd
.multilpa_occasion_names <- function(object, positions) {
  make.unique(make.names(paste0(.multilpa_time_name(object), "_",
                                as.character(positions))))
}

#' Profile assignments in occasion order
#'
#' Documented on `?get_results`, which is where a caller reaches this table from.
#'
#' @param x A fitted model of this package, fitted with `time =`.
#' @param format `"long"` for one row per observation, `"wide"` for one row per
#'   group with one column per occasion.
#' @return A base `data.frame`.
#' @noRd
.multilpa_sequences <- function(x, format = c("long", "wide")) {
  time_values <- .multilpa_require_time(x)
  format <- match.arg(format)
  long <- data.frame(
    group = x$group_values[x$group_index],
    group_class = x$group_classes[x$group_index],
    time = time_values,
    profile = x$subject_profiles,
    row.names = NULL, stringsAsFactors = FALSE
  )
  long <- long[order(long$group, long$time), , drop = FALSE]
  row.names(long) <- NULL
  if (identical(format, "long")) return(long)

  cells <- .multilpa_sequence_matrix(x, long)
  groups <- unique(long$group)
  occasions <- as.data.frame(lapply(seq_len(ncol(cells)), function(column) {
    factor(cells[, column], levels = seq_len(x$n_profiles))
  }), optional = TRUE)
  names(occasions) <- colnames(cells)
  # The group identifier is a column, not a row name: row names coerce to
  # character and would silently merge two groups that print alike.
  data.frame(
    group = groups,
    group_class = x$group_classes[match(groups, x$group_values)],
    occasions,
    row.names = NULL, stringsAsFactors = FALSE, check.names = FALSE)
}

#' How much data each group class contributes
#'
#' Documented on `?get_results`, which is where a caller reaches this table from.
#'
#' @param x A fitted model of this package, fitted with `time =`.
#' @return A base `data.frame`, one row per group class, including any class no
#'   group was assigned to.
#' @noRd
.multilpa_sequence_summary <- function(x) {
  time_values <- .multilpa_require_time(x)
  positions <- sort(unique(time_values))
  n_positions <- length(positions)
  # Native integer indices avoid unused factor levels and character rounding of
  # distinct numeric group identifiers in tapply().
  lengths_by_group <- tabulate(x$group_index, nbins = x$n_groups)
  position_index <- match(time_values, positions)
  first <- as.integer(tapply(position_index, x$group_index, min))
  last <- as.integer(tapply(position_index, x$group_index, max))
  # A group whose observations fall short of its own span skipped a position
  # inside it; a group that simply stopped early does not count as a gap.
  has_gap <- lengths_by_group < (last - first + 1L)
  class_by_group <- x$group_classes
  classes <- seq_len(x$n_group_classes)
  summarize <- function(class) {
    in_class <- class_by_group == class
    taken <- lengths_by_group[in_class]
    data.frame(
      group_class = class,
      groups = length(taken),
      observations = sum(taken),
      mean_length = if (length(taken)) mean(taken) else NA_real_,
      median_length = if (length(taken)) stats::median(taken) else NA_real_,
      shortest = if (length(taken)) min(taken) else NA_integer_,
      longest = if (length(taken)) max(taken) else NA_integer_,
      complete = sum(taken == n_positions),
      gaps = sum(has_gap[in_class]),
      row.names = NULL
    )
  }
  do.call(rbind, lapply(classes, summarize))
}
