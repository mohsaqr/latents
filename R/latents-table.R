# A tidy table that prints like a report.
#
# Result tables stay base data frames: numeric columns, one row per unit,
# usable directly in computation. Only printing differs. A table may carry a
# display plan (which columns to show, under which labels, in which format);
# the printer draws those columns aligned under a ruled header, with human
# labels, a fixed number of decimals per column, intervals as "[low, high]"
# and APA-style p-values. Long tables print their first rows.

#' Mark a data frame as a latents result table
#' @param table A data frame.
#' @param level The confidence level of its intervals, if any.
#' @param title A title printed above the table.
#' @param display A display plan (a list of `.latents_column()`), or `NULL`
#'   to show every column in a default format.
#' @return `table` with class `latents_table` prepended.
#' @noRd
.latents_table <- function(table, level = NULL, title = NULL, display = NULL) {
  stopifnot("a result table must be a data frame" = is.data.frame(table))
  row.names(table) <- NULL
  attr(table, "level") <- level
  attr(table, "title") <- title
  attr(table, "display") <- display
  class(table) <- c("latents_table", setdiff(class(table), "latents_table"))
  table
}

#' One displayed column of a result table
#'
#' @param label The header.
#' @param column The source column, or for `"ci"` the lower and upper bounds.
#' @param format `"text"`, `"label"` (class and term names made readable),
#'   `"number"`, `"share"` (three decimals), `"integer"`, `"ci"`, `"p"`,
#'   `"percent"` or `"logical"`.
#' @param decimals_from For `"ci"`, the column whose decimals the bounds use.
#' @noRd
.latents_column <- function(label, column, format = "number", decimals_from = NULL) {
  list(label = label, column = column, format = format, decimals_from = decimals_from)
}

#' Decimals for a column: two, three or four, so its smallest value keeps a
#' significant digit
#' @noRd
.latents_column_decimals <- function(values) {
  finite <- abs(values[is.finite(values) & values != 0])
  smallest <- if (length(finite) > 0L) min(finite) else NA_real_
  if (!is.finite(smallest) || smallest >= 0.1) return(2L)
  if (smallest >= 0.01) return(3L)
  4L
}

#' Readable class and term names
#' @noRd
.latents_label <- function(values) {
  values <- as.character(values)
  values <- sub("^class_([0-9]+)$", "Class \\1", values)
  values <- sub("^group_class_([0-9]+)$", "Group class \\1", values)
  values <- sub("^profile_([0-9]+)$", "Profile \\1", values)
  values <- gsub("(Intercept)", "Intercept", values, fixed = TRUE)
  values <- gsub(":", " x ", values, fixed = TRUE)
  values[values %in% c("common", "all")] <- "All classes"
  words <- c(variance = "Variance", sd = "SD", correlation = "Correlation",
             covariance = "Covariance")
  known <- values %in% names(words)
  values[known] <- words[values[known]]
  values
}

#' APA-style p-values: no leading zero, "<.001" below a thousandth
#' @noRd
.latents_format_p <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "<.001", sub("^0", "", sprintf("%.3f", p))))
}

#' The cells of one displayed column, as text
#' @noRd
.latents_cells <- function(x, spec) {
  values <- x[[spec$column[1L]]]
  fixed <- function(v, decimals) {
    ifelse(is.na(v), "", formatC(v, format = "f", digits = decimals))
  }
  switch(spec$format,
    text = ifelse(is.na(values), "", as.character(values)),
    label = .latents_label(values),
    integer = ifelse(is.na(values), "", format(round(values), trim = TRUE)),
    number = {
      finite <- values[is.finite(values)]
      if (length(finite) > 0L && all(abs(finite - round(finite)) < 1e-9)) {
        ifelse(is.na(values), "", format(round(values), trim = TRUE))
      } else {
        large <- length(finite) > 0L && max(abs(finite)) >= 100
        fixed(values, if (large) 2L else .latents_column_decimals(values))
      }
    },
    share = fixed(values, 3L),
    percent = ifelse(is.na(values), "", sprintf("%.1f%%", 100 * values)),
    p = .latents_format_p(values),
    logical = ifelse(is.na(values), "", ifelse(values, "yes", "no")),
    ci = {
      decimals <- .latents_column_decimals(x[[spec$decimals_from %||% spec$column[1L]]])
      low <- fixed(x[[spec$column[1L]]], decimals)
      high <- fixed(x[[spec$column[2L]]], decimals)
      low <- formatC(low, width = max(nchar(low)))
      high <- formatC(high, width = max(nchar(high)))
      ifelse(is.na(x[[spec$column[1L]]]), "", sprintf("[%s, %s]", low, high))
    })
}

#' The display plan of a table without one: every column, by type
#' @noRd
.latents_default_display <- function(x) {
  lapply(names(x), function(name) {
    values <- x[[name]]
    format <- if (is.logical(values)) "logical" else if (is.numeric(values)) {
      if (grepl("(^|_)p_(value|adjusted)$", name)) "p" else
        if (all(is.na(values) | abs(values - round(values)) < 1e-9) &&
            !grepl("share|weight|probability|entropy|posterior", name)) "integer" else
          "number"
    } else "text"
    .latents_column(name, name, format)
  })
}

#' Draw a one-row table as label-value lines
#' @noRd
.latents_render_card <- function(x, display) {
  frame <- as.data.frame(x)
  labels <- vapply(display, `[[`, character(1), "label")
  values <- vapply(display, function(spec) .latents_cells(frame, spec)[1L], character(1))
  # An entry the fit does not have (for example groups of a single-level
  # fit) is left out rather than printed empty.
  labels <- labels[nzchar(values)]
  values <- values[nzchar(values)]
  width <- max(nchar(labels, type = "width"))
  c(if (!is.null(attr(x, "title"))) c(attr(x, "title"), ""),
    sprintf("%s  %s", formatC(labels, width = -width), values))
}

#' Draw a result table as aligned text lines
#' @param x A `latents_table`.
#' @param n Most rows drawn.
#' @return A character vector of lines.
#' @noRd
.latents_render <- function(x, n = 20L) {
  display <- attr(x, "display") %||% .latents_default_display(x)
  if (isTRUE(attr(x, "card")) && nrow(x) == 1L) return(.latents_render_card(x, display))
  marks <- attr(x, "marks")
  title <- attr(x, "title")
  frame <- as.data.frame(x)
  shown <- utils::head(frame, n)
  cells <- lapply(display, function(spec) .latents_cells(shown, spec))
  labels <- vapply(display, `[[`, character(1), "label")
  right <- vapply(display, function(spec) !spec$format %in% c("text", "label"), logical(1))
  widths <- vapply(seq_along(cells), function(j) {
    as.integer(max(nchar(labels[j], type = "width"), nchar(cells[[j]], type = "width"), 1L))
  }, integer(1))
  pad <- function(text, width, align_right) {
    gap <- strrep(" ", pmax(width - nchar(text, type = "width"), 0L))
    if (align_right) paste0(gap, text) else paste0(text, gap)
  }
  line_of <- function(values) {
    paste(vapply(seq_along(values), function(j) pad(values[j], widths[j], right[j]),
                 character(1)), collapse = "  ")
  }
  body <- vapply(seq_len(nrow(shown)), function(i) {
    line_of(vapply(cells, `[`, character(1), i))
  }, character(1))
  if (length(marks) > 0L) body <- paste0(body, utils::head(marks, nrow(shown)))
  c(if (!is.null(title)) c(title, ""),
    line_of(labels),
    paste(vapply(widths, function(w) strrep("-", w), character(1)), collapse = "  "),
    body,
    if (nrow(frame) > n) sprintf("... %d more rows", nrow(frame) - n))
}

#' Print a latents result table
#'
#' Prints the columns people report, with readable labels and numbers; the
#' object itself is unchanged, a base data frame with every numeric column.
#'
#' @param x A `latents_table`.
#' @param n The most rows printed.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.latents_table <- function(x, n = 20L, ...) {
  if (nrow(x) == 0L) {
    cat(c(attr(x, "title"), "<no rows>"), sep = "\n")
    return(invisible(x))
  }
  cat(.latents_render(x, n), sep = "\n")
  invisible(x)
}

#' @rdname print.latents_table
#' @param row.names,optional Passed to the data frame method.
#' @return For `as.data.frame()`, the same table as a plain base data frame.
#' @export
as.data.frame.latents_table <- function(x, row.names = NULL, optional = FALSE, ...) {
  attributes(x)[c("level", "title", "display", "marks", "card")] <- NULL
  class(x) <- setdiff(class(x), "latents_table")
  x
}
