# Naming the candidates of an enumeration grid, and finding one by its name.
#
# One name per candidate, used everywhere: the grid's `candidate` column,
# `grid$<name>`, `get_candidate(grid, "<name>")` and the rows of
# `compare_models(grid, ...)`. Names are derived from the grid's table and
# fits each time, so a grid saved by an earlier version is named the same way.

#' The name of every candidate of a grid
#'
#' `profiles3_groups4` for 3 profiles and 4 group classes; `classes` instead of
#' `profiles` when every indicator is categorical; no `_groups` part on a
#' single-level grid; `_<structure>` appended when the grid crosses several
#' covariance structures; `<model>_groups2` on a group-class family grid,
#' from its `model` column (`additive_groups2`, `additive_equal_groups2`).
#'
#' @param x A `multilpa_enumeration`, `latents_transition_enumeration` or
#'   `latents_family_enumeration`.
#' @return A character vector, one unique name per row of the grid's table.
#' @noRd
.latents_candidate_names <- function(x) {
  table <- .subset2(x, "table")
  if (inherits(x, "latents_family_enumeration")) {
    names <- paste0(table$model, "_groups", table$n_group_classes)
  } else {
    fitted <- Filter(Negate(is.null), .subset2(x, "fits"))
    single_level <- length(fitted) > 0L && isTRUE(fitted[[1L]]$single_level)
    names <- paste0(.latents_class_word(x), table$n_profiles,
                    if (single_level) "" else paste0("_groups", table$n_group_classes))
    models <- table$model
    if (length(unique(models[!is.na(models)])) > 1L) {
      names <- ifelse(is.na(models), names, paste(names, models, sep = "_"))
    }
  }
  stopifnot("every candidate of a grid has its own name" = !anyDuplicated(names))
  names
}

#' Name the caller's data in a grid made through a wrapper
#'
#' `enumerate_lpa()` and `enumerate_lca()` call [enumerate_classes()] with
#' their own `data` argument, so the grid and its candidates would record
#' `data = data`. The wrapper passes its caller's expression here instead.
#'
#' @param grid A grid returned by [enumerate_classes()].
#' @param data_expression The wrapper caller's expression for `data`.
#' @return The grid, its own call and every candidate's naming that data.
#' @noRd
.latents_candidate_data <- function(grid, data_expression) {
  fits <- lapply(.subset2(grid, "fits"), function(fit) {
    if (is.null(fit) || is.null(fit$call)) return(fit)
    fit$call$data <- data_expression
    fit
  })
  grid[["fits"]] <- fits
  call <- .subset2(grid, "call")
  if (is.call(call)) {
    call$data <- data_expression
    grid[["call"]] <- call
  }
  grid
}

#' Is this an enumeration grid of this package?
#' @param x Any object.
#' @return A single logical.
#' @noRd
.latents_is_grid <- function(x) {
  inherits(x, c("multilpa_enumeration", "latents_transition_enumeration",
                "latents_family_enumeration"))
}

#' The row of a grid one candidate names
#'
#' A candidate is named by its name (`"classes3_groups4"`), by its class
#' counts as a pair (`c(3, 4)`), or by one count, which means one group class.
#'
#' @param x An enumeration grid.
#' @param candidate One name, or one or two class counts.
#' @return The row index.
#' @noRd
.latents_candidate_row <- function(x, candidate) {
  names <- .latents_candidate_names(x)
  table <- .subset2(x, "table")
  available <- function() {
    paste(sprintf("\"%s\"", names), collapse = ", ")
  }
  if (is.character(candidate) && length(candidate) == 1L && !is.na(candidate)) {
    at <- match(candidate, names)
    if (is.na(at)) {
      stop(errorCondition(sprintf(
        "\"%s\" is not a candidate of this grid. Its candidates are %s.",
        candidate, available()), class = "latents_unknown_candidate", call = NULL))
    }
    return(at)
  }
  counts_given <- is.numeric(candidate) && length(candidate) %in% 1:2 &&
    all(is.finite(candidate)) && all(candidate >= 1) && all(candidate == floor(candidate))
  if (!counts_given || inherits(x, "latents_family_enumeration")) {
    stop(errorCondition(sprintf(paste(
      "%s does not name a candidate. Use its name, one of %s%s."),
      deparse1(candidate), available(),
      if (inherits(x, "latents_family_enumeration")) "" else
        ", or its class counts, c(3, 4)"),
      class = "latents_bad_argument", call = NULL))
  }
  group_classes <- if (length(candidate) == 2L) candidate[2L] else 1
  at <- which(table$n_profiles == candidate[1L] &
                table$n_group_classes == group_classes)
  if (length(at) == 0L) {
    stop(errorCondition(sprintf(
      "No candidate has %d %s and %d group classes. Its candidates are %s.",
      as.integer(candidate[1L]), .latents_class_word(x), as.integer(group_classes),
      available()),
      class = "latents_unknown_candidate", call = NULL))
  }
  if (length(at) > 1L) {
    # A grid crossed with covariance structures has several candidates per
    # pair of class counts; only the name tells them apart.
    stop(errorCondition(sprintf(
      "%d candidates have these class counts: %s. Name one.", length(at),
      paste(sprintf("\"%s\"", names[at]), collapse = ", ")),
      class = "latents_unknown_candidate", call = NULL))
  }
  at
}

#' The fitted model in one row of a grid
#' @param x An enumeration grid.
#' @param at A row index.
#' @return The fit; a failed candidate raises `latents_failed_candidate`.
#' @noRd
.latents_candidate_at <- function(x, at) {
  fit <- .subset2(x, "fits")[[at]]
  if (is.null(fit)) {
    stop(errorCondition(sprintf(
      "Candidate \"%s\" could not be fitted: %s",
      .latents_candidate_names(x)[at], .subset2(x, "table")$error[at]),
      class = "latents_failed_candidate", call = NULL))
  }
  fit
}

#' A grid's table with each candidate's name first
#' @param x An enumeration grid.
#' @param table The table to name, by default the stored one.
#' @return `table` with a leading `candidate` column.
#' @noRd
.latents_named_table <- function(x, table = .subset2(x, "table")) {
  table <- table[setdiff(names(table), "candidate")]
  cbind(data.frame(candidate = .latents_candidate_names(x),
                   stringsAsFactors = FALSE), table)
}

#' `grid$<name>`: a component, or the candidate of that name
#'
#' The grid's own components (`table`, `fits`, `call`) come first, matched as
#' `$` matches on any list; any other name is looked up among the candidates.
#' A name that is neither returns `NULL`, as `$` does, and so does a candidate
#' whose fit failed; [get_candidate()] raises an error for both.
#'
#' @param x An enumeration grid.
#' @param name The name after `$`.
#' @return The component or the fitted model.
#' @noRd
.latents_grid_dollar <- function(x, name) {
  component <- .subset2(x, name, exact = FALSE)
  if (!is.null(component)) return(component)
  at <- match(name, .latents_candidate_names(x))
  if (is.na(at)) return(NULL)
  .subset2(x, "fits")[[at]]
}

#' @noRd
#' @export
`$.multilpa_enumeration` <- function(x, name) .latents_grid_dollar(x, name)

#' @noRd
#' @export
`$.latents_transition_enumeration` <- function(x, name) .latents_grid_dollar(x, name)

#' @noRd
#' @export
`$.latents_family_enumeration` <- function(x, name) .latents_grid_dollar(x, name)

#' Names an editor completes after `grid$`
#' @param x An enumeration grid.
#' @param pattern The partial name typed so far.
#' @return The components and candidate names matching it.
#' @noRd
.latents_grid_dollar_names <- function(x, pattern = "") {
  choices <- c(names(unclass(x)), .latents_candidate_names(x))
  grep(pattern, choices, value = TRUE)
}

#' @noRd
#' @exportS3Method utils::.DollarNames
.DollarNames.multilpa_enumeration <- function(x, pattern = "") {
  .latents_grid_dollar_names(x, pattern)
}

#' @noRd
#' @exportS3Method utils::.DollarNames
.DollarNames.latents_transition_enumeration <- function(x, pattern = "") {
  .latents_grid_dollar_names(x, pattern)
}

#' @noRd
#' @exportS3Method utils::.DollarNames
.DollarNames.latents_family_enumeration <- function(x, pattern = "") {
  .latents_grid_dollar_names(x, pattern)
}

#' Candidate and minimum-criterion tables of a family or transition grid
#'
#' The minimum is taken over the candidates that converged, as the profile
#' grid's criteria table takes it: an unconverged fit is reported in the
#' candidates table, never returned as the minimum.
#'
#' @param table The named candidates table.
#' @param what `"candidates"`, `"best"` or `"all"`.
#' @param criterion The criterion column.
#' @return The table, the minimum row, or a list of both; `"all"` without a
#'   converged candidate returns the candidates alone, and `"best"` raises
#'   `latents_failed_candidate`.
#' @noRd
.latents_grid_result_tables <- function(table, what, criterion) {
  if (identical(what, "candidates")) return(table)
  usable <- which(is.finite(table[[criterion]]) & table$converged %in% TRUE)
  if (length(usable) == 0L) {
    if (identical(what, "all")) return(list(candidates = table))
    stop(errorCondition(
      "No candidate was fitted and converged, so no criterion has a minimum.",
      class = "latents_failed_candidate", call = NULL))
  }
  minimum <- table[usable[which.min(table[[criterion]][usable])], , drop = FALSE]
  if (identical(what, "all")) return(list(candidates = table, best = minimum))
  minimum
}
