# The start search: which starts a fit tries, and how it chooses among them.
#
# `search = "kmeans"` is the search of latents 0.96.1 and earlier: `n_starts`
# starts (Ward's clustering first, then k-means from random centres; group
# classes from k-means of the group compositions or from random compositions),
# each run to the end. `search = "emEM"` screens a larger pool by short EM
# runs and runs only the best `n_starts` to the end (Biernacki, Celeux &
# Govaert, 2003); the pool begins with the kmeans search's own starts, drawn
# in the same order, followed by random partitions and random centres.

#' Settings of the emEM search that are not arguments
#'
#' Chosen by the start-search benchmark (SEARCH-PLAN.md, section 6): ten
#' screening iterations found a lower maximum than twenty on 9 of 200 EEE
#' datasets and a higher one on none.
#'
#' @return A list of `pool_factor` (the pool is this many times `n_starts`
#'   when `n_screen` is `NULL`), `iterations` (the screening budget) and
#'   `protect_first` (whether start 1, Ward's, always runs to the end).
#' @noRd
.latents_search_defaults <- function() {
  list(pool_factor = 5L, iterations = 20L, protect_first = TRUE)
}

#' Validate `n_screen` and resolve the screening plan of a fit
#'
#' @param search `"kmeans"` or `"emEM"`.
#' @param n_starts Starts run to the end.
#' @param n_screen `NULL` or the pool size.
#' @param max_iter The fit's iteration budget.
#' @param supplied Whether the caller supplied a `start`.
#' @return `NULL` when nothing is screened (`"kmeans"`, or a supplied start
#'   that runs alone), otherwise the `screen` list for
#'   `.latents_run_starts()`.
#' @noRd
.latents_search_screen <- function(search, n_starts, n_screen, max_iter,
                                   supplied = FALSE) {
  if (!is.null(n_screen)) {
    if (!is.numeric(n_screen) || length(n_screen) != 1L || !is.finite(n_screen) ||
        n_screen != round(n_screen) || n_screen < n_starts ||
        n_screen > .Machine$integer.max) {
      stop(errorCondition(sprintf(paste(
        "`n_screen` must be a single whole number no smaller than `n_starts`",
        "(%d): it is the number of starts screened, of which `n_starts` are run",
        "to the end."), as.integer(n_starts)),
        class = "latents_bad_argument", call = NULL))
    }
    if (identical(search, "kmeans")) {
      stop(errorCondition(paste(
        "`n_screen` sets the pool of the emEM search; `search = \"kmeans\"`",
        "screens nothing. Drop `n_screen`, or use `search = \"emEM\"`."),
        class = "latents_bad_argument", call = NULL))
    }
  }
  if (identical(search, "kmeans")) return(NULL)
  # A supplied start that is the only start is evaluated as it is.
  if (supplied && n_starts == 1L) return(NULL)
  defaults <- .latents_search_defaults()
  n_pool <- as.integer(n_screen %||% (defaults$pool_factor * n_starts))
  list(n_pool = n_pool, iterations = min(defaults$iterations, max_iter),
       protect = if (supplied || defaults$protect_first) 1L else integer())
}

#' The kind of start at a position of the pool
#'
#' Positions up to `n_starts` are the kmeans search's starts; after them the
#' pool alternates random partitions and random centres.
#'
#' @param start_index Position in the pool.
#' @param n_starts Starts run to the end.
#' @return `"kmeans"`, `"partition"` or `"centres"` (`"kmeans"` also for
#'   start 1, which the caller initializes by Ward's clustering).
#' @noRd
.latents_start_kind <- function(start_index, n_starts) {
  if (start_index <= n_starts) return("kmeans")
  if ((start_index - n_starts) %% 2L == 1L) "partition" else "centres"
}

#' The label of a start's kind in the starts table
#' @inheritParams .latents_start_kind
#' @return `"ward"`, `"kmeans"`, `"partition"` or `"centres"`.
#' @noRd
.latents_start_kind_label <- function(start_index, n_starts) {
  if (start_index == 1L) "ward" else .latents_start_kind(start_index, n_starts)
}

#' The search controls a fit records
#'
#' Kept separately from what the search found, so a refit (bootstrap,
#' sensitivity) repeats the search that was requested.
#'
#' @param search `"kmeans"` or `"emEM"`.
#' @param n_starts Starts run to the end.
#' @param screen The resolved screening plan, or `NULL`.
#' @param select_start The requested selection rule.
#' @param boundary_policy The requested `boundary` (the fit's own `boundary`
#'   field is whether the selected solution is at the bound).
#' @return A list of `search`, `n_starts`, `n_screen` (`NULL` when nothing
#'   was screened), `select_start` and `boundary_policy`.
#' @noRd
.latents_search_record <- function(search, n_starts, screen, select_start,
                                   boundary_policy) {
  list(search = search, n_starts = as.integer(n_starts),
       n_screen = screen$n_pool, select_start = select_start,
       boundary_policy = boundary_policy)
}

#' The fitting controls a bootstrap refit of a model uses
#'
#' A refit repeats the search and the stopping rule the model was fitted
#' with, so the statistic and its reference distribution come from the same
#' procedure; an explicit `n_starts`, `max_iter` or `tol` overrides them for
#' every model. A fit made before these controls were recorded falls back to
#' the values `bootstrap_lrt()` used to apply to every refit (10 starts, 1000
#' iterations, `1e-8`, the default search, selection and acceleration).
#'
#' @param model A fitted model.
#' @param n_starts,max_iter,tol `NULL` to inherit, or an override.
#' @return A named list of `search`, `n_starts`, `n_screen`, `max_iter`,
#'   `tol`, `acceleration`, `select_start` and `boundary`.
#' @noRd
.latents_refit_controls <- function(model, n_starts = NULL, max_iter = NULL,
                                    tol = NULL) {
  recorded <- !is.null(model$search)
  starts <- as.integer(n_starts %||% (if (recorded) model$n_starts) %||% 10L)
  search <- model$search %||% "kmeans"
  # A changed number of starts takes the default pool for that number.
  n_screen <- if (identical(search, "emEM") &&
                  identical(starts, as.integer(model$n_starts))) model$n_screen
  list(search = search, n_starts = starts, n_screen = n_screen,
       max_iter = as.integer(max_iter %||% (if (recorded) model$max_iter) %||% 1000L),
       tol = tol %||% (if (recorded) model$tol) %||% 1e-8,
       acceleration = model$acceleration %||% "squarem",
       select_start = model$select_start %||% "likelihood",
       boundary = model$boundary_policy %||% "allow")
}

#' The refit controls of a bootstrap, one row per model
#' @param controls A named list (`null`, `alternative`) from
#'   `.latents_refit_controls()`.
#' @return A data frame, one row per model, with the role and every control
#'   (`n_screen` is `NA` when nothing is screened).
#' @noRd
.latents_refit_controls_frame <- function(controls) {
  rows <- lapply(names(controls), function(role) {
    control <- controls[[role]]
    data.frame(role = role, search = control$search, n_starts = control$n_starts,
               n_screen = control$n_screen %||% NA_integer_,
               max_iter = control$max_iter, tol = control$tol,
               acceleration = control$acceleration,
               select_start = control$select_start, boundary = control$boundary,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Add the emEM columns to a starts table
#'
#' Under emEM every start of the pool has a row: its kind, whether it was run
#' to the end, and its likelihood after screening (for a start run to the end,
#' `log_likelihood` is its final one). A kmeans table is returned unchanged.
#'
#' @param table The engine's starts table, one row per start.
#' @param found The result of `.latents_run_starts()`.
#' @param n_starts Starts run to the end.
#' @param supplied Whether start 1 was supplied by the caller.
#' @return The table, with `kind`, `stage` and `screened_log_likelihood`
#'   under emEM.
#' @noRd
.latents_search_columns <- function(table, found, n_starts, supplied = FALSE) {
  if (is.null(found$stage)) return(table)
  table$kind <- vapply(seq_len(nrow(table)), function(start_index) {
    if (start_index == 1L && supplied) "supplied" else
      .latents_start_kind_label(start_index, n_starts)
  }, character(1))
  table$stage <- found$stage
  table$screened_log_likelihood <- found$screened_scores
  table
}

#' Whether each start was run to the end
#' @param found The result of `.latents_run_starts()`.
#' @return A logical vector, one per start; all `TRUE` for the kmeans search.
#' @noRd
.latents_final_stage <- function(found) {
  if (is.null(found$stage)) rep(TRUE, length(found$attempts)) else
    found$stage == "final"
}

