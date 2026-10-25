#' Execute a self-contained start task without creating workers
#' @param task A start index, isolated runner, and initial values or checkpoint.
#' @return A full or screened start attempt.
#' @noRd
.latents_execute_start <- function(task) {
  if (!is.null(task$screen)) {
    .latents_screen_attempt(task$run, task$start_index, task$initial, task$screen)
  } else {
    .latents_start_attempt(task$run, task$start_index, task$initial,
                           resume = task$resume)
  }
}

#' Build flat tasks while retaining candidate and original start order
#' @param candidates Collected setup results, one per grid row.
#' @param phase Screening or final runs.
#' @param screened,ranked Per-candidate screening results and ranking records.
#' @return Self-contained worker tasks; no finaliser or whole grid is captured.
#' @noRd
.latents_start_tasks <- function(candidates, phase = c("screen", "final"),
                             screened = NULL, ranked = NULL) {
  phase <- match.arg(phase)
  lapply(seq_along(candidates), function(i) {
    plan <- candidates[[i]]$value
    if (inherits(plan, "error")) return(list())
    screening <- !is.null(plan$screen)
    if (phase == "screen" && !screening) return(list())
    if (phase == "final" && screening && inherits(ranked[[i]]$value, "error")) {
      return(list())
    }
    starts <- if (phase == "screen") seq_along(plan$initial) else
      if (screening) ranked[[i]]$value$finalists else seq_along(plan$initial)
    lapply(starts, function(start_index) {
      list(candidate = i, start_index = start_index, run = plan$run,
           initial = if (phase == "final" && screening) NULL else
             plan$initial[[start_index]],
           resume = if (phase == "final" && screening)
             screened[[i]][[start_index]]$value$checkpoint else NULL,
           screen = if (phase == "screen") plan$screen else NULL)
    })
  }) |> unlist(recursive = FALSE, use.names = FALSE)
}

#' Run one flat search phase and group its collected results in grid order
#' @param tasks Self-contained start tasks.
#' @param costs Estimated candidate costs in grid order.
#' @param cores Total worker budget.
#' @param n_candidates Number of grid rows, including failed setups.
#' @param execute Function running one self-contained task.
#' @param on_lost Optional terminal process-loss callback for the mapper.
#' @return Collected attempts grouped by candidate, then original task order.
#' @noRd
.latents_pool_phase <- function(tasks, costs, cores, n_candidates,
                                  execute, on_lost = NULL) {
  candidate <- vapply(tasks, `[[`, integer(1), "candidate")
  results <- if (length(tasks)) .latents_map_candidates(
    tasks, execute, cores, costs = costs[candidate], signal = FALSE,
    on_lost = on_lost) else list()
  lapply(seq_len(n_candidates), function(i) results[candidate == i])
}


#' Prepare every candidate of a pooled grid
#'
#' A seeded candidate sets its own seed before drawing its starts
#' (`.latents_local_seed()`), so preparing it on a worker gives the plan this
#' process would. Unseeded candidates draw from the caller's random stream,
#' which only this process can continue in grid order; they, and every grid
#' on one core, are prepared here in grid order, as they always have been.
#' @param tasks Self-contained setup tasks in grid order.
#' @param setup A function of one task returning its plan, or the candidate's
#'   error condition.
#' @param costs Estimated candidate costs in grid order.
#' @param cores Total worker budget.
#' @param distribute Whether the setups may run on workers: every candidate
#'   is seeded and at least one has a setup worth distributing.
#' @return Collected setups (`value` and `conditions`) in grid order.
#' @noRd
.latents_pool_setup <- function(tasks, setup, costs, cores, distribute) {
  stopifnot("`setup` must be a function" = is.function(setup),
            "one cost per setup task" = length(costs) == length(tasks),
            "`distribute` must be TRUE or FALSE" =
              isTRUE(distribute) || isFALSE(distribute))
  if (!distribute || cores <= 1L || length(tasks) <= 1L) {
    return(lapply(tasks, function(task) {
      .latents_collect_conditions(function() setup(task))
    }))
  }
  .latents_map_candidates(tasks, setup, cores, costs = costs, signal = FALSE,
    on_lost = function(task, message) {
      errorCondition(message, class = "latents_candidate_failed", call = NULL)
    })
}

#' Report one phase of a pooled grid, when progress was asked for
#'
#' The message has its own class, `latents_progress`, so it can be told apart
#' from a fit's qualifications, which do not depend on the number of cores.
#' @param verbose Whether to report.
#' @param format,... A `sprintf()` format and its values.
#' @return NULL, invisibly.
#' @noRd
.latents_pool_progress <- function(verbose, format, ...) {
  if (isTRUE(verbose)) {
    # Base R has no messageCondition(); this is the same structure.
    message(structure(class = c("latents_progress", "message", "condition"),
                      list(message = paste0(sprintf(format, ...), "\n"), call = NULL)))
  }
  invisible(NULL)
}

#' A count with its noun, singular for one
#' @param n A whole number.
#' @param noun The singular noun; the plural adds an "s".
#' @return One string such as `"1 worker"` or `"12 workers"`.
#' @noRd
.latents_count <- function(n, noun) {
  n <- as.integer(n)
  sprintf("%d %s%s", n, noun, if (identical(n, 1L)) "" else "s")
}

#' A duration in words for progress messages
#' @param seconds Elapsed seconds.
#' @return One string such as `"under a minute"`, `"42 min"` or `"1 h 12 min"`.
#' @noRd
.latents_format_duration <- function(seconds) {
  stopifnot("`seconds` must be one non-negative number" =
              is.numeric(seconds) && length(seconds) == 1L && !is.na(seconds) &&
              seconds >= 0)
  minutes <- as.integer(round(seconds / 60))
  if (minutes < 1L) return("under a minute")
  if (minutes < 60L) return(sprintf("%d min", minutes))
  sprintf("%d h %d min", minutes %/% 60L, minutes %% 60L)
}
