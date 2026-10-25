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

