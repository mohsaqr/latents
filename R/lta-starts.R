#' Bind only forced fitting inputs into a runner or finaliser
#' @param fun A namespace function accepting a context argument.
#' @param context A list of evaluated fitting inputs.
#' @param record Whether to retain covariance M-step counts on the attempt.
#' @return A closure with no promises referring to the public caller.
#' @noRd
.lta_bind_context <- function(fun, context, record = FALSE) {
  force(fun)
  force(context)
  force(record)
  function(...) {
    if (!record) return(fun(..., context = context))
    .multilpa_reset_structure_log()
    attempt <- tryCatch(fun(..., context = context), error = function(error) {
      if (inherits(error, "latents_parallel_failure")) stop(error)
      list(error = conditionMessage(error))
    })
    attr(attempt, "lta_structure_log") <- as.list(.multilpa_structure_log)
    attempt
  }
}

#' Keep setup and finalisation in the parent and isolate worker inputs
#' @param context Evaluated setup inputs.
#' @param initial Ordered, wrapped starting values drawn during setup.
#' @param runner,finaliser Namespace functions for one start and finalisation.
#' @param run_vars Names of the setup inputs needed on workers.
#' @param n_observed The observed values the EM stopping rule counts: with
#'   `context$tol` it gives the start-selection tie (see
#'   `.multilpa_select_start()`).
#' @return A prepared fitting plan shared by lta and enumeration.
#' @noRd
.lta_plan <- function(context, initial, runner, finaliser, run_vars, n_observed) {
  stopifnot("runner inputs are present" = all(run_vars %in% names(context)))
  list(initial = initial, n_starts = context$n_starts,
       screen = context$screen, select_start = context$select_start,
       tie = context$tol * max(1, n_observed),
       run = .lta_bind_context(runner, context[run_vars], record = TRUE),
       finalise = .lta_bind_context(finaliser, context),
       structure_log = if (!is.null(context$structure))
         as.list(.multilpa_structure_log) else list(capped = 0L),
       report_structure = !is.null(context$structure))
}

#' Run the starts of a single prepared fit using the shared search
#' @param plan A prepared fitting plan.
#' @param cores Core budget for a direct fit.
#' @return A fitted transition model.
#' @noRd
.lta_run_plan <- function(plan, cores) {
  starts <- tryCatch(.latents_run_starts(
    plan$n_starts, plan$run, plan$select_start, initial = plan$initial,
    cores = cores, screen = plan$screen, tie = plan$tie), error = function(error) {
      .lta_report_failed_plan(plan, error)
      stop(error)
    })
  .lta_finish_plan(plan, starts)
}

#' Restore a fit's covariance counts from its parent and worker phases
#' @param plan A prepared fitting plan.
#' @param logs Covariance count records from start attempts.
#' @return NULL, invisibly.
#' @noRd
.lta_restore_structure_log <- function(plan, logs) {
  logs <- c(list(plan$structure_log), logs)
  counts <- vapply(logs, function(log) log$capped %||% 0L, integer(1))
  .multilpa_reset_structure_log()
  if (any(counts > 0L)) {
    latest <- logs[[utils::tail(which(counts > 0L), 1L)]]
    latest$capped <- sum(counts)
    invisible(list2env(latest, envir = .multilpa_structure_log))
  }
  invisible(NULL)
}

#' Preserve covariance qualifications when all starts fail
#' @param plan A prepared fitting plan.
#' @param error The selection error carrying its attempts' count records.
#' @return NULL, invisibly.
#' @noRd
.lta_report_failed_plan <- function(plan, error) {
  logs <- attr(error, "lta_structure_logs")
  if (plan$report_structure && !is.null(logs)) {
    .lta_restore_structure_log(plan, logs)
    .multilpa_report_structure_log()
  }
  invisible(NULL)
}

#' Finalise in the parent and report accumulated covariance qualifications
#' @param plan A prepared fitting plan.
#' @param starts The shared start-selection result.
#' @return A fitted transition model.
#' @noRd
.lta_finish_plan <- function(plan, starts) {
  if (plan$report_structure) {
    logs <- attr(starts, "lta_structure_logs") %||%
      lapply(starts$attempts, attr, which = "lta_structure_log")
    .lta_restore_structure_log(plan, logs)
    on.exit(.multilpa_report_structure_log(), add = TRUE)
  }
  plan$finalise(starts)
}

#' Execute one transition start through the shared task runner
#' @param task A self-contained start task.
#' @return A full or screened start attempt.
#' @noRd
.lta_execute_start <- function(task) .latents_execute_start(task)

#' Build transition tasks using the shared plan layout
#' @param candidates Collected setup results in grid order.
#' @param phase Screening or final runs.
#' @param screened,ranked Per-candidate screening results and ranking records.
#' @return Self-contained start tasks.
#' @noRd
.lta_start_tasks <- function(candidates, phase = c("screen", "final"),
                             screened = NULL, ranked = NULL) {
  .latents_start_tasks(candidates, match.arg(phase), screened, ranked)
}

#' Prepare one transition candidate for the shared pool
#' @param task The candidate's forced `.lta_setup()` arguments.
#' @return Its plan, or the candidate's error condition.
#' @noRd
.lta_setup_task <- function(task) {
  tryCatch(do.call(.lta_setup, task$arguments, quote = TRUE),
           error = .lta_candidate_error)
}

#' Run a transition phase in the shared worker pool
#' @param tasks Self-contained start tasks.
#' @param costs Estimated candidate costs in grid order.
#' @param cores Total worker budget.
#' @param n_candidates Number of grid rows, including failed setups.
#' @return Collected attempts grouped by candidate and original task order.
#' @noRd
.lta_pool_phase <- function(tasks, costs, cores, n_candidates) {
  .latents_pool_phase(tasks, costs, cores, n_candidates, .lta_execute_start)
}

#' Keep ordinary candidate errors while propagating worker loss
#' @param error A caught error condition.
#' @return The condition for the enumeration's error column.
#' @noRd
.lta_candidate_error <- function(error) {
  if (inherits(error, "latents_parallel_failure")) stop(error)
  error
}
