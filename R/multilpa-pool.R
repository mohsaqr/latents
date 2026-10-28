#' Whether a candidate uses a separate fitting engine
#' @param arguments Forced arguments of one profile candidate.
#' @return TRUE for membership covariates or joined categorical pairs.
#' @noRd
.multilpa_pool_fallback <- function(arguments) {
  if (length(arguments$profile_covariates) || length(arguments$group_covariates)) {
    return(TRUE)
  }
  length(arguments$local_dependence) > 0L &&
    any(.latents_parse_pairs(arguments$local_dependence)$item_1 %in%
          arguments$categorical)
}

#' Estimate measurement work per profile for start scheduling
#' @param data Observed data frame.
#' @param vars Indicator names.
#' @param extra Additional candidate arguments.
#' @return A positive approximate parameter width per profile.
#' @noRd
.multilpa_pool_width <- function(data, vars, extra) {
  stopifnot("data is a data frame" = is.data.frame(data),
            "vars are indicator names" = is.character(vars))
  discrete <- intersect(intersect(vars, names(data)), c(extra$categorical, extra$ordinal))
  category_width <- sum(vapply(discrete, function(name) {
    max(1, length(unique(data[[name]][!is.na(data[[name]])])) - 1)
  }, numeric(1)))
  # Means and variances (or a count mean and dispersion) give two entries;
  # categorical probabilities contribute one fewer than their observed levels.
  max(1, category_width + 2 * length(setdiff(vars, discrete)))
}

#' Execute a flat start or one candidate using a separate fitting engine
#' @param task Self-contained start inputs, or forced candidate arguments.
#' @return A start attempt, a complete fallback fit, or its error condition.
#' @noRd
.multilpa_execute_start <- function(task) {
  if (is.null(task$arguments)) return(.latents_execute_start(task))
  tryCatch({
    fit <- do.call(multilpa, task$arguments, quote = TRUE)
    fit$call <- task$call
    fit
  }, error = .lta_candidate_error)
}

#' Represent a terminal worker loss in the candidate's existing failure form
#' @param task One start or fallback candidate task.
#' @param message Worker-loss explanation.
#' @return A failed-start entry or a failed-candidate condition.
#' @noRd
.multilpa_pool_loss <- function(task, message) {
  if (is.null(task$arguments)) list(error = message) else
    errorCondition(message, class = "latents_candidate_failed", call = NULL)
}

#' Prepare one profile candidate for the shared pool
#' @param task The candidate's forced arguments and call.
#' @return Its plan; the arguments alone for a candidate using a separate
#'   engine; or the candidate's error condition.
#' @noRd
.multilpa_setup_task <- function(task) {
  tryCatch({
    if (.multilpa_pool_fallback(task$arguments)) {
      # These engines retain serial starts. Their whole candidate occupies
      # one pool worker, alongside the ordinary candidates' individual starts.
      list(arguments = task$arguments, call = task$call)
    } else {
      do.call(.multilpa_setup, c(task$arguments, list(.call = task$call)),
              quote = TRUE)
    }
  }, error = .lta_candidate_error)
}

#' Build profile tasks, including serial-engine candidate fallbacks
#' @param candidates Collected plans in grid order.
#' @param phase Screening or final runs.
#' @param screened,ranked Per-candidate screening results and rankings.
#' @return Self-contained tasks for one shared pool.
#' @noRd
.multilpa_start_tasks <- function(candidates, phase, screened = NULL, ranked = NULL) {
  tasks <- .latents_start_tasks(candidates, phase, screened, ranked)
  if (phase == "screen") return(tasks)
  fallback <- lapply(seq_along(candidates), function(i) {
    plan <- candidates[[i]]$value
    if (inherits(plan, "error") || is.null(plan$arguments)) return(list())
    list(list(candidate = i, arguments = plan$arguments, call = plan$call))
  }) |> unlist(recursive = FALSE, use.names = FALSE)
  c(tasks, fallback)
}

#' Prepare, pool and finalise a profile grid without nested workers
#' @param arguments Forced argument lists in grid order.
#' @param calls Candidate calls to retain on the fits.
#' @param costs Positive estimated costs in grid order.
#' @param cores Total worker budget.
#' @param verbose Whether to report each phase with a message.
#' @return Collected fits or error conditions, one per grid row.
#' @noRd
.multilpa_enumerate_fits <- function(arguments, calls, costs, cores,
                                     verbose = FALSE) {
  force(arguments)
  force(calls)
  force(costs)
  force(cores)
  indices <- seq_along(arguments)
  started <- proc.time()[["elapsed"]]
  worker_report <- .latents_worker_reporter()
  on.exit(worker_report$signal(), add = TRUE)
  .latents_pool_progress(verbose, "Preparing %s on %s.",
                         .latents_count(length(indices), "candidate"),
                         .latents_count(min(cores, length(indices)), "worker"))
  candidates <- worker_report$collect(function() .latents_pool_setup(
    lapply(indices, function(i) {
      list(candidate = i, arguments = arguments[[i]], call = calls[[i]])
    }),
    .multilpa_setup_task, costs, cores,
    # Candidates for a separate engine have no setup to distribute.
    distribute = all(vapply(arguments, function(current) !is.null(current$seed),
                            logical(1))) &&
      !all(vapply(arguments, .multilpa_pool_fallback, logical(1)))))
  screen_tasks <- .multilpa_start_tasks(candidates, "screen")
  if (length(screen_tasks)) {
    .latents_pool_progress(verbose, "Screening %s on %s.",
                           .latents_count(length(screen_tasks), "start"),
                           .latents_count(min(cores, length(screen_tasks)), "worker"))
  }
  screened <- worker_report$collect(function() .latents_pool_phase(
    screen_tasks, costs, cores, length(indices),
    .multilpa_execute_start, .multilpa_pool_loss))
  ranked <- lapply(indices, function(i) {
    plan <- candidates[[i]]$value
    .latents_collect_conditions(function() tryCatch({
      if (inherits(plan, "error") || is.null(plan$screen)) NULL else
        .latents_rank_screened(lapply(screened[[i]], `[[`, "value"),
                                plan$n_starts, plan$screen)
    }, error = .lta_candidate_error))
  })
  final_tasks <- .multilpa_start_tasks(candidates, "final", screened, ranked)
  .latents_pool_progress(verbose, "Fitting %s on %s.",
                         .latents_count(length(final_tasks), "start"),
                         .latents_count(min(cores, max(1L, length(final_tasks))), "worker"))
  finished <- worker_report$collect(function() .latents_pool_phase(
    final_tasks, costs, cores, length(indices), .multilpa_execute_start,
    .multilpa_pool_loss))
  on.exit(.latents_pool_progress(verbose, "Done in %s.", .latents_format_duration(
    proc.time()[["elapsed"]] - started)), add = TRUE)
  lapply(indices, function(i) {
    plan <- candidates[[i]]$value
    final <- .latents_collect_conditions(function() tryCatch({
      if (inherits(plan, "error")) stop(plan)
      if (!is.null(plan$arguments)) return(finished[[i]][[1L]]$value)
      if (inherits(ranked[[i]]$value, "error")) stop(ranked[[i]]$value)
      attempts <- lapply(finished[[i]], `[[`, "value")
      found <- if (is.null(plan$screen)) {
        .latents_select_starts(attempts, plan$select_start, interior = plan$interior,
                               tie = plan$tie)
      } else {
        .latents_select_screened(lapply(screened[[i]], `[[`, "value"),
                                  attempts, ranked[[i]]$value$finalists,
                                  ranked[[i]]$value$screened_scores,
                                  plan$select_start, interior = plan$interior,
                                  tie = plan$tie)
      }
      .multilpa_finish_plan(plan, found)
    }, error = function(error) {
      if (inherits(plan, "latents_multilpa_plan")) {
        .multilpa_report_failed_plan(plan, error)
      }
      .lta_candidate_error(error)
    }))
    final$conditions <- c(candidates[[i]]$conditions,
      unlist(lapply(screened[[i]], `[[`, "conditions"), recursive = FALSE),
      ranked[[i]]$conditions,
      unlist(lapply(finished[[i]], `[[`, "conditions"), recursive = FALSE),
      final$conditions)
    final
  })
}
