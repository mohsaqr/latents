# The estimation kernel: one EM driver for every model in the package.
#
# Every estimator here maximizes a mixture likelihood by EM from several
# starting points. The engines differ in what they compute -- the E-step
# (posteriors and the log likelihood) and the M-step (the parameters that
# maximize the posterior-weighted complete-data likelihood) -- but not in how
# they iterate: step, re-evaluate, guard against a decrease beyond rounding,
# stop when the relative gain falls below the tolerance, record the history.
# That loop lives here once. An engine hands the driver its two functions
# and its own settings, and the arithmetic of the iteration is the same for
# all of them.
#
# The settings that differ between engines are their contract, not
# accidents, and stay theirs: how much rounding a decrease may show before it
# is a defect (the models with inner Newton solves tolerate more), which
# condition reports it, whether SQUAREM extrapolates, whether a non-finite
# likelihood ends the run, and whether an inner solve must also settle.

#' Settings of one engine's EM iteration
#'
#' @param decrease_tolerance Relative decrease of the log likelihood, beyond
#'   which a step is a defect: a decrease larger than
#'   `decrease_tolerance * (1 + |log L|)` raises the decrease condition.
#' @param decrease_condition `"plain"` raises the unclassed message
#'   "EM likelihood decreased beyond numerical roundoff."; `"classed"` raises
#'   `latents_em_decrease` stating the size of the decrease; a string other
#'   than those is used verbatim as an unclassed message.
#' @param guard_decrease Whether to check for a decrease at all (a fit under a
#'   prior climbs the posterior, so its likelihood may fall).
#' @param stop_on_nonfinite Whether a non-finite likelihood after a step ends
#'   the run, keeping the step's parameters but the previous E-step.
#' @param accelerate `NULL` for plain EM, or a list of `flatten`, `restore` and
#'   `valid` functions over the parameter list for SQUAREM.
#' @param max_stalled With a maximization that reports `solved`, how many
#'   steps in a row may settle the likelihood without settling the inner solve
#'   before the run stops unconverged; `Inf` when the maximization has no inner
#'   solve.
#' @param stop_scale The number of observed values the likelihood sums over
#'   (rows of positive weight times their observed indicators), or `NULL`.
#'   With it, a step or cycle has settled when its gain is at most
#'   `tol * max(1, stop_scale)`: a change in log likelihood does not depend on
#'   the units of the data, so neither does this rule. Without it, the older
#'   rule `tol * (1 + |log L|)` applies, which a change of units moves, because
#'   rescaling the data by `c` shifts `log L` by `-n d log c` (review
#'   2026-10-26 round 2, K2-6).
#' @return A list of settings.
#' @noRd
.latents_em_settings <- function(decrease_tolerance = 1e-10,
                                 decrease_condition = "plain",
                                 guard_decrease = TRUE,
                                 stop_on_nonfinite = FALSE,
                                 accelerate = NULL,
                                 max_stalled = Inf,
                                 stop_scale = NULL) {
  stopifnot(
    "`decrease_tolerance` must be a single positive number" =
      is.numeric(decrease_tolerance) && length(decrease_tolerance) == 1L &&
      decrease_tolerance > 0,
    "`decrease_condition` must be a single string" =
      is.character(decrease_condition) && length(decrease_condition) == 1L,
    "`accelerate` must be NULL or a list of flatten, restore and valid" =
      is.null(accelerate) ||
      (is.list(accelerate) && all(c("flatten", "restore", "valid") %in%
                                    names(accelerate))),
    "`stop_scale` must be NULL or a single non-negative number" =
      is.null(stop_scale) || (is.numeric(stop_scale) && length(stop_scale) == 1L &&
                                is.finite(stop_scale) && stop_scale >= 0)
  )
  list(decrease_tolerance = decrease_tolerance,
       decrease_condition = decrease_condition,
       guard_decrease = isTRUE(guard_decrease),
       stop_on_nonfinite = isTRUE(stop_on_nonfinite),
       accelerate = accelerate, max_stalled = max_stalled,
       stop_scale = stop_scale)
}

#' The number of observed values a likelihood sums over
#'
#' The scale of the unit-free EM stopping rule (`stop_scale` of
#' `.latents_em_settings()`): every non-missing cell of the indicator blocks,
#' in rows of positive sampling weight.
#'
#' @param ... Matrices with one row per row of data (continuous indicators,
#'   categorical codes, ordinal and count data); `NULL` entries are skipped.
#' @param row_weights `NULL`, or one sampling weight per row.
#' @return A single non-negative number.
#' @noRd
.latents_observed_values <- function(..., row_weights = NULL) {
  blocks <- Filter(Negate(is.null), list(...))
  if (length(blocks) == 0L) return(0)
  counts <- Reduce(`+`, lapply(blocks, function(block) {
    rowSums(!is.na(as.matrix(block)))
  }))
  if (!is.null(row_weights)) counts <- counts[which(row_weights > 0)]
  as.numeric(sum(counts))
}

#' How close to the best log likelihood a start must end to replicate it
#'
#' `1e-6` per observed value. A change of units shifts every start's log
#' likelihood by the same amount, so a tolerance on their gaps must not scale
#' with `|log L|`, as `1e-6 * (1 + |log L|)` did; one per observed value is
#' unit-free and about as wide for data in ordinary units. It is wider than
#' the EM stopping step, because a start stopped by that rule can still lie
#' several steps below a maximum it approaches slowly.
#'
#' @param n_observed The observed values the likelihood sums over (the
#'   `stop_scale` of the family's EM).
#' @return A single positive number.
#' @noRd
.latents_replication_tolerance <- function(n_observed) {
  1e-6 * max(1, n_observed)
}

#' The likelihood change below which an EM step has settled
#' @param settings From `.latents_em_settings()`.
#' @param tol Relative convergence tolerance.
#' @param log_likelihood The log likelihood the step started from, used only
#'   by the older rule when the engine supplied no `stop_scale`.
#' @return A single non-negative number.
#' @noRd
.latents_em_threshold <- function(settings, tol, log_likelihood) {
  if (is.null(settings$stop_scale)) return(tol * (1 + abs(log_likelihood)))
  tol * max(1, settings$stop_scale)
}

#' Raise an engine's decrease condition
#' @param settings From `.latents_em_settings()`.
#' @param decrease The size of the decrease, positive.
#' @noRd
.latents_em_decrease <- function(settings, decrease) {
  switch(settings$decrease_condition,
    plain = stop("EM likelihood decreased beyond numerical roundoff."),
    classed = stop(errorCondition(sprintf(paste(
      "EM likelihood decreased by %.3g, beyond numerical roundoff. This is a",
      "defect; please report it with the call that produced it."),
      decrease), class = "latents_em_decrease", call = NULL)),
    stop(settings$decrease_condition))
}

#' One plain EM step with the engine's decrease guard
#'
#' @param state The current parameter state.
#' @param expectation Its E-step, carrying `log_likelihood`.
#' @param evaluate E-step: a function of a state.
#' @param maximize M-step: a function of a state and its E-step, returning
#'   the next state.
#' @param settings From `.latents_em_settings()`.
#' @return A list with the next `state` and its `expectation`.
#' @noRd
.latents_em_step <- function(state, expectation, evaluate, maximize, settings) {
  updated_state <- maximize(state, expectation)
  updated <- evaluate(updated_state)
  if (settings$stop_on_nonfinite && !is.finite(updated$log_likelihood)) {
    return(list(state = updated_state, expectation = updated, finite = FALSE))
  }
  gain <- updated$log_likelihood - expectation$log_likelihood
  if (settings$guard_decrease &&
      gain < -settings$decrease_tolerance * (1 + abs(expectation$log_likelihood))) {
    .latents_em_decrease(settings, -gain)
  }
  list(state = updated_state, expectation = updated, finite = TRUE)
}

#' Run EM from one starting state
#'
#' Plain EM steps, or SQUAREM cycles when the settings carry an accelerator
#' and at least three steps of the budget remain. Convergence is a small gain
#' over one step or cycle, `|gain| <= tol * max(1, stop_scale)` (or the older
#' `tol * (1 + |log L|)` for an engine without `stop_scale`; see
#' `.latents_em_settings()`); a maximization that reports
#' `attr(state, "solved")` must also have settled its inner solve. The
#' iteration count is the number of M-steps whose result the run kept: two
#' for a SQUAREM cycle that falls back to plain EM, three for one that
#' extrapolates.
#'
#' @param state Starting parameter state.
#' @param evaluate E-step: a function of a state returning a list with
#'   `log_likelihood`.
#' @param maximize M-step: a function of a state and its E-step returning the
#'   next state.
#' @param max_iter Maximum number of EM steps.
#' @param tol Relative convergence tolerance.
#' @param settings From `.latents_em_settings()`.
#' @param pause Iteration count after which the run stops at the end of the
#'   current step or cycle, resumable; `Inf` runs to the end. Whether a cycle
#'   is accelerated is decided against `max_iter`, not `pause`, so a run
#'   paused and resumed takes exactly the steps of one that was not.
#' @param resume `NULL`, or the `checkpoint` of a paused run, which the run
#'   continues (its E-step is recomputed from the state, which is what the run
#'   computed there, except after a non-finite step, when it is carried).
#' @return A list with `state`, `expectation`, `converged`, `iterations`,
#'   `history` (the log likelihood at the start and after every step) and
#'   `settled` (whether the last step settled the likelihood); with a finite
#'   `pause`, also `checkpoint`, everything a resumed run needs.
#' @noRd
.latents_em <- function(state, evaluate, maximize, max_iter, tol,
                        settings = .latents_em_settings(), pause = Inf,
                        resume = NULL) {
  if (is.null(resume)) {
    expectation <- evaluate(state)
    history <- expectation$log_likelihood
    converged <- FALSE
    settled <- FALSE
    stalled <- 0L
    iteration <- 0L
    ended <- FALSE
  } else {
    state <- resume$state
    expectation <- resume$expectation %||% evaluate(state)
    history <- resume$history
    converged <- resume$converged
    settled <- resume$settled
    stalled <- resume$stalled
    iteration <- resume$iterations
    ended <- resume$ended
  }
  step <- function(point, point_expectation) {
    .latents_em_step(point, point_expectation, evaluate, maximize, settings)
  }
  # EM is a fixed-point iteration: every step starts from the last point.
  while (!ended && iteration < max_iter && iteration < pause && !converged &&
         stalled < settings$max_stalled) {
    accelerate <- !is.null(settings$accelerate) && max_iter - iteration >= 3L
    if (accelerate) {
      cycle <- .latents_squarem_cycle(step, state, expectation, evaluate,
                                      settings$accelerate)
      # The M-steps whose result the path kept (review 2026-10-26, K2-7).
      evaluations <- length(cycle$history)
    } else {
      plain <- step(state, expectation)
      if (!plain$finite) {
        state <- plain$state
        ended <- TRUE
        break
      }
      cycle <- list(state = plain$state, expectation = plain$expectation,
                    history = plain$expectation$log_likelihood)
      evaluations <- 1L
    }
    # Over a SQUAREM cycle the gain is at least one plain step's gain, so an
    # accelerated fit never stops earlier than plain EM would.
    gain <- cycle$expectation$log_likelihood - expectation$log_likelihood
    settled <- abs(gain) <= .latents_em_threshold(settings, tol,
                                                  expectation$log_likelihood)
    solved <- attr(cycle$state, "solved") %||% TRUE
    converged <- settled && isTRUE(solved)
    stalled <- if (settled && !isTRUE(solved)) stalled + 1L else 0L
    iteration <- iteration + evaluations
    history <- c(history, cycle$history)
    state <- cycle$state
    expectation <- cycle$expectation
  }
  fit <- list(state = state, expectation = expectation, converged = converged,
              iterations = iteration, history = history, settled = settled)
  if (is.finite(pause) || !is.null(resume)) {
    fit$checkpoint <- list(state = state, history = history, converged = converged,
                           settled = settled, stalled = stalled,
                           iterations = iteration, ended = ended)
    # A run ended by a non-finite step keeps the E-step before it, which
    # evaluating its state again would not give.
    if (ended) fit$checkpoint$expectation <- expectation
  }
  fit
}

#' One SQUAREM cycle over an engine's parameter state
#'
#' Varadhan, R., & Roland, C. (2008). Simple and globally convergent methods
#' for accelerating the convergence of any EM algorithm. Scandinavian Journal
#' of Statistics, 35, 335-353.
#'
#' Two EM steps from theta0 give the first and second differences
#' r = theta1 - theta0 and v = theta2 - theta1 - r; the extrapolation
#' theta' = theta0 - 2 alpha r + alpha^2 v uses alpha = -|r| / |v| (scheme
#' S3), backtracked towards -1 (plain EM) until the state is valid, and one
#' EM step from theta' returns the iterate to the constraint set.
#'
#' @param step A function of a state and its E-step giving the next
#'   `list(state, expectation)`.
#' @param state Current state.
#' @param expectation Its E-step.
#' @param evaluate E-step function.
#' @param accelerator A list of `flatten`, `restore` and `valid`.
#' @return A list with `state`, `expectation`, `history` (the likelihood after
#'   each M-step the cycle kept: two when it falls back to plain EM, three when
#'   the extrapolated step is taken) and `accelerated`.
#' @noRd
.latents_squarem_cycle <- function(step, state, expectation, evaluate,
                                   accelerator) {
  first <- step(state, expectation)
  second <- step(first$state, first$expectation)
  theta0 <- accelerator$flatten(state)
  r <- accelerator$flatten(first$state) - theta0
  v <- accelerator$flatten(second$state) - theta0 - 2 * r
  plain <- list(state = second$state, expectation = second$expectation,
                history = c(first$expectation$log_likelihood,
                            second$expectation$log_likelihood),
                accelerated = FALSE)
  norm_v <- sqrt(sum(v^2))
  if (!is.finite(norm_v) || norm_v <= 0) return(plain)
  alpha <- min(-1, -sqrt(sum(r^2)) / norm_v)
  # Halving the distance to alpha = -1, which is always valid, reaches a
  # valid point in a bounded number of tries.
  candidate <- NULL
  tries <- 0L
  while (is.null(candidate) && alpha < -1 && tries < 30L) {
    theta <- theta0 - 2 * alpha * r + alpha^2 * v
    proposal <- accelerator$restore(state, theta)
    if (accelerator$valid(proposal)) candidate <- proposal else
      alpha <- (alpha - 1) / 2
    tries <- tries + 1L
  }
  if (is.null(candidate)) return(plain)
  # The E-step at the extrapolated point can still fail numerically (every
  # component's density underflows); that point is then simply not taken.
  landed <- tryCatch(step(candidate, evaluate(candidate)),
                     error = function(error) error)
  if (inherits(landed, "error") ||
      !is.finite(landed$expectation$log_likelihood) ||
      landed$expectation$log_likelihood < second$expectation$log_likelihood) {
    return(plain)
  }
  list(state = landed$state, expectation = landed$expectation,
       history = c(first$expectation$log_likelihood,
                   second$expectation$log_likelihood,
                   landed$expectation$log_likelihood),
       accelerated = TRUE)
}

#' Run every start and choose the best
#'
#' Each start runs to completion or fails; a failure is recorded with its
#' message (it appears in the starts table) rather than ending the fit, and
#' only when every start fails is the fit refused, naming the distinct
#' reasons. The best start is the highest likelihood under the selection
#' rule, with ties resolved as `.multilpa_select_start()` resolves them.
#'
#' @param n_starts Number of starts.
#' @param run A function of the start index returning a list with
#'   `expectation$log_likelihood` and `converged`.
#' @param select_start `"likelihood"` or `"converged"`.
#' @param score A function of a successful attempt giving its log
#'   likelihood; the default reads `attempt$expectation$log_likelihood`.
#' @param interior `NULL`, or a function of a successful attempt returning
#'   `TRUE` when no variance is at its lower bound; when supplied, starts at
#'   the bound are excluded from the choice (see `.multilpa_select_start()`).
#' @param prepare `NULL`, or a function of the start index returning that
#'   start's initial values; `run` is then a function of the start index and
#'   those values. Every start is prepared first, in order and on this
#'   process, so the random numbers are drawn exactly as when each start is
#'   prepared and run in turn (the iterations draw none); only the
#'   iterations are spread over `cores`.
#' @param cores Number of cores for the iterations, from
#'   `.latents_check_cores()`; used with `prepare` or supplied `initial` values.
#' @param screen `NULL` for the search of latents 0.96.1 and earlier, every
#'   start run to the end; or, for the emEM search, a list of `n_pool`,
#'   `iterations` and `protect` (see `.latents_run_screened()`), which needs
#'   `prepare` or already prepared `initial` values.
#' @param initial `NULL`, or ordered wrapped values from
#'   `.latents_prepare_starts()`, prepared in the caller's seed scope.
#' @param tie `NULL`, or the fit's EM stopping threshold (`tol` per observed
#'   value): starts closer than this are tied (see `.multilpa_select_start()`).
#' @return A list with `attempts`, `valid`, `scores`, `converged` and
#'   `best_start`; with `screen`, also `stage` and `screened_scores`.
#' @noRd
.latents_run_starts <- function(n_starts, run, select_start,
                                score = function(attempt) {
                                  attempt$expectation$log_likelihood
                                }, interior = NULL, prepare = NULL,
                                cores = 1L, screen = NULL, initial = NULL,
                                tie = NULL) {
  if (!is.null(screen)) {
    return(.latents_run_screened(n_starts, run, select_start, score, interior,
                                 prepare, cores, screen, initial, tie))
  }
  force(run)
  attempts <- if (is.null(prepare) && is.null(initial)) {
    lapply(seq_len(n_starts), function(start_index) {
      tryCatch(run(start_index),
               error = function(error) list(error = conditionMessage(error)))
    })
  } else {
    if (is.null(initial)) initial <- .latents_prepare_starts(n_starts, prepare)
    .latents_map_candidates(seq_len(n_starts), function(start_index) {
      .latents_start_attempt(run, start_index, initial[[start_index]])
    }, cores)
  }
  .latents_select_starts(attempts, select_start, score, interior, tie)
}

#' Prepare starting values in their original random-number order
#' @param n_starts Number of starts to prepare.
#' @param prepare Function of a start index returning its initial values.
#' @return One wrapped `values` or `error` entry per start, in start order.
#' @noRd
.latents_prepare_starts <- function(n_starts, prepare) {
  stopifnot("`prepare` must be a function" = is.function(prepare))
  lapply(seq_len(n_starts), function(start_index) {
    tryCatch(list(values = prepare(start_index)),
             error = function(error) list(error = conditionMessage(error)))
  })
}

#' Run one prepared start or continue its checkpoint
#' @param run Function of a start index and initial values; screening also
#'   requires `pause` and `resume` arguments.
#' @param start_index The start's original index.
#' @param initial Wrapped initial values from `.latents_prepare_starts()`;
#'   `NULL` when continuing a checkpoint.
#' @param pause Screening iteration budget, or `Inf` for a complete run.
#' @param resume Checkpoint to continue, or `NULL` for a new start.
#' @return The attempt, or an `error` entry for an ordinary start failure.
#' @noRd
.latents_start_attempt <- function(run, start_index, initial, pause = Inf,
                                   resume = NULL) {
  stopifnot("`run` must be a function" = is.function(run))
  if (!is.null(initial$error)) return(initial)
  tryCatch({
    if (!is.null(resume)) {
      run(start_index, NULL, resume = resume)
    } else if (is.finite(pause)) {
      run(start_index, initial$values, pause = pause)
    } else {
      # Ordinary prepared runners need only accept the original two arguments.
      run(start_index, initial$values)
    }
  }, error = function(error) {
    if (inherits(error, "latents_parallel_failure")) stop(error)
    list(error = conditionMessage(error))
  })
}

#' Retain accumulated covariance qualifications when every start fails
#' @param message The unchanged all-starts failure message.
#' @param attempts Start attempts, including every completed search phase.
#' @return A `latents_all_starts_failed` error condition, with structure logs
#'   only when an attempt recorded them.
#' @noRd
.latents_start_failure <- function(message, attempts) {
  condition <- errorCondition(message, class = "latents_all_starts_failed",
                               call = NULL)
  logs <- lapply(attempts, attr, which = "lta_structure_log")
  if (any(!vapply(logs, is.null, logical(1)))) {
    attr(condition, "lta_structure_logs") <- logs
  }
  condition
}

#' Select among complete start attempts
#' @param attempts Complete attempts in their original start order.
#' @param select_start,score,interior,tie As for `.latents_run_starts()`.
#' @return As `.latents_run_starts()` without screening fields.
#' @noRd
.latents_select_starts <- function(attempts, select_start,
                                   score = function(attempt) {
                                     attempt$expectation$log_likelihood
                                   }, interior = NULL, tie = NULL) {
  valid <- vapply(attempts, function(attempt) is.null(attempt$error), logical(1))
  if (!any(valid)) {
    stop(.latents_start_failure(sprintf(
      "All %d starts failed: %s", length(attempts),
      paste(unique(vapply(attempts, `[[`, character(1), "error")), collapse = "; ")),
      attempts))
  }
  scores <- vapply(attempts, function(attempt) {
    if (is.null(attempt$error)) score(attempt) else -Inf
  }, numeric(1))
  converged <- vapply(attempts, function(attempt) isTRUE(attempt$converged),
                      logical(1))
  flags <- if (!is.null(interior)) {
    vapply(attempts, function(attempt) is.null(attempt$error) && isTRUE(interior(attempt)),
           logical(1))
  }
  list(attempts = attempts, valid = valid, scores = scores,
       converged = converged,
       best_start = .multilpa_select_start(scores, converged, select_start, flags,
                                           tie = tie))
}

#' Screen a pool of starts and run the best to the end (emEM)
#'
#' Biernacki, C., Celeux, G., & Govaert, G. (2003). Choosing starting values
#' for the EM algorithm for getting the highest likelihood in multivariate
#' Gaussian mixture models. Computational Statistics & Data Analysis, 41,
#' 561-575.
#'
#' Every start of the pool is prepared in order on this process (so the
#' random numbers are drawn as when they are prepared one by one), then run
#' for `screen$iterations` EM iterations, paused at the end of the step or
#' cycle that reaches them. The protected starts and then those with the
#' highest screened likelihood (ties by pool order) are resumed to the end,
#' `n_starts` in all; a resumed start takes exactly the steps it would have
#' taken unscreened. The best is chosen among them by
#' `.multilpa_select_start()`; the other starts stay in the record with their
#' screened likelihood.
#'
#' @param n_starts Number of starts run to the end.
#' @param run A function of the start index, its initial values, `pause` and
#'   `resume`, returning an attempt that carries a `checkpoint` when paused.
#' @param select_start,score,interior,prepare,cores,initial,tie As for
#'   `.latents_run_starts()`.
#' @param screen A list of `n_pool` (starts screened, at least `n_starts`),
#'   `iterations` (the screening budget) and `protect` (pool indices that go
#'   to the end whatever their screened likelihood, such as a supplied start).
#' @return As `.latents_run_starts()`, with `attempts` over the whole pool, a
#'   screened-only attempt carrying `screened = TRUE`, its `state`,
#'   `converged`, `iterations` and `expectation$log_likelihood`; `stage`
#'   (`"final"` or `"screened"`) and `screened_scores`.
#' @noRd
.latents_run_screened <- function(n_starts, run, select_start, score, interior,
                                  prepare, cores, screen, initial = NULL,
                                  tie = NULL) {
  stopifnot(
    "screening needs `prepare` or `initial`" =
      is.function(prepare) || is.list(initial),
    "`screen$n_pool` must be at least `n_starts`" = screen$n_pool >= n_starts
  )
  force(run)
  force(score)
  force(screen)
  worker_report <- .latents_worker_reporter()
  on.exit(worker_report$signal(), add = TRUE)
  pool <- seq_len(screen$n_pool)
  if (is.null(initial)) {
    initial <- .latents_prepare_starts(screen$n_pool, prepare)
  }
  screened <- worker_report$collect(function() {
    .latents_map_candidates(pool, function(start_index) {
      .latents_screen_attempt(run, start_index, initial[[start_index]], screen,
                              score)
    }, cores)
  })
  ranked <- .latents_rank_screened(screened, n_starts, screen)
  finalists <- ranked$finalists
  continued <- worker_report$collect(function() {
    .latents_map_candidates(finalists, function(start_index) {
      .latents_start_attempt(run, start_index, NULL,
                             resume = screened[[start_index]]$checkpoint)
    }, cores)
  })
  .latents_select_screened(screened, continued, finalists, ranked$screened_scores,
                            select_start, score, interior, tie)
}

#' Screen one prepared start and retain only its resumable state
#' @param run,start_index,initial As for `.latents_start_attempt()`.
#' @param screen Screening controls from `.latents_run_screened()`.
#' @param score Function giving the likelihood of a successful attempt.
#' @return A reduced screened attempt, or an `error` entry.
#' @noRd
.latents_screen_attempt <- function(run, start_index, initial, screen,
                                    score = function(attempt) {
                                      attempt$expectation$log_likelihood
                                    }) {
  attempt <- .latents_start_attempt(run, start_index, initial,
                                    pause = screen$iterations)
  if (!is.null(attempt$error)) return(attempt)
  # Keep only what ranking and resuming need, rather than every posterior
  # in a large screening pool.
  result <- list(screened = TRUE, state = attempt$checkpoint$state,
       converged = isTRUE(attempt$converged), iterations = attempt$iterations,
       expectation = list(log_likelihood = score(attempt)),
       checkpoint = attempt$checkpoint)
  attr(result, "lta_structure_log") <- attr(attempt, "lta_structure_log")
  result
}

#' Rank screened starts while retaining protected starts
#' @param screened Reduced attempts from `.latents_screen_attempt()`.
#' @param n_starts Number of finalists.
#' @param screen Screening controls from `.latents_run_screened()`.
#' @param score Function reading a reduced attempt's screened likelihood.
#' @return The ordered `finalists` and all `screened_scores`.
#' @noRd
.latents_rank_screened <- function(screened, n_starts, screen,
                                   score = function(attempt) {
                                     attempt$expectation$log_likelihood
                                   }) {
  pool <- seq_along(screened)
  screened_valid <- vapply(screened, function(attempt) is.null(attempt$error),
                           logical(1))
  screened_scores <- vapply(screened, function(attempt) {
    if (is.null(attempt$error)) score(attempt) else -Inf
  }, numeric(1))
  usable <- screened_valid & is.finite(screened_scores)
  if (!any(usable)) {
    stop(.latents_start_failure(sprintf(
      "All %d starts failed: %s", screen$n_pool,
      paste(unique(vapply(screened[!screened_valid], `[[`, character(1), "error")),
            collapse = "; ")),
      screened))
  }
  protected <- intersect(screen$protect %||% integer(), pool[usable])
  ranked <- setdiff(pool[usable][order(-screened_scores[usable], pool[usable])],
                    protected)
  finalists <- utils::head(c(protected, ranked), n_starts)
  list(finalists = finalists, screened_scores = screened_scores)
}

#' Merge finalist continuations and select the best complete start
#' @param screened Reduced attempts in original pool order.
#' @param continued Complete attempts in `finalists` order.
#' @param finalists Indices of starts carried to completion.
#' @param screened_scores Screened likelihoods in original pool order.
#' @param select_start,score,interior,tie As for `.latents_run_starts()`.
#' @return As `.latents_run_starts()` with screening fields.
#' @noRd
.latents_select_screened <- function(screened, continued, finalists,
                                     screened_scores, select_start,
                                     score = function(attempt) {
                                       attempt$expectation$log_likelihood
                                     }, interior = NULL, tie = NULL) {
  pool <- seq_along(screened)
  attempts <- lapply(screened, function(attempt) {
    attempt$checkpoint <- NULL
    attempt
  })
  attempts[finalists] <- lapply(continued, function(attempt) {
    attempt$checkpoint <- NULL
    attempt
  })
  stage <- ifelse(pool %in% finalists, "final", "screened")
  valid <- vapply(attempts, function(attempt) is.null(attempt$error), logical(1))
  final_valid <- valid & stage == "final"
  if (!any(final_valid)) {
    stop(.latents_start_failure(sprintf(
      "All %d starts carried to the end failed: %s", length(finalists),
      paste(unique(vapply(attempts[finalists], `[[`, character(1), "error")),
            collapse = "; ")),
      c(screened, continued)))
  }
  scores <- vapply(seq_along(attempts), function(start_index) {
    attempt <- attempts[[start_index]]
    if (!valid[start_index]) -Inf else if (final_valid[start_index]) score(attempt) else
      attempt$expectation$log_likelihood
  }, numeric(1))
  converged <- vapply(attempts, function(attempt) isTRUE(attempt$converged),
                      logical(1))
  # Only a start run to the end can be reported.
  eligible_scores <- ifelse(final_valid, scores, -Inf)
  flags <- if (!is.null(interior)) {
    vapply(seq_along(attempts), function(start_index) {
      final_valid[start_index] && isTRUE(interior(attempts[[start_index]]))
    }, logical(1))
  }
  result <- list(attempts = attempts, valid = valid, scores = scores,
                 converged = converged, stage = stage, screened_scores = screened_scores,
                 best_start = .multilpa_select_start(eligible_scores, converged & final_valid,
                                                     select_start, flags, tie = tie))
  # A finalist replaces its screened attempt, but both phases contributed
  # covariance M-steps to the fit's qualification. Keep those counts separately.
  logs <- lapply(c(screened, continued), attr, which = "lta_structure_log")
  if (any(!vapply(logs, is.null, logical(1)))) attr(result, "lta_structure_logs") <- logs
  result
}

#' Quasi-Newton finish from an EM point
#'
#' EM crawls near a maximum when classes overlap; from its point L-BFGS-B
#' maximizes the log likelihood with the analytic gradient, inside lower
#' bounds on the coordinates that have them. The point is kept only if the
#' likelihood does not fall. Convergence is L-BFGS-B's own, except that its
#' line-search stop (code 52), which rounding can trigger at the maximum on
#' some platforms, is checked by one restart from where it stopped: if the
#' likelihood cannot be raised by more than the tolerance, the relative-change
#' rule EM uses is met.
#'
#' @param theta Packed starting coordinates.
#' @param lower Lower bounds, `-Inf` where there is none; `theta` is raised to
#'   them first.
#' @param value Negative log likelihood as a function of the coordinates.
#' @param gradient Its gradient.
#' @param tol Convergence tolerance for the code-52 check, per observed value
#'   as in the EM stopping rule.
#' @param control `optim()` control list.
#' @param scale The number of observed values the likelihood sums over (the
#'   family's `stop_scale`).
#' @param parscale `NULL`, or `optim()`'s `parscale`: one typical magnitude per
#'   coordinate, the reference standard deviation for a coordinate in the
#'   units of the data.
#' @return A list with `par`, `value` (the negative log likelihood there),
#'   `converged`, `improved` (whether the finish was kept) and `start_value`.
#' @noRd
.latents_quasi_newton <- function(theta, lower, value, gradient, tol, control,
                                  scale, parscale = NULL) {
  theta <- pmax(theta, lower)
  start_value <- value(theta)
  # L-BFGS-B stops on a change relative to |f|, and a change of units moves f
  # by a constant (the density Jacobian). Minimizing the change from the start
  # plus `scale` measures the stop against the observed values instead, which
  # is unit-free and about |log L| for data in ordinary units.
  shifted <- function(v) value(v) - start_value + scale
  if (!is.null(parscale)) control$parscale <- parscale
  search <- function(from) {
    found <- stats::optim(from, shifted, gradient, method = "L-BFGS-B",
                          lower = lower, control = control)
    # Evaluated again, not shifted back: (f - f0 + s) + f0 - s can round one
    # ulp above f0 at a point that did not move, which would read as a loss.
    found$value <- value(found$par)
    found
  }
  found <- search(theta)
  if (!is.finite(found$value) || found$value > start_value) {
    return(list(par = theta, value = start_value, converged = FALSE,
                improved = FALSE, start_value = start_value))
  }
  converged <- found$convergence == 0L
  if (identical(found$convergence, 52L)) {
    again <- search(found$par)
    if (is.finite(again$value) && again$value <= found$value) {
      converged <- again$convergence == 0L ||
        found$value - again$value <= tol * max(1, scale)
      found <- again
    }
  }
  list(par = found$par, value = found$value, converged = converged,
       improved = TRUE, start_value = start_value)
}
