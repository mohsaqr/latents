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
#' @return A list of settings.
#' @noRd
.latents_em_settings <- function(decrease_tolerance = 1e-10,
                                 decrease_condition = "plain",
                                 guard_decrease = TRUE,
                                 stop_on_nonfinite = FALSE,
                                 accelerate = NULL,
                                 max_stalled = Inf) {
  stopifnot(
    "`decrease_tolerance` must be a single positive number" =
      is.numeric(decrease_tolerance) && length(decrease_tolerance) == 1L &&
      decrease_tolerance > 0,
    "`decrease_condition` must be a single string" =
      is.character(decrease_condition) && length(decrease_condition) == 1L,
    "`accelerate` must be NULL or a list of flatten, restore and valid" =
      is.null(accelerate) ||
      (is.list(accelerate) && all(c("flatten", "restore", "valid") %in%
                                    names(accelerate)))
  )
  list(decrease_tolerance = decrease_tolerance,
       decrease_condition = decrease_condition,
       guard_decrease = isTRUE(guard_decrease),
       stop_on_nonfinite = isTRUE(stop_on_nonfinite),
       accelerate = accelerate, max_stalled = max_stalled)
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
#' and at least three steps of the budget remain. Convergence is the relative
#' gain over one step or cycle, `|gain| <= tol * (1 + |log L|)`; a
#' maximization that reports `attr(state, "solved")` must also have settled
#' its inner solve.
#'
#' @param state Starting parameter state.
#' @param evaluate E-step: a function of a state returning a list with
#'   `log_likelihood`.
#' @param maximize M-step: a function of a state and its E-step returning the
#'   next state.
#' @param max_iter Maximum number of EM steps.
#' @param tol Relative convergence tolerance.
#' @param settings From `.latents_em_settings()`.
#' @return A list with `state`, `expectation`, `converged`, `iterations`,
#'   `history` (the log likelihood at the start and after every step) and
#'   `settled` (whether the last step settled the likelihood).
#' @noRd
.latents_em <- function(state, evaluate, maximize, max_iter, tol,
                        settings = .latents_em_settings()) {
  expectation <- evaluate(state)
  history <- expectation$log_likelihood
  converged <- FALSE
  settled <- FALSE
  stalled <- 0L
  iteration <- 0L
  step <- function(point, point_expectation) {
    .latents_em_step(point, point_expectation, evaluate, maximize, settings)
  }
  # EM is a fixed-point iteration: every step starts from the last point.
  while (iteration < max_iter && !converged && stalled < settings$max_stalled) {
    accelerate <- !is.null(settings$accelerate) && max_iter - iteration >= 3L
    if (accelerate) {
      cycle <- .latents_squarem_cycle(step, state, expectation, evaluate,
                                      settings$accelerate)
      evaluations <- 3L
    } else {
      plain <- step(state, expectation)
      if (!plain$finite) {
        state <- plain$state
        break
      }
      cycle <- list(state = plain$state, expectation = plain$expectation,
                    history = plain$expectation$log_likelihood)
      evaluations <- 1L
    }
    # Over a SQUAREM cycle the gain is at least one plain step's gain, so an
    # accelerated fit never stops earlier than plain EM would.
    gain <- cycle$expectation$log_likelihood - expectation$log_likelihood
    settled <- abs(gain) <= tol * (1 + abs(expectation$log_likelihood))
    solved <- attr(cycle$state, "solved") %||% TRUE
    converged <- settled && isTRUE(solved)
    stalled <- if (settled && !isTRUE(solved)) stalled + 1L else 0L
    iteration <- iteration + evaluations
    history <- c(history, cycle$history)
    state <- cycle$state
    expectation <- cycle$expectation
  }
  list(state = state, expectation = expectation, converged = converged,
       iterations = iteration, history = history, settled = settled)
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
#' @return A list with `state`, `expectation`, `history` (three likelihoods)
#'   and `accelerated`.
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
                            second$expectation$log_likelihood,
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
                     error = function(error) NULL)
  if (is.null(landed) ||
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
#' @return A list with `attempts`, `valid`, `scores`, `converged` and
#'   `best_start`.
#' @noRd
.latents_run_starts <- function(n_starts, run, select_start,
                                score = function(attempt) {
                                  attempt$expectation$log_likelihood
                                }) {
  attempts <- lapply(seq_len(n_starts), function(start_index) {
    tryCatch(run(start_index),
             error = function(error) list(error = conditionMessage(error)))
  })
  valid <- vapply(attempts, function(attempt) is.null(attempt$error), logical(1))
  if (!any(valid)) {
    stop(errorCondition(sprintf(
      "All %d starts failed: %s", n_starts,
      paste(unique(vapply(attempts, `[[`, character(1), "error")), collapse = "; ")),
      class = "latents_all_starts_failed", call = NULL))
  }
  scores <- vapply(attempts, function(attempt) {
    if (is.null(attempt$error)) score(attempt) else -Inf
  }, numeric(1))
  converged <- vapply(attempts, function(attempt) isTRUE(attempt$converged),
                      logical(1))
  list(attempts = attempts, valid = valid, scores = scores,
       converged = converged,
       best_start = .multilpa_select_start(scores, converged, select_start))
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
#' @param tol Relative convergence tolerance for the code-52 check.
#' @param control `optim()` control list.
#' @return A list with `par`, `value` (the negative log likelihood there),
#'   `converged`, `improved` (whether the finish was kept) and `start_value`.
#' @noRd
.latents_quasi_newton <- function(theta, lower, value, gradient, tol, control) {
  theta <- pmax(theta, lower)
  start_value <- value(theta)
  found <- stats::optim(theta, value, gradient, method = "L-BFGS-B",
                        lower = lower, control = control)
  if (!is.finite(found$value) || found$value > start_value) {
    return(list(par = theta, value = start_value, converged = FALSE,
                improved = FALSE, start_value = start_value))
  }
  converged <- found$convergence == 0L
  if (identical(found$convergence, 52L)) {
    again <- stats::optim(found$par, value, gradient, method = "L-BFGS-B",
                          lower = lower, control = control)
    if (is.finite(again$value) && again$value <= found$value) {
      converged <- again$convergence == 0L ||
        found$value - again$value <= tol * (1 + abs(again$value))
      found <- again
    }
  }
  list(par = found$par, value = found$value, converged = converged,
       improved = TRUE, start_value = start_value)
}
