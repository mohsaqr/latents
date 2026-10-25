test_that("preparing starts preserves draw order and ordinary runner contracts", {
  set.seed(931)
  draws <- stats::rnorm(4L)
  set.seed(931)
  prepared <- .latents_prepare_starts(4L, function(index) {
    value <- stats::rnorm(1L)
    if (index == 2L) stop("failed preparation")
    value
  })
  expect_identical(prepared[[2L]], list(error = "failed preparation"))
  expect_identical(vapply(prepared[c(1L, 3L, 4L)], `[[`, numeric(1), "values"),
                   draws[c(1L, 3L, 4L)])
  run <- function(index, initial) {
    list(expectation = list(log_likelihood = initial), converged = TRUE)
  }
  attempts <- lapply(seq_along(prepared), function(index) {
    .latents_start_attempt(run, index, prepared[[index]])
  })
  found <- .latents_select_starts(attempts, "likelihood")
  set.seed(931)
  direct <- .latents_run_starts(4L, run, "likelihood", prepare = function(index) {
    value <- stats::rnorm(1L)
    if (index == 2L) stop("failed preparation")
    value
  })
  expect_identical(found, direct)
  expect_identical(found, .latents_run_starts(4L, run, "likelihood",
                                             initial = prepared))
  expect_identical(found$valid, c(TRUE, FALSE, TRUE, TRUE))
})

test_that("start attempts retain ordinary errors and rethrow pool failures", {
  failed <- .latents_start_attempt(function(index, initial) stop("ordinary"),
                                    1L, list(values = 0))
  expect_identical(failed, list(error = "ordinary"))
  expect_error(.latents_start_attempt(function(index, initial) {
    stop(errorCondition("worker lost", class = "latents_parallel_failure"))
  }, 1L, list(values = 0)), class = "latents_parallel_failure")
  expect_error(.latents_select_starts(list(failed), "likelihood"),
               class = "latents_all_starts_failed")
})

test_that("separate screening and selection preserve protected starts and raw scoring", {
  run <- function(index, initial, pause = Inf, resume = NULL) {
    if (index == 2L) stop("failed run")
    final <- !is.null(resume)
    list(raw_score = if (final) -index else c(-9, -8, -3, -1)[index],
         converged = final, iterations = if (final) 20L else pause,
         checkpoint = list(state = index), state = index)
  }
  score <- function(attempt) attempt$raw_score
  screen <- list(n_pool = 4L, iterations = 3L, protect = 1L)
  initial <- .latents_prepare_starts(4L, identity)
  screened <- lapply(seq_along(initial), function(index) {
    .latents_screen_attempt(run, index, initial[[index]], screen, score)
  })
  ranked <- .latents_rank_screened(screened, 2L, screen)
  expect_identical(ranked$finalists, c(1L, 4L))
  expect_identical(ranked$screened_scores, c(-9, -Inf, -3, -1))
  continued <- lapply(ranked$finalists, function(index) {
    .latents_start_attempt(run, index, NULL, resume = screened[[index]]$checkpoint)
  })
  found <- .latents_select_screened(screened, continued, ranked$finalists,
                                     ranked$screened_scores, "likelihood", score)
  direct <- .latents_run_starts(2L, run, "likelihood", score = score,
                                 prepare = identity, screen = screen)
  expect_identical(found, direct)
  expect_identical(found, .latents_run_starts(2L, run, "likelihood", score = score,
                                             initial = initial, screen = screen))
  expect_identical(found$best_start, 1L)
  expect_true(all(vapply(found$attempts, function(attempt) {
    is.null(attempt$checkpoint)
  }, logical(1))))
  expect_error(.latents_rank_screened(list(list(error = "failed")), 1L,
                                      list(n_pool = 1L)),
               class = "latents_all_starts_failed")
  expect_error(.latents_select_screened(screened, list(list(error = "failed")),
                                        1L, ranked$screened_scores, "likelihood"),
               class = "latents_all_starts_failed")
})

test_that("prepared start worker closures survive socket serialization replay", {
  replay <- function(indices, fun, cores, ...) {
    worker <- unserialize(serialize(fun, NULL))
    lapply(indices, worker)
  }
  local_mocked_bindings(.latents_map_candidates = replay)
  fit <- function(caller_values) {
    .latents_run_starts(3L, function(index, initial, pause = Inf, resume = NULL) {
      current <- if (is.null(resume)) initial else resume$state
      list(expectation = list(log_likelihood = current),
           converged = TRUE, iterations = if (is.finite(pause)) pause else 4L,
           checkpoint = list(state = current))
    }, "likelihood", prepare = function(index) caller_values[index], cores = 2L)
  }
  found <- fit(c(1, 3, 2))
  expect_identical(found$scores, c(1, 3, 2))
  expect_identical(found$best_start, 2L)
  screened <- .latents_run_starts(2L,
    function(index, initial, pause = Inf, resume = NULL) {
      current <- if (is.null(resume)) initial else resume$state
      list(expectation = list(log_likelihood = current),
           converged = TRUE, iterations = if (is.finite(pause)) pause else 4L,
           checkpoint = list(state = current))
    }, "likelihood", prepare = function(index) index, cores = 2L,
    screen = list(n_pool = 4L, iterations = 2L, protect = 1L))
  expect_identical(screened$best_start, 4L)
})

test_that("all-start failure conditions retain every phase's covariance counts", {
  failed <- structure(list(error = "failed"), lta_structure_log = list(capped = 1L))
  screened <- structure(list(expectation = list(log_likelihood = -1)),
                          lta_structure_log = list(capped = 1L))
  conditions <- list(
    tryCatch(.latents_select_starts(rep(list(failed), 2L), "likelihood"),
              error = identity),
    tryCatch(.latents_rank_screened(rep(list(failed), 3L), 1L,
                                    list(n_pool = 3L)), error = identity),
    tryCatch(.latents_select_screened(rep(list(screened), 3L),
                                      rep(list(failed), 2L), c(1L, 3L),
                                      rep(-1, 3L), "likelihood"), error = identity))
  expect_identical(vapply(conditions, conditionMessage, character(1)),
                   c("All 2 starts failed: failed", "All 3 starts failed: failed",
                     "All 2 starts carried to the end failed: failed"))
  expect_true(all(vapply(conditions, inherits, logical(1),
                         "latents_all_starts_failed")))
  expect_identical(vapply(conditions, function(condition) {
    sum(vapply(attr(condition, "lta_structure_logs"), `[[`, integer(1), "capped"))
  }, integer(1)), c(2L, 3L, 5L))
  generic <- .latents_start_failure("ordinary failure", list(list(error = "failed")))
  expect_identical(generic, errorCondition("ordinary failure",
                                            class = "latents_all_starts_failed",
                                            call = NULL))
})
