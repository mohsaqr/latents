test_that("worker-loss phase reports combine once without hiding task warnings", {
  exercise <- function(fail = FALSE) {
    report <- .latents_worker_reporter()
    on.exit(report$signal(), add = TRUE)
    report$collect(function() {
      .latents_warn_worker_lost(2L, 1L, 2L)
      warning(warningCondition("task qualification", class = "report_test_warning"))
    })
    report$collect(function() .latents_warn_worker_lost(1L, 1L, 1L))
    if (fail) stop(errorCondition("fitting failed", class = "report_test_error"))
    42L
  }
  invisible(lapply(c(FALSE, TRUE), function(fail) {
    seen <- list()
    value <- withCallingHandlers(tryCatch(exercise(fail), error = identity),
      warning = function(condition) {
        seen[[length(seen) + 1L]] <<- condition
        invokeRestart("muffleWarning")
      })
    expect_length(seen, 2L)
    expect_s3_class(seen[[1L]], "report_test_warning")
    report <- seen[[2L]]
    expect_s3_class(report, "latents_worker_lost")
    expect_identical(report$lost, 3L)
    expect_identical(report$recovered, 2L)
    expect_identical(report$remaining, 1L)
    expect_identical(report$retry_cores, 2L)
    if (fail) expect_s3_class(value, "report_test_error") else expect_identical(value, 42L)
  }))
})

test_that("a direct screened search combines its two mapper warnings", {
  mapper <- .latents_map_candidates
  local_mocked_bindings(.latents_map_candidates = function(indices, fun, cores, ...) {
    .latents_warn_worker_lost(1L, 1L, 1L)
    mapper(indices, fun, cores = 1L, ...)
  })
  seen <- list()
  value <- withCallingHandlers(.latents_run_starts(1L,
    run = function(index, initial, pause = Inf, resume = NULL) {
      state <- if (is.null(resume)) initial else resume$state
      list(expectation = list(log_likelihood = state), converged = TRUE,
           iterations = 1L, checkpoint = list(state = state))
    }, select_start = "likelihood", prepare = identity,
    screen = list(n_pool = 2L, iterations = 1L, protect = 1L), cores = 2L),
    latents_worker_lost = function(condition) {
      seen[[length(seen) + 1L]] <<- condition
      invokeRestart("muffleWarning")
    })
  expect_identical(value$best_start, 1L)
  expect_length(seen, 1L)
  expect_identical(seen[[1L]]$lost, 2L)
  expect_identical(seen[[1L]]$recovered, 2L)
})
