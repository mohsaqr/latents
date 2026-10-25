# Process-loss tests use marker files because forked workers cannot update the
# parent's memory. Guard every kill so a broken serial retry cannot kill R.
recovery_skip_forks <- function() {
  skip_on_cran()
  skip_on_os("windows")
  skip_if(!.latents_fork_safe(), "This BLAS cannot run safely in forked workers")
}

recovery_capture <- function(expr) {
  conditions <- list()
  value <- withCallingHandlers(expr,
    warning = function(condition) {
      conditions[[length(conditions) + 1L]] <<- condition
      invokeRestart("muffleWarning")
    },
    message = function(condition) invokeRestart("muffleMessage"))
  list(value = value, conditions = conditions)
}

recovery_worker_warnings <- function(result) {
  Filter(\(condition) inherits(condition, "latents_worker_lost"),
         result$conditions)
}

recovery_kill <- function(parent_pid) {
  stopifnot("only an isolated worker can be killed" = Sys.getpid() != parent_pid)
  tools::pskill(Sys.getpid(), tools::SIGKILL)
}

recovery_panel <- function() {
  .latents_local_seed(8401L)
  panel <- data.frame(person = rep(seq_len(12L), each = 5L),
                      beep = rep(seq_len(5L), 12L))
  state <- (panel$person + panel$beep) %% 2L
  panel$a <- 2 * state + stats::rnorm(nrow(panel), sd = 0.6)
  panel$b <- -1.5 * state + stats::rnorm(nrow(panel), sd = 0.7)
  panel
}

recovery_grid <- function(panel, cores, search = "kmeans") {
  enumerate_classes(panel, c("a", "b"), "person", time = "beep",
    n_profiles = 2:3, n_group_classes = 1L, model = "VVI", n_starts = 2L,
    max_iter = 10L, search = search,
    n_screen = if (search == "emEM") 4L else NULL, seed = 12L, cores = cores)
}

recovery_without_calls <- function(value) {
  value$call <- NULL
  value$fits <- lapply(value$fits, function(fit) {
    if (!is.null(fit)) fit$call <- NULL
    fit
  })
  value
}

test_that("a killed fork task is retried once and recovers its exact value", {
  recovery_skip_forks()
  marker <- tempfile("latents-recovered-task-")
  on.exit(unlink(marker), add = TRUE)
  parent_pid <- Sys.getpid()
  values <- list(c(0.25, 0.75), c(0.4, 0.6), c(0.2, 0.8))
  fun <- function(i) {
    if (i == 2L && !file.exists(marker)) {
      writeLines("first attempt", marker)
      recovery_kill(parent_pid)
    }
    list(start = i, initial = values[[i]])
  }
  expected <- lapply(seq_along(values), \(i) list(start = i, initial = values[[i]]))
  result <- recovery_capture(.latents_map_candidates(
    seq_along(values), fun, cores = 2L, backend = "fork"))
  expect_identical(result$value, expected)
  expect_true(file.exists(marker))
  lost <- recovery_worker_warnings(result)
  expect_length(lost, 1L)
  expect_equal(lost[[1L]]$lost, 1L)
  expect_equal(lost[[1L]]$recovered, 1L)
  expect_equal(lost[[1L]]$remaining, 0L)
  expect_match(conditionMessage(lost[[1L]]), "cores")
})

test_that("a permanently killed fork task leaves other task results intact", {
  recovery_skip_forks()
  parent_pid <- Sys.getpid()
  result <- recovery_capture(.latents_map_candidates(seq_len(3L), function(i) {
    if (i == 2L) recovery_kill(parent_pid)
    i * 10L
  }, cores = 2L, backend = "fork"))
  expect_identical(result$value[c(1L, 3L)], list(10L, 30L))
  expect_match(result$value[[2L]]$error,
               "worker process ended without returning a result", ignore.case = TRUE)
  expect_match(result$value[[2L]]$error, "memory", ignore.case = TRUE)
  lost <- recovery_worker_warnings(result)
  expect_length(lost, 1L)
  expect_equal(lost[[1L]]$lost, 1L)
  expect_equal(lost[[1L]]$recovered, 0L)
  expect_equal(lost[[1L]]$remaining, 1L)
})

test_that("a socket retry loses only the task that kills its worker", {
  skip_on_cran()
  skip_on_os("windows")
  # Use the same port probe as the other live socket tests in this suite.
  port <- get("defaultClusterOptions", envir = asNamespace("parallel"))$port
  probe <- tryCatch(serverSocket(port), error = identity)
  if (inherits(probe, "error")) {
    skip(paste("Socket ports unavailable:", conditionMessage(probe)))
  }
  close(probe)
  parent_pid <- Sys.getpid()
  # Task 4 kills every worker it runs on. Its chunk mates in the retry must
  # still come back, which they did not when the retry ran in whole chunks.
  result <- recovery_capture(.latents_map_candidates(seq_len(6L), function(i) {
    if (i == 4L) recovery_kill(parent_pid)
    i * 10L
  }, cores = 3L, backend = "socket"))
  expect_identical(result$value[-4L], list(10L, 20L, 30L, 50L, 60L))
  expect_match(result$value[[4L]]$error,
               "worker process ended without returning a result", ignore.case = TRUE)
  lost <- recovery_worker_warnings(result)
  expect_length(lost, 1L)
  expect_equal(lost[[1L]]$remaining, 1L)
})

test_that("ordinary errors in start code are recorded without a process-loss retry", {
  recovery_skip_forks()
  directory <- tempfile("latents-start-errors-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  runner <- function(start_index, initial) {
    cat("attempt\n", file = file.path(directory, as.character(start_index)),
        append = TRUE)
    if (start_index == 2L) {
      stop(errorCondition("deliberate start error", class = "recovery_test_error"))
    }
    list(value = initial)
  }
  result <- recovery_capture(.latents_map_candidates(seq_len(3L), function(i) {
    .latents_start_attempt(runner, i, list(values = i * 10L))
  }, cores = 2L, backend = "fork"))
  expect_identical(result$value,
                   list(list(value = 10L), list(error = "deliberate start error"),
                        list(value = 30L)))
  expect_length(recovery_worker_warnings(result), 0L)
  expect_identical(vapply(seq_len(3L), function(i) {
    length(readLines(file.path(directory, as.character(i))))
  }, integer(1)), rep(1L, 3L))
})

test_that("recovered flat starts preserve the complete transition grid", {
  recovery_skip_forks()
  panel <- recovery_panel()
  expected <- recovery_capture(recovery_grid(panel, 1L))
  marker <- tempfile("latents-recovered-start-")
  on.exit(unlink(marker), add = TRUE)
  parent_pid <- Sys.getpid()
  execute <- .lta_execute_start
  local_mocked_bindings(.lta_execute_start = function(task) {
    if (task$candidate == 1L && task$start_index == 1L && !file.exists(marker)) {
      writeLines("first attempt", marker)
      recovery_kill(parent_pid)
    }
    execute(task)
  })
  result <- recovery_capture(recovery_grid(panel, 2L))
  expect_identical(recovery_without_calls(result$value),
                   recovery_without_calls(expected$value))
  lost <- recovery_worker_warnings(result)
  expect_length(lost, 1L)
  expect_equal(lost[[1L]]$lost, 1L)
  expect_equal(lost[[1L]]$recovered, 1L)
})

test_that("a recovered direct lta start preserves the undisturbed fitted object", {
  recovery_skip_forks()
  panel <- recovery_panel()
  fit <- function(cores) lta(panel, c("a", "b"), "person", 2L, time = "beep",
    n_group_classes = 1L, model = "VVI", n_starts = 2L, max_iter = 10L,
    seed = 12L, cores = cores)
  expected <- recovery_capture(fit(1L))
  marker <- tempfile("latents-recovered-direct-start-")
  on.exit(unlink(marker), add = TRUE)
  parent_pid <- Sys.getpid()
  run <- .lta_run_start
  local_mocked_bindings(.lta_run_start = function(start_index, ...) {
    if (start_index == 1L && !file.exists(marker)) {
      writeLines("first attempt", marker)
      recovery_kill(parent_pid)
    }
    run(start_index, ...)
  })
  result <- recovery_capture(fit(2L))
  result$value$call <- NULL
  expected$value$call <- NULL
  expect_identical(result$value, expected$value)
  lost <- recovery_worker_warnings(result)
  expect_length(lost, 1L)
  expect_equal(lost[[1L]]$recovered, 1L)
})

test_that("a genuine transition start error is neither retried nor called worker loss", {
  recovery_skip_forks()
  panel <- recovery_panel()
  marker <- tempfile("latents-genuine-start-error-")
  on.exit(unlink(marker), add = TRUE)
  run <- .lta_run_start
  local_mocked_bindings(.lta_run_start = function(start_index, ...) {
    if (start_index == 1L) {
      cat("attempt\n", file = marker, append = TRUE)
      stop(errorCondition("deliberate transition start error",
                           class = "recovery_test_error"))
    }
    run(start_index, ...)
  })
  result <- recovery_capture(recovery_grid(panel, 2L))
  expect_true(all(is.na(result$value$table$error)))
  expect_length(readLines(marker), 2L)
  invisible(lapply(result$value$fits, function(fit) {
    starts <- get_results(fit, "starts")
    expect_identical(starts$error[1L], "deliberate transition start error")
    expect_true(is.na(starts$error[2L]))
  }))
  expect_length(recovery_worker_warnings(result), 0L)
})

test_that("a permanently lost flat start is reported alongside fitted candidates", {
  recovery_skip_forks()
  panel <- recovery_panel()
  expected <- recovery_capture(recovery_grid(panel, 1L))
  parent_pid <- Sys.getpid()
  execute <- .lta_execute_start
  local_mocked_bindings(.lta_execute_start = function(task) {
    if (task$candidate == 1L && task$start_index == 1L) recovery_kill(parent_pid)
    execute(task)
  })
  result <- recovery_capture(recovery_grid(panel, 2L))
  expect_true(all(is.na(result$value$table$error)))
  expect_true(all(is.finite(result$value$table$log_likelihood)))
  starts <- get_results(result$value$fits[[1L]], "starts")
  expect_match(starts$error[1L], "worker process ended", ignore.case = TRUE)
  expect_match(starts$error[1L], "memory", ignore.case = TRUE)
  expect_true(is.na(starts$error[2L]))
  expect_identical(starts$log_likelihood[2L],
                   get_results(expected$value$fits[[1L]], "starts")$log_likelihood[2L])
  expect_identical(result$value$fits[[1L]]$n_failed_starts, 1L)
  expect_identical(result$value$fits[[2L]], expected$value$fits[[2L]])
  failed <- Filter(\(condition) inherits(condition, "latents_failed_starts"),
                   result$conditions)
  expect_length(failed, 1L)
  expect_match(conditionMessage(failed[[1L]]), "1 of 2 starts failed", fixed = TRUE)
  lost <- recovery_worker_warnings(result)
  expect_length(lost, 1L)
  expect_equal(lost[[1L]]$remaining, 1L)
})

test_that("losing all starts fails only their candidate in a transition grid", {
  recovery_skip_forks()
  panel <- recovery_panel()
  expected <- recovery_capture(recovery_grid(panel, 1L))
  parent_pid <- Sys.getpid()
  execute <- .lta_execute_start
  local_mocked_bindings(.lta_execute_start = function(task) {
    if (task$candidate == 1L) recovery_kill(parent_pid)
    execute(task)
  })
  result <- recovery_capture(recovery_grid(panel, 2L))
  expect_null(result$value$fits[[1L]])
  expect_match(result$value$table$error[1L], "All 2 starts failed", fixed = TRUE)
  expect_match(result$value$table$error[1L], "worker process ended", ignore.case = TRUE)
  expect_true(is.na(result$value$table$log_likelihood[1L]))
  expect_identical(result$value$fits[[2L]], expected$value$fits[[2L]])
  expect_identical(result$value$table[2L, ], expected$value$table[2L, ])
  lost <- recovery_worker_warnings(result)
  expect_length(lost, 1L)
  expect_equal(lost[[1L]]$lost, 2L)
  expect_equal(lost[[1L]]$remaining, 2L)
})

test_that("emEM aggregates recovered screening and finalist losses into one warning", {
  recovery_skip_forks()
  panel <- recovery_panel()
  expected <- recovery_capture(recovery_grid(panel, 1L, "emEM"))
  directory <- tempfile("latents-recovered-emem-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  parent_pid <- Sys.getpid()
  execute <- .lta_execute_start
  local_mocked_bindings(.lta_execute_start = function(task) {
    phase <- if (is.null(task$screen)) "final" else "screen"
    marker <- file.path(directory, phase)
    if (task$candidate == 1L && task$start_index == 1L && !file.exists(marker)) {
      writeLines("first attempt", marker)
      recovery_kill(parent_pid)
    }
    execute(task)
  })
  result <- recovery_capture(recovery_grid(panel, 2L, "emEM"))
  expect_identical(recovery_without_calls(result$value),
                   recovery_without_calls(expected$value))
  expect_true(all(file.exists(file.path(directory, c("screen", "final")))))
  lost <- recovery_worker_warnings(result)
  expect_length(lost, 1L)
  expect_equal(lost[[1L]]$lost, 2L)
  expect_equal(lost[[1L]]$recovered, 2L)
  expect_equal(lost[[1L]]$remaining, 0L)
})
