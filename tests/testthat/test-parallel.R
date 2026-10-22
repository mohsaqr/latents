# Fitting candidates on several cores changes the time an enumeration takes and
# nothing it returns: every candidate sets and restores its own seed.

strip_calls <- function(grid) {
  grid$call <- NULL
  grid$fits <- lapply(grid$fits, function(fit) {
    if (!is.null(fit)) fit$call <- NULL
    fit
  })
  grid
}

collect_conditions <- function(expr) {
  seen <- character()
  value <- withCallingHandlers(expr,
    warning = function(condition) {
      seen <<- c(seen, paste("warning:", conditionMessage(condition)))
      invokeRestart("muffleWarning")
    },
    message = function(condition) {
      seen <<- c(seen, paste("message:", conditionMessage(condition)))
      invokeRestart("muffleMessage")
    })
  list(value = value, conditions = seen)
}

test_that("a profile grid on several cores is identical to one core", {
  skip_on_cran()
  skip_on_os("windows")
  data <- subset(course_engagement, student <= 40)
  vars <- c("browse", "lectures", "forum_read")
  one <- enumerate_classes(data, vars, "student", n_profiles = 1:3,
                           n_group_classes = 1:2, model = "VVI", n_starts = 2,
                           seed = 1)
  two <- enumerate_classes(data, vars, "student", n_profiles = 1:3,
                           n_group_classes = 1:2, model = "VVI", n_starts = 2,
                           seed = 1, cores = 2)
  expect_identical(strip_calls(two), strip_calls(one))
  expect_identical(as.data.frame(two), as.data.frame(one))
})

test_that("transition and family grids on several cores are identical to one core", {
  skip_on_cran()
  skip_on_os("windows")
  data <- subset(course_engagement, student <= 30)
  vars <- c("browse", "lectures")
  transition <- lapply(1:2, \(cores) enumerate_classes(
    data, vars, "student", time = "sequence", n_profiles = 2:3, n_starts = 2,
    seed = 1, cores = cores))
  expect_identical(strip_calls(transition[[2]]), strip_calls(transition[[1]]))
  family <- lapply(1:2, \(cores) collect_conditions(enumerate_classes(
    data, vars, "student", family = c("additive", "dispersion"),
    n_group_classes = 1:2, seed = 1, cores = cores)))
  expect_identical(strip_calls(family[[2]]$value), strip_calls(family[[1]]$value))
})

test_that("candidates' warnings reach the caller in grid order, as on one core", {
  skip_on_cran()
  skip_on_os("windows")
  vars <- c("cognitive_strategies", "intrinsic_value", "self_efficacy",
            "self_regulation", "test_anxiety")
  runs <- lapply(1:2, \(cores) collect_conditions(enumerate_lpa(
    srl, vars, n_profiles = 6:8, model = "VVI", n_starts = 3, seed = 1,
    cores = cores, search = "kmeans")))
  expect_true(length(runs[[1]]$conditions) > 0L)
  expect_identical(runs[[2]]$conditions, runs[[1]]$conditions)
  expect_identical(as.data.frame(runs[[2]]$value)$warnings,
                   as.data.frame(runs[[1]]$value)$warnings)
})

test_that("the socket backend returns values and conditions in order", {
  skip_on_cran()
  skip_if_not(requireNamespace("latents", quietly = TRUE))
  fun <- function(i) {
    if (i == 2L) warning(warningCondition(paste("warned", i), class = "test_warning"))
    message(paste("message", i))
    i * 10L
  }
  seen <- list()
  values <- withCallingHandlers(
    .latents_map_candidates(1:3, fun, cores = 2L, backend = "socket"),
    warning = function(condition) {
      seen[[length(seen) + 1L]] <<- condition
      invokeRestart("muffleWarning")
    },
    message = function(condition) {
      seen[[length(seen) + 1L]] <<- condition
      invokeRestart("muffleMessage")
    })
  expect_identical(values, list(10L, 20L, 30L))
  expect_identical(vapply(seen, conditionMessage, character(1)),
                   c("message 1\n", "warned 2", "message 2\n", "message 3\n"))
  expect_true(inherits(seen[[2L]], "test_warning"))
})

test_that("cores is validated by class, and one core is lapply", {
  data <- subset(course_engagement, student <= 10)
  invisible(lapply(list(0, 1.5, NA, "2", c(1, 2)), function(cores) {
    expect_error(enumerate_classes(data, "browse", "student", n_profiles = 1,
                                   n_group_classes = 1, cores = cores),
                 class = "latents_bad_argument")
  }))
  expect_identical(.latents_map_candidates(1:3, \(i) i^2, cores = 1L), list(1, 4, 9))
})

test_that("forked workers are refused only on macOS with Accelerate's BLAS", {
  veclib <- paste0("/System/Library/Frameworks/Accelerate.framework/Versions/A/",
                   "Frameworks/vecLib.framework/Versions/A/libBLAS.dylib")
  reference <- "/Library/Frameworks/R.framework/Resources/lib/libRblas.0.dylib"
  expect_false(.latents_fork_safe(veclib, "Darwin"))
  expect_true(.latents_fork_safe(reference, "Darwin"))
  expect_true(.latents_fork_safe(veclib, "Linux"))
  expect_true(.latents_fork_safe("/usr/lib/libopenblas.so", "Linux"))
})

test_that("a forked worker killed by a signal is a parallel failure, not a NULL candidate", {
  skip_on_cran()
  skip_on_os("windows")
  # The worker for candidate 2 ends itself as a segmentation fault would.
  fun <- function(i) {
    if (i == 2L) tools::pskill(Sys.getpid(), tools::SIGKILL)
    i
  }
  expect_error(suppressWarnings(
    .latents_map_candidates(1:3, fun, cores = 2L, backend = "fork")),
    class = "latents_parallel_failure")
})

test_that("a socket worker that dies is a parallel failure", {
  skip_on_cran()
  skip_if_not(requireNamespace("latents", quietly = TRUE))
  fun <- function(i) {
    # SIGTERM is defined on Windows and Unix; SIGKILL is NA on Windows.
    if (i == 2L) tools::pskill(Sys.getpid(), tools::SIGTERM)
    i
  }
  expect_error(.latents_map_candidates(1:3, fun, cores = 2L, backend = "socket"),
               class = "latents_parallel_failure")
})

test_that("one fit's starts on several cores give the same fit as one core", {
  skip_on_cran()
  skip_on_os("windows")
  items <- c("time_with_friends", "on_social_media", "sports")
  sample_data <- subset(student_esm, day <= 3)
  serial <- multilca(sample_data, vars = items, id = "student", n_profiles = 3,
                     n_group_classes = 2, n_starts = 6, seed = 1)
  spread <- multilca(sample_data, vars = items, id = "student", n_profiles = 3,
                     n_group_classes = 2, n_starts = 6, seed = 1, cores = 2)
  expect_identical(spread$log_likelihood, serial$log_likelihood)
  expect_identical(get_results(spread, "starts"), get_results(serial, "starts"))
  profiles <- c("happy", "relaxed")
  serial_lpa <- lpa(student_esm, profiles, n_profiles = 2, n_starts = 5,
                    seed = 2)
  # Two workers respect check limits; five starts keep the split uneven.
  spread_lpa <- lpa(student_esm, profiles, n_profiles = 2, n_starts = 5,
                    seed = 2, cores = 2)
  expect_identical(spread_lpa$log_likelihood, serial_lpa$log_likelihood)
  expect_identical(as.data.frame(spread_lpa), as.data.frame(serial_lpa))
})

test_that("a fit refuses a `cores` that is not a whole number of at least 1", {
  items <- c("time_with_friends", "on_social_media", "sports")
  expect_error(lca(student_esm, items, n_classes = 2, n_starts = 1, cores = 0),
               class = "latents_bad_argument")
  expect_error(lca(student_esm, items, n_classes = 2, n_starts = 1, cores = 1.5),
               class = "latents_bad_argument")
})

test_that("socket cleanup cannot replace a worker failure or leak remaining connections", {
  exercise <- function() {
    connections <- lapply(1:2, function(i) rawConnection(raw(), "w"))
    on.exit(invisible(lapply(connections, function(connection) {
      tryCatch(close(connection), error = function(condition) NULL)
    })), add = TRUE)
    nodes <- structure(lapply(connections, function(connection) list(con = connection)),
                       class = c("SOCKcluster", "cluster"))
    testthat::local_mocked_bindings(
      makePSOCKcluster = function(...) nodes,
      parLapply = function(...) stop("error reading from connection"),
      stopCluster = function(...) stop("error writing to connection during cleanup"),
      .package = "parallel")
    warnings <- character()
    failure <- tryCatch(withCallingHandlers(
      .latents_map_candidates(1:3, identity, cores = 2, backend = "socket"),
      warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }), error = identity)
    expect_s3_class(failure, "latents_parallel_failure")
    expect_match(conditionMessage(failure), "error reading from connection")
    expect_length(warnings, 0L)
    invisible(lapply(connections, function(connection) {
      expect_error(isOpen(connection), "invalid connection")
    }))
  }
  exercise()
})

test_that("cleanup closes a live socket even when another node is already closed", {
  skip_on_cran()
  cluster <- parallel::makePSOCKcluster(2)
  close(cluster[[1L]]$con)
  expect_silent(.latents_stop_cluster(cluster))
  invisible(lapply(cluster, function(node) {
    expect_error(isOpen(node$con), "invalid connection")
  }))
  # Repeated cleanup of already closed nodes also remains harmless.
  expect_silent(.latents_stop_cluster(cluster))
})
