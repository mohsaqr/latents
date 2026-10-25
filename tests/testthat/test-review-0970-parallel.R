# Review regressions use serialized socket payloads when ports are unavailable.
review_parallel_data <- function() {
  data.frame(person = rep(seq_len(12L), each = 6L),
    y = sin(seq_len(72L)) + rep(seq_len(12L) / 8, each = 6L),
    z = cos(seq_len(72L)) + rep(seq_len(12L) / 10, each = 6L))
}

review_parallel_rng <- function(fun) {
  kind <- RNGkind()
  had_seed <- exists(".Random.seed", globalenv(), inherits = FALSE)
  seed <- if (had_seed) get(".Random.seed", globalenv())
  on.exit({
    suppressWarnings(do.call(RNGkind, as.list(kind)))
    if (had_seed) assign(".Random.seed", seed, globalenv()) else
      if (exists(".Random.seed", globalenv(), inherits = FALSE)) {
        rm(".Random.seed", envir = globalenv())
      }
  })
  fun()
}

review_parallel_worker <- function(fun) {
  review_parallel_rng(function() {
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(998L)
    old_seed <- .Random.seed
    old_kind <- RNGkind()
    old_warn <- getOption("warn")
    on.exit(options(warn = old_warn))
    options(warn = 0L)
    result <- fun()
    expect_identical(RNGkind(), old_kind)
    expect_identical(.Random.seed, old_seed)
    expect_identical(getOption("warn"), 0L)
    result
  })
}

review_parallel_mock_socket <- function(.env = parent.frame()) {
  local_mocked_bindings(
    .latents_make_cluster = function(cores) list(),
    .latents_stop_cluster = function(cluster) invisible(NULL),
    .env = .env)
  local_mocked_bindings(parLapplyLB = function(cl, X, fun, chunk.size) {
    shipped <- unserialize(serialize(fun, NULL))
    review_parallel_worker(function() lapply(X, shipped))
  }, .package = "parallel", .env = .env)
  local_mocked_bindings(.latents_fork_safe = function(...) FALSE, .env = .env)
}

review_parallel_strip_calls <- function(value) {
  value$call <- NULL
  value$fits <- lapply(value$fits, function(fit) {
    if (!is.null(fit)) fit$call <- NULL
    fit
  })
  value
}

test_that("C1 seeded socket tasks inherit all three parent RNG kinds", {
  review_parallel_mock_socket()
  review_parallel_rng(function() {
    kinds <- expand.grid(kind = c("Wichmann-Hill", "Marsaglia-Multicarry",
      "Super-Duper", "Mersenne-Twister", "Knuth-TAOCP", "Knuth-TAOCP-2002",
      "L'Ecuyer-CMRG"), normal.kind = c("Buggy Kinderman-Ramage",
      "Ahrens-Dieter", "Box-Muller", "Inversion", "Kinderman-Ramage"),
      sample.kind = c("Rounding", "Rejection"), stringsAsFactors = FALSE)
    invisible(lapply(seq_len(nrow(kinds)), function(row) {
      kind <- as.list(kinds[row, ])
      suppressWarnings(do.call(RNGkind, kind))
      set.seed(12L)
      state <- .Random.seed
      draw <- function(i) {
        .latents_local_seed(i)
        c(stats::rnorm(3L), sample.int(100L, 5L))
      }
      expected <- .latents_map_candidates(1:2, draw, 1L)
      actual <- .latents_map_candidates(1:2, draw, 2L, backend = "socket")
      expect_identical(actual, expected, info = paste(unlist(kind), collapse = "/"))
      expect_identical(.Random.seed, state)
      expect_identical(RNGkind(), unname(unlist(kind)))
    }))
  })
})

test_that("C1 family and whole-profile fallback candidates preserve seeded fits", {
  review_parallel_mock_socket()
  local_mocked_bindings(.multilpa_pool_fallback = function(arguments) TRUE)
  review_parallel_rng(function() {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(11L)
    invisible(lapply(c("additive", "profiles"), function(family) {
      fit <- function(cores) suppressWarnings(do.call(enumerate_classes, c(
        list(data = review_parallel_data(), vars = c("y", "z"), id = "person",
          family = family, n_group_classes = 1:2, n_starts = 2L,
          max_iter = 15L, seed = 19L, cores = cores),
        if (family == "profiles") list(n_profiles = 2L, model = "VVI"))))
      expected <- fit(1L)
      actual <- fit(2L)
      expect_identical(review_parallel_strip_calls(actual),
                       review_parallel_strip_calls(expected), info = family)
    }))
  })
})

test_that("C1 live socket family enumeration agrees under L'Ecuyer-CMRG", {
  skip_on_cran()
  port <- get("defaultClusterOptions", asNamespace("parallel"))$port
  probe <- tryCatch(serverSocket(port), error = identity)
  if (inherits(probe, "error")) skip(paste("Socket ports unavailable:",
                                          conditionMessage(probe)))
  close(probe)
  local_mocked_bindings(.latents_fork_safe = function(...) FALSE)
  review_parallel_rng(function() {
    RNGkind("L'Ecuyer-CMRG")
    fit <- function(cores) suppressWarnings(enumerate_classes(
      review_parallel_data(), c("y", "z"), "person", family = "additive",
      n_group_classes = 1:2, n_starts = 2L, max_iter = 15L,
      seed = 19L, cores = cores))
    expect_identical(review_parallel_strip_calls(fit(2L)),
                     review_parallel_strip_calls(fit(1L)))
  })
})

test_that("C2 Windows retirement accepts exited workers and rejects live survivors", {
  # Clone the platform-dependent functions so Windows behavior is testable
  # on Unix without modifying base's platform binding.
  context <- new.env(parent = environment(.latents_abort_cluster))
  context$.Platform <- list(OS.type = "windows")
  context$.latents_stop_cluster <- function(cluster) invisible(NULL)
  invisible(lapply(c(".latents_abort_cluster", ".latents_windows_alive"), function(name) {
    if (!exists(name, envir = context, inherits = TRUE)) return(NULL)
    fun <- get(name, envir = context)
    environment(fun) <- context
    assign(name, fun, context)
  }))
  signals <- integer()
  local_mocked_bindings(pskill = function(pid, signal) {
    signals <<- c(signals, signal)
    rep(FALSE, length(pid))
  }, .package = "tools")
  listing <- '"other.exe","222222","Console","1","1,024 K"'
  context$system2 <- function(command, args, stdout, stderr) {
    expect_identical(command, "tasklist")
    listing
  }
  cluster <- structure(list(), latents_pids = 111111L)
  expect_true(context$.latents_abort_cluster(cluster))
  listing <- '"Rscript.exe","111111","Console","1","1,024 K"'
  expect_false(context$.latents_abort_cluster(cluster))
  listing <- structure(character(), status = 1L)
  expect_false(context$.latents_abort_cluster(cluster))
  expect_false(0L %in% signals)
})

test_that("C3 socket abnormal exits abort workers while normal exits stop them", {
  aborted <- stopped <- 0L
  local_mocked_bindings(
    .latents_make_cluster = function(cores) list(),
    .latents_stop_cluster = function(cluster) { stopped <<- stopped + 1L },
    .latents_abort_cluster = function(cluster) { aborted <<- aborted + 1L; TRUE })
  failure <- NULL
  local_mocked_bindings(parLapplyLB = function(cl, X, fun, chunk.size) {
    if (!is.null(failure)) stop(failure)
    lapply(X, fun)
  }, .package = "parallel")
  expect_identical(.latents_socket_pass(1:2, identity, 2L), list(1L, 2L))
  expect_identical(c(aborted, stopped), c(0L, 1L))
  failure <- structure(list(message = "interrupted", call = NULL),
                       class = c("interrupt", "condition"))
  expect_s3_class(tryCatch(.latents_socket_pass(1:2, identity, 2L),
                           interrupt = identity), "interrupt")
  expect_identical(c(aborted, stopped), c(1L, 1L))
  failure <- simpleError("ordinary infrastructure failure")
  expect_error(.latents_socket_pass(1:2, identity, 2L), "infrastructure")
  expect_identical(c(aborted, stopped), c(2L, 1L))
})

test_that("C3 isolated fork abnormal exits kill only their owned child", {
  skip_on_os("windows")
  killed <- integer()
  failure <- NULL
  local_mocked_bindings(mcparallel = function(expr, ...) list(pid = 111111L),
    mccollect = function(job, ...) {
      if (!is.null(failure)) stop(failure)
      list(10L)
    }, .package = "parallel")
  local_mocked_bindings(pskill = function(pid, signal) {
    expect_identical(signal, tools::SIGKILL)
    killed <<- c(killed, pid)
    TRUE
  }, .package = "tools")
  expect_identical(.latents_fork_pass(1L, identity, 1L), list(10L))
  expect_length(killed, 0L)
  failure <- structure(list(message = "interrupted", call = NULL),
                       class = c("interrupt", "condition"))
  expect_s3_class(tryCatch(.latents_fork_pass(1L, identity, 1L),
                           interrupt = identity), "interrupt")
  expect_identical(killed, 111111L)
  failure <- simpleError("collection failed")
  expect_error(.latents_fork_pass(1L, identity, 1L), "collection failed")
  expect_identical(killed, rep(111111L, 2L))
})

test_that("C4 socket recovery attempts every task at most twice", {
  attempts <- integer(6L)
  local_mocked_bindings(
    .latents_make_cluster = function(cores) list(),
    .latents_stop_cluster = function(cluster) invisible(NULL),
    .latents_abort_cluster = function(cluster) TRUE)
  local_mocked_bindings(parLapplyLB = function(cl, X, fun, chunk.size) {
    # Model all tasks of the broken public API call having started.
    attempts[X] <<- attempts[X] + 1L
    if (4L %in% X) stop("error reading from connection")
    lapply(X, fun)
  }, .package = "parallel")
  expect_warning(result <- .latents_map_candidates(1:6, identity, 3L,
    backend = "socket"), class = "latents_worker_lost")
  expect_identical(result[-4L], as.list(c(1L, 2L, 3L, 5L, 6L)))
  expect_match(result[[4L]]$error, "worker process ended")
  expect_true(all(attempts <= 2L))
  expect_identical(attempts[4L], 2L)
})

test_that("C5a converted family warnings remain failed candidate rows", {
  review_parallel_mock_socket()
  local_mocked_bindings(multilpa = function(...) {
    warning(warningCondition("candidate qualification", class = "review_warning"))
    stop("fitting should already have failed")
  })
  # add = TRUE: a bare on.exit() would replace the deferred restoration of
  # the mocks above, leaving multilpa() mocked for every later test file.
  old <- options(warn = 2L)
  on.exit(options(old), add = TRUE)
  fit <- function(cores) .additive_enumerate(review_parallel_data(), c("y", "z"),
    "person", "additive", 1:2, "varying", 19L, list(),
    quote(enumerate_classes(data)), cores)
  expected <- fit(1L)
  actual <- fit(2L)
  expect_identical(actual, expected)
  expect_identical(expected$table$error,
    rep("(converted from warning) candidate qualification", 2L))
  expect_identical(expected$table$warnings, rep("candidate qualification", 2L))
})

test_that("C5b collected conditions precede uncaught errors in grid order", {
  local_mocked_bindings(mclapply = function(X, FUN, mc.cores, mc.preschedule) {
    lapply(X, FUN)
  }, .package = "parallel")
  callback <- function(i) {
    message(paste("message", i))
    warning(warningCondition(paste("warning", i), class = "review_warning"))
    if (i == 3L) stop(errorCondition("task failed", class = "review_error"))
    i
  }
  settings <- expand.grid(signal = c(TRUE, FALSE), cores = 1:2)
  invisible(lapply(seq_len(nrow(settings)), function(row) {
    seen <- character()
    result <- tryCatch(withCallingHandlers(.latents_map_candidates(1:3,
      callback, settings$cores[row], costs = c(1, 3, 2),
      signal = settings$signal[row], backend = "fork"),
      message = function(condition) {
        seen <<- c(seen, trimws(conditionMessage(condition)))
        invokeRestart("muffleMessage")
      }, warning = function(condition) {
        seen <<- c(seen, conditionMessage(condition))
        expect_s3_class(condition, "review_warning")
        invokeRestart("muffleWarning")
      }), error = identity)
    expect_s3_class(result, "review_error")
    expect_identical(seen, as.vector(rbind(paste("message", 1:3),
                                          paste("warning", 1:3))))
  }))
})

test_that("worker execution settings restore on errors and interrupts without a seed", {
  review_parallel_rng(function() {
    RNGkind("L'Ecuyer-CMRG")
    failures <- list(simpleError("task failed"),
      structure(list(message = "interrupted", call = NULL),
                class = c("interrupt", "condition")))
    tasks <- lapply(failures, function(failure) .latents_worker_context(function(i) {
      set.seed(i)
      stop(failure)
    }))
    RNGkind("Mersenne-Twister")
    kind <- RNGkind()
    rm(".Random.seed", envir = globalenv())
    invisible(lapply(seq_along(tasks), function(i) {
      result <- tryCatch(tasks[[i]](19L), error = identity, interrupt = identity)
      expect_identical(result, failures[[i]])
      expect_identical(RNGkind(), kind)
      expect_false(exists(".Random.seed", globalenv(), inherits = FALSE))
    }))
  })
})
