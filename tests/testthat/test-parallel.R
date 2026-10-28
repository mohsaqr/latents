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
      # Progress messages name the worker count, which is what differs.
      if (!inherits(condition, "latents_progress")) {
        seen <<- c(seen, paste("message:", conditionMessage(condition)))
      }
      invokeRestart("muffleMessage")
    })
  list(value = value, conditions = seen)
}

skip_if_socket_ports_unavailable <- function() {
  # Probe the port PSOCK uses, so sandbox restrictions skip the live-socket
  # tests while a worker startup regression still fails.
  port <- get("defaultClusterOptions", envir = asNamespace("parallel"))$port
  socket <- tryCatch(serverSocket(port), error = identity)
  if (inherits(socket, "error")) {
    skip(paste("Socket ports unavailable:", conditionMessage(socket)))
  }
  close(socket)
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
  skip_if_socket_ports_unavailable()
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

test_that("a forked worker killed by a signal becomes a failed candidate", {
  skip_on_cran()
  skip_on_os("windows")
  # The worker for candidate 2 ends itself as a segmentation fault would.
  fun <- function(i) {
    if (i == 2L) tools::pskill(Sys.getpid(), tools::SIGKILL)
    i
  }
  expect_warning(
    result <- .latents_map_candidates(1:3, fun, cores = 2L, backend = "fork"),
    class = "latents_worker_lost")
  expect_identical(result[c(1L, 3L)], list(1L, 3L))
  expect_match(result[[2L]]$error, "likely killed for memory", fixed = TRUE)
})

test_that("a socket worker that dies becomes a failed candidate", {
  skip_on_cran()
  skip_if_socket_ports_unavailable()
  skip_if_not(requireNamespace("latents", quietly = TRUE))
  fun <- function(i) {
    # SIGTERM is defined on Windows and Unix; SIGKILL is NA on Windows.
    if (i == 2L) tools::pskill(Sys.getpid(), tools::SIGTERM)
    i
  }
  expect_warning(
    result <- .latents_map_candidates(1:3, fun, cores = 2L, backend = "socket"),
    class = "latents_worker_lost")
  expect_identical(result[[3L]], 3L)
  expect_match(result[[2L]]$error, "likely killed for memory", fixed = TRUE)
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
      clusterCall = function(cl, fun) rep(list(Sys.getpid()), length(cl)),
    makePSOCKcluster = function(...) nodes,
      parLapplyLB = function(...) stop("error reading from connection"),
      stopCluster = function(...) stop("error writing to connection during cleanup"),
      .package = "parallel")
    warnings <- character()
    failure <- tryCatch(withCallingHandlers(
      .latents_map_candidates(1:3, identity, cores = 2, backend = "socket"),
      warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }), error = identity)
    expect_true(all(vapply(failure, function(result) {
      grepl("likely killed for memory", result$error, fixed = TRUE)
    }, logical(1))))
    expect_length(warnings, 1L)
    expect_match(warnings, "3 worker tasks were lost", fixed = TRUE)
    invisible(lapply(connections, function(connection) {
      expect_error(isOpen(connection), "invalid connection")
    }))
  }
  exercise()
})

test_that("cleanup closes a live socket even when another node is already closed", {
  skip_on_cran()
  skip_if_socket_ports_unavailable()
  cluster <- parallel::makePSOCKcluster(2)
  close(cluster[[1L]]$con)
  expect_silent(.latents_stop_cluster(cluster))
  invisible(lapply(cluster, function(node) {
    expect_error(isOpen(node$con), "invalid connection")
  }))
  # Repeated cleanup of already closed nodes also remains harmless.
  expect_silent(.latents_stop_cluster(cluster))
})

test_that("every enumeration path sends its data, not the caller's variable, to workers", {
  # A socket worker unserializes the candidates' closure in a session that
  # lacks the caller's global variables; this replays exactly that.
  original_mapper <- .latents_map_candidates
  worker_view <- function(indices, fun, cores, signal = TRUE, ...) {
    shipped <- unserialize(serialize(list(indices = indices, fun = fun), NULL))
    rm(list = "worker_test_data", envir = globalenv())
    original_mapper(shipped$indices, shipped$fun, cores = 1L, signal = signal)
  }
  local_mocked_bindings(.latents_map_candidates = worker_view)
  grids <- list(
    transition = quote(enumerate_classes(worker_test_data, c("browse", "lectures"),
      "student", time = "sequence", n_profiles = 2, n_starts = 1, seed = 1, cores = 2)),
    family = quote(enumerate_classes(worker_test_data, c("browse", "lectures"),
      "student", family = "additive", n_group_classes = 1, n_starts = 1, seed = 1,
      cores = 2)),
    profiles = quote(enumerate_classes(worker_test_data, c("browse", "lectures"),
      "student", n_profiles = 2, n_group_classes = 1, n_starts = 1, seed = 1, cores = 2)))
  errors <- vapply(grids, function(grid) {
    assign("worker_test_data", subset(course_engagement, student <= 30),
           envir = globalenv())
    on.exit(suppressWarnings(rm(list = "worker_test_data", envir = globalenv())))
    table <- as.data.frame(suppressWarnings(eval(grid, globalenv())))
    sum(!is.na(table$error))
  }, integer(1))
  expect_identical(errors, c(transition = 0L, family = 0L, profiles = 0L))
})

test_that("the shared pool dispatches longest tasks first within its worker budget", {
  dispatched <- NULL
  workers <- NULL
  prescheduled <- NULL
  local_mocked_bindings(
    mclapply = function(X, FUN, mc.cores, mc.preschedule) {
      dispatched <<- X
      workers <<- mc.cores
      prescheduled <<- mc.preschedule
      lapply(X, FUN)
    }, .package = "parallel")
  result <- .latents_map_candidates(seq_len(5L), identity, cores = 3L,
    backend = "fork", costs = c(1, 4, 2, 4, 3))
  expect_identical(dispatched, c(2L, 4L, 5L, 3L, 1L))
  expect_identical(workers, 3L)
  expect_false(prescheduled)
  expect_identical(result, as.list(seq_len(5L)))
  .latents_map_candidates(seq_len(2L), identity, cores = 3L,
                         backend = "fork", costs = c(1, 2))
  expect_identical(workers, 2L)
})

test_that("collected task conditions retain their classes and original task order", {
  callback <- function(i) {
    warning(warningCondition(paste("warning", i), class = "pool_test_warning"))
    message(structure(list(message = paste("message", i), call = NULL),
      class = c("pool_test_message", "message", "condition")))
    i
  }
  collected <- .latents_map_candidates(seq_len(3L), callback, cores = 1L,
    costs = c(1, 3, 2), signal = FALSE)
  expect_identical(lapply(collected, `[[`, "value"), as.list(seq_len(3L)))
  expect_identical(vapply(collected, function(task) {
    conditionMessage(task$conditions[[1L]])
  }, character(1)), paste("warning", seq_len(3L)))
  expect_s3_class(collected[[1L]]$conditions[[1L]], "pool_test_warning")
  expect_s3_class(collected[[1L]]$conditions[[2L]], "pool_test_message")
  replay <- collect_conditions(.latents_signal_conditions(
    collected[[1L]]$conditions))
  expect_identical(replay$conditions, c("warning: warning 1", "message: message 1"))
  serial <- collect_conditions(.latents_map_candidates(seq_len(3L), callback,
    cores = 1L, costs = c(1, 3, 2)))
  expect_identical(serial$value, as.list(seq_len(3L)))
  expect_identical(serial$conditions, as.vector(rbind(
    paste("warning: warning", seq_len(3L)), paste("message: message", seq_len(3L)))))
})

test_that("the socket pool sends forced task inputs and dynamically balances them", {
  global_names <- c("pool_indices", "pool_callback", "pool_cores", "pool_costs")
  assign("pool_indices", seq_len(4L), envir = globalenv())
  assign("pool_callback", function(i) i * 3L, envir = globalenv())
  assign("pool_cores", 3L, envir = globalenv())
  assign("pool_costs", c(1, 4, 2, 3), envir = globalenv())
  on.exit(rm(list = intersect(global_names, ls(envir = globalenv())),
             envir = globalenv()))
  dispatched <- integer()
  workers <- NULL
  chunk <- NULL
  local_mocked_bindings(
    clusterCall = function(cl, fun) rep(list(Sys.getpid()), length(cl)),
    makePSOCKcluster = function(cores) {
      workers <<- cores
      list()
    },
    parLapplyLB = function(cl, X, fun, chunk.size) {
      dispatched <<- c(dispatched, X)
      chunk <<- chunk.size
      shipped <- unserialize(serialize(fun, NULL))
      rm(list = intersect(global_names, ls(envir = globalenv())), envir = globalenv())
      lapply(X, shipped)
    }, .package = "parallel")
  values <- .latents_map_candidates(pool_indices, pool_callback, pool_cores,
    backend = "socket", costs = pool_costs, signal = FALSE)
  expect_identical(dispatched, c(2L, 4L, 3L, 1L))
  expect_identical(workers, 3L)
  expect_identical(chunk, 1L)
  expect_identical(lapply(values, `[[`, "value"), list(3L, 6L, 9L, 12L))
  expect_true(all(vapply(values, function(task) length(task$conditions) == 0L,
                        logical(1))))
})


test_that("a serialized task collector captures only its forced task function", {
  collector <- .latents_task_collector(identity)
  expect_identical(ls(envir = environment(collector), all.names = TRUE), "fun")
  shipped <- unserialize(serialize(collector, NULL))
  expect_identical(ls(envir = environment(shipped), all.names = TRUE), "fun")
  expect_identical(shipped(list(value = 3L)),
                   list(value = list(value = 3L), conditions = list()))
})

test_that("socket restart and retry startup failures retain completed chunks", {
  creations <- 0L
  dispatched <- list()
  local_mocked_bindings(
    clusterCall = function(cl, fun) rep(list(Sys.getpid()), length(cl)),
    makePSOCKcluster = function(cores) {
      creations <<- creations + 1L
      if (creations > 1L) stop("cannot create another worker")
      list()
    },
    parLapplyLB = function(cl, X, fun, chunk.size) {
      dispatched[[length(dispatched) + 1L]] <<- X
      if (3L %in% X) stop("error reading from connection")
      lapply(X, fun)
    }, .package = "parallel")
  seen <- NULL
  result <- withCallingHandlers(
    .latents_map_candidates(seq_len(6L), identity, cores = 2L,
                             backend = "socket"),
    latents_worker_lost = function(condition) {
      seen <<- condition
      invokeRestart("muffleWarning")
    })
  expect_identical(result[seq_len(2L)], list(1L, 2L))
  expect_true(all(vapply(result[3:6], function(value) {
    grepl("likely killed for memory", value$error, fixed = TRUE)
  }, logical(1))))
  expect_identical(dispatched, list(1:2, 3:4))
  expect_identical(creations, 3L)
  expect_identical(seen$lost, 4L)
  expect_identical(seen$recovered, 0L)
})

test_that("initial socket startup failure remains a classed infrastructure error", {
  local_mocked_bindings(
    clusterCall = function(cl, fun) rep(list(Sys.getpid()), length(cl)),
    makePSOCKcluster = function(...) stop("cannot create the first worker"),
    .package = "parallel")
  expect_error(.latents_map_candidates(1:2, identity, cores = 2L,
                                        backend = "socket"),
               class = "latents_parallel_failure")
})

test_that("a fork retry that cannot spawn retains earlier recovered tasks", {
  # Windows has no fork, so `parallel` has no mcparallel() to mock.
  skip_on_os("windows")
  creations <- 0L
  local_mocked_bindings(
    mcparallel = function(expr, ...) {
      creations <<- creations + 1L
      if (creations > 1L) stop("unable to fork")
      list(value = force(expr))
    },
    mccollect = function(job) list(job$value),
    .package = "parallel")
  result <- .latents_fork_pass(seq_len(3L), .latents_task_collector(identity),
                               cores = 1L, retry = TRUE)
  expect_identical(result[[1L]], list(value = 1L, conditions = list()))
  expect_identical(result[-1L], list(NULL, NULL))
  expect_identical(creations, 2L)
})

test_that("a terminal loss callback error cannot hide the worker warning", {
  # Windows has no fork, so `parallel` has no mcparallel() to mock.
  skip_on_os("windows")
  local_mocked_bindings(
    mclapply = function(X, ...) rep(list(NULL), length(X)),
    mcparallel = function(...) stop("unable to fork"),
    .package = "parallel")
  seen <- NULL
  result <- tryCatch(withCallingHandlers(
    .latents_map_candidates(1:2, identity, cores = 2L, backend = "fork",
      on_lost = function(task, message) {
        stop(errorCondition("terminal callback failed", class = "callback_error"))
      }),
    latents_worker_lost = function(condition) {
      seen <<- condition
      invokeRestart("muffleWarning")
    }), error = identity)
  expect_s3_class(result, "callback_error")
  expect_s3_class(seen, "latents_worker_lost")
  expect_identical(seen$remaining, 2L)
})

test_that("socket chunk recovery retries only unfinished tasks on a fresh smaller cluster", {
  creations <- integer()
  dispatched <- list()
  failed <- FALSE
  local_mocked_bindings(
    clusterCall = function(cl, fun) rep(list(Sys.getpid()), length(cl)),
    makePSOCKcluster = function(cores) {
      creations <<- c(creations, cores)
      list()
    },
    parLapplyLB = function(cl, X, fun, chunk.size) {
      dispatched[[length(dispatched) + 1L]] <<- X
      if (3L %in% X && !failed) {
        failed <<- TRUE
        stop("error reading from connection")
      }
      lapply(X, unserialize(serialize(fun, NULL)))
    }, .package = "parallel")
  seen <- NULL
  result <- withCallingHandlers(
    .latents_map_candidates(seq_len(5L), identity, cores = 2L,
                             backend = "socket"),
    latents_worker_lost = function(condition) {
      seen <<- condition
      invokeRestart("muffleWarning")
    })
  expect_identical(result, as.list(seq_len(5L)))
  expect_identical(dispatched, list(1:2, 3:4, 5L, 3L, 4L))
  expect_identical(creations, c(2L, 2L, 1L))
  expect_identical(seen$lost, 2L)
  expect_identical(seen$recovered, 2L)
  expect_identical(seen$remaining, 0L)
})

test_that("an uncaught task R error retains its class and is never retried", {
  attempts <- integer(3L)
  worker_requests <- integer()
  failure <- errorCondition("a genuine task error", class = "test_task_error")
  local_mocked_bindings(
    mclapply = function(X, FUN, mc.cores, mc.preschedule) {
      worker_requests <<- c(worker_requests, mc.cores)
      lapply(X, FUN)
    }, .package = "parallel")
  result <- tryCatch(.latents_map_candidates(seq_len(3L), function(task) {
    attempts[task] <<- attempts[task] + 1L
    if (task == 2L) stop(failure)
    task
  }, cores = 2L, backend = "fork"), error = identity)
  expect_identical(result, failure)
  expect_s3_class(result, "test_task_error")
  expect_identical(attempts, rep(1L, 3L))
  expect_identical(worker_requests, 2L)
})

test_that("broken socket workers are retired before replacement pools start", {
  skip_on_os("windows")
  active <- integer()
  workers <- integer()
  signals <- integer()
  overlap <- FALSE
  busy <- FALSE
  failed <- FALSE
  next_pid <- 100000L
  local_mocked_bindings(.latents_stop_cluster = function(cluster) {
    if (!busy) active <<- setdiff(active, attr(cluster, "latents_pids"))
    invisible(NULL)
  })
  local_mocked_bindings(
    makePSOCKcluster = function(cores) {
      overlap <<- overlap || length(active) > 0L
      workers <<- c(workers, cores)
      pids <- next_pid + seq_len(cores)
      next_pid <<- next_pid + cores
      active <<- c(active, pids)
      structure(rep(list(list()), cores), fake_pids = pids)
    },
    clusterCall = function(cl, fun) as.list(attr(cl, "fake_pids")),
    parLapplyLB = function(cl, X, fun, chunk.size) {
      if (!failed) {
        failed <<- TRUE
        busy <<- TRUE
        stop("error reading from connection")
      }
      busy <<- FALSE
      lapply(X, fun)
    }, .package = "parallel")
  local_mocked_bindings(pskill = function(pid, signal) {
    signals <<- c(signals, signal)
    expect_false(any(pid == Sys.getpid()))
    if (signal == 0L) return(pid %in% active)
    if (signal == tools::SIGKILL) active <<- setdiff(active, pid)
    rep(TRUE, length(pid))
  }, .package = "tools")
  expect_warning(
    result <- .latents_map_candidates(seq_len(3L), identity, cores = 2L,
                                       backend = "socket"),
    class = "latents_worker_lost")
  expect_identical(result, as.list(seq_len(3L)))
  expect_false(overlap)
  expect_identical(workers, c(2L, 2L, 1L))
  expect_true(all(workers > 0L & workers <= 2L))
  expect_true(tools::SIGTERM %in% signals)
  expect_true(tools::SIGKILL %in% signals)
  expect_length(active, 0L)
})

test_that("socket PID registration failures close idle workers before failing", {
  cleanups <- 0L
  local_mocked_bindings(.latents_stop_cluster = function(cluster) {
    cleanups <<- cleanups + 1L
    invisible(NULL)
  })
  local_mocked_bindings(
    makePSOCKcluster = function(cores) rep(list(list()), cores),
    clusterCall = function(...) stop("cannot read worker PIDs"),
    .package = "parallel")
  expect_error(.latents_map_candidates(1:2, identity, cores = 2L,
                                        backend = "socket"),
               class = "latents_parallel_failure")
  expect_identical(cleanups, 1L)
})

test_that("unretired socket workers prevent replacement and give a classed failure", {
  creations <- 0L
  seen <- NULL
  local_mocked_bindings(.latents_abort_cluster = function(cluster) FALSE)
  local_mocked_bindings(
    makePSOCKcluster = function(cores) {
      creations <<- creations + 1L
      list()
    },
    clusterCall = function(cl, fun) list(),
    parLapplyLB = function(...) stop("error reading from connection"),
    .package = "parallel")
  result <- tryCatch(withCallingHandlers(
    .latents_map_candidates(1:3, identity, cores = 2L, backend = "socket"),
    latents_worker_lost = function(condition) {
      seen <<- condition
      invokeRestart("muffleWarning")
    }), error = identity)
  expect_s3_class(result, "latents_parallel_failure")
  expect_s3_class(seen, "latents_worker_lost")
  expect_identical(creations, 1L)
  expect_match(conditionMessage(result), "could not be stopped safely", fixed = TRUE)
})

test_that("a seeded multilevel LCA grid prepared on workers is identical to one core", {
  skip_on_cran()
  skip_on_os("windows")
  data <- subset(student_esm, student <= 40)
  items <- c("time_with_friends", "on_social_media", "sports", "walking", "reading")
  lca_grid <- function(...) {
    enumerate_classes(data, items, "student", n_profiles = 2:3,
                      n_group_classes = 1:2, type = "lca", n_starts = 3,
                      seed = 1, ...)
  }
  one <- lca_grid()
  progress <- character()
  two <- withCallingHandlers(lca_grid(cores = 2),
    latents_progress = function(condition) {
      progress <<- c(progress, conditionMessage(condition))
      invokeRestart("muffleMessage")
    })
  expect_identical(strip_calls(two), strip_calls(one))
  expect_identical(as.data.frame(two), as.data.frame(one))
  # Several cores report each phase once by default; nothing per start.
  expect_length(progress, 3L)
  expect_match(progress[1L], "^Preparing 4 candidates on 2 workers\\.")
  expect_match(progress[2L], "^Fitting 12 starts on 2 workers\\.")
  expect_match(progress[3L], "^Done in ")
  expect_no_message(lca_grid(cores = 2, verbose = FALSE))
})
