# Fitting the candidates of an enumeration on several cores.
#
# Each candidate is an independent fit that sets and restores its own seed
# (.latents_local_seed()), so it gives the same result whichever process fits
# it. Running candidates in parallel therefore changes the time an
# enumeration takes and nothing it returns.

#' Validate `cores`
#' @param cores The caller's `cores`.
#' @return The number of cores as an integer.
#' @noRd
.latents_check_cores <- function(cores) {
  if (!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) ||
      cores < 1 || cores != as.integer(cores)) {
    stop(errorCondition("`cores` must be a single whole number of at least 1.",
                        class = "latents_bad_argument", call = NULL))
  }
  cores <- as.integer(cores)
  # More workers than the machine has cores gains nothing and, as socket
  # workers (Windows, or macOS with Accelerate), each is a full R process
  # holding the package and a copy of the data: `cores = 24` on a ten-core
  # Mac started 24 of them and exhausted memory (2026-10-26).
  available <- tryCatch(parallel::detectCores(), error = function(condition) NA_integer_)
  if (is.finite(available) && available >= 1L && cores > available) {
    warning(warningCondition(sprintf(paste(
      "`cores = %d` exceeds the %d cores this machine has; %d workers are",
      "used. Each socket worker is a separate R process."),
      cores, available, available), class = "latents_too_many_cores",
      call = NULL))
    cores <- as.integer(available)
  }
  cores
}

#' Whether forked workers can use this R's linear algebra
#'
#' Apple's Accelerate (vecLib) BLAS is not safe in a forked process once the
#' parent has used it: the children die with a segmentation fault in the
#' first BLAS call (`backsolve()`, `crossprod()`). R on macOS can be switched
#' to it from the reference BLAS, so the fork is refused only where it is.
#'
#' @param blas Path of the BLAS library, as `extSoftVersion()` reports it.
#' @param sysname Operating system name, as `Sys.info()` reports it.
#' @return `TRUE` when forked workers are safe.
#' @noRd
.latents_fork_safe <- function(blas = extSoftVersion()[["BLAS"]],
                               sysname = Sys.info()[["sysname"]]) {
  !(identical(sysname, "Darwin") && grepl("Accelerate|vecLib", blas))
}

#' Stop every socket worker and close every owned connection
#'
#' A dead worker can make stopCluster fail while sending DONE. Clean nodes
#' independently and always close their connections so a cleanup failure
#' cannot replace the original worker condition or leave sockets for GC.
#' @param cluster A socket cluster created by makePSOCKcluster().
#' @return NULL, invisibly.
#' @noRd
.latents_stop_cluster <- function(cluster) {
  stopifnot("cluster is a list" = is.list(cluster))
  invisible(lapply(cluster, function(node) {
    single <- structure(list(node), class = class(cluster))
    tryCatch(parallel::stopCluster(single),
      error = function(condition) NULL,
      finally = tryCatch(close(node$con), error = function(condition) NULL))
  }))
  invisible(NULL)
}

#' Collect warnings and messages without changing their classes
#'
#' @param fun A zero-argument callback whose conditions should be collected.
#' @return A list with the callback's `value` and its ordered `conditions`.
#' @noRd
.latents_collect_conditions <- function(fun) {
  stopifnot("fun is a function" = is.function(fun))
  force(fun)
  conditions <- list()
  value <- withCallingHandlers(fun(),
    warning = function(condition) {
      # Muffling must not disable warning-to-error conversion inside the
      # candidate's own error handler (including on a socket worker).
      # Returning lets R perform its normal conversion at the signal site;
      # throwing here would bypass handlers nested inside this collector.
      if (isTRUE(getOption("warn") >= 2L)) return(invisible(NULL))
      conditions[[length(conditions) + 1L]] <<- condition
      invokeRestart("muffleWarning")
    },
    message = function(condition) {
      conditions[[length(conditions) + 1L]] <<- condition
      invokeRestart("muffleMessage")
    })
  list(value = value, conditions = conditions)
}

#' Give each worker task the parent's RNG kinds and warning policy
#'
#' Whole-candidate tasks seed themselves inside the fitting engine. Fresh
#' socket workers otherwise seed with their own default RNG kinds. Restore
#' the kinds and .Random.seed after each task, including errors and interrupts.
#' As with .latents_local_seed(), base R does not expose the Box-Muller cache
#' for restoration; set.seed() resets it before each seeded candidate.
#' @param fun A function of one task.
#' @return A task function with the parent's execution settings captured.
#' @noRd
.latents_worker_context <- function(fun) {
  stopifnot("fun is a function" = is.function(fun))
  force(fun)
  rng_kind <- RNGkind()
  warn <- getOption("warn")
  function(task) {
    old_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = globalenv())
    old_warn <- getOption("warn")
    on.exit({
      suppressWarnings(do.call(RNGkind, as.list(old_kind))) # nolint: undesirable_function_linter, object_name_linter.
      if (had_seed) assign(".Random.seed", old_seed, envir = globalenv()) else
        if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
          rm(".Random.seed", envir = globalenv())
        }
      options(warn = old_warn)
    }, add = TRUE)
    # Rounding is a supported sample kind but RNGkind warns when restoring
    # it. The parent's already-selected kind needs no new qualification.
    suppressWarnings(do.call(RNGkind, as.list(rng_kind))) # nolint: undesirable_function_linter.
    options(warn = warn)
    fun(task)
  }
}

#' Wrap a task function without capturing its scheduler
#'
#' @param fun A function of one task.
#' @return A function returning each task's value and collected conditions.
#' @noRd
.latents_task_collector <- function(fun) {
  stopifnot("fun is a function" = is.function(fun))
  force(fun)
  function(i) .latents_collect_conditions(function() {
    tryCatch(fun(i), error = function(condition) {
      structure(list(condition = .latents_bounded_condition(condition)),
                class = "latents_task_error")
    })
  })
}

#' Signal collected warnings and messages in their original order
#'
#' @param conditions A list of warning and message conditions.
#' @return NULL, invisibly.
#' @noRd
.latents_signal_conditions <- function(conditions) {
  stopifnot("conditions is a list" = is.list(conditions))
  invisible(lapply(conditions, function(condition) {
    if (inherits(condition, "warning")) warning(condition) else message(condition)
  }))
  invisible(NULL)
}

#' Run one fork pass, including isolated one-worker retries
#' @param indices Forced task payloads in dispatch order.
#' @param collected The task collector.
#' @param cores Worker limit for this pass.
#' @param retry Whether this is the retry, so a spawn failure remains a loss.
#' @return Collected task results; missing worker results remain NULL/try-error.
#' @noRd
.latents_fork_pass <- function(indices, collected, cores, retry = FALSE) {
  withCallingHandlers({
    if (cores > 1L) {
      tryCatch(parallel::mclapply(indices, collected, mc.cores = cores,
                                  mc.preschedule = FALSE),
        error = function(condition) {
          if (!retry) stop(condition)
          rep(list(NULL), length(indices))
        })
    } else {
      # mclapply(mc.cores = 1) runs in the parent. A retry must stay isolated
      # even with one worker, because the task may kill its process again.
      unavailable <- FALSE
      run_task <- function(task) {
        job <- NULL
        finished <- FALSE
        on.exit(if (!finished && !is.null(job)) suspendInterrupts({
          tryCatch(tools::pskill(job$pid, tools::SIGKILL),
            error = function(condition) NULL)
          # A killed child need not deliver a result during cleanup.
          suppressWarnings(tryCatch(parallel::mccollect(job, wait = FALSE), # nolint: undesirable_function_linter.
            error = function(condition) NULL,
            interrupt = function(condition) NULL))
        }), add = TRUE)
        job <- parallel::mcparallel(collected(task), mc.set.seed = TRUE)
        result <- parallel::mccollect(job)
        finished <- TRUE
        if (is.null(result)) NULL else result[[1L]]
      }
      lapply(indices, function(task) {
        if (unavailable) return(NULL)
        tryCatch(run_task(task), error = function(condition) {
          if (!retry) stop(condition)
          unavailable <<- TRUE
          NULL
        })
      })
    }
  }, warning = function(condition) {
    if (grepl("did not deliver a result|scheduled cores.*error|parallel function calls?.*error",
              conditionMessage(condition))) invokeRestart("muffleWarning")
  })
}

#' Check Windows worker liveness without sending a signal
#'
#' Signal zero also calls TerminateProcess on Windows. A successful tasklist
#' snapshot instead identifies exited workers when TerminateProcess fails.
#' Failed or malformed snapshots leave liveness unknown, never retired.
#' @param pids Owned worker process IDs.
#' @return A logical vector, with NA when liveness cannot be checked.
#' @noRd
.latents_windows_alive <- function(pids) {
  stopifnot("pids are numeric" = is.numeric(pids))
  unknown <- rep(NA, length(pids))
  # Nonzero exit status and malformed PIDs are diagnosed as unknown below.
  listing <- suppressWarnings(tryCatch( # nolint: undesirable_function_linter.
    system2("tasklist", c("/FO", "CSV", "/NH"), stdout = TRUE, stderr = FALSE),
    error = function(condition) NULL))
  if (!length(listing) || !is.null(attr(listing, "status")) &&
      attr(listing, "status") != 0L) return(unknown)
  processes <- tryCatch(utils::read.csv(text = listing, header = FALSE,
    colClasses = "character", stringsAsFactors = FALSE),
    error = function(condition) NULL)
  if (is.null(processes) || ncol(processes) < 2L) return(unknown)
  running <- suppressWarnings(as.integer(processes[[2L]])) # nolint: undesirable_function_linter.
  if (anyNA(running) || !length(running)) return(unknown)
  pids %in% running
}

#' Start socket workers and record the processes owned by this pool
#' @param cores Number of workers.
#' @return A socket cluster with its worker process IDs attached.
#' @noRd
.latents_make_cluster <- function(cores) {
  cluster <- parallel::makePSOCKcluster(cores)
  registered <- FALSE
  on.exit(if (!registered) .latents_stop_cluster(cluster), add = TRUE)
  pids <- unlist(parallel::clusterCall(cluster, Sys.getpid), use.names = FALSE)
  if (is.null(pids)) pids <- integer()
  stopifnot("one positive worker PID per socket node" =
              is.numeric(pids) && length(pids) == length(cluster) &&
              all(is.finite(pids) & pids > 0 & pids == as.integer(pids)))
  attr(cluster, "latents_pids") <- as.integer(pids[pids != Sys.getpid()])
  registered <- TRUE
  cluster
}

#' Retire owned socket workers before replacing a broken cluster
#'
#' stopCluster queues DONE behind an EM task already running on a surviving
#' worker. Terminate those processes first, so replacement workers cannot
#' overlap the abandoned computations. A zombie may still have a PID after
#' SIGKILL, but cannot execute another task. On Windows pskill uses the
#' TerminateProcess API for every signal. If termination cannot be confirmed,
#' callers raise a classed infrastructure error without starting replacements.
#' @param cluster A cluster whose owned process IDs have been registered.
#' @return Whether every owned worker has exited or received a fatal signal.
#' @noRd
.latents_abort_cluster <- function(cluster) {
  pids <- attr(cluster, "latents_pids")
  if (is.null(pids)) {
    .latents_stop_cluster(cluster)
    return(FALSE)
  }
  pids <- pids[pids != Sys.getpid()]
  if (!length(pids)) {
    .latents_stop_cluster(cluster)
    return(TRUE)
  }
  signal_workers <- function(pids, signal) {
    tryCatch(tools::pskill(pids, signal),
      error = function(condition) rep(FALSE, length(pids)))
  }
  terminated <- signal_workers(pids, tools::SIGTERM)
  .latents_stop_cluster(cluster)
  if (identical(.Platform$OS.type, "windows")) {
    if (all(terminated)) return(TRUE)
    alive <- .latents_windows_alive(pids[!terminated])
    return(all(!is.na(alive) & !alive))
  }
  # The short pause lets graceful termination release memory first. SIGKILL
  # then retires any busy survivor; a zombie must not block recovery forever.
  Sys.sleep(0.05)
  alive <- tryCatch(tools::pskill(pids, 0L),
    error = function(condition) rep(NA, length(pids)))
  pending <- is.na(alive) | alive
  if (!any(pending)) return(TRUE)
  killed <- signal_workers(pids[pending], tools::SIGKILL)
  if (all(killed)) return(TRUE)
  remaining <- pids[pending][!killed]
  alive <- tryCatch(tools::pskill(remaining, 0L),
    error = function(condition) rep(NA, length(remaining)))
  all(!is.na(alive) & !alive)
}

#' Run bounded socket chunks while retaining completed chunks
#'
#' The exported parallel API returns no partial results when parLapplyLB loses
#' a connection. Bound each call to one task per worker: a death loses at most
#' that chunk, while completed chunks remain available. Restart the cluster
#' before the next chunk; even a repeatedly dying task cannot discard later
#' chunks. A retry uses one-task calls from the outset: every task is attempted
#' at most twice, and a repeatedly dying task cannot lose its retry's mates.
#' @param indices Forced task payloads in dispatch order.
#' @param collected The task collector.
#' @param cores Worker limit for this pass.
#' @param retry Whether this is a retry whose startup failure remains a loss.
#' @return Collected task results; an unfinished chunk contains NULL results.
#' @noRd
.latents_socket_pass <- function(indices, collected, cores, retry = FALSE) {
  if (retry) cores <- 1L
  cluster <- NULL
  attempted <- FALSE
  unavailable <- FALSE
  finished <- FALSE
  on.exit(if (!is.null(cluster)) suspendInterrupts({
    if (finished) .latents_stop_cluster(cluster) else {
      tryCatch(.latents_abort_cluster(cluster), error = function(condition) NULL)
    }
  }), add = TRUE)
  chunks <- split(seq_along(indices), ceiling(seq_along(indices) / cores))
  results <- lapply(chunks, function(positions) {
    if (unavailable) return(rep(list(NULL), length(positions)))
    if (is.null(cluster)) {
      cluster <<- tryCatch(.latents_make_cluster(cores),
        error = function(condition) {
          if (retry || attempted) {
            unavailable <<- TRUE
            return(NULL)
          }
          stop(errorCondition(sprintf("Could not start the worker cluster (%s).",
            conditionMessage(condition)), class = "latents_parallel_failure",
            call = NULL))
        })
      if (is.null(cluster)) return(rep(list(NULL), length(positions)))
    }
    attempted <<- TRUE
    tryCatch(parallel::parLapplyLB(cluster, indices[positions], collected,
                                  chunk.size = 1L),
      error = function(condition) {
        if (!grepl("connection|unserialize|serialize|socket|node.*(died|terminated)",
                   conditionMessage(condition), ignore.case = TRUE)) stop(condition)
        retired <- .latents_abort_cluster(cluster)
        cluster <<- NULL
        if (!retired) {
          .latents_warn_worker_lost(length(positions), 0L, 1L)
          stop(errorCondition(paste(
            "Owned worker processes could not be stopped safely.",
            "No replacement workers were started; run again with fewer `cores`."),
            class = "latents_parallel_failure", call = NULL))
        }
        rep(list(NULL), length(positions))
      })
  })
  results <- unname(do.call(c, results))
  finished <- TRUE
  results
}

#' Report worker losses from one pool or several phases
#' @param lost Number of tasks whose first worker result was missing.
#' @param recovered Number recovered by the retry.
#' @param retry_cores Worker limit used for the retry.
#' @return NULL, invisibly.
#' @noRd
.latents_warn_worker_lost <- function(lost, recovered, retry_cores) {
  remaining <- lost - recovered
  warning(warningCondition(sprintf(paste(
    "%d worker tasks were lost; the retry recovered %d and %d remain failed.",
    "Worker processes may have been killed for memory. Try fewer `cores`."),
    lost, recovered, remaining), class = "latents_worker_lost",
    lost = as.integer(lost), recovered = as.integer(recovered),
    remaining = as.integer(remaining), retry_cores = as.integer(retry_cores),
    call = NULL))
  invisible(NULL)
}

#' Apply `fun` to candidates or starts in a shared worker pool
#'
#' With several cores, tasks use dynamically assigned forked processes where
#' the platform allows them, or a load-balanced socket cluster on Windows and
#' on macOS with Accelerate's BLAS (`.latents_fork_safe()`). Optional costs send
#' the longest estimated tasks first. Results and conditions retain the order
#' of `indices`, independently of dispatch or completion order. Conditions
#' retain their classes and can instead be returned for later signalling.
#'
#' Missing worker results are retried once using half as many fork workers,
#' rounded up, or one socket worker to isolate repeated losses. A classed
#' warning reports process losses and recovery. Socket calls use
#' bounded chunks because the exported parallel API cannot retain results from
#' within a call whose connection breaks; previously completed chunks survive.
#' @param indices Task identifiers.
#' @param fun A function of one identifier returning that task's result.
#'   It catches each candidate's own errors.
#' @param cores Number of cores, from `.latents_check_cores()`.
#' @param backend `"auto"`; `"socket"` or `"fork"` to force one backend
#'   (used by the tests).
#' @param costs NULL or positive finite estimated costs aligned with `indices`.
#' @param signal Whether to signal collected conditions in the parent.
#' @param on_lost NULL or a function of a task and a message, returning its
#'   terminal failure value when its retry also loses its worker.
#' @return A list in the order of `indices`. With `signal = FALSE`, each element
#'   contains `value` and `conditions`; otherwise each is the task's value.
#' @noRd
.latents_map_candidates <- function(indices, fun, cores,
                                    backend = c("auto", "socket", "fork"),
                                    costs = NULL, signal = TRUE,
                                    on_lost = NULL) {
  # Socket workers must receive values, never promises naming caller globals.
  force(indices)
  force(fun)
  force(cores)
  force(costs)
  force(signal)
  force(on_lost)
  backend <- match.arg(backend)
  stopifnot("signal is one non-missing logical value" =
              is.logical(signal) && length(signal) == 1L && !is.na(signal),
            "on_lost is NULL or a function" =
              is.null(on_lost) || is.function(on_lost))
  if (!is.null(costs)) {
    stopifnot("one positive finite cost per task" =
                is.numeric(costs) && length(costs) == length(indices) &&
                all(is.finite(costs) & costs > 0))
  }
  if (length(indices) == 0L) return(list())
  if (signal && is.null(costs) && cores <= 1L) {
    return(lapply(indices, fun))
  }
  positions <- if (is.null(costs)) seq_along(indices) else
    order(-costs, seq_along(indices))
  collected <- .latents_task_collector(if (cores > 1L)
    .latents_worker_context(fun) else fun)
  socket <- identical(backend, "socket") || (identical(backend, "auto") &&
    (identical(.Platform$OS.type, "windows") || !.latents_fork_safe()))
  workers <- min(cores, length(indices))
  run_pass <- if (socket) .latents_socket_pass else .latents_fork_pass
  results <- if (cores <= 1L) {
    lapply(indices[positions], collected)
  } else {
    run_pass(indices[positions], collected, workers)
  }
  # The collector always returns a list, including when fun returns NULL.
  # Only missing processes (NULL/try-error) are eligible for another attempt.
  is_lost <- function(result) is.null(result) || inherits(result, "try-error")
  lost <- vapply(results, is_lost, logical(1))
  if (any(lost)) {
    retry_cores <- if (socket) 1L else as.integer(ceiling(workers / 2))
    retried <- run_pass(indices[positions[lost]], collected,
                        min(retry_cores, sum(lost)), retry = TRUE)
    results[lost] <- retried
    still_lost <- vapply(results, is_lost, logical(1))
    recovered <- sum(lost) - sum(still_lost)
    .latents_warn_worker_lost(sum(lost), recovered, retry_cores)
    if (any(still_lost)) {
      loss_message <- paste("The worker process ended without returning a result,",
                            "likely killed for memory. Try fewer `cores`.")
      results[still_lost] <- lapply(indices[positions[still_lost]], function(task) {
        value <- if (is.null(on_lost)) list(error = loss_message) else
          on_lost(task, loss_message)
        list(value = value, conditions = list())
      })
    }
  }
  results <- results[order(positions)]
  task_errors <- vapply(results, function(result) {
    inherits(result$value, "latents_task_error")
  }, logical(1))
  if (any(task_errors)) {
    # There will be no return value for deferred signalling. Deliver the
    # completed tasks' qualifications before propagating the original error.
    invisible(lapply(results, function(result) {
      .latents_signal_conditions(result$conditions)
    }))
    stop(results[[which(task_errors)[1L]]]$value$condition)
  }
  if (!signal) return(results)
  invisible(lapply(results, function(result) {
    .latents_signal_conditions(result$conditions)
  }))
  lapply(results, `[[`, "value")
}
