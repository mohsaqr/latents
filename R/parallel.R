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
  as.integer(cores)
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

#' Apply `fun` to every candidate, on one core or several
#'
#' With one core this is `lapply()`. With several, candidates are fitted by
#' forked processes (`parallel::mclapply()`) where the platform allows it,
#' and by a socket cluster on Windows and on macOS with Accelerate's BLAS,
#' which forked processes cannot use (`.latents_fork_safe()`). A worker's warnings and messages would
#' otherwise be lost or interleaved, so they are collected in the worker and
#' signalled again here, candidate by candidate in the order of the grid, as
#' the conditions themselves, classes included.
#'
#' @param indices The candidates' row numbers.
#' @param fun A function of one row number returning that candidate's result.
#'   It catches each candidate's own errors.
#' @param cores Number of cores, from `.latents_check_cores()`.
#' @param backend `"auto"`; `"socket"` or `"fork"` to force one backend
#'   (used by the tests).
#' @return A list, one element per candidate, in the order of `indices`.
#' @noRd
.latents_map_candidates <- function(indices, fun, cores,
                                    backend = c("auto", "socket", "fork")) {
  backend <- match.arg(backend)
  if (cores <= 1L || length(indices) < 2L) return(lapply(indices, fun))
  cores <- min(cores, length(indices))
  collected <- function(i) {
    conditions <- list()
    value <- withCallingHandlers(fun(i),
      warning = function(condition) {
        conditions[[length(conditions) + 1L]] <<- condition
        invokeRestart("muffleWarning")
      },
      message = function(condition) {
        conditions[[length(conditions) + 1L]] <<- condition
        invokeRestart("muffleMessage")
      })
    list(value = value, conditions = conditions)
  }
  socket <- identical(backend, "socket") || (identical(backend, "auto") &&
    (identical(.Platform$OS.type, "windows") || !.latents_fork_safe()))
  results <- if (socket) {
    cluster <- parallel::makePSOCKcluster(cores)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    # A socket worker that dies breaks its connection, and parLapply() stops
    # with that error rather than returning the other candidates.
    tryCatch(parallel::parLapply(cluster, indices, collected),
      error = function(condition) {
        if (!grepl("connection", conditionMessage(condition), fixed = TRUE)) stop(condition)
        stop(errorCondition(sprintf(paste(
          "A worker process failed and the candidates were lost (%s).",
          "Run again with fewer `cores`, or `cores = 1`."),
          conditionMessage(condition)),
          class = "latents_parallel_failure", call = NULL))
      })
  } else {
    parallel::mclapply(indices, collected, mc.cores = cores,
                       mc.preschedule = FALSE)
  }
  # A worker that failed returns an error object, and one killed by a signal
  # (out of memory, a segmentation fault) returns NULL with only a warning;
  # `collected()` itself never returns NULL. Either is a failure of the run,
  # not of the candidate.
  lost <- vapply(results, is.null, logical(1))
  failed <- lost | vapply(results, inherits, logical(1), "try-error")
  if (any(failed)) {
    first <- which(failed)[1L]
    reason <- if (lost[first]) "the process ended without returning a result" else
      trimws(as.character(results[[first]]))
    stop(errorCondition(sprintf(paste(
      "%d of %d candidates were lost because their worker process failed",
      "(%s). Run again with fewer `cores`, or `cores = 1`."),
      sum(failed), length(results), reason),
      class = "latents_parallel_failure", call = NULL))
  }
  invisible(lapply(results, function(result) {
    lapply(result$conditions, function(condition) {
      if (inherits(condition, "warning")) warning(condition) else message(condition)
    })
  }))
  lapply(results, `[[`, "value")
}
