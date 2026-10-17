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

#' Apply `fun` to every candidate, on one core or several
#'
#' With one core this is `lapply()`. With several, candidates are fitted by
#' forked processes (`parallel::mclapply()`) where the platform allows it,
#' and by a socket cluster on Windows. A worker's warnings and messages would
#' otherwise be lost or interleaved, so they are collected in the worker and
#' signalled again here, candidate by candidate in the order of the grid, as
#' the conditions themselves, classes included.
#'
#' @param indices The candidates' row numbers.
#' @param fun A function of one row number returning that candidate's result.
#'   It catches each candidate's own errors.
#' @param cores Number of cores, from `.latents_check_cores()`.
#' @param backend `"auto"`, or `"socket"` to use a socket cluster on any
#'   platform (used by the tests).
#' @return A list, one element per candidate, in the order of `indices`.
#' @noRd
.latents_map_candidates <- function(indices, fun, cores, backend = c("auto", "socket")) {
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
  socket <- identical(backend, "socket") || identical(.Platform$OS.type, "windows")
  results <- if (socket) {
    cluster <- parallel::makePSOCKcluster(cores)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::parLapply(cluster, indices, collected)
  } else {
    parallel::mclapply(indices, collected, mc.cores = cores,
                       mc.preschedule = FALSE)
  }
  # A worker that died (out of memory, killed) returns an error object rather
  # than a candidate; that is a failure of the run, not of the candidate.
  failed <- vapply(results, inherits, logical(1), "try-error")
  if (any(failed)) {
    stop(errorCondition(sprintf(paste(
      "%d of %d candidates were lost because their worker process failed",
      "(%s). Run again with fewer `cores`, or `cores = 1`."),
      sum(failed), length(results),
      trimws(as.character(results[[which(failed)[1L]]]))),
      class = "latents_parallel_failure", call = NULL))
  }
  invisible(lapply(results, function(result) {
    lapply(result$conditions, function(condition) {
      if (inherits(condition, "warning")) warning(condition) else message(condition)
    })
  }))
  lapply(results, `[[`, "value")
}
