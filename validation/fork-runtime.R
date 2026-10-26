# Runtime preflight for the retained fork-based validation protocols.
# This guard changes neither a requested core count nor the RNG schedule.

#' Check a validation study's requested fork runtime
#' @param cores Requested positive whole-number core count.
#' @param os_type,sysname,blas Runtime identifiers, injectable for pure checks.
#' @return The requested core count, unchanged; refuses unsafe multi-core runs.
.validation_fork_cores <- function(cores, os_type = .Platform$OS.type,
                                    sysname = Sys.info()[["sysname"]],
                                    blas = extSoftVersion()[["BLAS"]]) {
  stopifnot("cores must be a single positive whole number" =
              is.numeric(cores) && length(cores) == 1L && is.finite(cores) &&
                cores >= 1 && cores == floor(cores))
  unsafe <- identical(os_type, "windows") ||
    (identical(sysname, "Darwin") &&
       grepl("Accelerate|vecLib", paste(blas, collapse = " "), ignore.case = TRUE))
  if (cores > 1 && unsafe) {
    stop(errorCondition(paste(
      "This validation protocol uses forked mclapply workers, which are unsafe",
      "with Apple Accelerate/vecLib and unavailable on Windows.",
      "Use an R runtime with a fork-safe reference BLAS, or explicitly configure",
      "a one-core study. Treat a serial configuration as a separate run because",
      "its RNG schedule can differ from the retained fork-based results."),
      class = "latents_validation_unsafe_fork", call = NULL))
  }
  cores
}
