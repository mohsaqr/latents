# `R CMD check --as-cran` sets _R_CHECK_LIMIT_CORES_, and then parallel refuses
# to start more than two workers. Tests that start a real cluster of three
# skip there; they run in the ordinary suite.
skip_if_cores_limited <- function() {
  limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  skip_if(nzchar(limit) && !identical(tolower(limit), "false"),
          "_R_CHECK_LIMIT_CORES_ allows at most two workers")
}

# Live multicore assertions cannot verify dispatch when the supported fallback
# deliberately serializes every request on a platform without a core report.
skip_if_unknown_cores <- function() {
  available <- .latents_available_cores()
  skip_if(!is.finite(available) || available < 2L,
          "This platform cannot report at least two available processors")
}
