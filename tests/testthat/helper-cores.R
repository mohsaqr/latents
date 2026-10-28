# `R CMD check --as-cran` sets _R_CHECK_LIMIT_CORES_, and then parallel refuses
# to start more than two workers. Tests that start a real cluster of three
# skip there; they run in the ordinary suite.
skip_if_cores_limited <- function() {
  limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  skip_if(nzchar(limit) && !identical(tolower(limit), "false"),
          "_R_CHECK_LIMIT_CORES_ allows at most two workers")
}
