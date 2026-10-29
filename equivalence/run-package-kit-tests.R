# Run one complete comparison kit in CI, preserving outcomes even on failure.
# From the repository root:
#   Rscript equivalence/run-package-kit-tests.R pool-packages

#' Require completed comparisons, without silently accepting skipped tests
#' @param results The data frame returned by testthat.
#' @return Invisibly TRUE, or an error when any test failed, errored or skipped.
require_completed_kit_tests <- function(results) {
  fields <- c("passed", "failed", "error", "skipped")
  if (!is.data.frame(results) || !nrow(results) ||
      !all(fields %in% names(results)) || anyNA(results[fields])) {
    stop("The kit produced no complete test outcome table.", call. = FALSE)
  }
  if (any(results$failed > 0 | results$error | results$skipped) ||
      sum(results$passed) == 0L) {
    stop(sprintf("Kit did not complete: %d failed expectations, %d errors, %d skips.",
                 sum(results$failed), sum(results$error), sum(results$skipped)),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Run a named kit and retain its test outcomes and source identity
#' @param kit One of the comparison kits listed in the CI matrix.
#' @param output Directory for fresh evidence.
#' @return Invisibly the test outcome table; errors for incomplete comparisons.
run_package_kit_tests <- function(kit, output = file.path("tmp", "equivalence-kits", kit)) {
  kits <- c("lca-packages", "lpa-packages", "regression-packages", "lta-packages",
            "pool-packages", "enumeration-packages", "growth-published",
            "stata-published", "mplus-growth", "latentgold-coverage")
  stopifnot("select one registered comparison kit" =
              is.character(kit) && length(kit) == 1L && kit %in% kits)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  source("equivalence/verification-matrix/provenance.R", local = TRUE)
  snapshot <- matrix_capture_source(root = ".")
  Sys.setenv(NOT_CRAN = "true")
  pkgload::load_all(".", quiet = TRUE)
  results <- as.data.frame(testthat::test_dir(
    "equivalence", filter = paste0("^", kit, "$"), package = "latents",
    load_package = "none", reporter = "summary", stop_on_failure = FALSE))
  path <- file.path(output, "test-results.csv")
  utils::write.csv(results[setdiff(names(results), "result")], path, row.names = FALSE)
  matrix_record_evidence(path, snapshot, root = ".")
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(output, "SESSION.txt"))
  require_completed_kit_tests(results)
  invisible(results)
}

if (sys.nframe() == 0L) {
  arguments <- commandArgs(trailingOnly = TRUE)
  if (length(arguments) != 1L) stop("Supply exactly one comparison kit name.", call. = FALSE)
  run_package_kit_tests(arguments[[1L]])
}
