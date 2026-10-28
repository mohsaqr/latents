Sys.setenv(NOT_CRAN = "true", VECLIB_MAXIMUM_THREADS = "1",
           OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
output <- file.path("validation", "spike-fixes-2026-10-22")
pkgload::load_all(quiet = TRUE)
suite_paths <- c(unit = "tests/testthat", equivalence = "tests/equivalence")
invisible(lapply(names(suite_paths), function(suite) {
  cat("\nRUNNING", suite, format(Sys.time()), "\n")
  elapsed <- system.time({
    results <- testthat::test_dir(suite_paths[[suite]], load_package = "none",
                                  reporter = "summary", stop_on_failure = FALSE)
  })
  saveRDS(results, file.path(output, paste0(suite, "-results.rds")))
  frame <- as.data.frame(results)
  utils::write.csv(frame[, setdiff(names(frame), "result")],
                    file.path(output, paste0(suite, "-results.csv")),
                    row.names = FALSE)
  print(colSums(frame[, c("nb", "failed", "skipped", "error", "warning", "passed")]))
  print(elapsed)
  cat("FINISHED", suite, format(Sys.time()), "\n")
}))
