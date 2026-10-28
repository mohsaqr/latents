Sys.setenv(NOT_CRAN = "true", VECLIB_MAXIMUM_THREADS = "1",
           OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
.libPaths(c("/private/tmp/latents-fixes-20261022/rlib", .libPaths()))
pkgload::load_all(quiet = TRUE)
suite_paths <- c(unit = "tests/testthat", equivalence = "tests/equivalence")
invisible(lapply(names(suite_paths), function(suite) {
  cat("\nRUNNING", suite, format(Sys.time()), "\n")
  elapsed <- system.time({
    results <- testthat::test_dir(suite_paths[[suite]], load_package = "none",
                                 reporter = "summary", stop_on_failure = FALSE)
  })
  saveRDS(results, file.path("/private/tmp/latents-fixes-20261022",
                            paste0(suite, "-results.rds")))
  frame <- as.data.frame(results)
  write.csv(frame[, setdiff(names(frame), "result")],
            file.path("/private/tmp/latents-fixes-20261022",
                      paste0(suite, "-results.csv")), row.names = FALSE)
  print(colSums(frame[, c("nb", "failed", "skipped", "error", "warning", "passed")]))
  print(elapsed)
  cat("FINISHED", suite, format(Sys.time()), "\n")
}))
