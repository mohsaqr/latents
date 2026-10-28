Sys.setenv(NOT_CRAN = "true", VECLIB_MAXIMUM_THREADS = "1",
           OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
.libPaths(c("/private/tmp/latents-review-20261022/rlib", .libPaths()))
pkgload::load_all(quiet = TRUE)
results <- testthat::test_dir("tests/testthat", load_package = "none",
                              filter = "^(parallel|search)$", reporter = "summary",
                              stop_on_failure = FALSE)
saveRDS(results, "/private/tmp/latents-review-20261022/socket-results.rds")
frame <- as.data.frame(results)
print(colSums(frame[, c("nb", "failed", "error", "warning", "skipped", "passed")]))
