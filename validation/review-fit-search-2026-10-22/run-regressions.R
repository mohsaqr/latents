# Run from the repository root. Failures document the review snapshot.
Sys.setenv(VECLIB_MAXIMUM_THREADS = 1, OMP_NUM_THREADS = 1,
           OPENBLAS_NUM_THREADS = 1)
output <- file.path("validation", "review-fit-search-2026-10-22")
results <- testthat::test_file(file.path(output, "test-review-regressions.R"),
                              reporter = "summary", stop_on_failure = FALSE)
saveRDS(results, file.path(output, "regression-results.rds"))
table <- as.data.frame(results)
table$result <- NULL
utils::write.csv(table, file.path(output, "regression-results.csv"),
                 row.names = FALSE)
cat("\nReview totals:", nrow(table), "cases;", sum(table$passed),
    "supporting passes;", sum(table$failed), "expected failures;",
    sum(table$error), "errors;", sum(table$warning), "warnings.\n")
stopifnot(nrow(table) == 5L, sum(table$failed) == 6L,
          sum(table$passed) == 7L, !any(table$error), sum(table$warning) == 0L)
