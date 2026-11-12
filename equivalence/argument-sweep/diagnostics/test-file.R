# Run test files of tests/testthat against the working tree, from the project root:
#   Rscript equivalence/argument-sweep/diagnostics/test-file.R fixes-1111 three-step
Sys.setenv(NOT_CRAN = "true")
filter <- paste(commandArgs(trailingOnly = TRUE), collapse = "|")
results <- as.data.frame(testthat::test_local(".", filter = filter, reporter = "check",
                                              stop_on_failure = FALSE))
cat(sprintf("\n%d tests, %d expectations: %d failed, %d errors, %d skipped, %d warnings\n",
            nrow(results), sum(results$nb), sum(results$failed), sum(results$error),
            sum(results$skipped), sum(results$warning)))
if (sum(results$failed) + sum(results$error) > 0L) quit(status = 1L)
