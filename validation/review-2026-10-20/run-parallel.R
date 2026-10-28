pkgload::load_all(".", quiet = TRUE)
results <- testthat::test_dir("tests/testthat", filter = "^parallel$",
  reporter = "summary", load_package = "none", stop_on_failure = FALSE)
saveRDS(results, "validation/review-2026-10-20/parallel-tests.rds")
frame <- as.data.frame(results)
write.csv(frame[, setdiff(names(frame), "result")],
  "validation/review-2026-10-20/parallel-tests.csv", row.names = FALSE)
cat("\nTOTALS\n")
print(colSums(frame[, c("nb", "failed", "error", "warning", "skipped", "passed")]))
