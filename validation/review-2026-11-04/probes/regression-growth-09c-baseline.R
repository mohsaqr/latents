# In-memory mutations of table formulas; test files run from a scratchpad copy, so no snapshot is touched.
suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
# Baseline of 09b: the same files with no mutation.
files <- c("test-mixture-regression.R", "test-review-1026-mixreg.R", "test-growth-mixture.R", "test-review-1026-growth-r2.R")
invisible(lapply(files, function(f) {
  res <- testthat::test_file(file.path("rg-mut", "testthat", f), reporter = "silent", load_package = "none", stop_on_failure = FALSE)
  invisible(lapply(res, function(test) invisible(lapply(test$results, function(r) {
    if (inherits(r, c("expectation_failure", "expectation_error"))) cat(sprintf("[%s] %s\n    %s\n", f, test$test, substr(gsub("\\s+", " ", conditionMessage(r)), 1, 230)))
  }))))
  cat(f, "done\n")
}))
