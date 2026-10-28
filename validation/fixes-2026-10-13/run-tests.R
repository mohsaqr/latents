Sys.setenv(NOT_CRAN = "true")
out <- "validation/fixes-2026-10-13"
pkgload::load_all(quiet = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out, "SESSION.txt"))
run_suite <- function(path, label) {
  cat("Starting", label, format(Sys.time()), "\n")
  result <- testthat::test_dir(path, load_package = "none", reporter = "summary",
                               stop_on_failure = FALSE, stop_on_warning = FALSE)
  frame <- as.data.frame(result)
  saveRDS(result, file.path(out, paste0(label, ".rds")))
  frame$result <- vapply(frame$result, function(x) paste(capture.output(str(x)), collapse = " "), character(1))
  utils::write.csv(frame, file.path(out, paste0(label, ".csv")), row.names = FALSE)
  cat("Completed", label, format(Sys.time()), "\n")
  cat("Cases:", nrow(frame), "passed:", sum(frame$passed), "failed:", sum(frame$failed),
       "errors:", sum(frame$error), "warnings:", sum(frame$warning),
       "skips:", sum(frame$skipped), "\n")
  invisible(frame)
}
unit <- run_suite("tests/testthat", "unit-tests")
independent <- run_suite("tests/equivalence", "equivalence-tests")
stopifnot(sum(unit$failed) == 0, !any(unit$error), sum(unit$warning) == 0,
          sum(independent$failed) == 0, !any(independent$error), sum(independent$warning) == 0)
