arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L, requireNamespace("testthat", quietly = TRUE))
output <- arguments[[1L]]
raw <- readRDS(file.path(output, "results.rds"))
cases <- as.data.frame(raw)
cat("One complete run, no result replacement or aggregation:\n")
print(colSums(cases[c("nb", "failed", "error", "warning", "skipped", "passed")]))
cat("Cases:", nrow(cases), "\n")
cat("Test files:", length(unique(cases$file)), "\n")
cat("Snapshot commit:", readLines(file.path(output, "commit.txt")), "\n")
issues <- do.call(rbind, lapply(raw, function(case) {
  failures <- Filter(function(result) !inherits(result, "expectation_success"),
                     case$results)
  if (!length(failures)) return(NULL)
  data.frame(file = case$file, test = case$test,
    type = vapply(failures, function(issue) class(issue)[[1L]], character(1)),
    message = vapply(failures, conditionMessage, character(1)))
}))
if (is.null(issues)) issues <- data.frame(file = character(), test = character(),
                                         type = character(), message = character())
utils::write.csv(issues, file.path(output, "issues.csv"), row.names = FALSE)
stopifnot(!any(cases$failed | cases$error | cases$warning | cases$skipped))
cat("PASS: every case passed without failures, errors, warnings or skips.\n")
