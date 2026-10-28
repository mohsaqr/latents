arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L)
output <- arguments[[1L]]
stopifnot(requireNamespace("testthat", quietly = TRUE))
raw_initial <- readRDS(file.path(output, "results.rds"))
initial <- as.data.frame(raw_initial)
initial_issues <- do.call(rbind, lapply(seq_along(raw_initial), function(index) {
  results <- raw_initial[[index]]$results
  issues <- Filter(function(result) !inherits(result, "expectation_success"), results)
  if (!length(issues)) return(NULL)
  data.frame(file = raw_initial[[index]]$file, test = raw_initial[[index]]$test,
    type = vapply(issues, function(issue) class(issue)[[1L]], character(1)),
    message = vapply(issues, conditionMessage, character(1)))
}))
utils::write.csv(initial_issues, file.path(output, "initial-issues.csv"), row.names = FALSE)
replacements <- c("claims-corrected.rds", "enumeration-current-library.rds")
reruns <- do.call(rbind, lapply(file.path(output, replacements), function(path) {
  as.data.frame(readRDS(path))
}))
final <- rbind(initial[!initial$file %in% reruns$file, , drop = FALSE], reruns)
saveRDS(final, file.path(output, "final-cases.rds"))
utils::write.csv(final[setdiff(names(final), "result")],
                 file.path(output, "final-cases.csv"), row.names = FALSE)
cat("Original full-suite totals:\n")
print(colSums(initial[c("nb", "failed", "skipped", "error", "warning", "passed")]))
cat("\nFinal totals after whole-file replacements:\n")
print(colSums(final[c("nb", "failed", "skipped", "error", "warning", "passed")]))
cat("Cases:", nrow(final), "\n")
cat("Replaced files:", paste(unique(reruns$file), collapse = ", "), "\n")
print(final[final$failed > 0 | final$error | final$warning > 0 | final$skipped > 0,
            setdiff(names(final), "result")])
quit(status = as.integer(any(final$failed > 0 | final$error)))
