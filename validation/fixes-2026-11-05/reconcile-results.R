# Retain the original sweep and replace only complete-file outcomes from later
# reruns. Command-line paths are in chronological order; the newest file wins.
args <- commandArgs(TRUE)
if (length(args) == 1L && grepl("[.]txt$", args)) args <- readLines(args)
stopifnot(length(args) >= 2L, all(file.exists(args)))
library(testthat)
destination <- dirname(args[[1L]])
merged <- readRDS(args[[1L]])
raw <- as.data.frame(merged)
provenance <- setNames(rep(args[[1L]], length(unique(raw$file))), unique(raw$file))
invisible(lapply(args[-1L], function(path) {
  next_results <- readRDS(path)
  files <- unique(vapply(next_results, `[[`, "", "file"))
  stopifnot(all(files %in% names(provenance)))
  merged <<- c(Filter(function(case) !case$file %in% files, merged), next_results)
  provenance[files] <<- path
}))
class(merged) <- "testthat_results"
tab <- as.data.frame(merged)
saveRDS(merged, file.path(destination, "validated-results.rds"))
write.csv(tab[, setdiff(names(tab), "result"), drop = FALSE],
          file.path(destination, "validated-results.csv"), row.names = FALSE)
columns <- c("passed", "failed", "error", "warning", "skipped")
by_file <- aggregate(tab[, columns], list(file = tab$file), sum)
by_file$source_result <- unname(provenance[by_file$file])
by_file$cases <- as.integer(table(tab$file)[by_file$file])
write.csv(by_file, file.path(destination, "validated-by-file.csv"), row.names = FALSE)
cat("Original sweep:", length(unique(raw$file)), "files;", nrow(raw), "cases\n")
print(colSums(raw[, columns]))
cat("Complete sweep plus latest whole-file reruns:", nrow(by_file), "files;",
    nrow(tab), "cases\n")
print(colSums(tab[, columns]))
cat("\nWhole-file rerun provenance:\n")
print(by_file[by_file$source_result != args[[1L]], ], row.names = FALSE)
cat("\nRemaining conditions:\n")
invisible(lapply(merged, function(case) {
  conditions <- Filter(function(x) !inherits(x, "expectation_success"), case$results)
  if (length(conditions)) {
    cat("\n", case$file, ": ", case$test, "\n", sep = "")
    invisible(lapply(conditions, function(x) cat(paste(class(x), collapse = ", "),
      "\n", conditionMessage(x), "\n")))
  }
}))
stopifnot(all(colSums(tab[, c("failed", "error", "warning")]) == 0L))
