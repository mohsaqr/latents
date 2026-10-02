# Run from the repository root with NOT_CRAN=true Rscript validation/review-2026-10-02-followup/verify.R.
# --verify-saved verifies existing results against their original source manifest.
output <- file.path("validation", "review-2026-10-02-followup")
stopifnot(dir.exists(output), identical(Sys.getenv("NOT_CRAN"), "true"),
          requireNamespace("testthat", quietly = TRUE))
source_files <- sort(c(list.files("R", full.names = TRUE, pattern = "\\.R$"),
                       list.files("man", full.names = TRUE, pattern = "\\.Rd$"),
                       list.files("tests", full.names = TRUE, recursive = TRUE),
                       "DESCRIPTION", "NAMESPACE", "NEWS.md"))
# Plot devices may write this generated artifact while tests execute.
source_files <- source_files[basename(source_files) != "Rplots.pdf"]
source_hash <- tools::md5sum(source_files)
manifest_path <- file.path(output, "SOURCE.csv")
if ("--verify-saved" %in% commandArgs(trailingOnly = TRUE)) {
  manifest <- read.csv(manifest_path)
  manifest <- manifest[basename(manifest$path) != "Rplots.pdf", ]
  stopifnot(identical(manifest$path, source_files),
            identical(manifest$md5, unname(source_hash)))
  results <- readRDS(file.path(output, "tests-final.rds"))
} else {
  write.csv(data.frame(path = names(source_hash), md5 = unname(source_hash)),
            manifest_path, row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(output, "SESSION.txt"))
  results <- testthat::test_local(reporter = "summary", stop_on_failure = FALSE)
  saveRDS(results, file.path(output, "tests-final.rds"))
}
table <- as.data.frame(results)
write.csv(table[setdiff(names(table), "result")], file.path(output, "tests-final.csv"),
          row.names = FALSE)
summary <- colSums(table[c("failed", "warning", "skipped", "passed")])
print(summary)
cat("Test cases:", nrow(table), "\n")
cat("Errors:", sum(table$error), "\n")
stopifnot(all(summary[c("failed", "warning", "skipped")] == 0), !any(table$error))
stopifnot(identical(source_hash, tools::md5sum(source_files)))
# Strip historical generated-artifact entries only after checking all source hashes.
write.csv(data.frame(path = names(source_hash), md5 = unname(source_hash)),
          manifest_path, row.names = FALSE)
cat("Source identity verified:", length(source_hash), "files.\n")
