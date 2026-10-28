arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L,
          arguments %in% c("unit", "equivalence", "references"))
kind <- arguments[[1L]]
destination <- "validation/release-0.96.1-2026-10-20"
source_files <- list.files("R", "[.]R$", full.names = TRUE)
before <- tools::md5sum(source_files)
started <- Sys.time()
commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
pkgload::load_all(".", quiet = TRUE)
cat("Commit:", commit, "\nPackage:", as.character(packageVersion("latents")), "\n")
folder <- switch(kind, unit = "tests/testthat", equivalence = "tests/equivalence",
                 references = "equivalence")
filter <- if (kind == "references") paste0("^(claims|latentgold|mclust(|-prior|-noise)|",
  "missing-patterns|mplus-(inference|robust|mixture-weights|lmr|missing-full)|",
  "numerical-equivalence)$") else NULL
results <- testthat::test_dir(folder, filter = filter, reporter = "summary",
  load_package = "none", stop_on_failure = FALSE)
saveRDS(results, file.path(destination, paste0(kind, "-final.rds")))
table <- as.data.frame(results)
write.csv(table[, setdiff(names(table), "result")],
  file.path(destination, paste0(kind, "-final.csv")), row.names = FALSE)
totals <- colSums(table[, c("nb", "failed", "error", "warning", "skipped", "passed")])
cat("\nTOTALS\n")
print(totals)
cat("Cases:", nrow(table), "\nElapsed:", as.numeric(difftime(Sys.time(), started,
  units = "secs")), "seconds\n")
unchanged <- identical(before, tools::md5sum(source_files))
cat("Production source unchanged:", unchanged, "\n")
write.csv(data.frame(file = names(before), md5 = unname(before)),
  file.path(destination, paste0(kind, "-final-source.csv")), row.names = FALSE)
print(sessionInfo())
if (!unchanged || any(totals[c("failed", "error", "warning", "skipped")] > 0)) {
  quit(status = 1L)
}
