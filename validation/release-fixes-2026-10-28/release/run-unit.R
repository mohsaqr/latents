library(latents)
cat("Package:", as.character(packageVersion("latents")), "\n")
cat("Location:", find.package("latents"), "\n")
print(sessionInfo())
base <- "validation/release-fixes-2026-10-28/release"
started <- Sys.time()
results <- testthat::test_dir("tests/testthat", load_package = "installed",
                              package = "latents", reporter = "summary",
                              stop_on_failure = FALSE)
table <- as.data.frame(results)
saveRDS(results, file.path(base, "unit-results.rds"))
write.csv(table[setdiff(names(table), "result")], file.path(base, "unit-results.csv"),
          row.names = FALSE)
cat("\nElapsed:", as.numeric(difftime(Sys.time(), started, units = "secs")), "seconds\n")
print(colSums(table[c("nb", "failed", "skipped", "error", "warning", "passed")]))
cat("Files:", length(unique(table$file)), "test blocks:", nrow(table), "\n")
quit(status = as.integer(any(table$failed > 0 | table$error | table$warning > 0)))
