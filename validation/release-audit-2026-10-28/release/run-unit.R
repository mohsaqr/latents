Sys.setenv(NOT_CRAN = "true")
.libPaths(c("/private/tmp/latents-release-audit-20261028/lib", .libPaths()))
library(latents)
cat("Package:", as.character(packageVersion("latents")), "\n")
cat("Location:", find.package("latents"), "\n")
print(sessionInfo())
started <- Sys.time()
results <- testthat::test_dir("tests/testthat", load_package = "installed",
                              package = "latents", reporter = "summary",
                              stop_on_failure = FALSE)
table <- as.data.frame(results)
saveRDS(results, "validation/release-audit-2026-10-28/release/unit-results.rds")
write.csv(table[setdiff(names(table), "result")],
          "validation/release-audit-2026-10-28/release/unit-results.csv", row.names = FALSE)
cat("\nElapsed:", as.numeric(difftime(Sys.time(), started, units = "secs")), "seconds\n")
print(colSums(table[intersect(c("nb", "failed", "skipped", "error", "warning", "passed"), names(table))]))
quit(status = as.integer(any(table$failed > 0 | table$error)))
