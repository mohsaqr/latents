# Replace the corrected step-routing test file without double-counting its rerun.
library(testthat)
base <- "validation/release-fixes-2026-10-28/release"
initial <- as.data.frame(readRDS(file.path(base, "unit-results.rds")))
rerun <- as.data.frame(readRDS(file.path(base, "engine-fixes-results.rds")))
stopifnot(identical(unique(rerun$file), "test-engine-fixes.R"))
combined <- rbind(initial[!initial$file %in% unique(rerun$file), ], rerun)
stopifnot(setequal(unique(combined$file), unique(initial$file)))
cat("Consolidated package tests (replace corrected engine-fixes test file):\n")
print(colSums(combined[c("nb", "failed", "skipped", "error", "warning", "passed")]))
cat("Files:", length(unique(combined$file)), "test blocks:", nrow(combined), "\n")
write.csv(combined[setdiff(names(combined), "result")],
          file.path(base, "unit-consolidated.csv"), row.names = FALSE)
stopifnot(!any(combined$failed > 0 | combined$error | combined$warning > 0))
