library(testthat)
directory <- "validation/review-2026-10-20"
initial <- as.data.frame(readRDS(file.path(directory, "unit-tests.rds")))
development <- as.data.frame(readRDS(file.path(directory, "development-unit-tests.rds")))
combined <- rbind(initial[!initial$file %in% unique(development$file), ], development)
parallel_path <- file.path(directory, "parallel-tests.rds")
if (file.exists(parallel_path)) {
  parallel <- as.data.frame(readRDS(parallel_path))
  combined <- rbind(combined[!combined$file %in% unique(parallel$file), ], parallel)
}
equivalence <- as.data.frame(readRDS(file.path(directory, "equivalence-full.rds")))
columns <- c("nb", "failed", "error", "warning", "skipped", "passed")
summary <- rbind(
  data.frame(suite = "unit combined without duplicate files", cases = nrow(combined),
             as.list(colSums(combined[, columns]))),
  data.frame(suite = "development equivalence", cases = nrow(equivalence),
             as.list(colSums(equivalence[, columns]))))
print(summary, row.names = FALSE)
write.csv(summary, file.path(directory, "validation-summary.csv"), row.names = FALSE)
write.csv(combined[, setdiff(names(combined), "result")],
           file.path(directory, "unit-combined.csv"), row.names = FALSE)
manifest <- read.csv(file.path(directory, "source-manifest.csv"))
stopifnot(identical(unname(tools::md5sum(manifest$file)), manifest$md5))
cat("Production source fingerprints unchanged.\n")
cat("Equivalence failures:\n")
print(equivalence[equivalence$failed > 0 | equivalence$error, setdiff(names(equivalence), "result")])
