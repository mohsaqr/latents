library(testthat)
base <- "validation/release-fixes-2026-10-28/release"
read <- function(name) as.data.frame(readRDS(file.path(base, paste0(name, "-results.rds"))))
initial <- read("external-tests")
rerun <- rbind(read("nb-reference"), read("vignette-reference"))
stopifnot(all(unique(rerun$file) %in% unique(initial$file)))
combined <- rbind(initial[!initial$file %in% unique(rerun$file), ], rerun)
stopifnot(setequal(unique(combined$file), unique(initial$file)))
retained <- read("retained-references")
counts <- c("nb", "failed", "skipped", "error", "warning", "passed")
cat("Consolidated external tests (replace two corrected test files):\n")
print(colSums(combined[counts]))
cat("Retained references:\n")
print(colSums(retained[counts]))
write.csv(combined[setdiff(names(combined), "result")],
          file.path(base, "external-tests-consolidated.csv"), row.names = FALSE)
stopifnot(!any(combined$failed > 0 | combined$error),
           !any(retained$failed > 0 | retained$error))
