library(latents)
label <- commandArgs(trailingOnly = TRUE)[[1L]]
cat("Library:", find.package("latents"), "\n")
set.seed(1)
indicators <- c("cognitive_strategies", "intrinsic_value", "self_efficacy",
                "self_regulation", "test_anxiety")
str(srl[indicators])
grid <- enumerate_lpa(srl, indicators, n_profiles = 1:6)
table <- get_results(grid, "candidates")
print(table[c("model", "n_profiles", "bic_individual", "converged", "boundary")])
base <- "validation/release-fixes-2026-10-28/release"
write.csv(table, file.path(base, paste0("workflow-grid-", label, ".csv")), row.names = FALSE)
saveRDS(grid, file.path(base, paste0("workflow-grid-", label, ".rds")))
