# compare_models(test = "lrt") against the references of make-kit.R.
# Run from the project root: Rscript equivalence/lrt/compare.R
# Writes equivalence/lrt/comparison.csv and prints it.
devtools::load_all(".", quiet = TRUE)
kit <- file.path("equivalence", "lrt")
source(file.path(kit, "common.R"))
result <- lrt_comparison(kit)
utils::write.csv(result, file.path(kit, "comparison.csv"), row.names = FALSE)
print(result, digits = 10)
cat(sprintf("\n%d of %d quantities agree\n", sum(result$agrees), nrow(result)))
