# Compare three_step() on a categorical distal outcome with Latent GOLD 6.1's
# listings in returned/. Run from the project root after make-kit.R and
# run-kit.sh; writes comparison.csv.
pkgload::load_all(".", quiet = TRUE)
source(file.path("equivalence", "latentgold", "compare-functions.R"))
base <- file.path("equivalence", "latentgold-categorical-distal")
source(file.path(base, "common.R"))
references <- lapply(stats::setNames(distal_cases, distal_cases), function(case) {
  distal_reference(distal_case(case))
})
comparison <- distal_comparison(references, file.path(base, "returned"))
utils::write.csv(comparison, file.path(base, "comparison.csv"), row.names = FALSE)
print(table(comparison$status))
print(subset(comparison, status == "differ"), row.names = FALSE)
