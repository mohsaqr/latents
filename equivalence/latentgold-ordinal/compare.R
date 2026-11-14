# Compare latents with Latent GOLD 6.1 on the ordinal/count kit. Run from the
# package root after make-kit.R and the Latent GOLD batch runs (listings
# copied to equivalence/latentgold-ordinal/returned/):
#   Rscript equivalence/latentgold-ordinal/compare.R
# equivalence/test-latentgold-ordinal.R runs the same comparison under testthat.
devtools::load_all(".", quiet = TRUE)
base <- file.path("equivalence", "latentgold-ordinal")
source(file.path(base, "compare-functions.R"))
comparison <- compare_ordinal_kit(base)
print(comparison, digits = 10, row.names = FALSE)
utils::write.csv(comparison, file.path(base, "comparison.csv"), row.names = FALSE)
