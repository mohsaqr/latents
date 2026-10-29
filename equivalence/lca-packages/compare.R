# Compares lca() with five independent R implementations of the single-level
# latent class model on eight datasets. Run from the project root:
#   Rscript equivalence/lca-packages/compare.R
# Writes comparison.csv and RUN.txt beside this script. Exits non-zero when a
# quantity disagrees without a registered known difference, or when a
# registered difference no longer matches a disagreement.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
directory <- file.path("equivalence", "lca-packages")
source(file.path("equivalence", "package-kit-runner.R"))
source(file.path(directory, "compare-functions.R"))
if (!run_package_kit(directory, compare_all_lca, names(.fitters))) quit(status = 1L, save = "no")
