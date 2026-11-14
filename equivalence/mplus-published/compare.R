# lca(), lpa(), lta() and the VLMR statistics against genuine Mplus output
# printed on the UCLA OARC Mplus pages (README.md). Run from the project root:
#   Rscript equivalence/mplus-published/compare.R
# Writes comparison.csv and RUN.txt beside this script. Exits non-zero when a
# quantity disagrees without a registered known difference, or when a
# registered difference no longer matches a disagreement.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
directory <- file.path("equivalence", "mplus-published")
source(file.path("equivalence", "package-kit-runner.R"))
source(file.path(directory, "compare-functions.R"))
if (!run_package_kit(directory, function() compare_all_mplus(directory), character())) {
  quit(status = 1L, save = "no")
}
