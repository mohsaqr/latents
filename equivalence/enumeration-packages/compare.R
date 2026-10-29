# enumerate_lpa(), enumerate_lca() and bootstrap_lrt() against mclust and
# poLCA. Run from the project root:
#   Rscript equivalence/enumeration-packages/compare.R
# Writes comparison.csv, search.csv and RUN.txt beside this script. Exits
# non-zero when a quantity disagrees without a registered known difference.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
directory <- file.path("equivalence", "enumeration-packages")
source(file.path("equivalence", "package-kit-runner.R"))
source(file.path(directory, "compare-functions.R"))
if (!run_package_kit(directory, compare_all_enumeration, c("mclust", "poLCA"))) {
  quit(status = 1L, save = "no")
}
