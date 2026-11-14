# Recovery of known parameters from simulab data, for every public function
# (README.md). About an hour; run from the project root:
#   Rscript equivalence/simulab-recovery/compare.R
# Writes comparison.csv and RUN.txt beside this script. Exits non-zero when a
# quantity disagrees without a registered known difference, or when a
# registered difference no longer matches a disagreement.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
directory <- file.path("equivalence", "simulab-recovery")
source(file.path("equivalence", "package-kit-runner.R"))
source(file.path(directory, "compare-functions.R"))
if (!run_package_kit(directory, compare_all_recovery, c("simulab", "mice"))) {
  quit(status = 1L, save = "no")
}
