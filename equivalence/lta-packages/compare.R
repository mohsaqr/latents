# Compares lta() with LMest and seqHMM on RLMSdat, and with the published
# analysis of Bartolucci, Pandolfi and Pennoni (2017). Run from the project
# root:
#   Rscript equivalence/lta-packages/compare.R
# Writes comparison.csv and RUN.txt beside this script. Exits non-zero when a
# quantity disagrees without a registered known difference, or when a
# registered difference no longer matches a disagreement.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
directory <- file.path("equivalence", "lta-packages")
source(file.path("equivalence", "package-kit-runner.R"))
source(file.path(directory, "compare-functions.R"))
if (!run_package_kit(directory, compare_all_lta_packages, c("LMest", "seqHMM"))) {
  quit(status = 1L, save = "no")
}
