# Reproduces the paquid growth mixtures of Proust-Lima, Philipps and Liquet
# (2017) with mixture_regression(), against the published summary table and
# posterior classification and the current lcmm. Run from the project root:
#   Rscript equivalence/growth-published/compare.R
# Writes comparison.csv and RUN.txt beside this script. Exits non-zero when a
# quantity disagrees without a registered known difference, or when a
# registered difference no longer matches a disagreement.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
directory <- file.path("equivalence", "growth-published")
source(file.path("equivalence", "package-kit-runner.R"))
source(file.path(directory, "compare-functions.R"))
if (!run_package_kit(directory, function() compare_all_growth_published(directory), "lcmm")) {
  quit(status = 1L, save = "no")
}
