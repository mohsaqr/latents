# Growth mixtures in mixture_regression() against genuine Mplus 9 Demo runs
# (JStats fixtures, kept in returned/). Run from the project root:
#   Rscript equivalence/mplus-growth/compare.R
# Writes comparison.csv and RUN.txt beside this script. Exits non-zero when a
# quantity disagrees without a registered known difference.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
directory <- file.path("equivalence", "mplus-growth")
source(file.path("equivalence", "package-kit-runner.R"))
source(file.path(directory, "compare-functions.R"))
if (!run_package_kit(directory, function() compare_all_mplus_growth(directory), character())) {
  quit(status = 1L, save = "no")
}
