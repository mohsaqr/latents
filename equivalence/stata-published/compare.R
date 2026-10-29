# lpa() and lca() against genuine Stata 19 output printed in its manual
# (Examples 52g and 50g) and Masyn (2013) Table 25.12. Run from the project
# root:
#   Rscript equivalence/stata-published/compare.R
# Writes comparison.csv and RUN.txt beside this script. Exits non-zero when a
# quantity disagrees without a registered known difference, or when a
# registered difference no longer matches a disagreement.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
directory <- file.path("equivalence", "stata-published")
source(file.path("equivalence", "package-kit-runner.R"))
source(file.path(directory, "compare-functions.R"))
if (!run_package_kit(directory, function() compare_all_stata(directory), character())) {
  quit(status = 1L, save = "no")
}
