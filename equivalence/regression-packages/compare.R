# Compares mixture_regression() with independent R implementations (mixtools,
# MoEClust, flexmix) and with published results (SAS PROC FMM Example 40.3,
# Grun and Leisch 2008, Benaglia et al. 2009). Run from the project root:
#   Rscript equivalence/regression-packages/compare.R
# Writes comparison.csv and RUN.txt beside this script. Exits non-zero when a
# quantity disagrees without a registered known difference, or when a
# registered difference no longer matches a disagreement.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
directory <- file.path("equivalence", "regression-packages")
source(file.path("equivalence", "package-kit-runner.R"))
source(file.path(directory, "compare-functions.R"))
if (!run_package_kit(directory, compare_all_regression, names(.fitters))) quit(status = 1L, save = "no")
