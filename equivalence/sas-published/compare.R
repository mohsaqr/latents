# lca(), lpa(), enumerate_lpa() and mixture_regression() against SAS output
# printed by UCLA OARC (PROC LCA), in SAS and R (PROC LCA) and in the SAS/STAT
# 14.2 User's Guide (PROC FMM) (README.md). Run from the project root:
#   Rscript equivalence/sas-published/compare.R
# Writes comparison.csv and RUN.txt beside this script. Exits non-zero when a
# quantity disagrees without a registered known difference, or when a
# registered difference no longer matches a disagreement.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
directory <- file.path("equivalence", "sas-published")
source(file.path("equivalence", "package-kit-runner.R"))
source(file.path(directory, "compare-functions.R"))
if (!run_package_kit(directory, function() compare_all_sas(directory), character())) {
  quit(status = 1L, save = "no")
}
