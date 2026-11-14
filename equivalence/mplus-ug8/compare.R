# Growth and latent transition models against the Mplus 8.8 User's Guide,
# chapter 8 (README.md). Run from the project root:
#   Rscript equivalence/mplus-ug8/compare.R
# Writes comparison.csv and RUN.txt beside this script. Exits non-zero when a
# quantity disagrees without a registered known difference, or when a
# registered difference no longer matches a disagreement.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
directory <- file.path("equivalence", "mplus-ug8")
source(file.path("equivalence", "package-kit-runner.R"))
source(file.path(directory, "compare-functions.R"))
if (!run_package_kit(directory, function() compare_all_ug8(directory), character())) {
  quit(status = 1L, save = "no")
}
