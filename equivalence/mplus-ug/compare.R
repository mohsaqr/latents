# The Mplus 8.8 User's Guide examples latents fits (chapters 7 and 8),
# against Mplus's own output (README.md). Run from the project root:
#   Rscript equivalence/mplus-ug/compare.R
# Writes comparison.csv and RUN.txt beside this script. Exits non-zero when a
# quantity disagrees without a registered known difference, or when a
# registered difference no longer matches a disagreement.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
directory <- file.path("equivalence", "mplus-ug")
source(file.path("equivalence", "package-kit-runner.R"))
source(file.path(directory, "compare-functions.R"))
if (!run_package_kit(directory, function() compare_all_ug(directory), character())) {
  quit(status = 1L, save = "no")
}
