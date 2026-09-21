# Run the equivalence tests: every comparison of multilpa against an external
# implementation or a recorded reference (Mplus, mclust, depmixS4, tidySEM,
# the JStats latent transition fixtures, the in-house exhaustive enumeration).
#
# From the project root:
#   Rscript equivalence/run.R                 # every file
#   Rscript equivalence/run.R jstats mclust   # files whose name matches
#
# This folder is in .Rbuildignore and never ships with the package; the
# package's own tests/testthat/ holds unit and invariant tests only.
requested <- commandArgs(trailingOnly = TRUE)
filter <- if (length(requested) > 0L) paste(requested, collapse = "|") else NULL

# The package is loaded here rather than by `load_package = "source"`: that
# mode sources the helpers of tests/testthat/ and skips this folder's own.
pkgload::load_all(".", quiet = TRUE)
testthat::test_dir(file.path("equivalence"), filter = filter,
                   package = "multilpa", load_package = "none",
                   stop_on_failure = TRUE)
