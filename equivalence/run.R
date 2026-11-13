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

# This folder never goes to CRAN, so a `skip_on_cran()` here must not fire.
# Without this, testthat reads every plain Rscript run as a CRAN run and the
# mclust comparisons were skipped by `Rscript equivalence/run.R` itself.
Sys.setenv(NOT_CRAN = "true")

# The package is loaded here rather than by `load_package = "source"`: that
# mode sources the helpers of tests/testthat/ and skips this folder's own.
pkgload::load_all(".", quiet = TRUE)
# Every file reports: testthat otherwise quits after ten failed expectations,
# and one reference package's changed return value (six failures in a file)
# ended a run with thirteen files unread.
testthat::set_max_fails(Inf)
started <- Sys.time()
results <- testthat::test_dir(file.path("equivalence"), filter = filter,
                              package = "latents", load_package = "none",
                              stop_on_failure = FALSE)

# Every run reports to CVS, failing runs included (equivalence/cvs/cvs-report.R),
# labelled by its filter when it has one (the ci shards report disjoint files).
source(file.path("equivalence", "cvs", "cvs-report.R"))
filter_label <- gsub("^-+|-+$", "", substr(gsub("[^a-z0-9]+", "-", tolower(filter %||% "")), 1L, 40L))
label <- if (nzchar(filter_label)) paste0("equivalence-", filter_label) else "equivalence"
cvs_write_report(.cvs_from_testthat(results, "equivalence", "equiv"), label, started)

outcome <- as.data.frame(results)
if (any(outcome$failed > 0L | outcome$error)) {
  quit(status = 1L, save = "no")
}
