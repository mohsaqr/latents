# Run the equivalence tests: every comparison of multilpa against an external
# implementation or a recorded reference (Mplus, mclust, depmixS4, tidySEM,
# the JStats latent transition fixtures, the in-house exhaustive enumeration).
#
# From the project root:
#   Rscript equivalence/run.R                 # every file
#   Rscript equivalence/run.R jstats mclust   # files whose name matches
#   Rscript equivalence/run.R --dir=tests/equivalence   # the tests/equivalence/ files
#
# This folder is in .Rbuildignore and never ships with the package; the
# package's own tests/testthat/ holds unit and invariant tests only.
#
# Every run writes tmp/equivalence-tests/<label>-results.csv (one row per test,
# as.data.frame() of the testthat results) with its provenance sidecar (the
# runtime source captured before the run; verification-matrix/provenance.R),
# which the verification matrix reads as current evidence, and reports to CVS.
arguments <- commandArgs(trailingOnly = TRUE)
# --dir= runs another folder of equivalence tests (tests/equivalence/, the
# numerical references of the package's own models) the same way.
dir_argument <- startsWith(arguments, "--dir=")
test_directory <- if (any(dir_argument)) sub("^--dir=", "", arguments[dir_argument][[1L]]) else "equivalence"
stopifnot("--dir= names a folder of the project" = dir.exists(test_directory))
requested <- arguments[!dir_argument]
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
source(file.path("equivalence", "verification-matrix", "provenance.R"))
source_snapshot <- matrix_capture_source(".")
started <- Sys.time()
results <- testthat::test_dir(test_directory, filter = filter,
                              package = "latents", load_package = "none",
                              stop_on_failure = FALSE)

# Every run reports to CVS, failing runs included (equivalence/cvs/cvs-report.R),
# labelled by its filter when it has one (the ci shards report disjoint files).
source(file.path("equivalence", "cvs", "cvs-report.R"))
filter_label <- gsub("^-+|-+$", "", substr(gsub("[^a-z0-9]+", "-", tolower(filter %||% "")), 1L, 40L))
directory_label <- gsub("[^a-z0-9]+", "-", tolower(test_directory))
label <- if (nzchar(filter_label)) paste0(directory_label, "-", filter_label) else directory_label
cvs_write_report(.cvs_from_testthat(results, test_directory, "equiv"), label, started)

outcome <- as.data.frame(results)
results_dir <- file.path("tmp", "equivalence-tests")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
results_path <- file.path(results_dir, sprintf("%s-results.csv", label))
utils::write.csv(outcome[setdiff(names(outcome), "result")], results_path, row.names = FALSE)
matrix_record_evidence(results_path, source_snapshot, ".")
if (any(outcome$failed > 0L | outcome$error)) {
  quit(status = 1L, save = "no")
}
