# Run the package's own tests (tests/testthat) and report them to CVS.
#
# From the project root:
#   Rscript equivalence/cvs/test-local.R
#
# The same run as `testthat::test_local()`, which the ci job `test` ran before
# CVS reporting; the script exits non-zero when any test failed or errored.
started <- Sys.time()
results <- testthat::test_local(stop_on_failure = FALSE)
source(file.path("equivalence", "cvs", "cvs-report.R"))
cvs_write_report(.cvs_from_testthat(results, file.path("tests", "testthat"), "unit"), "unit",
                 started)
outcome <- as.data.frame(results)
if (any(outcome$failed > 0L | outcome$error)) {
  quit(status = 1L, save = "no")
}
