# The CVS reporter (equivalence/cvs/cvs-report.R): every run's outcome reaches
# CVS with the status it had, nothing is counted twice or dropped, and the file
# is one the CVS server attributes to latents.

source(testthat::test_path("cvs", "cvs-report.R"), local = TRUE)

cvs_statuses <- function(suites) {
  unlist(lapply(suites, function(suite) {
    vapply(suite$assertionResults, `[[`, character(1), "status")
  }))
}

cvs_titles <- function(suite) vapply(suite$assertionResults, `[[`, character(1), "title")

test_that("testthat results keep each block's outcome", {
  directory <- withr::local_tempdir()
  writeLines(c(
    'test_that("agrees", { expect_equal(1, 1); expect_true(TRUE) })',
    'test_that("disagrees", { expect_equal(1, 2) })',
    'test_that("no reference", { skip("reference not installed") })',
    'test_that("stops", { stop("boom") })',
    'test_that("passes then skips", { expect_true(TRUE); skip("platform") })',
    'test_that("agrees", { expect_true(TRUE) })'), file.path(directory, "test-one.R"))
  writeLines('test_that("other file", { expect_true(TRUE) })', file.path(directory, "test-two.R"))
  results <- testthat::test_dir(directory, stop_on_failure = FALSE, reporter = "silent")
  suites <- .cvs_from_testthat(results, "equivalence", "equiv")

  expect_identical(vapply(suites, `[[`, character(1), "name"),
                   c("equivalence/test-one.R", "equivalence/test-two.R"))
  expect_identical(cvs_titles(suites[[1L]]),
                   c("agrees", "disagrees", "no reference", "stops", "passes then skips", "agrees [2]"))
  expect_identical(cvs_statuses(suites[1L]),
                   c("passed", "failed", "skipped", "failed", "skipped", "passed"))
  expect_identical(suites[[1L]]$assertionResults[[1L]][["_cvs"]],
                   list(checkKind = "testthat", module = "one"))
  failure <- suites[[1L]]$assertionResults[[4L]]$failureMessages
  expect_match(unlist(failure), "boom")
  skipped <- suites[[1L]]$assertionResults[[3L]]$failureMessages
  expect_match(unlist(skipped), "reference not installed")
  # Every file declares its tier, since CVS cannot read one from an R file name.
  expect_identical(lapply(suites, `[[`, "_cvs"), rep(list(list(tier = "equiv")), 2L))
  unit <- .cvs_from_testthat(results, file.path("tests", "testthat"), "unit")
  expect_identical(unit[[1L]][["_cvs"]], list(tier = "unit"))
  expect_identical(unit[[1L]]$name, "tests/testthat/test-one.R")
  expect_error(.cvs_from_testthat(results, "equivalence", "harness"), "must be \"equiv\" or \"unit\"")
  expect_error(.cvs_from_testthat(results, "equivalence"), "tier")
  # Invariant: one assertion per test_that() block, none lost or added.
  expect_identical(length(cvs_statuses(suites)), nrow(as.data.frame(results)))
})

test_that("harness quantities keep their difference, tolerance and reference", {
  comparisons <- data.frame(
    suite = c("polca", "polca", "polca", "mclust"),
    source = c("poLCA 1.6.0.2", "poLCA 1.6.0.2", "poLCA 1.6.0.2", "mclust 6.1.3"),
    quantity = c("loglik", "loglik", "posteriors", "means"),
    difference = c(-2e-9, 3e-9, Inf, 0.5), tolerance = c(1e-6, 1e-6, 1e-6, 1e-6),
    agrees = c(TRUE, TRUE, FALSE, FALSE), output = c(NA, NA, "posteriors", "means"),
    n_values = c(1L, 1L, 100L, 6L), known_open = c(NA, NA, NA, "KO-1"),
    stringsAsFactors = FALSE)
  suites_table <- data.frame(suite = c("polca", "mclust", "lcmm", "glca"),
                             status = c("ran", "ran", "skipped", "failed"),
                             reason = c(NA, NA, "lcmm is not installed", "object not found"),
                             stringsAsFactors = FALSE)
  suites <- .cvs_from_harness(comparisons, suites_table)

  expect_identical(vapply(suites, `[[`, character(1), "name"),
                   file.path("equivalence", "harness", "suites",
                             c("polca.R", "mclust.R", "lcmm.R", "glca.R")))
  expect_identical(cvs_titles(suites[[1L]]), c("loglik", "loglik [2]", "posteriors"))
  expect_true(all(vapply(suites, function(suite) identical(suite[["_cvs"]]$tier, "equiv"), logical(1))))
  expect_identical(cvs_statuses(suites), c("passed", "passed", "failed", "failed", "skipped", "failed"))
  first <- suites[[1L]]$assertionResults[[1L]][["_cvs"]]
  expect_equal(first$delta, 2e-9)
  expect_equal(first$tolerance, 1e-6)
  expect_identical(first$rPackage, "poLCA 1.6.0.2")
  expect_null(first$rFunction)
  # JSON has no Inf: a non-finite difference is left out, the failure says why.
  third <- suites[[1L]]$assertionResults[[3L]]
  expect_null(third[["_cvs"]]$delta)
  expect_match(unlist(third$failureMessages), "exceeds tolerance")
  expect_match(unlist(suites[[2L]]$assertionResults[[1L]]$failureMessages), "known open defect KO-1")
  expect_match(unlist(suites[[4L]]$assertionResults[[1L]]$failureMessages), "object not found")
})

test_that("every other reason the harness fails a run reaches CVS as a failure", {
  comparisons <- data.frame(
    suite = c("lcmm", "lta-general"), source = "lcmm 2.2.2", quantity = c("loglik", "loglik"),
    difference = 0, tolerance = 1e-6, agrees = TRUE, output = NA_character_,
    n_values = 1L, known_open = NA_character_, stringsAsFactors = FALSE)
  suites_table <- data.frame(suite = c("lcmm", "lta-general"), status = "ran",
                             reason = NA_character_, stringsAsFactors = FALSE)
  coverage <- data.frame(
    suite = c("lcmm", "lta-general", "lta-general", "lta-general"),
    fit = c("two classes", "chain", "chain", "chain"),
    output = c("classes$share", "transitions$supported", "transitions$probability",
               "transitions$origin_expected_count"),
    n_produced = c(2L, 4L, 4L, 4L), n_compared = c(2L, 0L, 4L, 2L),
    status = c("compared", "not covered", "compared", "partial"),
    reason = c(NA, "no reference computes it", NA, NA), stringsAsFactors = FALSE)
  known_open <- data.frame(id = c("KO-1", "KO-2"), suite = c("lcmm", "lcmm"),
                           status = c("matched", "stale"), stringsAsFactors = FALSE)
  suites <- .cvs_from_harness(comparisons, suites_table, coverage, known_open)

  # Every quantity agrees, yet the run failed: CVS must show why.
  expect_identical(cvs_titles(suites[[1L]]), c("loglik", "known-open entry KO-2 matches"))
  expect_identical(cvs_titles(suites[[2L]]),
                   c("loglik", "output compared: chain: transitions$supported",
                     "output compared: chain: transitions$origin_expected_count"))
  expect_identical(cvs_statuses(suites), c("passed", "failed", "passed", "failed", "failed"))
  expect_identical(unlist(suites[[2L]]$assertionResults[[2L]]$failureMessages),
                   "not covered: 0 of 4 values compared; no reference computes it")
  expect_identical(unlist(suites[[2L]]$assertionResults[[3L]]$failureMessages),
                   "partial: 2 of 4 values compared")
  # Runs from before the coverage ledger have no coverage or register to pass.
  expect_identical(cvs_statuses(.cvs_from_harness(comparisons, suites_table)), c("passed", "passed"))
})

test_that("sweep statuses map to CVS statuses, and an unknown one stops", {
  results <- data.frame(
    module = c("gaussian", "gaussian", "gaussian", "lta", "lta"),
    case = c("a", "a", "b", "c", "d"),
    quantity = c("loglik", "posteriors", "loglik", "means", "loglik"),
    difference = c(0, 1e-9, 2, NA, Inf), tolerance = c(1e-8, 1e-6, 1e-8, 0, 0),
    status = c("agrees", "known difference", "disagrees", "no reference", "error"),
    note = c(NA, "k-means stops lower", NA, "no continuous indicator", "worker failed"),
    stringsAsFactors = FALSE)
  suites <- .cvs_from_sweep(results)
  expect_identical(vapply(suites, `[[`, character(1), "name"),
                   file.path("equivalence", "argument-sweep", "modules", c("gaussian.R", "lta.R")))
  expect_identical(cvs_statuses(suites), c("passed", "skipped", "failed", "skipped", "failed"))
  expect_identical(cvs_titles(suites[[1L]]), c("a: loglik", "a: posteriors", "b: loglik"))
  expect_identical(suites[[2L]][["_cvs"]], list(tier = "equiv"))
  expect_identical(unlist(suites[[1L]]$assertionResults[[2L]]$failureMessages),
                   "known difference: k-means stops lower")
  results$status[[1L]] <- "agreed"
  expect_error(.cvs_from_sweep(results), "unknown sweep status: agreed")
})

test_that("a report goes to the inbox and to tmp/cvs, or only to tmp/cvs when reporting is off", {
  suites <- .cvs_from_sweep(data.frame(
    module = "gaussian", case = c("a", "b", "c"), quantity = "loglik",
    difference = c(0, 1, 0), tolerance = 1e-8,
    status = c("agrees", "disagrees", "no reference"), note = NA_character_,
    stringsAsFactors = FALSE))
  started <- as.POSIXct("2026-11-13 10:00:00", tz = "UTC")

  local <- withr::local_tempdir()
  inbox <- withr::local_tempdir()
  withr::local_envvar(CVS_DISABLED = "", CVS_INBOX_DIR = inbox)
  expect_message(paths <- cvs_write_report(suites, "sweep", started, output = local),
                 "CVS report: .* [(]3 tests, 1 failed[)]")
  expect_length(paths, 2L)
  expect_identical(dirname(normalizePath(paths)), normalizePath(c(local, inbox)))
  # Named by the run's start, not by when the file was written.
  expect_identical(basename(paths[[2L]]),
                   sprintf("latents-sweep-%s.json", format(started, "%Y%m%dT%H%M%S")))
  expect_identical(list.files(inbox, all.files = TRUE, no.. = TRUE), basename(paths[[2L]]))

  report <- jsonlite::fromJSON(paths[[2L]], simplifyVector = FALSE)
  expect_identical(report[["_cvs"]]$target, "latents")
  expect_identical(report[["_cvs"]]$version, as.character(utils::packageVersion("latents")))
  expect_identical(report$startTime, as.numeric(started) * 1000)
  expect_identical(c(report$numTotalTests, report$numPassedTests, report$numFailedTests,
                     report$numPendingTests), c(3L, 1L, 1L, 1L))
  expect_identical(c(report$numPassedTestSuites, report$numFailedTestSuites), c(0L, 1L))
  expect_false(report$success)
  # Invariant: the totals are the assertions in the file, nothing more or less.
  assertions <- unlist(lapply(report$testResults, function(suite) {
    vapply(suite$assertionResults, `[[`, character(1), "status")
  }))
  expect_identical(length(assertions), report$numTotalTests)
  expect_identical(readLines(paths[[1L]]), readLines(paths[[2L]]))

  withr::local_envvar(CVS_DISABLED = "1")
  expect_message(off <- cvs_write_report(suites, "sweep", started, output = local),
                 "reporting off")
  expect_length(off, 1L)
  expect_identical(dirname(normalizePath(off)), normalizePath(local))
  expect_message(older <- cvs_write_report(suites, "sweep", started, output = local,
                                            version = "1.1.6"), "reporting off")
  expect_identical(jsonlite::fromJSON(older)[["_cvs"]]$version, "1.1.6")
  expect_error(cvs_write_report(suites, "Sweep Run", started, output = local), "lowercase")
  expect_error(cvs_write_report(list(), "sweep", started, output = local), "non-empty")
})
