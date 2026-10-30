# Fast safeguards: no models fitted and no repository reports regenerated.
# From the repository root: Rscript equivalence/verification-matrix/tests/run.R
library(testthat)
source("equivalence/verification-matrix/build.R")
source("equivalence/run-package-kit-tests.R")

test_that("shared kit helpers load from testthat's equivalence working directory", {
  environment <- new.env(parent = globalenv())
  withr::with_dir("equivalence", sys.source("package-kit-runner.R", envir = environment))
  expect_true(is.function(environment$apply_register))
  expect_true(is.function(environment$matrix_capture_source))
})

#' A miniature package tree for provenance mutations
#' @return Temporary directory, removed by each test.
source_fixture <- function() {
  root <- tempfile("matrix-source-")
  dir.create(file.path(root, "R"), recursive = TRUE)
  writeLines("export(example)", file.path(root, "NAMESPACE"))
  writeLines("example <- function() 1", file.path(root, "R", "example.R"))
  writeLines(c("Package: fixture", "Version: 1.0.0"), file.path(root, "DESCRIPTION"))
  utils::write.csv(data.frame(estimate = 1), file.path(root, "result.csv"), row.names = FALSE)
  root
}

test_that("source and result identities are both required for current evidence", {
  root <- source_fixture()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  result <- file.path(root, "result.csv")
  expect_false(matrix_evidence_current(result, root))
  snapshot <- matrix_capture_source(root)
  matrix_record_evidence(result, snapshot, root)
  expect_true(matrix_evidence_current(result, root))
  # Manifest row order can differ across tools or locale; identity is by path.
  sidecar <- paste0(result, ".provenance.csv")
  rows <- read.csv(sidecar)
  utils::write.csv(rows[rev(seq_len(nrow(rows))), ], sidecar, row.names = FALSE)
  expect_true(matrix_evidence_current(result, root))
  utils::write.csv(rbind(rows, rows[1L, ]), sidecar, row.names = FALSE)
  expect_false(matrix_evidence_current(result, root))
  utils::write.csv(rows, sidecar, row.names = FALSE)
  writeLines("example <- function() 2", file.path(root, "R", "example.R"))
  expect_false(matrix_evidence_current(result, root))
  expect_error(matrix_record_evidence(result, snapshot, root), "runtime source changed")
  writeLines("example <- function() 1", file.path(root, "R", "example.R"))
  expect_true(matrix_evidence_current(result, root))
  writeLines("additional <- TRUE", file.path(root, "R", "added.R"))
  expect_false(matrix_evidence_current(result, root))
  unlink(file.path(root, "R", "added.R"))
  writeLines("", file.path(root, "NAMESPACE"))
  expect_false(matrix_evidence_current(result, root))
  writeLines("export(example)", file.path(root, "NAMESPACE"))
  utils::write.csv(data.frame(estimate = 999), result, row.names = FALSE)
  expect_false(matrix_evidence_current(result, root))
})

test_that("retrospective attestations identify themselves and remain hash bound", {
  root <- source_fixture()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  result <- file.path(root, "result.csv")
  sidecar <- matrix_record_evidence(result, matrix_capture_source(root), root,
    basis = "retrospective source-equivalence attestation",
    note = "Synthetic test: historical source independently compared; not a rerun.")
  rows <- read.csv(sidecar)
  expect_true(all(rows$basis == "retrospective source-equivalence attestation"))
  expect_true(all(grepl("not a rerun", rows$note)))
  expect_true(matrix_evidence_current(result, root))
  rows$md5[rows$kind == "result"] <- "tampered"
  utils::write.csv(rows, sidecar, row.names = FALSE)
  expect_false(matrix_evidence_current(result, root))
})

test_that("actual testthat errors and skips cannot create source strength", {
  path <- tempfile(fileext = ".R")
  results_path <- tempfile(fileext = ".csv")
  on.exit(unlink(c(path, results_path)), add = TRUE)
  writeLines(c(
    'test_that("crashed reference", { expect_true(TRUE); stop("reference crashed") })',
    'test_that("skipped reference", { skip("optional package absent") })'), path)
  actual <- as.data.frame(testthat::test_file(path, reporter = "silent"))
  utils::write.csv(actual[setdiff(names(actual), "result")], results_path, row.names = FALSE)
  files <- data.frame(role = "tests", dir = "tests/equivalence", path = results_path,
                       run = "2026-10-29", code = "1.0.0", current = TRUE)
  imported <- .matrix_test_results(files)
  expect_identical(imported$nb, c(1L, 0L))
  expect_identical(imported$failed, c(1L, 0L))
  expect_identical(imported$incomplete, c(1L, 1L))
  evidence <- data.frame(evidence_id = c("crashed", "skipped"), origin = "testthat",
    location = "tests/equivalence/test-weighted-logits.R", item = imported$test,
    fn = "multilpa", reference = imported$test, kind = "R package", checks = imported$nb,
    values = NA_real_, failed = imported$failed, incomplete = imported$incomplete,
    run = imported$run, code = imported$code, current = TRUE, dataset = NA_character_,
    shared = FALSE, known = FALSE, one_class = FALSE)
  candidates <- data.frame(fn = character(), reference = character(), kind = character())
  cells <- .matrix_cells(evidence, candidates)
  expect_identical(cells$status[match(evidence$reference, cells$reference)], c("failed", "incomplete"))
  summary <- .matrix_summary(evidence, cells)
  expect_equal(summary$n_sources[summary$fn == "multilpa"], 0L)
  expect_equal(summary$failed_current[summary$fn == "multilpa"], 1L)
  expect_equal(summary$incomplete_current[summary$fn == "multilpa"], 2L)
})

test_that("internal checks are retained but excluded from external check totals", {
  evidence <- data.frame(evidence_id = c("external", "internal"), origin = "testthat",
    location = "synthetic", item = "two classes", fn = "lpa", reference = c("program", "self"),
    kind = c("R package", "internal"), checks = c(1L, 100L), values = c(3, 400),
    failed = 0L, incomplete = 0L, run = "2026-10-29", code = "1.0.0", current = TRUE,
    dataset = NA_character_, shared = FALSE, known = FALSE, one_class = FALSE)
  candidates <- data.frame(fn = character(), reference = character(), kind = character())
  cells <- .matrix_cells(evidence, candidates)
  summary <- .matrix_summary(evidence, cells)
  expect_equal(summary$checks_current[summary$fn == "lpa"], 1L)
  expect_equal(summary$internal_checks[summary$fn == "lpa"], 100L)
  report <- .matrix_markdown(evidence, cells, summary)
  expect_match(report[[3L]], "1 current checks.*3 values")
})

test_that("CI rejects empty, errored, skipped and partial outcome tables", {
  complete <- data.frame(passed = 3L, failed = 0L, error = FALSE, skipped = FALSE)
  expect_true(require_completed_kit_tests(complete))
  expect_error(require_completed_kit_tests(complete[FALSE, ]), "no complete")
  expect_error(require_completed_kit_tests(complete["passed"]), "no complete")
  invisible(lapply(c("error", "skipped", "failed"), function(column) {
    broken <- complete
    broken[[column]] <- 1L
    expect_error(require_completed_kit_tests(broken), "did not complete")
  }))
  complete$passed <- 0L
  expect_error(require_completed_kit_tests(complete), "did not complete")
})

source("equivalence/validation/check-package-kit-gates.R")

test_that("a 'sometimes' register entry explains a difference but is never stale", {
  environment <- new.env(parent = globalenv())
  withr::with_dir("equivalence", sys.source("package-kit-runner.R", envir = environment))
  comparison <- data.frame(case = c("a", "b"), reference = "pkg",
                           quantity = "maximized log likelihood",
                           agrees = c(FALSE, TRUE), difference = c(1, 0),
                           tolerance = 1e-6, n_values = 1)
  known <- data.frame(case = c("a", "b", "c"), reference = "pkg",
                      quantity = "^maximized log likelihood$", reason = "r",
                      occurs = c("sometimes", "sometimes", "always"))
  result <- environment$apply_register(comparison, known)
  expect_identical(result$comparison$status, c("known difference", "agrees"))
  # "b" did not occur but is optional; "c" did not occur and is required.
  expect_identical(result$stale$case, "c")
  # Without the column every entry is required, as before.
  result <- environment$apply_register(comparison, known[c("case", "reference", "quantity", "reason")])
  expect_identical(result$stale$case, c("b", "c"))
  known$occurs[1] <- "often"
  expect_error(environment$apply_register(comparison, known), "occurs")
})
