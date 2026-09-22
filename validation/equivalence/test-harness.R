# Tests for the harness itself.
#
# The equivalence suites test the package. Nothing tested the harness, which is
# how it came to sit broken for two releases while its saved report was still
# being cited. These run in seconds and fit no models.
#
# From the project root:
#   Rscript validation/equivalence/test-harness.R

suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("validation", "equivalence", "registry.R"))

if (!requireNamespace("testthat", quietly = TRUE)) {
  cat("SKIPPED: testthat is not installed, so the harness self-tests cannot run.\n")
  quit(status = 0L, save = "no")
}
library(testthat)

scratch <- tempfile("multilpa-harness-tests-")
dir.create(scratch)
on.exit(unlink(scratch, recursive = TRUE), add = TRUE, after = FALSE)

# --------------------------------------------------------------------------
test_that("check_validation_api names every retired argument, with its line", {
  stale <- file.path(scratch, "stale.R")
  writeLines(c(
    'demo <- function(students, columns) {',
    '  fit <- multilpa(students, indicators = columns, group = "school", n_profiles = 2)',
    '  three_step(fit, students, "y", level_ci = 0.9, p_adjust = "BH")',
    '  bootstrap_lrt(fit, fit, n_boot = 10)',
    '  enumerate_classes(students, vars = columns, id = "school", profiles = 2:3,',
    '                    group_classes = 1:2)',
    '  parameter_inference(object = fit)',
    '}'), stale)
  found <- check_validation_api(paths = stale)
  expect_equal(nrow(found), 8L)
  expect_setequal(found$argument,
                  c("indicators", "group", "level_ci", "p_adjust", "n_boot",
                    "profiles", "group_classes", "object"))
  expect_setequal(found$replacement,
                  c("vars", "id", "ci_level", "adjust", "iter", "n_profiles",
                    "n_group_classes", "x"))
  # The line is the line the argument is written on, not the line the call
  # starts on: `group_classes` is on the continuation line.
  expect_equal(found$line[found$argument == "group_classes"], 6L)
})

test_that("check_validation_api accepts current calls and dots-forwarded ones", {
  current <- file.path(scratch, "current.R")
  writeLines(c(
    'demo <- function(students, columns) {',
    '  multilpa(students, vars = columns, id = "school", n_profiles = 2)',
    '  # n_starts is not a formal of enumerate_classes(); it reaches multilpa()',
    '  # through the dots, and must not be reported as stale.',
    '  enumerate_classes(students, vars = columns, id = "school", n_starts = 6L)',
    '}'), current)
  expect_equal(nrow(check_validation_api(paths = current)), 0L)
})

test_that("check_validation_api names a call to a removed verb", {
  # get_data() was renamed get_results() in 0.4.0 and seven scripts kept the
  # old name, unseen, because only still-exported verbs were being checked.
  removed <- file.path(scratch, "removed.R")
  writeLines(c('x <- get_data(fit, "entropy")',
               'y <- fit_transitions(d, "y", "g", n_profiles = 2, time = "t")',
               'lmr_lrt(a, b)'), removed)
  found <- check_validation_api(paths = removed)
  expect_equal(found$line, 1:3)
  expect_equal(found$call, c("get_data", "fit_transitions", "lmr_lrt"))
  expect_equal(found$replacement[1:2], c("get_results", "lta"))
  expect_match(found$replacement[3L], "^deferred:")
  expect_error(.stop_if_stale_api(found), "get_data\\(\\) -> get_results",
               class = "multilpa_stale_validation_api")
})

test_that("a script defining a removed name itself is calling its own binding", {
  own <- file.path(scratch, "own.R")
  writeLines(c('lmr_lrt <- deferred_verb("lmr.R", "lmr_lrt")',
               'lmr_lrt(a, b)'), own)
  expect_equal(nrow(check_validation_api(paths = own)), 0L)
})

test_that("the whole validation tree matches the loaded package", {
  expect_equal(nrow(check_validation_api(root = ".")), 0L)
})

test_that(".stop_if_stale_api raises a catchable classed condition", {
  stale <- file.path(scratch, "stale2.R")
  writeLines('multilpa(d, indicators = v, id = "g", n_profiles = 2)', stale)
  expect_error(.stop_if_stale_api(check_validation_api(paths = stale)),
               class = "multilpa_stale_validation_api")
})

# --------------------------------------------------------------------------
test_that("a missing package is a skip with a reason, not an error", {
  expect_error(
    require_suite_packages("a.package.that.is.not.installed", reason = "a reference"),
    class = "multilpa_suite_skipped")
  expect_true(require_suite_packages("stats"))
  condition <- tryCatch(
    require_suite_packages("a.package.that.is.not.installed", reason = "a reference"),
    multilpa_suite_skipped = function(condition) condition)
  expect_match(conditionMessage(condition), "not installed")
  expect_match(conditionMessage(condition), "a reference")
})

test_that("a missing artifact is a skip with a reason", {
  expect_error(require_suite_files(file.path(scratch, "absent.out")),
               class = "multilpa_suite_skipped")
  expect_true(require_suite_files(file.path(scratch, "stale.R")))
})

test_that("a skipped suite is recorded, and does not stop the other suites", {
  suites <- file.path(scratch, "suites")
  dir.create(suites)
  writeLines(c(
    'suite_alpha <- function() {',
    '  require_suite_packages("a.package.that.is.not.installed",',
    '                         reason = "the reference this suite needs")',
    '}'), file.path(suites, "alpha.R"))
  writeLines(c(
    'suite_beta <- function() {',
    '  transform(compare_values("a quantity", 1, 1, tolerance = 1e-8),',
    '            source = "a reference")',
    '}'), file.path(suites, "beta.R"))
  writeLines(c(
    'suite_gamma <- function() stop("this suite is broken")'),
    file.path(suites, "gamma.R"))
  dir.create(file.path(scratch, "validation", "equivalence"), recursive = TRUE)
  file.rename(suites, file.path(scratch, "validation", "equivalence", "suites"))
  dir.create(file.path(scratch, "R"))
  writeLines("# no sources", file.path(scratch, "R", "empty.R"))

  run <- suppressMessages(run_equivalence(root = scratch, check_api = FALSE))
  status <- as.data.frame(run, what = "suites")
  expect_equal(status$status, c("skipped", "ran", "failed"))
  expect_match(status$reason[status$suite == "alpha"], "not installed")
  expect_match(status$reason[status$suite == "gamma"], "this suite is broken")
  expect_equal(nrow(as.data.frame(run, what = "comparisons")), 1L)
})

# --------------------------------------------------------------------------
test_that("compare_values refuses to guess a tolerance", {
  expect_error(compare_values("a quantity", 1, 1))
  expect_error(compare_values("a quantity", 1, 1, tolerance = 0))
  agreed <- compare_values("a quantity", 1, 1 + 1e-9, tolerance = 1e-8)
  expect_true(agreed$agrees)
  expect_false(compare_values("a quantity", 1, 1 + 1e-7, tolerance = 1e-8)$agrees)
})

test_that("compare_printed measures the gap against the printed precision", {
  # 15.3 states its quantity to within 0.05; 15.26 is 0.8 of that bound.
  printed <- compare_printed("G-squared", 15.3, 15.26, digits = 1)
  expect_true(printed$agrees)
  expect_equal(printed$difference, 0.8, tolerance = 1e-12)
  expect_false(compare_printed("G-squared", 15.3, 15.2, digits = 1)$agrees)
})

# --------------------------------------------------------------------------
test_that("seed_stability reports the spread and the number of optima", {
  # A deterministic stand-in for a fit: the same optimum under every seed.
  stable <- seed_stability("a stable model",
                           \(seed) structure(-100, class = "logLik"),
                           seeds = c(1L, 2L, 3L, 4L), tolerance = 1e-5)
  expect_equal(nrow(stable), 2L)
  expect_true(all(stable$agrees))
  expect_equal(stable$obtained, c(0, 1))

  # Two of four seeds land on an inferior mode: both rows must fail.
  unstable <- seed_stability("an unstable model",
                             \(seed) structure(if (seed %% 2L == 0L) -100 else -101,
                                               class = "logLik"),
                             seeds = c(1L, 2L, 3L, 4L), tolerance = 1e-5)
  expect_false(any(unstable$agrees))
  expect_equal(unstable$obtained, c(1, 2))
})

# --------------------------------------------------------------------------
test_that("a run records the version and fingerprint it was produced under", {
  suites <- file.path(scratch, "provenance", "validation", "equivalence", "suites")
  dir.create(suites, recursive = TRUE)
  writeLines(c(
    'suite_delta <- function() {',
    '  transform(compare_values("a quantity", 1, 1, tolerance = 1e-8),',
    '            source = "a reference")',
    '}'), file.path(suites, "delta.R"))
  root <- file.path(scratch, "provenance")
  dir.create(file.path(root, "R"))
  writeLines("# a source file", file.path(root, "R", "one.R"))

  run <- suppressMessages(run_equivalence(root = root, check_api = FALSE))
  session <- as.data.frame(run, what = "session")
  expect_true(all(c("multilpa", "multilpa at end", "R/ fingerprint",
                    "source unchanged during run") %in% session$component))
  expect_equal(equivalence_provenance(run, "multilpa"),
               as.character(utils::packageVersion("multilpa")))
  expect_equal(equivalence_provenance(run, "source unchanged during run"), "yes")
  expect_error(equivalence_provenance(run, "not a component"),
               class = "multilpa_bad_argument")
})

test_that("a source edit during a run is reported, not hidden", {
  root <- file.path(scratch, "drift")
  suites <- file.path(root, "validation", "equivalence", "suites")
  dir.create(suites, recursive = TRUE)
  dir.create(file.path(root, "R"))
  writeLines("# before", file.path(root, "R", "one.R"))
  # The suite edits R/ while it runs, which is what eight agents working in one
  # tree do to a long harness run.
  writeLines(c(
    'suite_epsilon <- function() {',
    sprintf('  writeLines("# after", %s)', deparse(file.path(root, "R", "one.R"))),
    '  transform(compare_values("a quantity", 1, 1, tolerance = 1e-8),',
    '            source = "a reference")',
    '}'), file.path(suites, "epsilon.R"))

  run <- suppressMessages(run_equivalence(root = root, check_api = FALSE))
  expect_equal(equivalence_provenance(run, "source unchanged during run"), "NO")
  expect_false(identical(equivalence_provenance(run, "R/ fingerprint"),
                         equivalence_provenance(run, "R/ fingerprint at end")))
  expect_match(paste(.equivalence_markdown(run), collapse = "\n"),
               "package source changed while this run was in progress")
})

cat("\nHarness self-tests complete.\n")
