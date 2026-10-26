# Tests for the harness itself.
#
# The equivalence suites test the package. Nothing tested the harness, which is
# how it came to sit broken for two releases while its saved report was still
# being cited. These run in seconds and fit no models.
#
# From the project root:
#   Rscript equivalence/harness/test-harness.R

suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("equivalence", "harness", "registry.R"))

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

test_that("check_validation_api skips R CMD check copies of the package tests", {
  # Archived *.Rcheck trees under validation/ hold the package's own tests,
  # which pass removed arguments on purpose; they stopped the harness.
  tree <- file.path(scratch, "api-root")
  check_copy <- file.path(tree, "validation", "old", "latents.Rcheck", "tests")
  dir.create(check_copy, recursive = TRUE)
  dir.create(file.path(tree, "equivalence"))
  stale_call <- 'multilpa(students, indicators = "x", id = "school", n_profiles = 2)'
  writeLines(stale_call, file.path(check_copy, "test-refusal.R"))
  expect_equal(nrow(check_validation_api(root = tree)), 0L)
  writeLines(stale_call, file.path(tree, "validation", "script.R"))
  expect_equal(check_validation_api(root = tree)$argument, "indicators")
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
  dir.create(file.path(scratch, "equivalence", "harness"), recursive = TRUE)
  file.rename(suites, file.path(scratch, "equivalence", "harness", "suites"))
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
                           \(seed) structure(-100, class = "logLik", converged = TRUE),
                           seeds = c(1L, 2L, 3L, 4L), tolerance = 1e-5)
  expect_equal(nrow(stable), 3L)
  expect_true(all(stable$agrees))
  expect_equal(stable$obtained, c(0, 1, 4))

  # Two of four seeds land on an inferior mode: both rows must fail.
  unstable <- seed_stability("an unstable model",
                             \(seed) structure(if (seed %% 2L == 0L) -100 else -101,
                                               class = "logLik", converged = TRUE),
                             seeds = c(1L, 2L, 3L, 4L), tolerance = 1e-5)
  expect_false(any(unstable$agrees[1:2]))
  expect_true(unstable$agrees[3L])
  expect_equal(unstable$obtained, c(1, 2, 4))
})

test_that("stability is independent of rounding bins and requires convergence", {
  stable <- seed_stability("bin edge", function(seed) {
    structure(-1e-5 * (100.5 + c(-0.01, 0.01, 0)[seed]),
              class = "logLik", converged = TRUE)
  }, 1:3, tolerance = 1e-5)
  expect_true(all(stable$agrees))
  expect_equal(stable$obtained, c(2e-7, 1, 3), tolerance = 1e-12)
  failed <- seed_stability("unconverged", function(seed) {
    structure(-100, class = "logLik", converged = FALSE)
  }, 1:3)
  expect_false(failed$agrees[3L])
  expect_equal(failed$obtained[3L], 0)
  fit <- function(seed) structure(-100, class = "logLik", converged = TRUE)
  expect_error(seed_stability("repeated", fit, c(1, 1, 1)), "distinct")
  expect_error(seed_stability("infinite", fit, c(1, 2, Inf)), "finite")
  expect_error(seed_stability("missing", fit, c(1, 2, NA_real_)), "finite")
  expect_error(seed_stability("fractional", fit, c(1, 2, 0.5)), "whole")
})

test_that("nonfinite tolerances and forged agreement cannot pass", {
  expect_error(compare_values("infinite", 0, Inf, tolerance = Inf))
  expect_false(compare_values("nonfinite", 0, Inf, tolerance = 1)$agrees)
  path <- file.path(scratch, "forged.R")
  writeLines(c('suite_forged <- function() {',
               ' x <- compare_values("wrong", 0, 10, tolerance=1e-5)',
               ' x$agrees <- TRUE; x$source <- "fake"; x', '}'), path)
  expect_false(.evaluate_suite("forged", path)$value$agrees)
})

# --------------------------------------------------------------------------
test_that("a run records the version and fingerprint it was produced under", {
  suites <- file.path(scratch, "provenance", "equivalence", "harness", "suites")
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
               as.character(utils::packageVersion("latents")))
  expect_equal(equivalence_provenance(run, "source unchanged during run"), "yes")
  expect_error(equivalence_provenance(run, "not a component"),
               class = "latents_bad_argument")
})

test_that("a source edit during a run is reported, not hidden", {
  root <- file.path(scratch, "drift")
  suites <- file.path(root, "equivalence", "harness", "suites")
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

# --------------------------------------------------------------------------
# Every output: compare_block, the ledger and the known-open register.

test_that("compare_block checks every cell and matches missing values exactly", {
  agreed <- compare_block("x", c(1, NA, Inf, 0.5), c(1 + 1e-9, NA, Inf, 0.5),
                          tolerance = 1e-6, output = "t$x")
  expect_true(agreed$agrees)
  expect_identical(agreed$n_values, 4L)
  expect_identical(agreed$output, "t$x")
  # One bad cell among many fails the whole column and is the one reported.
  bad <- compare_block("x", c(1, 2, 3), c(1, 2.5, 3), tolerance = 1e-6)
  expect_false(bad$agrees)
  expect_equal(c(bad$reference, bad$obtained), c(2, 2.5))
  # Missing in one program only, or infinite in one only, is a disagreement.
  expect_false(compare_block("x", c(1, NA), c(1, 1), tolerance = 1)$agrees)
  expect_false(compare_block("x", c(1, Inf), c(1, 1e300), tolerance = 1)$agrees)
  expect_false(compare_block("x", Inf, -Inf, tolerance = 1)$agrees)
  # Logical columns match exactly, whatever tolerance is asked for.
  expect_false(compare_block("x", c(TRUE, FALSE), c(TRUE, TRUE), tolerance = 10)$agrees)
  # A column missing in both survives the finiteness check of .evaluate_suite.
  path <- file.path(scratch, "missing-both.R")
  writeLines(c('suite_missing_both <- function() {',
               '  x <- compare_block("x", c(NA, NA), c(NA, NA), tolerance = 1e-6)',
               '  x$source <- "a reference"; x', '}'), path)
  expect_true(.evaluate_suite("missing_both", path)$value$agrees)
})

test_that("compare_tables finds a missing row and a wrong cell", {
  have <- list(t = data.frame(profile = 1:3, value = c(0.1, 0.2, 0.3)))
  same <- compare_tables("case", have, have, c(default = 1e-8))
  expect_true(all(same$agrees))
  lost <- compare_tables("case", list(t = have$t[1:2, ]), have, c(default = 1e-8))
  expect_false(all(lost$agrees))
  planted <- have
  planted$t$value[[3L]] <- 0.3 + 1e-4
  wrong <- compare_tables("case", planted, have, c(default = 1e-8))
  expect_false(wrong$agrees[wrong$output %in% "t$value"])
})

test_that("the ledger accounts for every output and refuses a stale claim", {
  inventory <- data.frame(output = c("t$a", "t$b", "t$c", "u$d"), table = c("t", "t", "t", "u"),
                          column = c("a", "b", "c", "d"), n_values = c(4L, 4L, 4L, 2L),
                          stringsAsFactors = FALSE)
  rows <- rbind(compare_block("a", 1:4, 1:4, 1e-8, output = "t$a"),
                compare_block("b", 1:2, 1:2, 1e-8, output = "t$b"))
  ledger <- output_ledger("fit", inventory, rows, exempt = c(u = "no reference"))
  status <- stats::setNames(ledger$status, ledger$output)
  expect_identical(unname(status[c("t$a", "t$b", "t$c", "u$d")]),
                   c("compared", "partial", "not covered", "exempt"))
  expect_error(output_ledger("fit", inventory, rows, exempt = c(`t$z` = "gone")),
               class = "multilpa_bad_ledger")
  expect_error(output_ledger("fit", inventory,
                             compare_block("z", 1, 1, 1e-8, output = "t$typo")),
               class = "multilpa_bad_ledger")
})

test_that("a known open defect is reported by id, and a stale entry is caught", {
  comparisons <- data.frame(suite = c("s", "s", "s"), quantity = c("a: x", "b: y", "c: z"),
                            agrees = c(FALSE, FALSE, TRUE), stringsAsFactors = FALSE)
  register <- data.frame(id = c("K-1", "K-2", "K-3"), suite = c("s", "s", "other"),
                         pattern = c("^a:", "^c:", "."), defect = "d", decision = "e",
                         stringsAsFactors = FALSE)
  marked <- .mark_known_open(comparisons, register, suites = "s")
  expect_identical(marked$comparisons$known_open, c("K-1", NA, NA))
  expect_identical(marked$register$status, c("open", "stale", "not run"))
  shared <- register[1L, ]
  shared$suite <- "other|s"
  both <- .mark_known_open(comparisons, shared, suites = c("s", "other"))
  expect_identical(both$comparisons$known_open, c("K-1", NA, NA))
  expect_identical(.mark_known_open(comparisons, shared, suites = "s")$register$status,
                   "not run")
  empty <- .mark_known_open(comparisons, register[0, ], suites = "s")
  expect_true(all(is.na(empty$comparisons$known_open)))
})

test_that("a column named like a key is a value where the table says so", {
  sequences <- data.frame(group = "a", time = 1:2, profile = 1L)
  expect_true("time" %in% .output_keys(sequences, "sequences"))
  effects <- data.frame(person = c("a", "b"), class = 1L, `(Intercept)` = c(0.1, 0.2),
                        time = c(0.3, 0.4), check.names = FALSE)
  expect_false("time" %in% .output_keys(effects, "random_effects"))
  expect_true("random_effects$time" %in% fit_outputs(list(random_effects = effects))$output)
})

test_that("a fit's inventory is read from the package's own catalogue", {
  data <- datasets::faithful[seq_len(60L), ]
  data$unit <- seq_len(nrow(data))
  fit <- suppressWarnings(multilpa(data, vars = c("eruptions", "waiting"), id = "unit",
                                   n_profiles = 2, n_group_classes = 1, n_starts = 2,
                                   seed = 1))
  inventory <- fit_outputs(fit)
  expect_true(all(c("posteriors$posterior", "bch_weights$weight",
                    "classification_errors$probability",
                    "profiles$mean_standard_error") %in% inventory$output))
  # The echoed input is not an output.
  expect_false(any(c("assignments$eruptions", "assignments$unit") %in% inventory$output))
})

test_that("exemptions drop what they name, and no exemption drops nothing", {
  expected <- list(t = data.frame(profile = 1:2, a = 1:2, b = 3:4), u = data.frame(k = 1, v = 2))
  expect_identical(.drop_exempt(expected, character()), expected)
  dropped <- .drop_exempt(expected, c(`t$b` = "reason", u = "reason"))
  expect_identical(names(dropped), "t")
  expect_identical(names(dropped$t), c("profile", "a"))
})

test_that("an emptied register is a register with no entries", {
  root <- file.path(scratch, "empty-register")
  dir.create(file.path(root, "equivalence", "harness"), recursive = TRUE)
  writeLines('"id","suite","pattern","defect","decision"',
             file.path(root, "equivalence", "harness", "known-open.csv"))
  register <- .known_open(root)
  expect_identical(nrow(register), 0L)
  expect_true(all(vapply(register, is.character, logical(1))))
  comparisons <- data.frame(suite = "s", quantity = "q", agrees = FALSE,
                            stringsAsFactors = FALSE)
  expect_true(is.na(.mark_known_open(comparisons, register, "s")$comparisons$known_open))
})

test_that("a failed suite leaves its known-open entries unjudged", {
  root <- file.path(scratch, "register")
  suites <- file.path(root, "equivalence", "harness", "suites")
  dir.create(suites, recursive = TRUE)
  dir.create(file.path(root, "R"))
  writeLines("# no sources", file.path(root, "R", "empty.R"))
  writeLines('suite_broken <- function() stop("broken")', file.path(suites, "broken.R"))
  writeLines(c('suite_fine <- function() {',
               '  transform(compare_values("q", 1, 1, tolerance = 1e-8), source = "r")',
               '}'), file.path(suites, "fine.R"))
  utils::write.csv(data.frame(id = c("A", "B"), suite = c("broken", "fine"),
                              pattern = c("x", "y"), defect = "d", decision = "e"),
                   file.path(root, "equivalence", "harness", "known-open.csv"),
                   row.names = FALSE)
  run <- suppressMessages(run_equivalence(root = root, check_api = FALSE))
  register <- as.data.frame(run, what = "known_open")
  expect_identical(register$status, c("not run", "stale"))
})

# The entropy charge of ICL, AWE and CLC (selection review round 2, S-1): the
# complete-data entropy H(group classes) + sum_i E_C[H(X_i | C, y_i)], the
# same under both conventions.
harness_entropy_case <- function(group_classes, row_entropy = NULL) {
  set.seed(3)
  n <- 12L
  groups <- rep(seq_len(4L), each = 3L)
  raw <- matrix(stats::runif(n * 2L), n)
  posteriors <- raw / rowSums(raw)
  group_posteriors <- if (group_classes > 1L) {
    zg <- matrix(stats::runif(4L * group_classes), 4L)
    zg <- zg / rowSums(zg)
    rownames(zg) <- as.character(seq_len(4L))
    zg
  }
  reference <- harness_reference(
    log_likelihood = -40, n_parameters = 5, posteriors = posteriors, groups = groups,
    profile_probabilities = matrix(0.5, group_classes, 2L),
    group_probabilities = rep(1 / group_classes, group_classes),
    group_posteriors = group_posteriors, row_conditional_entropy = row_entropy)
  fit <- structure(list(vars = character()), class = "multilpa")
  list(tables = reference_tables(reference, fit), posteriors = posteriors,
       group_posteriors = group_posteriors)
}

test_that("one group class: both conventions charge the row posterior entropy", {
  case <- harness_entropy_case(1L)
  criteria <- case$tables$information_criteria
  entropy <- -sum(case$posteriors * log(case$posteriors))
  expect_equal(criteria$icl_individual, 80 + 5 * log(12) + 2 * entropy)
  expect_equal(criteria$icl_groups, 80 + 5 * log(4) + 2 * entropy)
  expect_equal(criteria$clc_groups, criteria$clc_individual)
  expect_equal(criteria$awe_groups - criteria$awe_individual, 2 * 5 * (log(4) - log(12)))
  expect_identical(attr(case$tables, "exempt"), character())
  # Supplying the row entropies changes nothing when there is one group class.
  supplied <- harness_entropy_case(1L, -rowSums(case$posteriors * log(case$posteriors)))
  expect_identical(supplied$tables$information_criteria, criteria)
})

test_that("several group classes: the complete-data entropy, or an exemption by name", {
  row_entropy <- seq(0.1, 1.2, length.out = 12L)
  case <- harness_entropy_case(2L, row_entropy)
  criteria <- case$tables$information_criteria
  complete <- -sum(case$group_posteriors * log(case$group_posteriors)) + sum(row_entropy)
  expect_equal(criteria$clc_groups, 80 + 2 * complete)
  expect_equal(criteria$clc_individual, 80 + 2 * complete)
  expect_equal(criteria$icl_groups, 80 + 5 * log(4) + 2 * complete)
  expect_identical(attr(case$tables, "exempt"), character())
  # Without them the three criteria are left out, with a reason, never filled
  # in with another entropy.
  missing <- harness_entropy_case(2L)
  charged <- paste0("information_criteria$",
                    c(outer(c("awe", "icl", "clc"), c("groups", "individual"), paste,
                            sep = "_")))
  expect_setequal(names(attr(missing$tables, "exempt")), charged)
  expect_true(all(nzchar(attr(missing$tables, "exempt"))))
  left <- names(missing$tables$information_criteria)
  expect_false(any(sub("^information_criteria[$]", "", charged) %in% left))
  expect_true(all(c("bic_groups", "bic_individual", "aic") %in% left))
})

cat("\nHarness self-tests complete.\n")
