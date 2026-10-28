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

test_that("the runner rejects partial coverage and required-suite skips", {
  inventory <- data.frame(output = "t$value", table = "t", column = "value",
                          n_values = 3L)
  compared <- compare_block("one of three", 1, 1, 1e-8, output = "t$value")
  run <- list(suites = data.frame(status = "ran"),
              comparisons = transform(compared, known_open = NA_character_),
              coverage = output_ledger("fit", inventory, compared),
              known_open = data.frame(status = character()))
  expect_identical(run$coverage$status, "partial")
  expect_false(.equivalence_passes(run))
  all_values <- compare_block("all three", 1:3, 1:3, 1e-8, output = "t$value")
  run$comparisons <- transform(all_values, known_open = NA_character_)
  run$coverage <- output_ledger("fit", inventory, all_values)
  expect_true(.equivalence_passes(run, require_all = TRUE))
  run$suites$status <- "skipped"
  expect_true(.equivalence_passes(run))
  expect_false(.equivalence_passes(run, require_all = TRUE))
  run$suites$status <- "failed"
  expect_false(.equivalence_passes(run))
  run$suites$status <- "ran"
  run$comparisons$agrees <- FALSE
  expect_false(.equivalence_passes(run))
  run$comparisons$known_open <- "registered"
  expect_true(.equivalence_passes(run))
  run$known_open <- data.frame(status = "stale")
  expect_false(.equivalence_passes(run))
})

test_that("replication likelihood is independently checked, not exempted", {
  obtained <- list(model = data.frame(log_likelihood = -10),
                   replication = data.frame(n_starts = 4L, log_likelihood = -11))
  expected <- list(model = data.frame(log_likelihood = -10))
  compared <- compare_tables("replication", obtained, expected)
  likelihood <- subset(compared, output %in% "replication$log_likelihood")
  expect_equal(nrow(likelihood), 1L)
  expect_false(likelihood$agrees)
  ledger <- output_ledger("replication", fit_outputs(obtained), compared)
  expect_identical(ledger$status[ledger$output == "replication$n_starts"], "exempt")
  expect_identical(ledger$status[ledger$output == "replication$log_likelihood"], "compared")
  expect_false("replication$log_likelihood" %in% names(.output_exemptions))
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

test_that("mixed residual references cover every pair using weighted category codes", {
  data <- data.frame(category = factor(rep(c("z", "a", "b"), 4L),
                                        levels = c("z", "a", "b")),
                     y1 = as.numeric(seq_len(12L)),
                     y2 = c(3, 1, 5, 2, 8, 4, 7, 6, 9, 12, 10, 11))
  data$y1[[2L]] <- NA_real_
  data$y2[[5L]] <- NA_real_
  probability <- rep(c(0.25, 0.75), 6L)
  reference <- list(data = data, posteriors = cbind(probability, 1 - probability),
                     means = matrix(0, 2L, 2L, dimnames = list(NULL, c("y1", "y2"))),
                     variances = matrix(1, 2L, 2L, dimnames = list(NULL, c("y1", "y2"))))
  residuals <- .reference_residuals(reference, c("y1", "y2"), "category",
                                     list(category = levels(data$category)))
  expect_equal(nrow(residuals), 6L)
  mixed <- subset(residuals, kind == "mixed")
  expect_equal(nrow(mixed), 4L)
  expect_identical(unique(mixed$indicator_1), "category")
  expect_setequal(mixed$indicator_2, c("y1", "y2"))
  # An independent oracle: quarter-unit posterior weights become integer
  # replications, whose ordinary Pearson correlation must agree exactly.
  invisible(lapply(seq_len(nrow(mixed)), function(index) {
    row <- mixed[index, ]
    state <- as.integer(sub("profile_", "", row$profile))
    complete <- stats::complete.cases(data[c("category", row$indicator_2)])
    rows <- rep(which(complete), as.integer(4 * reference$posteriors[complete, state]))
    correlation <- stats::cor(as.integer(data$category[rows]), data[[row$indicator_2]][rows])
    expect_equal(row$observed, correlation, tolerance = 1e-14)
    expect_equal(row$residual, correlation, tolerance = 1e-14)
    expect_equal(row$effective_n, sum(reference$posteriors[complete, state]))
  }))
  expect_true(all(mixed$expected == 0))
  expect_true(all(is.na(as.matrix(mixed[c("statistic", "df", "p_value", "p_adjusted")]))))
  # No Gaussian pair exists here: the mixed reference must still return its
  # own rows, including an explicit unavailable correlation for empty data.
  reference$data$y1[] <- NA_real_
  empty <- .reference_residuals(reference, "y1", "category",
                                list(category = levels(data$category)))
  expect_equal(nrow(empty), 2L)
  expect_true(all(is.na(empty$residual)))
  expect_true(all(empty$effective_n == 0))
})

# --------------------------------------------------------------------------
# Verb tables (EQ-20): three_step(), r3step(), pool_imputations(), predict(),
# parameter_inference() and the enumeration's lmr table are not a fit's
# catalogue, so they reach the gate through table_ledger() / compare_verb().

# A run as run_equivalence() assembles it, for .equivalence_passes().
harness_run_of <- function(case) {
  list(suites = data.frame(status = "ran"),
       comparisons = transform(case$comparisons, known_open = NA_character_),
       coverage = case$ledger, known_open = data.frame(status = character()))
}
# A verb result shaped like three_step()'s, with an integer key column.
verb_obtained <- list(three_step = keyed_table(
  data.frame(level = "individuals", class = 1:3, estimate = c(0.5, 1.5, 2.5),
             standard_error = c(0.1, 0.2, 0.3)), c("level", "class")))

test_that("a verb column nobody compared fails the gate", {
  expected <- list(three_step = data.frame(level = "individuals", class = 1:3,
                                           estimate = c(0.5, 1.5, 2.5)))
  case <- compare_verb("three step", verb_obtained, expected)
  status <- stats::setNames(case$ledger$status, case$ledger$output)
  expect_identical(unname(status[c("three_step$estimate", "three_step$standard_error")]),
                   c("compared", "not covered"))
  expect_true(all(case$comparisons$agrees))
  expect_false(.equivalence_passes(harness_run_of(case)))
  # Exempted with a reason, the same column passes; every value is accounted.
  exempted <- compare_verb("three step", verb_obtained, expected,
                           exempt = c(`three_step$standard_error` = "no reference SE"))
  expect_identical(exempted$ledger$status[exempted$ledger$output == "three_step$standard_error"],
                   "exempt")
  expect_true(.equivalence_passes(harness_run_of(exempted), require_all = TRUE))
})

test_that("a partly compared verb column fails the gate", {
  # The reference covers two of the three classes: the third is unchecked.
  expected <- list(three_step = data.frame(level = "individuals", class = 1:2,
                                           estimate = c(0.5, 1.5),
                                           standard_error = c(0.1, 0.2)))
  case <- compare_verb("three step", verb_obtained, expected)
  expect_true(all(case$comparisons$agrees))
  expect_true(all(case$ledger$status == "partial"))
  expect_false(.equivalence_passes(harness_run_of(case)))
})

test_that("a sabotaged verb reference fails the gate", {
  expected <- list(three_step = data.frame(level = "individuals", class = 1:3,
                                           estimate = c(0.5, 1.5, 2.5),
                                           standard_error = c(0.1, 0.2, 0.3)))
  faithful_case <- compare_verb("three step", verb_obtained, expected)
  expect_true(.equivalence_passes(harness_run_of(faithful_case), require_all = TRUE))
  sabotaged <- expected
  sabotaged$three_step$standard_error[[3L]] <- 0.3 + 1e-4
  case <- compare_verb("three step", verb_obtained, sabotaged)
  expect_false(case$comparisons$agrees[case$comparisons$output %in% "three_step$standard_error"])
  expect_identical(unique(case$ledger$status), "compared")
  expect_false(.equivalence_passes(harness_run_of(case)))
  # A text value under declared keys is compared exactly.
  labelled <- list(predict = keyed_table(data.frame(row = 1:2, profile = c("profile_1", "profile_2")),
                                         "row"))
  wrong <- list(predict = data.frame(row = 1:2, profile = c("profile_1", "profile_1")))
  text <- compare_verb("predict", labelled, wrong)
  expect_false(text$comparisons$agrees[text$comparisons$output %in% "predict$profile"])
  expect_match(text$comparisons$quantity[text$comparisons$output %in% "predict$profile"],
               "first mismatch \"profile_1\" against \"profile_2\"", fixed = TRUE)
})

test_that("an uncovered verb column fails a whole run, not just its ledger", {
  root <- file.path(scratch, "verb-run")
  suites <- file.path(root, "equivalence", "harness", "suites")
  dir.create(suites, recursive = TRUE)
  dir.create(file.path(root, "R"))
  writeLines("# no sources", file.path(root, "R", "empty.R"))
  writeLines(c(
    'suite_verb <- function() {',
    '  obtained <- list(r3step = keyed_table(data.frame(term = c("(Intercept)", "x"),',
    '    estimate = c(-1, 2), standard_error = c(0.5, 0.4)), "term"))',
    '  expected <- list(r3step = data.frame(term = c("(Intercept)", "x"), estimate = c(-1, 2)))',
    '  case <- compare_verb("r3step", obtained, expected)',
    '  case$comparisons$source <- "a reference"',
    '  case',
    '}'), file.path(suites, "verb.R"))
  run <- suppressMessages(run_equivalence(root = root, check_api = FALSE))
  coverage <- as.data.frame(run, what = "coverage")
  expect_identical(coverage$status[coverage$output == "r3step$standard_error"], "not covered")
  expect_identical(as.data.frame(run, what = "suites")$outputs_not_covered, 1L)
  expect_false(.equivalence_passes(run))
})

test_that("declared keys make every other column a value, and must identify rows", {
  # predict()'s `profile` is the modal class: a value, which the catalogue
  # rules would take for a key.
  prediction <- data.frame(row = 1:3, unit = 4:6, profile = c(1L, 2L, 1L),
                           posterior = c(0.9, 0.8, 0.7))
  expect_false("predict$profile" %in% fit_outputs(list(predict = prediction))$output)
  keyed <- keyed_table(prediction, "row")
  expect_setequal(fit_outputs(list(predict = keyed))$output,
                  c("predict$unit", "predict$profile", "predict$posterior"))
  # Classes and attributes of a verb's result are dropped; factors become text.
  styled <- structure(transform(prediction, label = factor(c("a", "b", "a"))),
                      class = c("latents_table", "data.frame"), extra = "x")
  plain <- keyed_table(styled, "row")
  expect_identical(class(plain), "data.frame")
  expect_null(attr(plain, "extra"))
  expect_type(plain$label, "character")
  expect_error(keyed_table(prediction, "unit_typo"), class = "multilpa_bad_ledger")
  expect_error(keyed_table(rbind(prediction, prediction), "profile"),
               class = "multilpa_bad_ledger")
  # The expected side must be keyed as uniquely as the obtained side.
  doubled <- list(predict = rbind(prediction, prediction))
  expect_error(compare_verb("predict", list(predict = keyed), doubled),
               class = "multilpa_bad_ledger")
  expect_error(table_ledger("v", list(predict = keyed[0, ]), compare_values("q", 1, 1, 1)),
               class = "multilpa_bad_ledger")
  expect_error(table_ledger("v", keyed, compare_values("q", 1, 1, 1)),
               class = "multilpa_bad_ledger")
})

# --------------------------------------------------------------------------
# Reference tables for weights, coefficients, ordinal and count indicators.

test_that("the weights branch reproduces a weighted fit from mclust::me.weighted", {
  skip_if_not_installed("mclust")
  data <- datasets::faithful
  set.seed(20)
  data$w <- round(stats::runif(nrow(data), 0.2, 3), 2)
  # Three units of weight zero: in no count, entropy or sample size.
  data$w[c(5L, 50L, 150L)] <- 0
  data$unit <- seq_len(nrow(data))
  vars <- c("eruptions", "waiting")
  fit <- multilpa(data, vars = vars, id = "unit", n_profiles = 2, n_group_classes = 1,
                  weights = "w", n_starts = 4, seed = 1, tol = 1e-14, max_iter = 20000)
  x <- as.matrix(data[vars])
  # mclust's weighted EM (diagonal, varying: VVI), from a split of the data,
  # run to its fixed point.
  start <- cbind(as.numeric(data$eruptions <= 3), as.numeric(data$eruptions > 3))
  em <- mclust::me.weighted(x, modelName = "VVI", z = start, weights = data$w,
                            control = mclust::emControl(tol = c(1e-14, 1e-14),
                                                        itmax = c(1e5, 1e5)))
  means <- t(em$parameters$mean)
  variances <- t(apply(em$parameters$variance$sigma, 3L, diag))
  colnames(means) <- colnames(variances) <- vars
  # Densities and posteriors at mclust's estimates, in base R. me.weighted's
  # own `loglik` is on none of the package's scales, so the pseudo log
  # likelihood sum_j w_j log L_j is taken here with the weights rescaled to
  # sum to the number of units of positive weight (?multilpa, `weights`).
  joint <- vapply(1:2, function(k) {
    log(em$parameters$pro[k]) +
      rowSums(stats::dnorm(x, matrix(means[k, ], nrow(x), 2L, byrow = TRUE),
                           matrix(sqrt(variances[k, ]), nrow(x), 2L, byrow = TRUE),
                           log = TRUE))
  }, numeric(nrow(x)))
  top <- apply(joint, 1L, max)
  log_density <- top + log(rowSums(exp(joint - top)))
  rescaled <- data$w * sum(data$w > 0) / sum(data$w)
  build <- function(weights) {
    harness_reference(sum(rescaled * log_density), 9, exp(joint - log_density),
                      means = means, variances = variances,
                      profile_probabilities = matrix(em$parameters$pro, 1L),
                      group_probabilities = 1, data = data[vars], weights = weights)
  }
  no_se <- c(`profiles$mean_standard_error` = "no standard errors in this self-test",
             `profiles$variance_standard_error` = "no standard errors in this self-test")
  case <- compare_fit("faithful, weighted", fit, build(data$w), exempt = no_se)
  # Both EMs run to 1e-14: estimates agree to about 1e-8, inside the 1e-6
  # default and the 1e-4 OCC bound (outputs.R).
  expect_true(all(case$comparisons$agrees))
  expect_false(any(case$ledger$status %in% c("partial", "not covered")))
  expect_true(.equivalence_passes(harness_run_of(case), require_all = TRUE))
  counts <- reference_tables(align_reference(build(data$w), fit), fit)$counts
  expect_equal(sum(counts$effective_count[counts$level == "individuals"]), 269)
  expect_true(all(c("counts$effective_count", "classification$estimated_n",
                    "bch_weights$weight", "entropy$entropy_sum") %in%
                    case$ledger$output[case$ledger$status == "compared"]))
  # Sabotage: the reference ignoring the weights must fail, on counts too.
  unweighted <- compare_fit("faithful, weights ignored", fit,
                            build(rep(1, nrow(data))), exempt = no_se)
  expect_false(all(unweighted$comparisons$agrees[
    unweighted$comparisons$output %in% "counts$effective_count"]))
  expect_false(.equivalence_passes(harness_run_of(unweighted)))
  # A weighted fit with no reference weights is refused, not compared.
  unlabelled <- build(data$w)
  unlabelled$weights <- NULL
  expect_error(reference_tables(unlabelled, fit), class = "multilpa_bad_reference")
})

test_that("membership logits are re-expressed against the fit's last class", {
  # A reference whose classes are in another order and whose logits are
  # against its first class; terms name its own group-class numbering.
  logits <- cbind(profile_1 = c(0, 0, 0, 0), profile_2 = c(1, -2, 0.5, 0.25),
                  profile_3 = c(-1, 3, 0.75, -0.5))
  rownames(logits) <- c("group_class_1", "group_class_2", "x:group_class_1", "x:group_class_2")
  table <- .reference_coefficients(list(profile = logits))
  expect_identical(unique(table$outcome), c("profile_1", "profile_2"))
  expect_equal(table$estimate[table$outcome == "profile_1"], logits[, 1L] - logits[, 3L],
               ignore_attr = TRUE)
  # Invariance: adding a constant per term (another baseline) changes nothing.
  shifted <- .reference_coefficients(list(profile = logits + c(5, -1, 2, 7)))
  expect_equal(shifted$estimate, table$estimate, tolerance = 1e-14)
  # Group classes 1 and 2 swapped: the terms are renumbered with them.
  expect_identical(.renumber_group_classes(rownames(logits), c(2L, 1L)),
                   c("group_class_2", "group_class_1", "x:group_class_2", "x:group_class_1"))
  expect_identical(.reference_covariate_count(logits), 1L)
  expect_identical(.reference_covariate_count(
    matrix(c(0, 0, 1, 2), 2L, dimnames = list(c("(Intercept)", "z"), c("a", "b")))), 1L)
})

test_that("ordinal locations and count tables follow their definitions", {
  intercepts <- c(0.4, -0.3, 0.2)
  locations <- c(-1.2, 0.6, 0)
  probabilities <- .ordinal_probabilities(intercepts, locations)
  expect_equal(rowSums(probabilities), rep(1, 3L), tolerance = 1e-15)
  colnames(probabilities) <- c("low", "mid", "high", "top")
  count_means <- matrix(c(1.5, 4, 9), 3L, dimnames = list(NULL, "visits"))
  reference <- harness_reference(
    -10, 8, matrix(1 / 3, 4L, 3L), ordinal = list(rating = probabilities),
    count_means = count_means, count_dispersion = count_means / 10,
    standard_errors = list(location = matrix(c(0.1, 0.2, NA), 3L,
                                             dimnames = list(NULL, "rating")),
                           count_mean = count_means / 100))
  ordinal <- .reference_ordinal(reference, 3L)
  expect_equal(nrow(ordinal), 12L)
  # The location is recovered from the probabilities alone, the last at 0.
  expect_equal(ordinal$location[ordinal$category == "low"], locations, tolerance = 1e-12)
  expect_identical(ordinal$location[ordinal$profile == 3L], rep(0, 4L))
  expect_equal(ordinal$location_standard_error[ordinal$category == "top"], c(0.1, 0.2, NA))
  counts <- .reference_counts(reference, 3L)
  expect_identical(names(counts), c("profile", "indicator", "mean", "dispersion",
                                    "standard_error"))
  expect_equal(counts$standard_error, c(1.5, 4, 9) / 100)
})

test_that("one group class reports no group-level classification", {
  case <- harness_entropy_case(1L)
  expect_null(case$tables$group_posteriors)
  expect_identical(unique(case$tables$classification$level), "individuals")
  expect_identical(unique(case$tables$bch_weights$level), "individuals")
  expect_setequal(case$tables$entropy$level, c("individuals", "groups"))
  several <- harness_entropy_case(2L, seq(0.1, 1.2, length.out = 12L))
  expect_setequal(several$tables$classification$level, c("individuals", "groups"))
  expect_false(is.null(several$tables$group_posteriors))
})

# --------------------------------------------------------------------------
# The additive, cross-level and general transition families.

test_that("labelled-class families have their modal labels and bookkeeping accounted", {
  set.seed(5)
  teams <- 24L
  type <- rep(1:2, length.out = teams)
  data <- do.call(rbind, lapply(seq_len(teams), function(j) {
    k <- sample(1:2, 5L, replace = TRUE, prob = if (type[j] == 1L) c(0.85, 0.15) else c(0.15, 0.85))
    data.frame(team = sprintf("t%02d", j), y1 = stats::rnorm(5L, c(-1.5, 1.5)[k], 0.6),
               y2 = stats::rnorm(5L, c(-1, 1)[k], 0.6))
  }))
  fit <- multilpa(data, c("y1", "y2"), "team", n_profiles = 2, n_group_classes = 2,
                  family = "full_cross_level", seed = 1, n_starts = 1)
  tables <- fit_tables(fit)
  inventory <- fit_outputs(tables)
  # The modal classes are text labels: values, compared exactly.
  expect_true(all(c("assignments$profile", "assignments$group_class",
                    "groups$group_class") %in% inventory$output))
  # Bookkeeping in the one-row `fit` table is exempt by default.
  ledger <- output_ledger("cross-level", inventory, compare_values("q", 1, 1, 1))
  expect_identical(ledger$status[ledger$output %in% c("fit$converged", "fit$boundary")],
                   c("exempt", "exempt"))
  # A suite's own tables go through the same ledger, and a wrong label fails.
  expected <- list(assignments = tables$assignments)
  swapped <- expected
  swapped$assignments$profile <- ifelse(swapped$assignments$profile == "profile_1",
                                        "profile_2", "profile_1")
  same <- compare_fit_tables("cross-level", fit, expected)
  wrong <- compare_fit_tables("cross-level", fit, swapped)
  expect_true(same$comparisons$agrees[same$comparisons$output %in% "assignments$profile"])
  expect_false(wrong$comparisons$agrees[wrong$comparisons$output %in% "assignments$profile"])
  expect_identical(same$ledger$status[same$ledger$output == "composition$share"],
                   "not covered")
  expect_false(.equivalence_passes(harness_run_of(same)))
})

test_that("a row whose key is missing is found and compared, not counted absent", {
  # parameter_inference() gives group-class probabilities `term = NA`.
  obtained <- list(pi = keyed_table(
    data.frame(parameter = c("pi", "pi"), term = c("x", NA), estimate = c(1, 2)),
    c("parameter", "term")))
  close <- list(pi = data.frame(parameter = c("pi", "pi"), term = c("x", NA),
                                estimate = c(1, 2) + 1e-12))
  case <- compare_verb("na key", obtained, close, tolerance = c(default = 1e-9))
  expect_true(.equivalence_passes(harness_run_of(case)))
  # The NA row is still compared: a wrong value in it fails.
  wrong <- close
  wrong$pi$estimate[2L] <- 3
  expect_false(.equivalence_passes(harness_run_of(
    compare_verb("na key", obtained, wrong, tolerance = c(default = 1e-9)))))
  # And a reference row the fit lacks is still reported missing.
  extra <- list(pi = rbind(close$pi, data.frame(parameter = "pi", term = "z",
                                                estimate = 4)))
  expect_false(.equivalence_passes(harness_run_of(
    compare_verb("na key", obtained, extra, tolerance = c(default = 1e-9)))))
})

cat("\nHarness self-tests complete.\n")
