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

# Tagging: every assertion carries the reference it is evidence for, by the
# verification matrix's own rules (verification-matrix/attribution.R).

test_that("each kind of reference maps to its CVS axis and authority", {
  mapped <- .cvs_axis(
    c("poLCA", "Latent GOLD", "Mplus", "Mplus 9 Demo (JStats fixtures)", "Stata 19 [SEM] Example 52g (gsem)",
      "SAS PROC TRAJ (Andruff et al. 2009)", "Masyn (2013)", "numDeriv / dense likelihood", "internal identity", "tidySEM classification-error matrix",
      "kit provenance", "CVS reporter"),
    c("R package", "commercial program", "commercial program", "commercial program", "commercial program",
      "commercial program", "published", "independent code", "internal", "component check", "provenance", "infrastructure"))
  expect_identical(mapped$axis, c("r-equivalence", "latentgold", "mplus", "mplus", "stata", "sas", "paper",
                                  "independent", "internal", "internal", "internal", "internal"))
  expect_identical(mapped$authority, c("r", "latentgold-only", "mplus-only", "mplus-only", "stata-only",
                                       "sas-only", "paper", "independent", "internal", "internal", "internal", "internal"))
  expect_identical(mapped$software, c(NA, "latentgold", "Mplus", "Mplus", "stata", "sas", "paper", "independent",
                                      NA, NA, NA, NA))
  # A kind CVS has no axis for stays untagged rather than guessed.
  expect_true(is.na(.cvs_axis("StepMix", "Python package")$axis))
})

test_that("an unknown commercial program stops with a classed error", {
  expect_error(.cvs_axis("SPSS 29", "commercial program"), class = "latents_cvs_unknown_program")
  expect_error(.cvs_axis("Mplus and Latent GOLD", "commercial program"), class = "latents_cvs_unknown_program")
  expect_identical(.cvs_axis("IMMERSE lca-bookdown (Mplus 8.11)", "commercial program")$axis, "mplus")
  # Every commercial program latents' rule tables and matrix name as evidence has an
  # axis, so the error can only come from an edit of those tables, caught here and not
  # in a run. Candidates (references not yet wired in) are not evidence.
  rules <- cvs_rules()
  tests <- do.call(rbind, Map(matrix_rule_references, rules$tests$reference, rules$tests$kind))
  cells <- utils::read.csv(file.path("verification-matrix", "matrix.csv"), stringsAsFactors = FALSE)
  kits <- utils::read.csv(file.path("verification-matrix", "kit-evidence.csv"), stringsAsFactors = FALSE)
  every <- rbind(tests, rules$harness[c("reference", "kind")], cells[!cells$candidate, c("reference", "kind")],
                 kits[nzchar(kits$reference), c("reference", "kind")])
  commercial <- unique(every[every$kind == "commercial program", ])
  expect_gt(nrow(commercial), 3L)
  expect_false(anyNA(.cvs_axis(commercial$reference, commercial$kind)$axis))
})

cvs_tag_fixture <- function() {
  rules <- list(
    harness = data.frame(pattern = c("^poLCA", "^Mplus", "^invariance"),
                         reference = c("poLCA", "Mplus", "invariance"),
                         kind = c("R package", "commercial program", "internal"), stringsAsFactors = FALSE),
    tests = data.frame(
      dir = "equivalence", file = c("test-a.R", "test-a.R", "test-a.R", "test-b.R"),
      pattern = c("^the kit", "Stata", ".", "^published"),
      fn = c("none", "lpa;lca", "lpa", "lca"),
      reference = c("kit provenance", "Stata 19;Masyn (2013)", "glca;depmixS4", "Goodman (2002)"),
      kind = c("provenance", "commercial program;published", "R package", "published"),
      stringsAsFactors = FALSE))
  testthat_entry <- function(name, titles) {
    list(name = name, `_cvs` = list(tier = "equiv"), assertionResults = lapply(titles, function(title) {
      .cvs_assertion(basename(name), title, "passed", meta = list(checkKind = "testthat"))
    }))
  }
  harness <- .cvs_from_harness(
    data.frame(suite = "polca", source = c("poLCA 1.6.0.2", "Mplus 8.11", "invariance", "unmapped 1.0"),
               quantity = c("loglik", "loglik", "rows permuted", "x"), difference = 0, tolerance = 1e-6,
               agrees = TRUE, output = NA_character_, n_values = 1L, known_open = NA_character_,
               stringsAsFactors = FALSE),
    data.frame(suite = c("polca", "lcmm"), status = c("ran", "skipped"), reason = c(NA, "not installed"),
               stringsAsFactors = FALSE))
  sweep <- .cvs_from_sweep(data.frame(module = "gaussian", case = "a", quantity = "loglik", difference = 0,
                                      tolerance = 1e-8, status = "agrees", note = NA_character_,
                                      stringsAsFactors = FALSE))
  suites <- c(list(testthat_entry("equivalence/test-a.R",
                                  c("the kit is complete", "Stata and Masyn agree", "glca agrees", "glca agrees [2]")),
                   testthat_entry("equivalence/test-b.R", c("published table", "a test no rule names")),
                   testthat_entry("tests/testthat/test-a.R", "glca agrees")),
              harness, sweep)
  list(rules = rules, suites = suites)
}

cvs_tags <- function(suite) {
  keys <- c("axis", "authority", "reference")
  lapply(suite$assertionResults, function(a) stats::setNames(lapply(keys, function(k) a[["_cvs"]][[k]]), keys))
}

test_that("assertions are tagged by the matrix's rules, and the untaggable are left as they are", {
  fixture <- cvs_tag_fixture()
  tagged <- cvs_tag_suites(fixture$suites, fixture$rules)
  tags <- cvs_tags(tagged[[1L]])
  # Provenance is an internal check; a numbered repeat is attributed by its name.
  expect_identical(tags[[1L]], list(axis = "internal", authority = "internal", reference = "kit provenance"))
  expect_identical(tags[[3L]], list(axis = "r-equivalence", authority = "r", reference = "glca; depmixS4"))
  expect_identical(tags[[4L]], tags[[3L]])
  # References on two axes: the reference is kept, no axis is chosen.
  expect_identical(tags[[2L]], list(axis = NULL, authority = NULL, reference = "Stata 19; Masyn (2013)"))
  b <- cvs_tags(tagged[[2L]])
  expect_identical(b[[1L]], list(axis = "paper", authority = "paper", reference = "Goodman (2002)"))
  expect_identical(b[[2L]], list(axis = NULL, authority = NULL, reference = NULL))
  # The package's own tests have no attribution rule: untagged.
  expect_identical(tagged[[3L]], fixture$suites[[3L]])
  harness <- cvs_tags(tagged[[4L]])
  expect_identical(harness[[1L]], list(axis = "r-equivalence", authority = "r", reference = "poLCA"))
  expect_identical(harness[[2L]], list(axis = "mplus", authority = "mplus-only", reference = "Mplus"))
  expect_identical(harness[[3L]]$axis, "internal")
  expect_null(harness[[4L]]$reference)                     # a source no rule names
  expect_identical(tagged[[5L]], fixture$suites[[5L]])     # "(suite ran)" has no source
  # A sweep check carries the reference its run recorded (here, by the sweep's
  # rules: a plain log likelihood with no reference named is an internal check).
  expect_identical(cvs_tags(tagged[[6L]])[[1L]],
                   list(axis = "internal", authority = "internal", reference = "internal check"))
  expect_identical(cvs_axis_counts(tagged), data.frame(
    axis = c("(several axes)", "(unattributed)", "internal", "mplus", "paper", "r-equivalence"),
    assertions = c(1L, 4L, 3L, 1L, 1L, 3L), stringsAsFactors = FALSE))
})

test_that("tagging only adds fields: nothing is dropped, renamed or reordered", {
  fixture <- cvs_tag_fixture()
  tagged <- cvs_tag_suites(fixture$suites, fixture$rules)
  strip <- function(suites) lapply(suites, function(suite) {
    suite$assertionResults <- lapply(suite$assertionResults, function(a) {
      kept <- a[["_cvs"]][setdiff(names(a[["_cvs"]]), c("axis", "authority", "reference"))]
      a[["_cvs"]] <- if (length(kept)) kept else NULL
      a
    })
    suite
  })
  expect_identical(strip(tagged), strip(fixture$suites))
  expect_identical(sum(cvs_axis_counts(tagged)$assertions),
                   sum(vapply(fixture$suites, function(s) length(s$assertionResults), integer(1))))
  # Tagging is idempotent.
  expect_identical(cvs_tag_suites(tagged, fixture$rules), tagged)
})

test_that("a saved report is tagged again with its own name, run, start and version", {
  local <- withr::local_tempdir()
  withr::local_envvar(CVS_DISABLED = "1")
  suites <- .cvs_from_harness(
    data.frame(suite = "polca", source = "poLCA 1.6.0.2", quantity = "loglik", difference = 1e-9,
               tolerance = 1e-6, agrees = TRUE, output = NA_character_, n_values = 1L,
               known_open = NA_character_, stringsAsFactors = FALSE),
    data.frame(suite = "polca", status = "ran", reason = NA_character_, stringsAsFactors = FALSE))
  # A report as written before tagging.
  saved <- file.path(local, "latents-harness-20261112T192825-1.json")
  writeLines(jsonlite::toJSON(list(numTotalTests = 1L, startTime = 1.7629e12, testResults = suites,
                                   `_cvs` = list(target = "latents", run = "harness", version = "1.1.19")),
                              auto_unbox = TRUE, null = "null", digits = NA), saved)
  output <- file.path(local, "out")
  expect_message(path <- cvs_retag_report(saved, output = output), "axes: r-equivalence 1")
  expect_identical(basename(path), basename(saved))
  report <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  expect_identical(report[["_cvs"]], list(target = "latents", run = "harness", version = "1.1.19"))
  expect_identical(report$startTime, 1.7629e12)
  tag <- report$testResults[[1L]]$assertionResults[[1L]][["_cvs"]]
  expect_identical(tag[c("axis", "authority", "reference", "rPackage")],
                   list(axis = "r-equivalence", authority = "r", reference = "poLCA", rPackage = "poLCA 1.6.0.2"))
  expect_error(cvs_retag_report(file.path(local, "missing.json")), "existing file")
})

test_that("sweep checks reach CVS tagged with the reference their run recorded", {
  results <- data.frame(
    module = "gaussian", case = c("structure[data=iris,model=VVV,k=2]", "structure[data=iris,model=VVV,k=2]",
                                  "structure[data=iris,model=VVV,k=2]", "spelling[model=EEE]"),
    quantity = c("mclust's EM from the package's posteriors: means",
                 "log likelihood rebuilt from the result tables", "parameter count is mclust's",
                 "model code against volume/shape/orientation: log likelihood"),
    difference = 0, tolerance = 1e-6, n_values = 1L, status = "agrees", note = NA_character_,
    stringsAsFactors = FALSE)
  # Without reference columns (results saved before the sweep recorded them),
  # the matrix's sweep rules attribute each check.
  from_rules <- cvs_tags(cvs_tag_suites(.cvs_from_sweep(results))[[1L]])
  expect_identical(vapply(from_rules, `[[`, character(1), "axis"),
                   c("r-equivalence", "independent", "internal", "internal"))
  expect_identical(vapply(from_rules, `[[`, character(1), "reference"),
                   c("mclust", "likelihood rebuilt from the result tables", "mclust parameter counts",
                     "internal check"))
  # A run's own columns win over the rules.
  recorded <- cbind(results, fn = "lpa", reference = c("poLCA", "Rubin's rules written here",
                                                        "mclust parameter counts", "internal check"),
                    kind = c("R package", "independent code", "component check", "internal"),
                    stringsAsFactors = FALSE)
  from_run <- cvs_tags(cvs_tag_suites(.cvs_from_sweep(recorded))[[1L]])
  expect_identical(vapply(from_run, `[[`, character(1), "reference"), recorded$reference)
  expect_identical(from_run[[1L]]$axis, "r-equivalence")
})

test_that("the sweep's rules give every reference check a reference and stop on an unnamed check", {
  dir <- file.path("verification-matrix")
  calibration <- data.frame(
    module = c("gaussian", "categorical", "regression", "pooling", "pooling", "postfit", "postfit",
               "growth", "discrete", "families", "lta", "postfit", "plots", "covariates"),
    case = c("structure[data=srl,model=VVI,k=2]", "classes[data=carcinoma,k=1,complete]",
             "flexmix[data=NPreg,k=2]", "pooling[m=5]", "pooling[m=5]", "three-step[data=srl,k=2]",
             "parameter-inference[data=srl,k=2,model=VVI]", "one-class[data=scores,random=slope]",
             "negative-binomial[k=1]", "additive[g=2]", "transitions[data=x]",
             "bootstrap-and-sensitivity[data=srl]", "views", "class-regression-against-poLCA[k=2]"),
    quantity = c("optimum at least that of mclust's search", "one class: poLCA's log likelihood",
                 "flexmix's EM from the package's posteriors: log likelihood", "mitools::MIcombine: estimates",
                 "Rubin's rules: standard errors", "logits against stats::optim on the step-three likelihood",
                 "observed: standard errors of means and variances against numDeriv", "lme4's fixed effects",
                 "one profile: glm.nb's means", "group posteriors rebuilt from the result tables",
                 "log likelihood by the forward recursion from the result tables",
                 "sensitivity(): each seed's log likelihood is the direct fit's", "plot_profiles() draws, or is refused by class",
                 "poLCA's regression continued from the package's probabilities"),
    stringsAsFactors = FALSE)
  attributed <- matrix_sweep_attribution(calibration$module, calibration$case, calibration$quantity, dir)
  expect_identical(attributed$kind, c("R package", "component check", "R package", "R package",
                                      "independent code", "independent code", "independent code",
                                      "component check", "component check", "independent code",
                                      "independent code", "internal", "internal", "R package"))
  expect_identical(attributed$fn, c("lpa", "lca", "mixture_regression", "pool_imputations",
                                    "pool_imputations", "three_step", "lpa", "growth", "multilpa",
                                    "multilpa", "lta", "lpa", "", "lca"))
  # Every kind in the tables has a CVS axis.
  references <- matrix_read_rules("sweep-references.csv", dir)
  expect_false(anyNA(.cvs_axis(references$reference, references$kind)$axis))
  # A check no rule names stops the attribution, by class.
  local <- withr::local_tempdir()
  utils::write.csv(data.frame(module = "gaussian", pattern = "^structure", fn = "lpa"),
                   file.path(local, "sweep-functions.csv"), row.names = FALSE)
  utils::write.csv(references[references$pattern != ".", ], file.path(local, "sweep-references.csv"),
                   row.names = FALSE)
  expect_error(matrix_sweep_attribution("gaussian", "spelling[model=EEE]", "rows shuffled: log likelihood",
                                        local), class = "latents_sweep_unattributed")
})

test_that("every attribution rule of a reported test names references on one CVS axis", {
  # CVS takes one axis per assertion, so a block checking two references on
  # different axes (a package and a published table) counted for neither.
  # equivalence/ and tests/equivalence/ both report (jobs equivalence-*).
  rules <- cvs_rules()$tests
  expect_setequal(unique(rules$dir), c("equivalence", "tests/equivalence"))
  axes <- vapply(seq_len(nrow(rules)), function(r) {
    references <- matrix_rule_references(rules$reference[[r]], rules$kind[[r]])
    length(unique(.cvs_axis(references$reference, references$kind)$axis))
  }, integer(1))
  expect_identical(paste(rules$file, rules$pattern)[axes != 1L], character())
})
