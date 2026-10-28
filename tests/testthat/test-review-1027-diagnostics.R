# Review 2026-10-27: diagnostics() as a structured report. One row per check in
# five domains, each value read from a table the fit already computes, judged
# against a cited criterion. It used to print five unrelated numbers, printed
# NA for a one-class fit and for an empty group class, and refused the
# additive, cross-level and mixture-regression families.

diagnostic_data <- function() {
  withr::local_seed(11)
  n <- 600L
  state <- sample(3L, n, TRUE)
  data <- data.frame(id = rep(seq_len(60L), each = 10L), t = rep(1:10, 60L),
                     y1 = stats::rnorm(n, c(-2, 0, 2)[state]),
                     y2 = stats::rnorm(n, c(1, -1, 0)[state]),
                     y3 = stats::rnorm(n, c(0, 2, -2)[state]),
                     x = stats::rnorm(n), w = stats::runif(n, 0.5, 2))
  probabilities <- rbind(c(0.7, 0.2, 0.1), c(0.1, 0.7, 0.2), c(0.2, 0.1, 0.7))
  invisible(lapply(c("c1", "c2", "c3", "c4"), function(name) {
    data[[name]] <<- factor(vapply(state, function(k) sample(3L, 1L, prob = probabilities[k, ]),
                                   integer(1)))
  }))
  data
}

diagnostic_fits <- function() {
  data <- diagnostic_data()
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  continuous <- c("y1", "y2", "y3")
  categorical <- c("c1", "c2", "c3", "c4")
  list(
    lpa = quiet(lpa(data, continuous, n_profiles = 3, n_starts = 3, seed = 1)),
    lpa_one = quiet(lpa(data, continuous, n_profiles = 1, n_starts = 1, seed = 1)),
    lca = quiet(lca(data, categorical, n_classes = 3, n_starts = 3, seed = 1)),
    weighted = quiet(lpa(data, continuous, n_profiles = 3, weights = "w", n_starts = 3,
                         seed = 1)),
    noise = quiet(lpa(data, continuous, n_profiles = 3, noise = TRUE, n_starts = 3, seed = 1)),
    multilpa = quiet(multilpa(data, continuous, id = "id", n_profiles = 3,
                              n_group_classes = 2, n_starts = 3, seed = 1)),
    multilca = quiet(multilca(data, categorical, id = "id", n_profiles = 3,
                              n_group_classes = 2, n_starts = 3, seed = 1)),
    covariates = quiet(multilpa(data, continuous, id = "id", n_profiles = 3,
                                n_group_classes = 2, profile_covariates = "x", n_starts = 3,
                                seed = 1)),
    transitions = quiet(lta(data, continuous, id = "id", n_profiles = 3, time = "t",
                            n_starts = 3, seed = 1)),
    additive = quiet(multilpa(data, continuous, id = "id", family = "additive",
                              n_group_classes = 2, n_starts = 3, seed = 1)),
    cross_level = quiet(multilpa(data, continuous, id = "id", n_profiles = 3,
                                 n_group_classes = 2, family = "restricted_cross_level",
                                 n_starts = 2, seed = 1)),
    regression = quiet(mixture_regression(y1 ~ x, data, n_classes = 2, n_starts = 3,
                                          seed = 1)),
    growth = quiet(mixture_regression(y1 ~ t, data, n_classes = 2, id = "id",
                                      class_level = "group", random = "t", n_starts = 2,
                                      seed = 1)))
}

test_that("every family gets a report with the contract's columns and statuses", {
  skip_on_cran()
  fits <- diagnostic_fits()
  invisible(lapply(names(fits), function(name) {
    checked <- diagnostics(fits[[name]], plots = FALSE)
    checks <- get_results(checked, "checks")
    expect_named(checks, c("domain", "check", "level", "value", "display", "criterion",
                           "status", "detail", "reference"))
    expect_true(all(checks$status %in% c("ok", "check", "fail", "info", "n/a")), label = name)
    expect_true(all(c("Estimation", "Information") %in% checks$domain), label = name)
    # Nothing judged prints a missing value: a check that cannot be computed
    # is n/a, with the reason.
    judged <- checks[checks$status %in% c("ok", "check", "fail"), , drop = FALSE]
    expect_false(any(grepl("\\bNA\\b", judged$display)), label = name)
    expect_true(all(nzchar(checks$detail[checks$status == "n/a"])), label = name)
    shown <- utils::capture.output(print(checked))
    expect_false(any(grepl("\\bNA\\b", shown)), label = name)
    expect_match(shown[[1L]], "^Diagnostics: ", label = name)
  }))
})

test_that("every value is the value of the table it is read from", {
  skip_on_cran()
  fit <- diagnostic_fits()$multilca
  checked <- diagnostics(fit, plots = FALSE)
  checks <- get_results(checked, "checks")
  value <- function(check, level) checks$value[checks$check == check & checks$level == level]
  classification <- get_results(fit, "classification")
  invisible(lapply(c("individuals", "groups"), function(level) {
    rows <- classification[classification$level == level, , drop = FALSE]
    assigned <- rows[rows$n_modal > 0L, , drop = FALSE]
    expect_identical(value("lowest average posterior", level), min(assigned$average_posterior))
    expect_identical(value("lowest odds of correct classification", level),
                     min(assigned$odds_correct_classification))
    expect_identical(value("smallest class", level), min(rows$estimated_proportion))
    expect_identical(value("classes with no assigned unit", "individuals") +
                       0 * nrow(rows), as.numeric(sum(classification$n_modal[
                         classification$level == "individuals"] == 0L)))
    entropy <- get_results(fit, "entropy")
    expect_identical(value("relative entropy", level),
                     entropy$relative_entropy[entropy$level == level])
  }))
  expect_identical(value("best log likelihood reached by", "model"),
                   as.numeric(fit$n_best_replicated))
  residuals <- get_results(checked, "residuals")
  expect_identical(value("pairs dependent within a class", "model"),
                   as.numeric(sum(residuals$p_adjusted < 0.05)))
  expect_identical(value("log likelihood, parameters, BIC", "model"), fit$log_likelihood)
})

test_that("the posterior formulas match a hand computation", {
  posteriors <- rbind(c(0.9, 0.1), c(0.8, 0.2), c(0.3, 0.7), c(0.4, 0.6))
  table <- .latents_posterior_classification(posteriors, "individuals")
  expect_identical(table$n_assigned, c(2L, 2L))
  expect_equal(table$estimated_proportion, c(0.6, 0.4))
  expect_equal(table$average_posterior, c(0.85, 0.65))
  # Nagin's odds of correct classification: (AvePP / (1 - AvePP)) / (pi / (1 - pi)).
  expect_equal(table$odds_correct_classification,
               c((0.85 / 0.15) / (0.6 / 0.4), (0.65 / 0.35) / (0.4 / 0.6)))
  # A class no unit is assigned to has no average posterior.
  empty <- .latents_posterior_classification(rbind(c(0.9, 0.1), c(0.6, 0.4)), "groups")
  expect_identical(empty$n_assigned, c(2L, 0L))
  expect_true(is.na(empty$average_posterior[[2L]]))
})

test_that("rule-of-thumb quantities are reported without a verdict", {
  skip_on_cran()
  # AvePP, OCC, the smallest class, entropy and the replication count have
  # published cut-offs that are rules of thumb; the package reports them and
  # passes no judgement, whatever their value.
  checks <- get_results(diagnostics(diagnostic_fits()$multilca, plots = FALSE), "checks")
  reported <- c("best log likelihood reached by", "smallest class",
                "modal vs estimated proportions", "relative entropy",
                "lowest average posterior", "lowest odds of correct classification")
  rows <- checks[checks$check %in% reported, , drop = FALSE]
  expect_gt(nrow(rows), 0L)
  expect_true(all(rows$status == "info"))
  expect_true(all(rows$criterion == ""))
  expect_false("criteria" %in% names(formals(diagnostics)))
})

test_that("screened starts are reported apart, never as unfinished", {
  skip_on_cran()
  fit <- diagnostic_fits()$regression
  starts <- get_results(fit, "starts")
  skip_if_not("stage" %in% names(starts))
  table <- get_results(diagnostics(fit, plots = FALSE), "replication")
  expect_identical(table$n_screened, sum(starts$stage != "completed"))
  expect_identical(table$n_starts, sum(starts$stage == "completed"))
  checks <- get_results(diagnostics(fit, plots = FALSE), "checks")
  unfinished <- checks[checks$check == "starts failed or unfinished", , drop = FALSE]
  expect_identical(unfinished$value,
                   as.numeric(sum(starts$stage == "completed" & !starts$converged)))
})

test_that("a one-class fit, an empty class and a noise fit report n/a or check, not NA", {
  skip_on_cran()
  fits <- diagnostic_fits()
  one <- get_results(diagnostics(fits$lpa_one, plots = FALSE), "checks")
  expect_true(all(one$status[one$domain == "Classification"] == "n/a"))
  noise <- get_results(diagnostics(fits$noise, plots = FALSE), "checks")
  expect_identical(noise$status[noise$domain == "Local independence"], "n/a")
})

test_that("a class with no assigned unit is flagged and its AvePP left out of the minimum", {
  fit <- list()
  class(fit) <- "multilpa_cross_level"
  # Two groups, both most likely in class 1: class 2 has no member.
  fit$group_posteriors <- rbind(c(0.9, 0.1), c(0.6, 0.4))
  checks <- .latents_level_checks(fit, "groups", 2L, 0.4)
  empty <- checks[checks$check == "classes with no assigned group", , drop = FALSE]
  expect_identical(empty$status, "check")
  expect_identical(empty$value, 1)
  expect_match(empty$detail, "^class 2: no group")
  lowest <- checks[checks$check == "lowest average posterior", , drop = FALSE]
  expect_equal(lowest$value, 0.75)
  expect_identical(lowest$status, "info")
})

test_that("a diagnostics object saved before the checks still prints and answers", {
  skip_on_cran()
  checked <- diagnostics(diagnostic_fits()$lpa, plots = FALSE)
  checked$checks <- NULL
  checked$replication <- NULL
  expect_no_error(utils::capture.output(print(checked)))
  expect_identical(nrow(get_results(checked, "replication")), 1L)
  expect_false("checks" %in% names(get_results(checked, "all")))
})

test_that("residuals without a test are described, never reported as passing", {
  skip_on_cran()
  # Under sampling weights the bivariate residuals carry no p value; a count
  # of significant pairs would count tests that were never run.
  fit <- diagnostic_fits()$weighted
  residuals <- get_results(diagnostics(fit, plots = FALSE), "residuals")
  skip_if_not(all(is.na(residuals$p_adjusted)))
  checks <- get_results(diagnostics(fit, plots = FALSE), "checks")
  local <- checks[checks$domain == "Local independence", , drop = FALSE]
  expect_identical(local$status, "info")
  expect_identical(local$criterion, "not tested")
  expect_match(local$detail, "sampling weights")
  expect_equal(local$value, max(abs(residuals$residual)))
})

# Adversarial review of the redesign (Codex, 2026-10-27): each test failed on
# the code it was written against.

test_that("finished starts are every stage but screened, so emEM counts are right", {
  skip_on_cran()
  data <- diagnostic_data()
  fit <- suppressWarnings(lpa(data, c("y1", "y2", "y3"), n_profiles = 3, model = "EEE",
                              search = "emEM", n_starts = 3, seed = 1))
  starts <- get_results(fit, "starts")
  table <- get_results(diagnostics(fit, plots = FALSE), "replication")
  if ("stage" %in% names(starts)) {
    expect_identical(table$n_starts, sum(starts$stage != "screened"))
    expect_identical(table$n_screened, sum(starts$stage == "screened"))
  }
  expect_gt(table$n_starts, 0L)
  expect_lte(table$n_best_replicated, table$n_starts)
  expect_lte(table$n_converged + table$n_failed, table$n_starts)
})

test_that("a cross-level fit is checked at both levels, weighted where it can be", {
  skip_on_cran()
  data <- diagnostic_data()
  data$w <- rep(stats::runif(60L, 0.5, 2), each = 10L)
  fit <- suppressWarnings(multilpa(data, c("y1", "y2", "y3"), id = "id", n_profiles = 2,
                                   n_group_classes = 2, family = "restricted_cross_level",
                                   weights = "w", n_starts = 2, seed = 1))
  checks <- get_results(diagnostics(fit, plots = FALSE), "checks")
  expect_true(all(c("individuals", "groups") %in% checks$level))
  # The group level is weighted by the stored group weights.
  table <- .latents_level_classification(fit, "groups")
  weights <- fit$sampling_weights * nrow(fit$group_posteriors) / sum(fit$sampling_weights)
  expect_equal(table$estimated_proportion,
               unname(colSums(fit$group_posteriors * weights) / sum(weights)))
  # The sufficient statistics retain the row-to-group map for this level.
  individual <- .latents_level_classification(fit, "individuals")
  row_weights <- fit$sampling_weights[fit$sufficient_statistics$index]
  expect_equal(individual$estimated_n,
               unname(colSums(fit$subject_posteriors * row_weights)))
  expect_false(any(checks$status[checks$level == "individuals"] == "n/a"))
})

test_that("untested residual pairs never make the check pass", {
  residuals <- data.frame(profile = "profile_1", indicator_1 = c("a", "a", "b", "c"),
                          indicator_2 = c("b", "c", "c", "d"), residual = c(0.1, 0.2, 0.05, 0.3),
                          p_adjusted = c(0.4, 0.6, NA, NA))
  check <- .latents_residual_checks(residuals, "")
  expect_identical(check$status, "info")
  expect_match(check$display, "0 of 2 tested pairs \\(2 not tested\\)")
  empty <- residuals
  empty$residual <- NA_real_
  empty$p_adjusted <- NA_real_
  expect_identical(.latents_residual_checks(empty, "")$status,
                   "n/a")
})

test_that("a family without classification plots is diagnosed with plots = TRUE", {
  skip_on_cran()
  fit <- diagnostic_fits()$additive
  expect_message(checked <- diagnostics(fit), "no classification plots")
  expect_s3_class(checked, "multilpa_diagnostics")
})

test_that("an indicator named profile is printed with the others", {
  skip_on_cran()
  data <- diagnostic_data()
  data$profile <- data$y2
  fit <- suppressWarnings(lpa(data, c("y1", "profile"), n_profiles = 2, n_starts = 2, seed = 1))
  printed <- utils::capture.output(print(fit))
  expect_true(any(grepl("^ +profile +-?[0-9]", printed)))
  expect_true(any(grepl("^ +y1 +-?[0-9]", printed)))
})

test_that("diagnostics keeps 0.9.21's default table, the entropy", {
  skip_on_cran()
  checked <- diagnostics(diagnostic_fits()$lpa, plots = FALSE)
  expect_identical(as.data.frame(checked), get_results(checked, "entropy"))
  expect_identical(get_results(checked), get_results(checked, "entropy"))
  expect_identical(names(get_results(checked, "all"))[1:4],
                   c("entropy", "classification", "average_posteriors", "residuals"))
})

test_that("a failed seed keeps its place on the sensitivity plot and is named", {
  skip_if_not_installed("ggplot2")
  table <- data.frame(seed = 1:4, log_likelihood = c(-10, -10, NA, -12),
                      converged = c(TRUE, TRUE, NA, TRUE), iterations = c(5L, 6L, NA, 7L),
                      optimum = c(1L, 1L, NA, 2L), best = c(TRUE, TRUE, NA, FALSE),
                      agreement = c(1, 1, NA, 0.8))
  attr(table, "tolerance") <- 1e-3
  class(table) <- c("latents_sensitivity", "data.frame")
  plot <- plot(table)
  built <- ggplot2::ggplot_build(plot)
  expect_identical(built$layout$panel_params[[1L]]$x$get_limits(), as.character(1:4))
  expect_match(plot$labels$caption, "failed .* seed 3")
})
