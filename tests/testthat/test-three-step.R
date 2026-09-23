.step_data <- function(seed = 11L, n_groups = 50L, per = 10L, separation = 1.1,
                       difference = 10) {
  set.seed(seed)
  group <- rep(seq_len(n_groups), each = per)
  n <- length(group)
  truth <- 1L + as.integer(stats::runif(n) > 0.5)
  data.frame(g = group, truth = truth,
             a = stats::rnorm(n, ifelse(truth == 2L, separation, -separation)),
             b = stats::rnorm(n, ifelse(truth == 2L, separation, -separation)),
             y = stats::rnorm(n, ifelse(truth == 2L, difference, 0)))
}

.step_fit <- function(data, ...) {
  multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
           n_starts = 6, seed = 1, ...)
}

test_that("the error matrix conditions on the true class and is a distribution", {
  data <- .step_data()
  fit <- .step_fit(data)
  errors <- get_results(fit, "classification_errors", level = "individuals")

  expect_named(errors, c("level", "true_class", "assigned_class", "probability"))
  expect_equal(nrow(errors), 4L)
  # rows sum to one within each true class
  totals <- tapply(errors$probability, errors$true_class, sum)
  expect_equal(as.vector(totals), c(1, 1))
  expect_true(all(errors$probability >= 0 & errors$probability <= 1))
  # separated classes are assigned correctly most of the time
  correct <- subset(errors, true_class == assigned_class)
  expect_true(all(correct$probability > 0.85))
})

test_that("BCH weights sum to one within a unit and are signed", {
  data <- .step_data()
  fit <- .step_fit(data)
  weights <- get_results(fit, "bch_weights", level = "individuals")

  # one row per unit and class, not one column per class
  expect_named(weights, c("level", "unit", "assigned_class", "class", "weight"))
  expect_equal(nrow(weights), fit$n_observations * fit$n_profiles)
  expect_setequal(weights$class, seq_len(fit$n_profiles))
  within_unit <- tapply(weights$weight, weights$unit, sum)
  expect_equal(as.vector(within_unit), rep(1, fit$n_observations))
  # unlike posteriors, they are not probabilities
  expect_true(any(weights$weight < 0))
  expect_true(any(weights$weight > 1))
  # a unit's own class carries its largest weight
  heaviest <- tapply(seq_len(nrow(weights)), weights$unit,
                     function(rows) weights$class[rows][which.max(weights$weight[rows])])
  expect_equal(as.vector(heaviest),
               as.vector(tapply(weights$assigned_class, weights$unit, unique)))
})

test_that("BCH removes the attenuation that modal assignment creates", {
  data <- .step_data()
  fit <- .step_fit(data)
  # `contrast = "pairs"` reports class 2 minus class 1; orient it so the gap is
  # the high-outcome class minus the other, whichever way the labels fell.
  orientation <- if (which.max(fit$means[, 1L]) == 2L) 1 else -1
  gap <- function(method) {
    orientation * three_step(fit, data, "y", method = method,
                             contrast = "pairs")$estimate
  }

  # The generating difference is 10. Modal assignment shrinks it; so does
  # posterior weighting, more so; BCH does not.
  expect_lt(gap("modal"), 9.5)
  expect_lt(gap("proportional"), gap("modal"))
  expect_gt(gap("bch"), gap("modal"))
  expect_lt(abs(gap("bch") - 10), abs(gap("modal") - 10))
  expect_lt(abs(gap("bch") - 10), 1)
})

test_that("the result is tidy, with cluster-robust intervals", {
  data <- .step_data()
  fit <- .step_fit(data)
  result <- three_step(fit, data, "y")

  expect_named(result, c("level", "method", "class", "estimate",
                         "standard_error", "conf_low", "conf_high",
                         "effective_n"))
  expect_equal(nrow(result), fit$n_profiles)
  expect_true(all(result$standard_error > 0))
  expect_equal(result$conf_high - result$conf_low,
               2 * stats::qnorm(0.975) * result$standard_error)
  # the call is described by columns, not by attributes the caller must reach for
  expect_identical(unique(result$method), "bch")
  expect_identical(unique(result$level), "individuals")
  expect_identical(unique(three_step(fit, data, "y", method = "modal")$method),
                   "modal")
  # a wider level gives a wider interval
  wide <- three_step(fit, data, "y", ci_level = 0.99)
  expect_true(all(wide$conf_high - wide$conf_low >
                    result$conf_high - result$conf_low))
})

test_that("a group-level outcome is related to the group classes", {
  set.seed(3)
  n_groups <- 60L
  group <- rep(seq_len(n_groups), each = 10L)
  group_class <- rep(rep(c(1L, 2L), length.out = n_groups), each = 10L)
  n <- length(group)
  profile <- 1L + as.integer(stats::runif(n) <
                               ifelse(group_class == 2L, 0.85, 0.15))
  outcome <- rep(stats::rnorm(n_groups, rep(c(0, 5), length.out = n_groups)),
                 each = 10L)
  data <- data.frame(g = group, z = outcome,
                     a = stats::rnorm(n, ifelse(profile == 2L, 1.6, -1.6)),
                     b = stats::rnorm(n, ifelse(profile == 2L, 1.6, -1.6)))
  fit <- multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 2,
                  n_starts = 6, seed = 2)
  result <- three_step(fit, data, "z", level = "groups")

  expect_equal(nrow(result), fit$n_group_classes)
  expect_identical(unique(result$level), "groups")
  difference <- three_step(fit, data, "z", level = "groups", contrast = "pairs")
  expect_gt(abs(difference$estimate), 3)
  expect_lt(difference$p_value, 0.05)
  expect_equal(nrow(get_results(fit, "classification_errors", level = "groups")), 4L)
  expect_equal(nrow(get_results(fit, "bch_weights", level = "groups")),
               fit$n_groups * fit$n_group_classes)
})

test_that("a broken contract is refused", {
  data <- .step_data()
  fit <- .step_fit(data)

  expect_error(three_step(fit, data, "absent"), "must name a single column")
  # `y` varies within a group, so it cannot be a group-level outcome.
  expect_error(three_step(fit, data, "y", level = "groups"),
               class = "multilpa_bad_outcome")
  missing_outcome <- data
  missing_outcome$y[1L] <- NA
  expect_error(three_step(fit, missing_outcome, "y"), "must not be missing")
  expect_error(three_step(fit, data[1:10, ], "y"),
               class = "multilpa_bad_inference_data")
  expect_error(three_step(fit, data, "y", ci_level = 1), "`ci_level` must be")
})

test_that("inseparable classes are refused rather than inverted", {
  data <- .step_data()
  fit <- .step_fit(data)
  # A degenerate error matrix means the assignment carries no information.
  broken <- fit
  broken$subject_profiles <- rep(1L, fit$n_observations)
  expect_error(get_results(broken, "bch_weights", level = "individuals"), class = "multilpa_inseparable_classes")
  expect_error(three_step(broken, data, "y"),
               class = "multilpa_inseparable_classes")
})

test_that("classification errors remain available for covariate fits", {
  data <- .step_data()
  data$x <- stats::rnorm(nrow(data))
  fit <- multilpa(data, c("a", "b"), "g", n_profiles = 2,
                        n_group_classes = 2, profile_covariates = "x",
                        n_starts = 4, seed = 1)
  expect_equal(nrow(get_results(fit, "classification_errors", level = "individuals")), 4L)
  expect_error(three_step(fit, data, "y"),
               class = "multilpa_unsupported_three_step")
})

test_that("a one-class level is degenerate but does not error", {
  data <- .step_data()
  fit <- .step_fit(data)
  # The fit has a single group class, so its error matrix is 1 x 1. vapply()
  # drops such a result to a vector, which used to break the accessor.
  errors <- get_results(fit, "classification_errors", level = "groups")

  expect_equal(nrow(errors), 1L)
  expect_equal(errors$probability, 1)
  expect_equal(nrow(get_results(fit, "bch_weights", level = "groups")), fit$n_groups)
})

test_that("the weights satisfy the identities that define the BCH method", {
  data <- .step_data()
  fit <- .step_fit(data)
  pieces <- .multilpa_level_assignments(fit, "individuals")
  errors <- .multilpa_error_matrix(pieces)
  inverse <- solve(errors)
  weights <- get_results(fit, "bch_weights", level = "individuals")

  # Bolck, Croon & Hagenaars (2004): the weights are the inverse of the
  # misclassification matrix applied to assigned class membership.
  expect_equal(unname(rowSums(errors)), rep(1, nrow(errors)))
  expect_equal(unname(errors %*% inverse), diag(nrow(errors)))
  expect_equal(as.vector(tapply(weights$weight, weights$unit, sum)),
               rep(1, fit$n_observations))

  # The property that makes the third step unbiased: the weighted class sizes
  # reproduce the model's own estimated class sizes.
  expect_equal(as.vector(tapply(weights$weight, weights$class, sum)),
               unname(colSums(pieces$posteriors)))
})

test_that("the per-unit weights agree with the original table formulation", {
  data <- .step_data()
  fit <- .step_fit(data)
  pieces <- .multilpa_level_assignments(fit, "individuals")
  inverse <- solve(.multilpa_error_matrix(pieces))

  # Bolck et al. correct the assigned-by-external contingency table directly;
  # this package carries per-unit weights. They are the same method.
  bins <- cut(data$y, breaks = stats::quantile(data$y, c(0, 0.5, 1)),
              include.lowest = TRUE)
  weights <- merge(get_results(fit, "bch_weights", level = "individuals"),
                   data.frame(unit = seq_along(bins), bin = bins),
                   by = "unit")
  from_table <- t(inverse) %*% table(pieces$modal, bins)
  from_weights <- tapply(weights$weight, list(weights$class, weights$bin), sum)
  expect_equal(unname(from_table), unname(from_weights))
})

.r3step_data <- function(seed = 21L, slope = 1.2, intercept = -0.3) {
  set.seed(seed)
  g <- rep(seq_len(50), each = 12L)
  n <- length(g)
  x <- stats::rnorm(n)
  truth <- 1L + as.integer(stats::runif(n) < stats::plogis(intercept + slope * x))
  data.frame(g = g, x = x, truth = truth,
             a = stats::rnorm(n, ifelse(truth == 2L, 1.2, -1.2)),
             b = stats::rnorm(n, ifelse(truth == 2L, 1.2, -1.2)))
}

test_that("R3STEP recovers a covariate effect the naive regression attenuates", {
  data <- .r3step_data()
  fit <- multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                  n_starts = 5, seed = 1)
  corrected <- r3step(fit, data, "x")
  # r3step contrasts class 1 against the final class; orient to the generator.
  sign <- if (which.max(fit$means[, 1L]) != 1L) -1 else 1
  slope <- sign * corrected$estimate[corrected$term == "x"]
  assigned <- fit$subject_profiles
  naive <- sign * stats::coef(
    stats::glm(I(assigned == 1L) ~ x, data = data, family = stats::binomial()))[2L]

  expect_lt(naive, 1.2)
  expect_gt(slope, naive)
  expect_lt(abs(slope - 1.2), abs(naive - 1.2))
  expect_lt(abs(slope - 1.2), 0.25)
})

test_that("the R3STEP table is tidy and its intervals are consistent", {
  data <- .r3step_data()
  fit <- multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                  n_starts = 5, seed = 1)
  result <- r3step(fit, data, "x")

  expect_named(result, c("level", "outcome", "term", "estimate",
                         "standard_error", "statistic", "p_value",
                         "p_value_adjusted", "conf_low", "conf_high"))
  # one non-reference class, two terms
  expect_equal(nrow(result), 2L)
  expect_setequal(result$term, c("(Intercept)", "x"))
  expect_identical(attr(result, "reference_class"), 2L)
  expect_true(all(result$standard_error > 0))
  expect_equal(result$conf_high - result$conf_low,
               2 * stats::qnorm(0.975) * result$standard_error)
  expect_equal(result$statistic, result$estimate / result$standard_error)
})

test_that("the robust sandwich runs and more covariates are supported", {
  data <- .r3step_data()
  data$w <- stats::rnorm(nrow(data))
  fit <- multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                  n_starts = 5, seed = 1)
  observed <- r3step(fit, data, c("x", "w"))
  robust <- r3step(fit, data, c("x", "w"), vcov_type = "robust")

  expect_equal(nrow(observed), 3L)
  expect_setequal(observed$term, c("(Intercept)", "x", "w"))
  expect_equal(robust$estimate, observed$estimate)
  expect_true(all(robust$standard_error > 0))
  expect_identical(attr(robust, "vcov_type"), "robust")
  # the irrelevant covariate is not significant
  expect_gt(observed$p_value[observed$term == "w"], 0.05)
})

test_that("R3STEP works at the group level and refuses a varying covariate", {
  set.seed(5)
  n_groups <- 60L
  g <- rep(seq_len(n_groups), each = 10L)
  n <- length(g)
  w <- rep(stats::rnorm(n_groups), each = 10L)
  group_class <- 1L + as.integer(
    stats::runif(n_groups) < stats::plogis(1.5 * w[!duplicated(g)]))
  expanded <- rep(group_class, each = 10L)
  profile <- 1L + as.integer(stats::runif(n) <
                               ifelse(expanded == 2L, 0.85, 0.15))
  data <- data.frame(g = g, w = w, varying = stats::rnorm(n),
                     a = stats::rnorm(n, ifelse(profile == 2L, 1.6, -1.6)),
                     b = stats::rnorm(n, ifelse(profile == 2L, 1.6, -1.6)))
  fit <- multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 2,
                  n_starts = 6, seed = 2)

  result <- r3step(fit, data, "w", level = "groups")
  expect_equal(nrow(result), 2L)
  expect_identical(result$level[1L], "groups")
  expect_lt(result$p_value[result$term == "w"], 0.05)
  expect_error(r3step(fit, data, "varying", level = "groups"),
               class = "multilpa_bad_covariate")
})

test_that("a broken contract is refused", {
  data <- .r3step_data()
  fit <- multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                  n_starts = 4, seed = 1)

  expect_error(r3step(fit, data, "absent"), "must name columns")
  expect_error(r3step(fit, data, "truth", ci_level = 1), "`ci_level` must be")
  missing_covariate <- data
  missing_covariate$x[1L] <- NA
  expect_error(r3step(fit, missing_covariate, "x"), "must not be missing")
  data$copy <- data$x
  expect_error(r3step(fit, data, c("x", "copy")),
               class = "multilpa_bad_covariate")
  expect_error(r3step(fit, data, "x", level = "groups"),
               class = "multilpa_inseparable_classes")
})
