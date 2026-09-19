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
  errors <- classification_errors(fit)

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
  weights <- bch_weights(fit)
  columns <- grep("^weight_class_", names(weights), value = TRUE)

  expect_equal(nrow(weights), fit$n_observations)
  expect_equal(length(columns), fit$n_profiles)
  expect_equal(unname(rowSums(weights[columns])), rep(1, nrow(weights)))
  # unlike posteriors, they are not probabilities
  expect_true(any(unlist(weights[columns]) < 0))
  expect_true(any(unlist(weights[columns]) > 1))
  # a unit's own class carries its largest weight
  own <- vapply(seq_len(nrow(weights)), function(i)
    which.max(as.numeric(weights[i, columns])), integer(1))
  expect_equal(own, weights$assigned_class)
})

test_that("BCH removes the attenuation that modal assignment creates", {
  data <- .step_data()
  fit <- .step_fit(data)
  high <- which.max(fit$means[, 1L])
  gap <- function(method) {
    estimates <- three_step(fit, data, "y", method = method)$estimate
    estimates[high] - estimates[-high]
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

  expect_named(result, c("class", "estimate", "standard_error", "conf_low",
                         "conf_high", "effective_n"))
  expect_equal(nrow(result), fit$n_profiles)
  expect_true(all(result$standard_error > 0))
  expect_equal(result$conf_high - result$conf_low,
               2 * stats::qnorm(0.975) * result$standard_error)
  expect_identical(attr(result, "method"), "bch")
  expect_identical(attr(result, "level"), "individuals")
  # a wider level gives a wider interval
  wide <- three_step(fit, data, "y", level_ci = 0.99)
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
  expect_identical(attr(result, "level"), "groups")
  expect_gt(abs(diff(result$estimate)), 3)
  expect_equal(nrow(classification_errors(fit, level = "groups")), 4L)
  expect_equal(nrow(bch_weights(fit, level = "groups")), fit$n_groups)
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
               "one row per observation")
  expect_error(three_step(fit, data, "y", level_ci = 1), "`level_ci` must be")
})

test_that("inseparable classes are refused rather than inverted", {
  data <- .step_data()
  fit <- .step_fit(data)
  # A degenerate error matrix means the assignment carries no information.
  broken <- fit
  broken$subject_profiles <- rep(1L, fit$n_observations)
  expect_error(bch_weights(broken), class = "multilpa_inseparable_classes")
  expect_error(three_step(broken, data, "y"),
               class = "multilpa_inseparable_classes")
})

test_that("a covariate fit is supported", {
  data <- .step_data()
  data$x <- stats::rnorm(nrow(data))
  fit <- fit_covariates(data, c("a", "b"), "g", n_profiles = 2,
                        n_group_classes = 2, profile_covariates = "x",
                        n_starts = 4, seed = 1)
  result <- three_step(fit, data, "y")

  expect_equal(nrow(result), 2L)
  expect_true(all(result$standard_error > 0))
  expect_equal(nrow(classification_errors(fit)), 4L)
})

test_that("a one-class level is degenerate but does not error", {
  data <- .step_data()
  fit <- .step_fit(data)
  # The fit has a single group class, so its error matrix is 1 x 1. vapply()
  # drops such a result to a vector, which used to break the accessor.
  errors <- classification_errors(fit, level = "groups")

  expect_equal(nrow(errors), 1L)
  expect_equal(errors$probability, 1)
  expect_equal(nrow(bch_weights(fit, level = "groups")), fit$n_groups)
})
