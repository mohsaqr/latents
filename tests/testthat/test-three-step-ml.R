# The ML step-three model of a distal outcome (Vermunt 2010; Bakk, Tekle and
# Vermunt 2013): the assigned class is an indicator with the step-one
# classification errors held fixed, the outcome normal within each class.
# The comparisons with independent references are in tests/equivalence/test-three-step-ml.R.

# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

ml_data <- function(seed = 11, n_groups = 80, per = 10) {
  set.seed(seed)
  n <- n_groups * per
  truth <- 1L + as.integer(stats::runif(n) > 0.5)
  data.frame(g = rep(seq_len(n_groups), each = per),
             y1 = stats::rnorm(n, c(-1, 1)[truth]),
             y2 = stats::rnorm(n, c(-1, 1)[truth]),
             distal = stats::rnorm(n, c(0, 3)[truth], c(1, 2)[truth]),
             truth = truth)
}

ml_fit <- function(data) {
  quietly(multilpa(data, c("y1", "y2"), "g", n_profiles = 2, n_group_classes = 1,
                   n_starts = 3, seed = 1))
}

test_that("with no classification error the ML means are the class means", {
  # Hard posteriors make the error matrix the identity, so each unit belongs
  # to its assigned class and the step-three model is a normal model per class.
  set.seed(4)
  values <- c(stats::rnorm(30, 2), stats::rnorm(20, -1, 3))
  modal <- rep(1:2, c(30, 20))
  pieces <- list(posteriors = cbind(modal == 1, modal == 2) * 1, modal = modal,
                 n_classes = 2L)
  fit <- .multilpa_step_ml_fit(pieces, values, "varying")
  expect_true(fit$converged)
  expect_equal(fit$means, as.vector(tapply(values, modal, mean)), tolerance = 1e-10)
  expected_variances <- as.vector(tapply(values, modal, \(v) mean((v - mean(v))^2)))
  expect_equal(fit$variances, expected_variances, tolerance = 1e-10)
  pooled <- .multilpa_step_ml_fit(pieces, values, "equal")
  expect_equal(pooled$variances[1L], sum(30 * expected_variances[1L],
                                         20 * expected_variances[2L]) / 50,
               tolerance = 1e-10)
})

test_that("the analytic scores are the gradient of the step-three likelihood", {
  data <- ml_data()
  fit <- ml_fit(data)
  pieces <- .multilpa_level_assignments(fit, "individuals")
  invisible(lapply(c("varying", "equal"), \(variance) {
    step <- .multilpa_step_ml_fit(pieces, data$distal, variance)
    point <- .multilpa_step_ml_coordinates(step, variance) + c(0.1, -0.2, 0.05,
      if (identical(variance, "varying")) -0.1)
    analytic <- colSums(.multilpa_step_ml_scores(step, data$distal, variance, point))
    numeric <- vapply(seq_along(point), \(j) {
      h <- 1e-6
      up <- point; up[j] <- up[j] + h
      down <- point; down[j] <- down[j] - h
      (.multilpa_step_ml_log_likelihood(step, data$distal, variance, up) -
         .multilpa_step_ml_log_likelihood(step, data$distal, variance, down)) / (2 * h)
    }, numeric(1))
    expect_equal(analytic, numeric, tolerance = 1e-6)
    # At the estimate the score is zero.
    expect_lt(max(abs(colSums(.multilpa_step_ml_scores(step, data$distal, variance)))),
              1e-4)
  }))
})

test_that("ML recovers the class means and variances of a simulated outcome", {
  data <- ml_data()
  fit <- ml_fit(data)
  result <- three_step(fit, data, "distal", method = "ml")
  expect_s3_class(result, "data.frame")
  expect_named(result, c("level", "method", "class", "estimate", "standard_error",
                         "conf_low", "conf_high", "variance", "effective_n"))
  # The profiles are numbered by the fit; match them to the simulated classes.
  order <- order(result$estimate)
  expect_equal(result$estimate[order], c(0, 3), tolerance = 0.15)
  expect_equal(result$variance[order], c(1, 4), tolerance = 0.25)
  expect_true(all(result$conf_low < result$estimate & result$estimate < result$conf_high))
  expect_equal(sum(result$effective_n), nrow(data), tolerance = 1e-8)
  expect_identical(attr(result, "variance"), "varying")
})

test_that("shifting the outcome shifts the means; scaling scales means and errors", {
  data <- ml_data()
  fit <- ml_fit(data)
  base <- three_step(fit, data, "distal", method = "ml")
  shifted <- three_step(fit, transform(data, distal = distal + 100), "distal",
                        method = "ml")
  expect_equal(shifted$estimate, base$estimate + 100, tolerance = 1e-6)
  expect_equal(shifted$standard_error, base$standard_error, tolerance = 1e-6)
  scaled <- three_step(fit, transform(data, distal = distal * 10), "distal",
                       method = "ml")
  expect_equal(scaled$estimate, base$estimate * 10, tolerance = 1e-6)
  expect_equal(scaled$standard_error, base$standard_error * 10, tolerance = 1e-6)
  expect_equal(scaled$variance, base$variance * 100, tolerance = 1e-6)
})

test_that("ML standard errors do not depend on the outcome's units", {
  data <- ml_data()
  fit <- ml_fit(data)
  base <- three_step(fit, data, "distal", method = "ml")
  invisible(lapply(c(1e-4, 1e4, 1e8), \(scale) {
    scaled <- three_step(fit, transform(data, distal = distal * scale), "distal",
                         method = "ml")
    expect_equal(scaled$estimate / scale, base$estimate, tolerance = 1e-6)
    expect_equal(scaled$standard_error / scale, base$standard_error, tolerance = 1e-4)
  }))
})

test_that("ML pairs are the differences of its class means", {
  data <- ml_data()
  fit <- ml_fit(data)
  means <- three_step(fit, data, "distal", method = "ml")
  pairs <- three_step(fit, data, "distal", method = "ml", contrast = "pairs")
  expect_identical(nrow(pairs), 1L)
  expect_equal(pairs$estimate,
               means$estimate[means$class == pairs$class] -
                 means$estimate[means$class == pairs$reference_class],
               tolerance = 1e-10)
  expect_identical(unique(pairs$method), "ml")
})

test_that("ML arguments and outcomes it cannot take are refused by class", {
  data <- ml_data()
  fit <- ml_fit(data)
  expect_error(three_step(fit, data, "distal", method = "bch", variance = "equal"),
               class = "latents_bad_argument")
  expect_error(three_step(fit, transform(data, distal = 1), "distal", method = "ml"),
               class = "latents_bad_outcome")
  # The released methods are unchanged by the new argument and method,
  # abbreviations included: "m" is still "modal".
  expect_identical(three_step(fit, data, "distal"),
                   three_step(fit, data, "distal", method = "bch"))
  expect_identical(three_step(fit, data, "distal", method = "m"),
                   three_step(fit, data, "distal", method = "modal"))
  expect_identical(three_step(fit, data, "distal", method = "p"),
                   three_step(fit, data, "distal", method = "proportional"))
  expect_error(three_step(fit, data, "distal", method = "x"), "ml")
})

test_that("ML uncertainty at a large origin matches known normal means", {
  set.seed(4)
  values <- c(stats::rnorm(30, 2), stats::rnorm(20, -1, 3)) + 1e12
  modal <- rep(1:2, c(30, 20))
  pieces <- list(posteriors = cbind(modal == 1, modal == 2) * 1,
                 modal = modal, n_classes = 2L)
  centered <- values - mean(values)
  reference <- vapply(split(centered, modal), function(value) {
    sqrt(sum((value - mean(value))^2)) / length(value)
  }, numeric(1))
  invisible(lapply(c("varying", "equal"), function(variance) {
    fit <- .multilpa_step_ml_fit(pieces, values, variance)
    influence <- .multilpa_step_ml_influence(fit, values, variance)
    expect_equal(unname(sqrt(colSums(influence^2))), unname(reference),
                 tolerance = 1e-7)
  }))
})

test_that("public ML means and contrasts retain uncertainty at large outcome origins", {
  data <- ml_data(n_groups = 30, per = 8)
  fit <- ml_fit(data)
  # Compare with the actual representable centered data, avoiding attributing
  # the rounding of the supplied 1e12-offset values to the statistical method.
  shifted <- transform(data, distal = distal + 1e12)
  origin <- mean(shifted$distal)
  centered <- transform(shifted, distal = distal - origin)
  settings <- expand.grid(variance = c("equal", "varying"),
                          vcov_type = c("independent", "cluster"),
                          contrast = c("none", "pairs"), stringsAsFactors = FALSE)
  invisible(lapply(seq_len(nrow(settings)), function(index) {
    arguments <- as.list(settings[index, ])
    reference <- do.call(three_step, c(list(fit, centered, "distal", method = "ml"), arguments))
    result <- do.call(three_step, c(list(fit, shifted, "distal", method = "ml"), arguments))
    expect_equal(result$standard_error, reference$standard_error, tolerance = 1e-7)
    estimate <- if (arguments$contrast == "none") result$estimate - origin else result$estimate
    expect_equal(estimate, reference$estimate, tolerance = 2e-4)
    if (arguments$contrast == "none") {
      expect_equal(result$variance, reference$variance, tolerance = 1e-7)
      expect_equal(result$effective_n, reference$effective_n, tolerance = 1e-7)
    }
  }))
})
