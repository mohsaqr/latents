# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
# The comparisons with independent references are in tests/equivalence/test-review-1020-inference.R.
skip_on_cran()

test_that("empty groups do not evade robust validity or alter CR1", {
  set.seed(792)
  make_fit <- function(data, ...) quietly(multilpa(data, "y", "g", 1, 1,
    missing = "fiml", n_starts = 1, seed = 3, tol = 1e-12, ...))
  invisible(lapply(1:2, function(n_groups) {
    data <- data.frame(g = rep(seq_len(n_groups), each = 50),
                       y = stats::rnorm(50 * n_groups))
    padded <- rbind(data, data.frame(g = n_groups + seq_len(20), y = NA_real_))
    invisible(lapply(c("robust", "opg"), function(type) {
      expect_error(parameter_inference(make_fit(data), vcov_type = type),
                   class = "latents_too_few_groups")
      expect_error(parameter_inference(make_fit(padded), vcov_type = type),
                   class = "latents_too_few_groups")
      expect_error(parameter_inference(make_fit(padded, drop_empty = FALSE), vcov_type = type),
                   class = "latents_too_few_groups")
    }))
  }))
  data <- data.frame(g = rep(1:6, each = 30), y = stats::rnorm(180))
  padded <- rbind(data, data.frame(g = 7:26, y = NA_real_))
  invisible(lapply(c("robust", "opg"), function(type) {
    a <- parameter_inference(make_fit(data), vcov_type = type)
    b <- parameter_inference(make_fit(padded), vcov_type = type)
    expect_equal(attr(b, "covariance"), attr(a, "covariance"), tolerance = 1e-8)
    expect_equal(attr(b, "group_scores")[1:6, ], attr(a, "group_scores"),
                 tolerance = 1e-8)
    expect_lt(max(abs(attr(b, "group_scores")[-(1:6), ])), 1e-12)
  }))
  # Legacy counting intentionally uses every empty group in CR1.
  legacy <- parameter_inference(make_fit(padded, drop_empty = FALSE), vcov_type = "robust")
  current <- parameter_inference(make_fit(padded), vcov_type = "robust")
  expect_equal(attr(legacy, "covariance"),
               attr(current, "covariance") * (26 / 25) / (6 / 5), tolerance = 1e-8)
  # Normalizing weights changes the pseudo likelihood, but cancels in the sandwich.
  data$weight <- rep(c(1, 2, 1, 3, 2, 1), each = 30)
  padded <- rbind(data, data.frame(g = 7:26, y = NA_real_, weight = 2))
  zero <- rbind(data, data.frame(g = 7:26, y = stats::rnorm(20), weight = 0))
  weighted <- function(d) parameter_inference(make_fit(d, weights = "weight"))
  expect_equal(attr(weighted(padded), "covariance"), attr(weighted(data), "covariance"),
               tolerance = 1e-8)
  expect_equal(attr(weighted(zero), "covariance"), attr(weighted(data), "covariance"),
               tolerance = 1e-8)
  # One-dimensional scores have full numerical rank even when the lone
  # informative unit's score is only roundoff. Legacy CR1 cannot authorize
  # an otherwise invalid covariance, including with fixed measurement.
  one <- subset(data, g == 1)
  start <- starting_values(make_fit(one))
  one_padded <- rbind(one, data.frame(g = 7:26, y = NA_real_, weight = 1))
  fixed <- make_fit(one_padded, start = start, fixed = "means", drop_empty = FALSE)
  expect_error(parameter_inference(fixed, vcov_type = "robust"),
               class = "latents_too_few_groups")
})

test_that("direct-effect translation stays regular with missing indicators and sampling weights", {
  set.seed(215)
  n <- 400
  data <- data.frame(g = rep(1:40, each = 10), z = stats::rnorm(n),
                     w = stats::rnorm(n), weight = rep(c(1, 2), each = 200))
  class <- sample(c(-3, 3), n, TRUE)
  data$a <- class + 0.6 * data$z + stats::rnorm(n)
  data$b <- class - 0.4 * data$z + 0.3 * data$w + stats::rnorm(n)
  data$a[seq(5, n, 11)] <- NA_real_
  data$b[seq(8, n, 13)] <- NA_real_
  fit <- function(d) quietly(multilpa(d, c("a", "b"), "g", 2, 1,
    direct_effects = c("a ~ z", "b ~ z + w"), weights = "weight", missing = "fiml",
    n_starts = 1, seed = 1, tol = 1e-12))
  a <- fit(data)
  b <- fit(transform(data, z = z + 2026, w = w - 1000))
  expect_equal(a$log_likelihood, b$log_likelihood, tolerance = 1e-8)
  ia <- parameter_inference(a)
  ib <- parameter_inference(b)
  effects <- ia$parameter == "direct_effect"
  expect_equal(ib$estimate[effects], ia$estimate[effects], tolerance = 2e-6)
  expect_equal(ib$standard_error[effects], ia$standard_error[effects], tolerance = 2e-6)
  expect_equal(ib$p_value[effects], ia$p_value[effects], tolerance = 2e-6)
  expect_gt(attr(ib, "condition_ratio"), 1e-5)
})

test_that("irrelevant direct covariates cannot destabilize fitting or inference", {
  set.seed(803)
  data <- data.frame(g = rep(1:20, each = 20), z = stats::rnorm(400),
                     w = stats::rnorm(400), weight = 1)
  data$y <- 0.8 + 0.7 * data$z - 0.4 * data$w + stats::rnorm(400)
  fit <- function(d, ...) quietly(multilpa(d, "y", "g", 1, 1,
    direct_effects = "y ~ z + w", missing = "fiml", n_starts = 1, seed = 1,
    tol = 1e-12, ...))
  base <- fit(data)
  padding <- data.frame(g = 21:30, z = 1e5, w = -1e5, weight = 1, y = NA_real_)
  padded <- fit(rbind(data, padding))
  expect_true(padded$converged)
  expect_equal(padded$log_likelihood, base$log_likelihood, tolerance = 1e-9)
  expect_equal(parameter_inference(padded, vcov_type = "robust")$standard_error,
               parameter_inference(base, vcov_type = "robust")$standard_error,
               tolerance = 2e-7)
  padding$y <- stats::rnorm(10)
  padding$weight <- 0
  weighted <- fit(data, weights = "weight")
  zero_padded <- fit(rbind(data, padding), weights = "weight")
  expect_true(zero_padded$converged)
  expect_equal(zero_padded$direct_effects, weighted$direct_effects, tolerance = 2e-7)
  expect_equal(parameter_inference(zero_padded)$standard_error,
               parameter_inference(weighted)$standard_error, tolerance = 2e-7)
})

test_that("covariate and both transition engines preserve robust inference under empty padding", {
  set.seed(34)
  data <- data.frame(g = rep(1:40, each = 5), time = rep(1:5, 40),
                     z = stats::rnorm(200))
  data$y <- stats::rnorm(200, sample(c(-3, 3), 200, TRUE))
  padding <- data.frame(g = rep(41:60, each = 5), time = rep(1:5, 20),
                        z = stats::rnorm(100), y = NA_real_)
  fitters <- list(
    covariate = function(d, ...) quietly(multilpa(d, "y", "g", 2, 1,
      profile_covariates = "z", missing = "fiml", n_starts = 1, seed = 1,
      tol = 1e-12, ...)),
    homogeneous = function(d, ...) quietly(lta(d, "y", "g", 2, time = "time",
      missing = "fiml", n_starts = 1, seed = 1, tol = 1e-12, ...)),
    general = function(d, ...) quietly(lta(d, "y", "g", 2, time = "time",
      transitions = "occasion", missing = "fiml", n_starts = 1, seed = 1,
      tol = 1e-12, ...)))
  covariance <- function(fit, type) {
    if (inherits(fit, "multilpa_lta")) vcov(fit, type = type) else
      vcov(fit, vcov_type = type)
  }
  invisible(lapply(fitters, function(fitter) {
    base <- fitter(data)
    padded <- fitter(rbind(data, padding))
    expect_true(base$converged && padded$converged)
    expect_equal(base$log_likelihood, padded$log_likelihood, tolerance = 1e-9)
    invisible(lapply(c("robust", "opg"), function(type) {
      expect_equal(covariance(padded, type), covariance(base, type), tolerance = 3e-6)
    }))
    legacy <- fitter(rbind(data, padding), drop_empty = FALSE)
    expect_equal(covariance(legacy, "robust"), covariance(padded, "robust") *
                   (60 / 59) / (40 / 39), tolerance = 3e-6)
    data$weight <- rep(c(1, 2), each = 100)
    padding$weight <- 2
    weighted_base <- fitter(data, weights = "weight")
    weighted_pad <- fitter(rbind(data, padding), weights = "weight")
    expect_equal(covariance(weighted_pad, "robust"),
                 covariance(weighted_base, "robust"), tolerance = 3e-6)
  }))
})
