test_that("normal quadrature integrates polynomial moments", {
  rule <- .ri_quadrature(15L)
  expect_equal(sum(rule$weights), 1, tolerance = 1e-14)
  expect_equal(sum(rule$weights * rule$nodes), 0, tolerance = 1e-13)
  expect_equal(sum(rule$weights * rule$nodes^2), 1, tolerance = 1e-13)
  expect_equal(sum(rule$weights * rule$nodes^4), 3, tolerance = 1e-12)
  expect_error(.ri_quadrature(2L))
})

test_that("one-profile integration agrees with dense multivariate normal algebra", {
  x <- matrix(c(1, 0, 2, 3, 1, -1, 2, 1, 4, 2), 5L, 2L)
  groups <- c(1L, 1L, 2L, 2L, 2L)
  parameters <- list(means = matrix(c(0.4, 1.2), 1L),
    variances = matrix(c(0.7, 1.4), 1L), random_sd = 0.65,
    profile_probabilities = 1)
  actual <- .ri_evaluate(x, groups, parameters, .ri_quadrature(11L), TRUE)
  expected <- vapply(split(seq_len(nrow(x)), groups), function(rows) {
    residual <- as.vector(sweep(x[rows, , drop = FALSE], 2L, parameters$means[1L, ], "-"))
    covariance <- diag(rep(parameters$variances[1L, ], each = length(rows))) + parameters$random_sd^2
    -0.5 * (length(residual) * log(2 * pi) + as.numeric(determinant(covariance, logarithm = TRUE)$modulus) +
      drop(crossprod(residual, solve(covariance, residual))))
  }, numeric(1))
  expect_equal(actual$group_log_likelihood, expected, tolerance = 1e-12)
  expect_equal(actual$log_likelihood, sum(expected), tolerance = 1e-12)
  expect_equal(actual$subject_posteriors, matrix(1, 5L, 1L))
  expect_true(all(actual$random_intercept_sd < parameters$random_sd))
})

test_that("random-intercept mixture approaches ordinary mixture as variance vanishes", {
  x <- matrix(c(-1, 0, 2, 3, -2, 1, 2, 1, 4, 2), 5L, 2L)
  groups <- c(1L, 1L, 2L, 2L, 2L)
  parameters <- list(means = matrix(c(-1, 2, -0.5, 2), 2L),
    variances = matrix(c(0.7, 1, 0.8, 1.2), 2L), random_sd = 1e-10,
    profile_probabilities = c(0.4, 0.6))
  actual <- .ri_evaluate(x, groups, parameters, .ri_quadrature(21L), TRUE)
  component_density <- vapply(seq_len(2L), function(profile) {
    parameters$profile_probabilities[profile] *
      stats::dnorm(x[, 1L], parameters$means[profile, 1L], sqrt(parameters$variances[profile, 1L])) *
      stats::dnorm(x[, 2L], parameters$means[profile, 2L], sqrt(parameters$variances[profile, 2L]))
  }, numeric(5L))
  expect_equal(actual$log_likelihood, sum(log(rowSums(component_density))), tolerance = 1e-12)
  expect_equal(actual$subject_posteriors, component_density / rowSums(component_density), tolerance = 1e-12)
})

test_that("random-intercept estimates recover separated profiles and preserve RNG", {
  set.seed(50)
  group <- rep(seq_len(30L), each = 5L)
  intercept <- stats::rnorm(30L, sd = 0.45)
  profile <- sample(seq_len(2L), 150L, replace = TRUE)
  data <- data.frame(group = group,
    y1 = c(-2, 2)[profile] + intercept[group] + stats::rnorm(150L, sd = 0.55),
    y2 = c(-1.5, 1.5)[profile] + intercept[group] + stats::rnorm(150L, sd = 0.65))
  rng_before <- .Random.seed
  fit <- fit_ml_lpa_random_intercept(data, c("y1", "y2"), "group", 2L, n_starts = 2L, seed = 3L)
  expect_identical(.Random.seed, rng_before)
  expect_s3_class(fit, "ml_lpa_random_intercept")
  expect_false(inherits(fit, "ml_lpa"))
  expect_true(fit$converged)
  expect_false(fit$boundary)
  expect_false(any(fit$boundary_flags))
  expect_equal(fit$optimization_bounds$profile_logit, c(-25, 25))
  expect_gt(fit$random_sd, fit$optimization_bounds$random_sd[1L])
  expect_lt(fit$random_sd, fit$optimization_bounds$random_sd[2L])
  expect_true(fit$quadrature_check_passed)
  expect_equal(unname(fit$means[order(fit$means[, "y1"]), ]), matrix(c(-2, 2, -1.5, 1.5), 2L), tolerance = 0.2)
  expect_equal(fit$random_sd, 0.45, tolerance = 0.15)
  expect_equal(rowSums(fit$subject_posteriors), rep(1, 150L), tolerance = 1e-12)
  expect_length(fit$random_intercept_mean, 30L)
  expect_equal(nobs(fit), 30L)
  expect_equal(as.numeric(logLik(fit)), fit$log_likelihood)
  expect_equal(stats::AIC(fit), fit$aic)
  expect_equal(stats::BIC(fit), fit$bic)
  expect_output(print(fit), "Random-intercept LPA")
  expect_equal(fit$n_parameters, 10L)
  shared <- fit_ml_lpa_random_intercept(data, c("y1", "y2"), "group", 2L,
    variance_model = "equal", n_starts = 1L, seed = 3L)
  expect_equal(unname(shared$variances[1L, ]), unname(shared$variances[2L, ]))
  expect_equal(shared$n_parameters, 8L)
  expect_lte(shared$log_likelihood, fit$log_likelihood + 1e-6)
  expect_warning(coarse <- fit_ml_lpa_random_intercept(data, c("y1", "y2"), "group", 2L,
    quadrature_nodes = 3L, quadrature_check_nodes = 61L, n_starts = 1L, seed = 3L), "Quadrature likelihood check failed")
  expect_false(coarse$quadrature_check_passed)
})

test_that("one-profile univariate optimizer recovers a Gaussian random intercept", {
  set.seed(98)
  group <- rep(seq_len(80L), each = 6L)
  data <- data.frame(group = group, y = 1.7 + stats::rnorm(80L, sd = 0.6)[group] + stats::rnorm(480L, sd = 0.4))
  fit <- fit_ml_lpa_random_intercept(data, "y", "group", 1L, n_starts = 1L, seed = 3L)
  expect_true(fit$converged)
  expect_identical(fit$integration, "analytic")
  expect_equal(fit$quadrature_log_likelihood_difference, 0)
  expect_equal(unname(fit$means[1L, 1L]), 1.7, tolerance = 0.2)
  expect_equal(fit$random_sd, 0.6, tolerance = 0.15)
  expect_equal(unname(fit$variances[1L, 1L]), 0.16, tolerance = 0.05)
  expect_error(fit_ml_lpa_random_intercept(transform(data, y = NA_real_), "y", "group", 1L))
  expect_error(fit_ml_lpa_random_intercept(transform(data, group = 1L), "y", "group", 1L))
  expect_error(fit_ml_lpa_random_intercept(data, "y", "group", 1L, quadrature_nodes = 21L, quadrature_check_nodes = 11L))
  expect_error(fit_ml_lpa_random_intercept(data, "y", "group", 1L, n_starts = 0L))
  expect_error(fit_ml_lpa_random_intercept(data, "y", "group", 1L, seed = -1))
  expect_warning(fit_ml_lpa_random_intercept(data, "y", "group", 1L, n_starts = 1L, max_iter = 1L), "did not converge")
  expect_warning(bounded <- fit_ml_lpa_random_intercept(data, "y", "group", 1L,
    min_variance = 0.5, n_starts = 1L, seed = 3L), "boundary")
  expect_true(bounded$boundary_flags[["residual_variance"]])
  expect_equal(unname(bounded$variances[1L, 1L]), 0.5, tolerance = 1e-6)
})

test_that("random-intercept quadrature agrees with analytic Gaussian integration", {
  x <- matrix(c(-0.3, 0.1, 0.4, 0.8, -0.1, 0.5), 3L, 2L)
  single <- list(means = matrix(c(0.2, 0.3), 1L), variances = matrix(c(0.7, 0.9), 1L),
    random_sd = 0.4, profile_probabilities = 1)
  duplicated <- list(means = single$means[c(1L, 1L), ], variances = single$variances[c(1L, 1L), ],
    random_sd = 0.4, profile_probabilities = c(0.3, 0.7))
  analytic <- .ri_evaluate(x, c(1L, 1L, 2L), single, .ri_quadrature(31L), TRUE)
  numerical <- .ri_evaluate(x, c(1L, 1L, 2L), duplicated, .ri_quadrature(31L), TRUE)
  expect_equal(numerical$log_likelihood, analytic$log_likelihood, tolerance = 1e-12)
  expect_equal(numerical$random_intercept_mean, analytic$random_intercept_mean, tolerance = 1e-12)
  expect_equal(unname(numerical$random_intercept_sd), unname(analytic$random_intercept_sd), tolerance = 1e-12)
  expect_equal(numerical$subject_posteriors, matrix(rep(c(0.3, 0.7), each = 3L), 3L), tolerance = 1e-12)
})

test_that("one-profile random intercept matches genuine Mplus output", {
  fixture <- readRDS(test_path("..", "fixtures", "mplus", "random-intercept-one-profile.rds"))
  fit <- fit_ml_lpa_random_intercept(fixture$data, "y1", "group", 1L,
    n_starts = 1L, tol = 1e-11, seed = 983L)
  actual <- c(fit$variances, fit$means, fit$random_sd^2)
  expect_lt(max(abs(actual - fixture$mplus$parameters)), 1e-5)
  expect_lt(abs(fit$log_likelihood - fixture$mplus$log_likelihood), 5e-5)
  expect_lt(abs(fit$aic - fixture$mplus$aic), 1e-4)
  expect_lt(abs(fit$bic_individual - fixture$mplus$bic), 1e-4)
  expect_equal(fit$n_parameters, fixture$mplus$n_parameters)
})

test_that("mixture quadrature agrees with independent adaptive integration", {
  x <- matrix(c(-0.5, 1, 0.2, 0.7, 0.9, -0.3), 3L, 2L)
  parameters <- list(means = matrix(c(-0.6, 0.8, -0.2, 0.6), 2L),
    variances = matrix(c(0.5, 0.8, 0.6, 0.9), 2L),
    random_sd = 0.35, profile_probabilities = c(0.3, 0.7))
  value <- .ri_evaluate(x, c(1L, 1L, 2L), parameters, .ri_quadrature(61L))
  independent <- vapply(list(1:2, 3L), function(rows) {
    integrand <- function(intercepts) {
      vapply(intercepts, function(intercept) {
        densities <- vapply(seq_len(2L), function(profile) {
          parameters$profile_probabilities[profile] *
            stats::dnorm(x[rows, 1L], parameters$means[profile, 1L] + intercept, sqrt(parameters$variances[profile, 1L])) *
            stats::dnorm(x[rows, 2L], parameters$means[profile, 2L] + intercept, sqrt(parameters$variances[profile, 2L]))
        }, numeric(length(rows)))
        densities <- matrix(densities, length(rows), 2L)
        prod(rowSums(densities)) * stats::dnorm(intercept, sd = parameters$random_sd)
      }, numeric(1))
    }
    log(stats::integrate(integrand, -Inf, Inf, rel.tol = 1e-11)$value)
  }, numeric(1))
  expect_equal(unname(value$group_log_likelihood), independent, tolerance = 1e-10)
})
