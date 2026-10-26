# Fixed-posterior regression for a regular equal-variance local maximum.
# The small fixture retains simulated data/posteriors and an independently
# optimized reference point, never read by the production optimizer.

test_that("equal-variance step three preserves the observed-class Gaussian solution", {
  pieces <- list(posteriors = rbind(cbind(rep(1, 30), 0), cbind(rep(0, 30), 1)),
                 modal = rep(1:2, each = 30), n_classes = 2L)
  values <- c(rep(c(-2, -1, 0), 10), rep(c(3, 4, 5), 10))
  fit <- .multilpa_step_ml_fit(pieces, values, "equal")
  expect_equal(fit$means, c(-1, 4), tolerance = 1e-12)
  expect_equal(fit$variances, rep(2 / 3, 2), tolerance = 1e-12)
  expected <- sum(stats::dnorm(values, rep(c(-1, 4), each = 30), sqrt(2 / 3),
                                log = TRUE)) - 60 * log(2)
  expect_equal(fit$log_likelihood, expected, tolerance = 1e-12)
  expect_true(fit$converged)
  expect_false(fit$at_floor)
})

test_that("equal-variance step three reaches the retained higher regular maximum", {
  skip_on_cran()
  case <- readRDS(test_path("fixtures", "ml3-equal-rare-class.rds"))
  fit <- .multilpa_step_ml_fit(case$pieces, case$values, "equal")
  expect_true(fit$converged)
  expect_false(fit$at_floor)
  expect_lt(abs(fit$log_likelihood - case$known_log_likelihood), 1e-7)
  expect_lt(max(abs((fit$means - case$known_means) / case$known_standard_error)), 1e-6)
  influence <- .multilpa_step_ml_influence(fit, case$values, "equal")
  covariance <- crossprod(influence) * length(case$values) / (length(case$values) - 1)
  expect_equal(covariance, case$known_robust_covariance, tolerance = 1e-6,
               ignore_attr = TRUE)

  # Rebuild the conditional mixture likelihood directly from P, modal
  # assignments and normal densities. No production likelihood or score is
  # used by this numerical information/sandwich oracle.
  skip_if_not_installed("numDeriv")
  posterior <- case$pieces$posteriors
  modal <- case$pieces$modal
  k <- ncol(posterior)
  errors <- crossprod(posterior, outer(modal, seq_len(k), "==") * 1) / colSums(posterior)
  prior <- colMeans(posterior)
  values <- (case$values - mean(case$values)) / stats::sd(case$values)
  unit_loglik <- function(theta) {
    joint <- vapply(seq_len(k), function(j) {
      log(prior[j]) + log(errors[j, modal]) +
        stats::dnorm(values, theta[j], exp(theta[k + 1]), log = TRUE)
    }, numeric(length(values)))
    peak <- apply(joint, 1L, max)
    peak + log(rowSums(exp(joint - peak)))
  }
  point <- c(fit$outcome_coordinates$means,
              log(sqrt(fit$variances[1]) / fit$outcome_coordinates$scale))
  # The default ten-percent Hessian step is too wide for the rare class's
  # curved mixture contribution. Check covariance stability under a further
  # step refinement rather than treating the default as exact arithmetic.
  information <- -numDeriv::hessian(function(theta) sum(unit_loglik(theta)), point,
                                     method.args = list(d = .03))
  refined <- -numDeriv::hessian(function(theta) sum(unit_loglik(theta)), point,
                                 method.args = list(d = .01))
  expect_gt(min(eigen(information, symmetric = TRUE, only.values = TRUE)$values), 0)
  expect_lt(max(abs(numDeriv::grad(function(theta) sum(unit_loglik(theta)), point))), 1e-5)
  independent_scores <- numDeriv::jacobian(unit_loglik, point)
  independent_influence <- (independent_scores %*% solve(information))[, seq_len(k)] *
    stats::sd(case$values)
  independent_covariance <- crossprod(independent_influence) *
    length(values) / (length(values) - 1)
  refined_influence <- (independent_scores %*% solve(refined))[, seq_len(k)] *
    stats::sd(case$values)
  refined_covariance <- crossprod(refined_influence) * length(values) / (length(values) - 1)
  expect_equal(independent_covariance, refined_covariance,
               tolerance = 1e-5, ignore_attr = TRUE)
  expect_equal(covariance, independent_covariance, tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("equal-variance optimizer basins respect outcome units and origin", {
  skip_on_cran()
  case <- readRDS(test_path("fixtures", "ml3-equal-rare-class.rds"))
  baseline <- .multilpa_step_ml_fit(case$pieces, case$values, "equal")
  invisible(lapply(list(c(scale = 1e-4, offset = 0), c(scale = 1e4, offset = 0),
                       c(scale = 1, offset = 1e8)), function(change) {
    values <- change[["offset"]] + change[["scale"]] * case$values
    fit <- .multilpa_step_ml_fit(case$pieces, values, "equal")
    expect_lt(max(abs((fit$means - change[["offset"]]) / change[["scale"]] -
                       baseline$means) / case$known_standard_error), 1e-6)
    expect_equal(fit$variances / change[["scale"]]^2, baseline$variances,
                 tolerance = 1e-7)
    expect_lt(abs(fit$log_likelihood + length(values) * log(change[["scale"]]) -
                    baseline$log_likelihood), 1e-6)
    expect_true(fit$converged)
    expect_false(fit$at_floor)
  }))
})
