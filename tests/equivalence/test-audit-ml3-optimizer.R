# Comparisons with independent references moved from tests/testthat/test-audit-ml3-optimizer.R,
# so that they run in CI's equivalence job rather than under R CMD check.

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
