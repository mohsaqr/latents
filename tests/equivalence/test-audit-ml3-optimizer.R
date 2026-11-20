# Comparisons with independent references moved from tests/testthat/test-audit-ml3-optimizer.R,
# so that they run in CI's equivalence job rather than under R CMD check.

test_that("equal-variance step three reaches the retained higher regular maximum", {
  skip_on_cran()
  case <- readRDS(test_path("fixtures", "ml3-equal-rare-class.rds"))
  fit <- latents:::.multilpa_step_ml_fit(case$pieces, case$values, "equal")
  expect_true(fit$converged)
  expect_false(fit$at_floor)
  # The fixture's known_* values: an independent BFGS/Newton maximisation of
  # the step-three likelihood at the retained posteriors (case$source).
  retained <- "independent: ML three-step likelihood by BFGS/Newton at fixed posteriors (ml3 fixture)"
  note <- "Package side: the internal ML step three (.multilpa_step_ml_fit, .multilpa_step_ml_influence) behind three_step(method = \"ml\"); reference values stored in the fixture without a provenance header."
  expect_reference(fit$log_likelihood, case$known_log_likelihood, reference = retained,
                   covers = "three_step()", quantity = "step-three log likelihood, equal variances, rare class",
                   tolerance = 1e-7, measure = "absolute", note = note)
  # The means in units of their standard errors, as the original bound.
  classes <- sprintf("class_%d", seq_along(case$known_means))
  expect_reference(fit$means / case$known_standard_error,
                   case$known_means / case$known_standard_error, reference = retained,
                   covers = "three_step()", quantity = "class means in standard-error units",
                   tolerance = 1e-6, measure = "absolute", key = classes, note = note)
  influence <- latents:::.multilpa_step_ml_influence(fit, case$values, "equal")
  covariance <- crossprod(influence) * length(case$values) / (length(case$values) - 1)
  covariance <- matrix(covariance, length(classes), dimnames = list(classes, classes))
  expect_reference(covariance,
                   matrix(case$known_robust_covariance, length(classes), dimnames = list(classes, classes)),
                   reference = retained, covers = "three_step()",
                   quantity = "robust covariance of the class means",
                   tolerance = 1e-6, measure = "relative", note = note)

  # Rebuild the conditional mixture likelihood directly from P, modal
  # assignments and normal densities. No production likelihood or score is
  # used by this numerical information/sandwich oracle.
  skip_without_reference("numDeriv")
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
  # Unit checks on the reference itself: positive definite, and the package's
  # point is its stationary point.
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
  # The reference is stable under the step refinement (a check of the reference).
  expect_equal(independent_covariance, refined_covariance,
               tolerance = 1e-5, ignore_attr = TRUE)
  # 1e-5 relative, the original tolerance: numDeriv's Hessian at a step of
  # 0.03 is accurate to about that on this curved rare-class contribution.
  # Every entry is larger than 1 in magnitude, so relative and scaled agree.
  expect_reference(covariance,
                   matrix(independent_covariance, k, dimnames = list(classes, classes)),
                   reference = "numDeriv", covers = "three_step()",
                   quantity = "robust covariance of the class means from numerical scores and information",
                   tolerance = 1e-5, measure = "relative",
                   note = "Package side: the internal influence functions (.multilpa_step_ml_influence) behind three_step(method = \"ml\"); reference: numDeriv::jacobian() and numDeriv::hessian() of a step-three likelihood written in the test.")
})
