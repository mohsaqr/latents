test_that("one-step covariate model agrees with genuine multilevel Mplus ML", {
  fixture <- readRDS(test_path("..", "fixtures", "mplus", "twolevel-covariates.rds"))
  target <- fixture$expected
  fit <- fit_multilpa_covariates(fixture$data, c("y1", "y2"), "g", 2, 2, "z", "w",
                               n_starts = 5, seed = 812, tol = 1e-12)
  p <- order(fit$means[, 1])
  beta <- unname(fit$profile_coefficients)
  if (p[1] == 2L) beta <- -beta
  h <- order(beta[1:2, 1], decreasing = TRUE)
  beta[1:2, ] <- beta[h, ]
  gamma <- unname(fit$group_coefficients) * if (h[1] == 2L) -1 else 1
  expect_true(fit$converged)
  expect_false(fit$boundary)
  expect_equal(fit$n_parameters, target$n_parameters)
  expect_lt(abs(fit$log_likelihood - target$log_likelihood), 5e-5)
  expect_lt(max(abs(fit$means[p, ] - target$means)), 1e-5)
  expect_lt(max(abs(fit$variances[p, ] - target$variances)), 1e-5)
  expect_lt(max(abs(beta - target$profile_coefficients)), 1e-5)
  expect_lt(max(abs(gamma - target$group_coefficients)), 1e-5)
  expect_lt(max(abs(fit$subject_posteriors[, p] - target$subject_posteriors)), 1e-5)
  expect_lt(max(abs(fit$group_posteriors[, h] - target$group_posteriors)), 1e-5)
  expect_lt(abs(fit$aic - target$aic), 5e-4)
  expect_lt(abs(fit$bic_individual - target$bic_individual), 5e-4)
})
