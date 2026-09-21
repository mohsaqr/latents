test_that("one-profile random intercept matches genuine Mplus output", {
  fixture <- readRDS(equivalence_fixture("mplus", "random-intercept-one-profile.rds"))
  fit <- fit_random_intercept(fixture$data, "y1", "group", 1L,
    n_starts = 1L, tol = 1e-11, seed = 983L)
  actual <- c(fit$variances, fit$means, fit$random_sd^2)
  expect_lt(max(abs(actual - fixture$mplus$parameters)), 1e-5)
  expect_lt(abs(fit$log_likelihood - fixture$mplus$log_likelihood), 5e-5)
  expect_lt(abs(fit$aic - fixture$mplus$aic), 1e-4)
  expect_lt(abs(fit$bic_individual - fixture$mplus$bic), 1e-4)
  expect_equal(fit$n_parameters, fixture$mplus$n_parameters)
})

