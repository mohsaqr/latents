test_that("one-step covariate model agrees with genuine multilevel Mplus ML", {
  fixture <- readRDS(equivalence_fixture("mplus", "twolevel-covariates.rds"))
  target <- fixture$expected
  fit <- multilpa(fixture$data, c("y1", "y2"), "g", 2, 2, "z", "w",
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

.covariate_fit <- function(data, ...) {
  multilpa(data, c("y1", "y2"), "g", 2, 2, "z", "w",
           n_starts = 10, seed = 812, tol = 1e-12, ...)
}

test_that("standard errors reproduce a genuine Mplus covariate run", {
  fixture <- readRDS(equivalence_fixture("mplus", "twolevel-covariates.rds"))
  fit <- .covariate_fit(fixture$data)
  inference <- parameter_inference(fit, fixture$data)
  membership <- subset(inference, parameter == "coefficient")
  value <- function(term) membership$estimate[membership$term == term]
  error <- function(term) membership$standard_error[membership$term == term]

  # Mplus 9 covariates.out, MODEL RESULTS. Its group-class reference is the
  # opposite of this package's, so the group-level signs are mirrored.
  expect_equal(value("z"), -1.008, tolerance = 5e-4)
  expect_equal(error("z"), 0.117, tolerance = 5e-3)
  # Both group-class labellings are the same maximum: every start reaches it to
  # within 3e-11, so which label carries which intercept is decided by the last
  # bits and is not a property of the model. Assert the pair, not the labelling
  # -- the same reason `w` and `(Intercept)` below are wrapped in `abs()`.
  group_estimates <- sort(membership$estimate[grepl("^group_class",
                                                    membership$term)])
  expect_equal(group_estimates, sort(c(1.662, 1.662 - 3.233)),
               tolerance = 5e-4)
  # The intercept Mplus reports carries its standard error whichever label it
  # lands on, so pair the error with the estimate rather than with the name.
  reported <- which.min(abs(membership$estimate - 1.662))
  expect_equal(membership$standard_error[reported], 0.152, tolerance = 5e-3)
  expect_equal(abs(value("w")), 0.740, tolerance = 5e-4)
  expect_equal(error("w"), 0.311, tolerance = 5e-3)
  expect_equal(abs(value("(Intercept)")), 0.214, tolerance = 5e-3)
  expect_equal(error("(Intercept)"), 0.278, tolerance = 5e-3)
})

test_that("the group-class contrast matches Mplus through vcov()", {
  fixture <- readRDS(equivalence_fixture("mplus", "twolevel-covariates.rds"))
  fit <- .covariate_fit(fixture$data)
  covariance <- vcov(fit, fixture$data)
  rows <- grep("^profile[.]coefficient[.]profile_1[.]group_class", rownames(covariance))
  contrast <- c(-1, 1)
  inference <- parameter_inference(fit, fixture$data)
  estimates <- inference$estimate[grepl("^group_class", inference$term) &
                                    inference$level == "profile"]

  # Mplus reports CW#1 ON CB#1 = -3.233 (S.E. 0.226); this package carries one
  # intercept per group class, whose difference is the same quantity.
  # Signed only up to the group-class labelling; its magnitude is the quantity.
  expect_equal(abs(sum(contrast * estimates)), 3.233, tolerance = 1e-3)
  expect_equal(sqrt(drop(contrast %*% covariance[rows, rows] %*% contrast)),
               0.226, tolerance = 5e-3)
})

