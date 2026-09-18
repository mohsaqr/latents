test_that("weighted multinomial optimizer reproduces fractional counts", {
  set.seed(64)
  design <- cbind(1, rnorm(100))
  coefficients <- matrix(c(.2, .7, -.4, -.5), 2, 2)
  probabilities <- .ml_lpa_softmax(design, coefficients)
  fitted <- .ml_lpa_weighted_logits(design, probabilities * 3,
                                     matrix(0, 2, 2))
  expect_true(fitted$converged)
  expect_equal(fitted$coefficients, coefficients, tolerance = 1e-5)
  expect_equal(rowSums(probabilities), rep(1, 100), tolerance = 1e-14)
  expect_equal(.ml_lpa_softmax(matrix(1, 3, 1), matrix(numeric(), 1, 0)), matrix(1, 3, 1))
  expect_equal(.ml_lpa_weighted_logits(matrix(1, 3, 1), matrix(1, 3, 1),
                                         matrix(numeric(), 1, 0))$coefficients,
               matrix(numeric(), 1, 0))
})

test_that("covariate expectation reduces to the independently verified original", {
  set.seed(65)
  x <- matrix(rnorm(24), 12, 2)
  index <- rep(1:4, each = 3)
  parameters <- list(means = rbind(c(-1, -1), c(1, 1)),
                     variances = matrix(1, 2, 2),
                     profile_probabilities = rbind(c(.8, .2), c(.3, .7)),
                     group_probabilities = c(.6, .4))
  designs <- list(matrix(rep(c(1, 0), 12), 12, 2, byrow = TRUE),
                  matrix(rep(c(0, 1), 12), 12, 2, byrow = TRUE))
  beta <- matrix(qlogis(c(.8, .3)), 2, 1)
  gamma <- matrix(qlogis(.6), 1, 1)
  original <- .ml_lpa_expectation(x, index, parameters)
  result <- .ml_lpa_cov_expectation(x, index, parameters, designs, matrix(1, 4, 1), beta, gamma)
  expect_equal(result$log_likelihood, original$log_likelihood, tolerance = 1e-12)
  expect_equal(result$subject_posteriors, original$subject_posteriors, tolerance = 1e-12)
  expect_equal(result$group_posteriors, original$group_posteriors, tolerance = 1e-12)
})

test_that("one-step covariates recover slopes and have consistent diagnostics", {
  set.seed(223)
  group <- rep(1:60, each = 12)
  z <- rnorm(length(group))
  group_w <- rnorm(60)
  w <- group_w[group]
  h <- rbinom(60, 1, plogis(-.2 + .7 * group_w))[group]
  k <- rbinom(length(group), 1, plogis(-1.5 + 3 * h + .9 * z))
  d <- data.frame(group = group, z = z, w = w,
                  y1 = rnorm(length(group), ifelse(k == 1, -3, 3), .5),
                  y2 = rnorm(length(group), ifelse(k == 1, -2, 2), .6))
  rng <- .Random.seed
  fit <- fit_ml_lpa_covariates(d, c("y1", "y2"), "group", 2, 2,
                               "z", "w", n_starts = 3, seed = 23, tol = 1e-10)
  expect_identical(.Random.seed, rng)
  expect_true(fit$converged)
  expect_true(all(diff(fit$log_likelihood_history) > -1e-6))
  expect_equal(fit$n_parameters, 13)
  expect_equal(dim(fit$profile_coefficients), c(3L, 1L))
  expect_equal(dim(fit$group_coefficients), c(2L, 1L))
  expect_lt(abs(abs(fit$profile_coefficients["z", ]) - .9), .2)
  expect_lt(abs(abs(fit$group_coefficients["w", ]) - .7), .2)
  expect_equal(unname(rowSums(fit$group_priors)), rep(1, 60), tolerance = 1e-12)
  expect_equal(unname(rowSums(fit$subject_posteriors)), rep(1, nrow(d)), tolerance = 1e-12)
  expect_equal(unname(rowSums(fit$group_posteriors)), rep(1, 60), tolerance = 1e-12)
  expect_equal(as.numeric(logLik(fit)), fit$log_likelihood)
  expect_equal(nobs(fit), 60L)
  expect_equal(AIC(fit), fit$aic)
  expect_equal(BIC(fit), fit$bic)
  expect_output(print(fit), "with covariates")
  # Profiles are practically observed here: independent binomial ML checks slope.
  target <- as.numeric(k == if (fit$means[1, 1] < 0) 1 else 0)
  independent <- glm(target ~ factor(h) + z, family = binomial())
  expect_lt(abs(fit$profile_coefficients["z", ] - coef(independent)["z"]), .025)
  bad <- d
  bad$w[1] <- 20
  expect_error(fit_ml_lpa_covariates(bad, c("y1", "y2"), "group", 2, 2, "z", "w", n_starts = 1), "constant within")
  bad$z[1] <- NA_real_
  expect_error(fit_ml_lpa_covariates(bad, c("y1", "y2"), "group", 2, 2, "z", "w"), "finite numeric")
  bad$z <- I(cbind(d$z, d$z))
  expect_error(fit_ml_lpa_covariates(bad, c("y1", "y2"), "group", 2, 2, "z", "w"), "finite numeric")
  expect_error(fit_ml_lpa_covariates(d, c("y1", "y2"), "group", 2, 2, "y1"), "distinct")
})

test_that("no predictors agrees with the base model and single-level works", {
  set.seed(25)
  d <- data.frame(g = rep(1:30, each = 10), y = c(rnorm(150, -3), rnorm(150, 3)), z = rnorm(300))
  base <- fit_ml_lpa(d, "y", "g", 2, 1, variance_model = "equal", n_starts = 2, seed = 42, tol = 1e-10)
  fit <- fit_ml_lpa_covariates(d, "y", "g", 2, 1, variance_model = "equal", n_starts = 2, seed = 42, tol = 1e-10)
  expect_equal(fit$log_likelihood, base$log_likelihood, tolerance = 1e-7)
  expect_equal(fit$n_parameters, base$n_parameters)
  expect_equal(ncol(fit$group_coefficients), 0L)
  single <- fit_ml_lpa_covariates(d, "y", "g", 1, 1, n_starts = 1, seed = 42)
  expect_equal(as.numeric(single$means), mean(d$y), tolerance = 1e-12)
  expect_equal(single$n_parameters, 2L)
  expect_error(fit_ml_lpa_covariates(d, "y", "g", 1, 1, "z", n_starts = 1), "at least two profiles")
  d$constant <- 1
  expect_error(fit_ml_lpa_covariates(d, "y", "g", 2, 1, "constant", n_starts = 1), "rank deficient")
})
