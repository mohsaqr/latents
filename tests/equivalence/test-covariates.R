
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
  fit <- multilpa(d, c("y1", "y2"), "group", 2, 2,
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
  expect_error(multilpa(bad, c("y1", "y2"), "group", 2, 2, "z", "w", n_starts = 1), "constant within")
  bad$z[1] <- NA_real_
  expect_error(multilpa(bad, c("y1", "y2"), "group", 2, 2, "z", "w"), "finite numeric")
  bad$z <- I(cbind(d$z, d$z))
  expect_error(multilpa(bad, c("y1", "y2"), "group", 2, 2, "z", "w"), "finite numeric")
  expect_error(multilpa(d, c("y1", "y2"), "group", 2, 2, "y1"), "distinct")
})
