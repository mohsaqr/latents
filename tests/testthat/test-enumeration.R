test_that("enumeration retains failures and both BIC conventions", {
  set.seed(12)
  d <- data.frame(g = rep(1:30, each = 10), y = c(rnorm(150, -2), rnorm(150, 2)))
  result <- enumerate_ml_lpa(d, "y", "g", profiles = 1:2, group_classes = 1:2,
                             n_starts = 2, seed = 12)
  expect_equal(nrow(result$table), 4L)
  expect_length(result$fits, 4L)
  expect_null(result$fits[[3]])
  expect_match(result$table$error[3], "profile")
  expect_true(all(is.finite(result$table$aic[c(1, 2, 4)])))
  expect_equal(result$table$bic_groups[1], result$fits[[1]]$bic)
  expect_equal(result$table$bic_individual[1], result$fits[[1]]$bic_individual)
  expect_true(is.na(result$table$profile_entropy[1]))
  expect_gt(result$table$profile_entropy[2], 0)
  expect_lte(result$table$profile_entropy[2], 1)
  expect_error(enumerate_ml_lpa(d, "y", "g", profiles = 0), "profiles")
  expect_error(enumerate_ml_lpa(d, "y", "g", n_profiles = 2), "model counts")
})

test_that("simulation preserves group layout and full covariance moments", {
  set.seed(9)
  noise <- matrix(rnorm(2000), 1000, 2) %*% chol(matrix(c(1, .6, .6, 2), 2))
  d <- data.frame(g = rep(1:100, each = 10), a = noise[, 1], b = noise[, 2])
  model <- fit_ml_lpa(d, c("a", "b"), "g", 1, 1, covariance_model = "full", n_starts = 1)
  simulated <- .ml_lpa_simulate(model)
  expect_identical(simulated$g, d$g)
  expect_identical(names(simulated), c("a", "b", "g"))
  expect_equal(cov(simulated[c("a", "b")]), cov(d[c("a", "b")]), tolerance = .2)
  expect_false(anyNA(simulated))
  one <- fit_ml_lpa(d, "a", "g", 1, 1, covariance_model = "full", n_starts = 1)
  single_indicator <- .ml_lpa_simulate(one)
  expect_equal(dim(single_indicator), c(1000L, 2L))
  expect_identical(single_indicator$g, d$g)
  expect_equal(mean(single_indicator$a), mean(d$a), tolerance = .1)
})

test_that("bootstrap withholds p-values if a replicate cannot be fitted", {
  set.seed(53)
  d <- data.frame(g = rep(1:20, each = 5), y = c(rnorm(50, -2), rnorm(50, 2)))
  small <- fit_ml_lpa(d, "y", "g", 1, 1, n_starts = 2, seed = 12)
  large <- fit_ml_lpa(d, "y", "g", 2, 1, n_starts = 2, seed = 12)
  testthat::local_mocked_bindings(fit_ml_lpa = function(...) stop("test optimization failure"))
  expect_warning(result <- bootstrap_lrt_ml_lpa(small, large, d, n_boot = 2, seed = 1), "p_value is NA")
  expect_true(is.na(result$p_value))
  expect_equal(result$n_valid, 0L)
  expect_true(all(result$replicates$error == "test optimization failure"))
})

test_that("bootstrap refits generated data and reports finite simulation correction", {
  set.seed(53)
  d <- data.frame(g = rep(1:30, each = 8), y = rnorm(240))
  small <- fit_ml_lpa(d, "y", "g", 1, 1, variance_model = "equal", n_starts = 2, seed = 8)
  large <- fit_ml_lpa(d, "y", "g", 2, 1, variance_model = "equal", n_starts = 3, seed = 8,
                      max_iter = 3000, tol = 1e-7)
  rng <- .Random.seed
  result <- bootstrap_lrt_ml_lpa(small, large, d, n_boot = 3, n_starts = 3,
                                  max_iter = 3000, tol = 1e-7, seed = 42)
  expect_identical(.Random.seed, rng)
  expect_equal(result$n_valid, 3L)
  expect_equal(result$statistic, max(0, 2 * (large$log_likelihood - small$log_likelihood)))
  expect_equal(result$p_value, (1 + sum(result$replicates$statistic >= result$statistic)) / 4)
  expect_equal(nrow(result$replicates), 3L)
  expect_true(all(is.na(result$replicates$error)))
  expect_true(is.finite(result$monte_carlo_se))
  expect_error(bootstrap_lrt_ml_lpa(small, small, d, n_boot = 2), "one class")
  changed <- d
  changed$y <- d$y + 1
  expect_error(bootstrap_lrt_ml_lpa(small, large, changed, n_boot = 2), "reproduce")
  changed$y[1] <- NA_real_
  expect_error(bootstrap_lrt_ml_lpa(small, large, changed, n_boot = 2), "complete")
})

test_that("full covariance print and summary expose residual matrices", {
  set.seed(3)
  d <- data.frame(g = rep(1:10, each = 10), a = rnorm(100), b = rnorm(100))
  fit <- fit_ml_lpa(d, c("a", "b"), "g", 1, 1, covariance_model = "full", n_starts = 1)
  expect_output(print(fit), "full residual covariance")
  expect_equal(summary(fit)$covariances, fit$covariances)
  expect_output(print(summary(fit)), "covariance matrices")
})
