test_that("enumeration retains failures and both BIC conventions", {
  set.seed(12)
  d <- data.frame(g = rep(1:30, each = 10), y = c(rnorm(150, -2), rnorm(150, 2)))
  result <- enumerate_classes(d, "y", "g", n_profiles = 1:2, n_group_classes = 1:2,
                             n_starts = 2, seed = 12)
  grid <- as.data.frame(result)
  expect_equal(nrow(grid), 4L)
  expect_match(grid$error[3], "profile")
  expect_error(candidate_fit(result, n_profiles = 1, n_group_classes = 2),
               class = "multilpa_failed_candidate")
  # Every criterion the grid names must actually carry numbers. `aic` was
  # silently an all-missing column named `NA.`, so a finiteness assertion on
  # `grid$aic` passed vacuously against a column that did not exist.
  expect_true("aic" %in% names(grid))
  expect_true(all(is.finite(grid$aic[c(1, 2, 4)])))
  expect_true(all(c("kic", "clc_groups", "clc_individual") %in% names(grid)))
  one_profile <- candidate_fit(result, n_profiles = 1, n_group_classes = 1)
  expect_equal(grid$bic_groups[1], one_profile$bic)
  expect_equal(grid$bic_individual[1], one_profile$bic_individual)
  expect_true(is.na(grid$profile_entropy[1]))
  expect_gt(grid$profile_entropy[2], 0)
  expect_lte(grid$profile_entropy[2], 1)
  expect_error(enumerate_classes(d, "y", "g", n_profiles = 0), "n_profiles")
  # `n_profiles` names the grid here and the count everywhere else, so a single
  # value is a one-row grid rather than the error the old `profiles =` spelling
  # had to raise when the same word arrived through `...`.
  single <- enumerate_classes(d, "y", "g", n_profiles = 2, n_group_classes = 1,
                              n_starts = 2, seed = 1)
  expect_equal(nrow(as.data.frame(single)), 1L)
  # A name the grid does not own is forwarded to multilpa(), which refuses it.
  # The enumerator records every candidate's failure rather than raising, so the
  # refusal has to be visible in the table -- not silently dropped.
  forwarded <- enumerate_classes(d, "y", "g", n_profiles = 2, n_group_classes = 1,
                                 nonsense = 1, seed = 1)
  expect_match(as.data.frame(forwarded)$error, "unused argument")
})

test_that("simulation preserves group layout and full covariance moments", {
  set.seed(9)
  noise <- matrix(rnorm(2000), 1000, 2) %*% chol(matrix(c(1, .6, .6, 2), 2))
  d <- data.frame(g = rep(1:100, each = 10), a = noise[, 1], b = noise[, 2])
  model <- multilpa(d, c("a", "b"), "g", 1, 1, covariance_model = "full", n_starts = 1)
  simulated <- .multilpa_simulate(model)
  expect_identical(simulated$g, d$g)
  expect_identical(names(simulated), c("a", "b", "g"))
  expect_equal(cov(simulated[c("a", "b")]), cov(d[c("a", "b")]), tolerance = .2)
  expect_false(anyNA(simulated))
  one <- multilpa(d, "a", "g", 1, 1, covariance_model = "full", n_starts = 1)
  single_indicator <- .multilpa_simulate(one)
  expect_equal(dim(single_indicator), c(1000L, 2L))
  expect_identical(single_indicator$g, d$g)
  expect_equal(mean(single_indicator$a), mean(d$a), tolerance = .1)
})

test_that("bootstrap withholds p-values if a replicate cannot be fitted", {
  set.seed(53)
  d <- data.frame(g = rep(1:20, each = 5), y = c(rnorm(50, -2), rnorm(50, 2)))
  small <- multilpa(d, "y", "g", 1, 1, n_starts = 2, seed = 12)
  large <- multilpa(d, "y", "g", 2, 1, n_starts = 2, seed = 12)
  testthat::local_mocked_bindings(multilpa = function(...) stop("test optimization failure"))
  expect_warning(result <- bootstrap_lrt(small, large, d, iter = 2, seed = 1),
                 class = "multilpa_failed_replicates")
  expect_s3_class(result, "multilpa_bootstrap_lrt")
  ## Methods added in this sweep reach the generic only after the NAMESPACE is
  ## regenerated, so they are exercised by explicit call here. Dispatch itself
  ## is asserted in test-tidy-classes.R.
  test <- as.data.frame.multilpa_bootstrap_lrt(result)
  expect_identical(nrow(test), 1L)
  expect_true(is.na(test$p_value))
  expect_equal(test$n_valid, 0L)
  replicates <- as.data.frame.multilpa_bootstrap_lrt(result, what = "replicates")
  expect_true(all(replicates$error == "test optimization failure"))
  expect_error(plot.multilpa_bootstrap_lrt(result),
               class = "multilpa_nothing_to_plot")
})

test_that("bootstrap refits generated data and reports finite simulation correction", {
  set.seed(53)
  d <- data.frame(g = rep(1:30, each = 8), y = rnorm(240))
  small <- multilpa(d, "y", "g", 1, 1, variance_model = "equal", n_starts = 2, seed = 8)
  large <- multilpa(d, "y", "g", 2, 1, variance_model = "equal", n_starts = 3, seed = 8,
                      max_iter = 3000, tol = 1e-7)
  rng <- .Random.seed
  result <- bootstrap_lrt(small, large, d, iter = 3, n_starts = 3,
                                  max_iter = 3000, tol = 1e-7, seed = 42)
  expect_identical(.Random.seed, rng)
  test <- as.data.frame.multilpa_bootstrap_lrt(result)
  replicates <- as.data.frame.multilpa_bootstrap_lrt(result, what = "replicates")
  expect_identical(nrow(test), 1L)
  expect_equal(test$n_valid, 3L)
  expect_equal(test$statistic, max(0, 2 * (large$log_likelihood - small$log_likelihood)))
  expect_equal(test$p_value, (1 + sum(replicates$statistic >= test$statistic)) / 4)
  expect_equal(nrow(replicates), 3L)
  expect_true(all(is.na(replicates$error)))
  expect_true(is.finite(test$monte_carlo_se))
  expect_output(print.multilpa_bootstrap_lrt(result), "Parametric bootstrap")
  expect_output(print.summary_multilpa_bootstrap_lrt(
    summary.multilpa_bootstrap_lrt(result)), "simulated, not chi-square")
  expect_error(bootstrap_lrt(small, small, d, iter = 2), "one class")
  changed <- d
  changed$y <- d$y + 1
  expect_error(bootstrap_lrt(small, large, changed, iter = 2), "reproduce")
  changed$y[1] <- NA_real_
  expect_error(bootstrap_lrt(small, large, changed, iter = 2), "complete")
})

test_that("full covariance print and summary expose residual matrices", {
  set.seed(3)
  d <- data.frame(g = rep(1:10, each = 10), a = rnorm(100), b = rnorm(100))
  fit <- multilpa(d, c("a", "b"), "g", 1, 1, covariance_model = "full", n_starts = 1)
  expect_output(print(fit), "full residual covariance")
  expect_equal(summary(fit)$covariances, fit$covariances)
  expect_output(print(summary(fit)), "covariance matrices")
})
