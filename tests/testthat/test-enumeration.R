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

test_that("simulation preserves categorical types and declared factor order", {
  set.seed(7)
  data <- data.frame(
    g = rep(seq_len(20L), each = 5L),
    y = stats::rnorm(100L),
    category = ordered(rep(c("10", "2"), 50L), levels = c("10", "2")),
    rating = rep(c(10L, 2L), each = 50L))
  fit <- quietly(multilpa(data, c("y", "category", "rating"), "g", 2L, 1L,
                           categorical = c("category", "rating"),
                           n_starts = 2L, max_iter = 80L, seed = 1L))
  simulated <- .multilpa_simulate(fit)
  expect_identical(class(simulated$category), c("ordered", "factor"))
  expect_identical(levels(simulated$category), c("10", "2"))
  expect_type(simulated$rating, "integer")
  expect_identical(names(simulated), c("y", "category", "rating", "g"))
  encoded <- .multilpa_encode_categorical(
    simulated[c("category", "rating")])
  expect_identical(encoded$levels$category, fit$categorical_levels$category)
  expect_identical(encoded$levels$rating, fit$categorical_levels$rating)
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
  test <- get_results(result, "test")
  expect_identical(nrow(test), 1L)
  expect_true(is.na(test$p_value))
  expect_equal(test$n_valid, 0L)
  replicates <- get_results(result, "replicates")
  expect_true(all(replicates$error == "test optimization failure"))
  expect_error(plot.multilpa_bootstrap_lrt(result),
               class = "multilpa_nothing_to_plot")
})

test_that("bootstrap refits generated data and reports finite simulation correction", {
  set.seed(53)
  d <- data.frame(g = rep(1:30, each = 8), y = rnorm(240))
  small <- multilpa(d, "y", "g", 1, 1, variance_model = "equal", n_starts = 2, seed = 8)
  large <- multilpa(d, "y", "g", 2, 1, variance_model = "equal", n_starts = 1, seed = 8,
                      max_iter = 3000, tol = 1e-7)
  rng <- .Random.seed
  result <- bootstrap_lrt(small, large, d, iter = 3, n_starts = 1,
                                  max_iter = 3000, tol = 1e-7, seed = 42)
  expect_identical(.Random.seed, rng)
  test <- get_results(result, "test")
  replicates <- get_results(result, "replicates")
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
  expect_output(print(summary(fit)), "-- covariances", fixed = TRUE)
})

test_that("boundary convergence noise is not mistaken for a reversed likelihood", {
  skip_on_cran()
  vars <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
  null_fit <- multilpa(engagement_small, vars, "student", n_profiles = 2,
                       n_group_classes = 1, n_starts = 3, seed = 1, tol = 1e-8)
  alt_fit <- multilpa(engagement_small, vars, "student", n_profiles = 2,
                      n_group_classes = 2, n_starts = 3, seed = 1, tol = 1e-8)
  # Under the null the alternative converges to the null solution, so every
  # replicate statistic sits at zero plus EM noise. EM stops on a RELATIVE
  # change, so that noise is about `2 * tol * |log_likelihood|`, which on a
  # likelihood of this size is several times the fixed -1e-5 window this check
  # used to apply. The result was
  # that healthy replicates were reported as "Nonconvergence or reversed
  # likelihood" and a perfectly ordinary comparison withheld its p-value.
  # A replicate refit on simulated data can land a small class; that warning is
  # recorded in the replicate table and also propagates, so muffle just it.
  # Eight replicates rather than the twelve this test used on the smaller
  # bundled dataset it was written against: each replicate now refits 1422
  # observations on five indicators. The claim is that EVERY replicate is
  # valid, and eight carry it as well as twelve did.
  result <- quietly(bootstrap_lrt(null_fit, alt_fit, iter = 8, n_starts = 2,
                                  max_iter = 2000, tol = 1e-8, seed = 7))
  test <- get_results(result, "test")
  expect_identical(test$n_valid, 8L)
  expect_false(is.na(test$p_value))

  replicates <- get_results(result, "replicates")
  expect_true(all(replicates$valid))
  expect_true(all(is.na(replicates$error)))
  # A replicate that raised nothing reports NA, not an empty string sitting
  # beside an NA-valued error column. Some replicates legitimately do warn, so
  # the contract is "never empty", not "always absent".
  expect_false(any(!is.na(replicates$warnings) & !nzchar(replicates$warnings)))

  # An unconverged model is refused before any statistic is formed, by class.
  # This guard fires ahead of the reversed-likelihood one, and correctly so: a
  # reversal means the alternative was badly optimised, which is what this
  # catches. That ordering makes `multilpa_reversed_likelihood` hard to reach
  # from here, so it is not asserted in this test.
  poor_start <- list(means = matrix(rep(c(0, 1), each = length(vars)),
                                    nrow = 2, byrow = TRUE),
                     variances = matrix(1, 2, length(vars)),
                     profile_probabilities = matrix(0.5, 2, 2),
                     group_probabilities = c(0.5, 0.5))
  unconverged <- quietly(multilpa(engagement_small, vars, "student",
                                  n_profiles = 2, n_group_classes = 2,
                                  n_starts = 1, seed = 1, max_iter = 0,
                                  start = poor_start, tol = 1e-8))
  expect_error(bootstrap_lrt(null_fit, unconverged, iter = 2, n_starts = 1,
                             seed = 1),
               class = "multilpa_no_converge")
})
