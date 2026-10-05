# Enumeration over group-class families and the bootstrap likelihood-ratio
# test between nested group-class fits.

compare_truth <- list(means = rbind(c(-1, -1), c(1, 1)),
                      between = rbind(c(0.25, 0.4), c(0.4, 0.25)),
                      within = c(1, 1), weights = c(0.5, 0.5))
compare_data <- additive_draw(compare_truth, rep(8L, 100L), seed = 21)
compare_vars <- c("y1", "y2")

test_that("enumeration crosses families, restrictions and class counts", {
  skip_on_cran()
  grid <- enumerate_classes(compare_data, compare_vars, "group",
                            family = c("additive", "dispersion",
                                       "additive_dispersion"),
                            n_group_classes = 1:2, seed = 1, n_starts = 3)
  expect_s3_class(grid, "latents_family_enumeration")
  table <- get_results(grid)
  # Two class counts x (additive x2 + dispersion x1 + additive-dispersion x2).
  expect_identical(nrow(table), 10L)
  expect_setequal(unique(table$model),
                  c("additive", "additive_equal", "dispersion",
                    "additive_dispersion", "additive_dispersion_equal"))
  # With one class every family is the same model.
  one <- subset(table, n_group_classes == 1)
  expect_equal(one$log_likelihood, rep(one$log_likelihood[1], nrow(one)),
               tolerance = 1e-6)
  # Each candidate's row reproduces its own fit.
  fit <- candidate_fit(grid, n_group_classes = 2, model = "additive")
  expect_s3_class(fit, "multilpa_additive")
  expect_equal(subset(table, n_group_classes == 2 & model == "additive")$bic,
               fit$bic)
  best <- get_results(grid, "best")
  expect_identical(nrow(best), 1L)
  expect_equal(best$bic, min(table$bic, na.rm = TRUE))
  expect_identical(as.data.frame(grid), table)
  expect_error(candidate_fit(grid, n_group_classes = 2),
               class = "latents_unknown_candidate")
  expect_error(enumerate_classes(compare_data, compare_vars, "group",
                                 family = "additive", model = "VVI"),
               class = "latents_bad_argument")
  expect_error(enumerate_classes(compare_data, compare_vars, "group",
                                 family = "additive", categorical = "y1"),
               class = "latents_bad_argument")
  skip_if_not_installed("ggplot2")
  expect_true(ggplot2::is_ggplot(plot(grid)))
})

test_that("the bootstrap LRT simulates the null and refits both models", {
  skip_on_cran()
  null <- multilpa(compare_data, compare_vars, "group", n_group_classes = 1,
                   family = "additive", seed = 1)
  alternative <- multilpa(compare_data, compare_vars, "group", n_group_classes = 2,
                          family = "additive", seed = 1)
  set.seed(5)
  before <- .Random.seed
  result <- bootstrap_lrt(null, alternative, iter = 9, n_starts = 2, seed = 3)
  expect_identical(.Random.seed, before)
  expect_s3_class(result, "multilpa_bootstrap_lrt")
  test <- get_results(result, "test")
  expect_identical(test$null_family, "additive")
  expect_identical(test$alternative_group_classes, 2L)
  expect_equal(test$statistic,
               2 * (alternative$log_likelihood - null$log_likelihood))
  expect_identical(test$n_valid, 9L)
  # Two well-separated classes: no null replicate reaches the observed value.
  expect_equal(test$p_value, 1 / 10)
  expect_output(print(result), "additive, 1 group class")
  again <- bootstrap_lrt(null, alternative, iter = 9, n_starts = 2, seed = 3)
  expect_identical(get_results(again, "replicates"), get_results(result, "replicates"))
})

test_that("simulation from a group-class fit has the fitted structure", {
  fit <- multilpa(compare_data, compare_vars, "group", n_group_classes = 2,
                  family = "additive", seed = 1)
  set.seed(9)
  simulated <- .additive_simulate(fit)
  expect_identical(names(simulated), c("group", "y1", "y2"))
  expect_identical(simulated$group, compare_data$group)
  # Large-sample check of the simulator against the fitted moments.
  set.seed(10)
  many <- do.call(rbind, lapply(1:20, \(i) .additive_simulate(fit)))
  expected_mean <- colSums(fit$group_probabilities * fit$means)
  # 2000 simulated groups: the Monte Carlo SE of the mean is about 0.026, so
  # four SEs bound the discrepancy.
  expect_lt(max(abs(colMeans(many[, compare_vars]) - expected_mean)), 0.1)
})

test_that("bootstrap refuses non-nested or mismatched group-class fits", {
  skip_on_cran()
  additive <- multilpa(compare_data, compare_vars, "group", n_group_classes = 2,
                       family = "additive", seed = 1)
  dispersion <- multilpa(compare_data, compare_vars, "group", n_group_classes = 2,
                         family = "dispersion", seed = 1)
  expect_error(bootstrap_lrt(additive, dispersion, iter = 3),
               class = "latents_bad_argument")
  expect_error(bootstrap_lrt(additive, additive, iter = 3),
               class = "latents_bad_argument")
  other <- multilpa(compare_data[compare_data$group != "g01", ], compare_vars,
                    "group", n_group_classes = 2, family = "additive_dispersion",
                    seed = 1)
  expect_error(bootstrap_lrt(additive, other, iter = 3),
               class = "latents_bad_inference_data")
  profile <- suppressMessages(lpa(compare_data, compare_vars, n_profiles = 2,
                                  n_starts = 2, seed = 1))
  expect_error(bootstrap_lrt(additive, profile, iter = 3),
               class = "latents_bad_argument")
})
