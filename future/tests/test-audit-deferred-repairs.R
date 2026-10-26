test_that("deferred reporting adapters return fitted values and quadrature diagnostics", {
  set.seed(108)
  data <- data.frame(id = rep(1:12, each = 5),
                     y = rep(rnorm(12), each = 5) + rnorm(60))
  fit <- fit_random_intercept(data, "y", "id", 1, n_starts = 1, seed = 1)
  profiles <- get_results(fit)
  expect_equal(profiles$mean, unname(fit$means[1, 1]))
  expect_equal(profiles$variance, unname(fit$variances[1, 1]))
  expect_equal(nrow(get_results(fit, "random_intercepts")), 12L)
  expect_equal(nrow(get_results(fit, "assignments")), nrow(data))
  expect_equal(get_results(fit, "model")$log_likelihood, as.numeric(logLik(fit)))
  expect_output(print(fit), "Random-intercept LPA")
  expect_output(print(summary(fit), rows = 1), "Integration:")
  expect_identical(as.data.frame(summary(fit)), profiles)
  expect_error(get_results(fit, "unknown"), class = "latents_bad_argument")
  styled <- draw(plot(fit, labels = FALSE, symbols = 17, palette = "#123456"))
  points <- ggplot2::ggplot_build(styled)$data[[1L]]
  expect_true(all(points$shape == 17))
  expect_true(all(points$colour == "#123456"))
  labelled <- draw(plot(fit, labels = TRUE))
  expect_equal(length(labelled$layers), length(styled$layers) + 1L)
})

test_that("the deferred adjusted LR rejects different count measurements", {
  smaller <- structure(list(n_profiles = 1L, n_group_classes = 1L,
    n_observations = 10L, n_groups = 10L, n_parameters = 1L,
    log_likelihood = -100, extra_data = list(count_data = matrix(20, 10, 1))),
    class = "multilpa")
  larger <- smaller
  larger$n_profiles <- 2L
  larger$n_parameters <- 3L
  larger$log_likelihood <- -90
  larger$extra_data$count_data[,] <- 2
  expect_error(lmr_lrt(smaller, larger), class = "multilpa_incomparable_models")
  larger$extra_data <- smaller$extra_data
  result <- lmr_lrt(smaller, larger)
  expect_equal(result$statistic, 20)
  expect_true(is.na(result$p_value))
  expect_equal(result$adjustment_factor, 1 + 1 / (2 * log(10)))
})
