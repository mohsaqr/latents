# Independent references and edge cases missed by the post-submission additions.
review_regression_data <- function() {
  withr::with_seed(943, {
    data <- expand.grid(time = 0:3, id = 1:32)
    data$school <- (data$id - 1L) %/% 4L + 1L
    data$off <- 4 + 2 * data$time
    data$y <- 3 + 0.6 * data$time + data$off +
      rep(stats::rnorm(32), each = 4) + stats::rnorm(nrow(data), sd = 0.7)
    data$binary <- stats::rbinom(nrow(data), 1, stats::plogis(-0.4 + 0.3 * data$time))
    data$weight <- rep(seq_len(32), each = 4)
    data
  })
}

test_that("one-class binary regression agrees with an ordinary logistic GLM", {
  data <- review_regression_data()
  reference <- stats::glm(binary ~ time, data, family = stats::binomial())
  fit <- mixture_regression(binary ~ time, data, n_classes = 1,
                            family = "binomial", n_starts = 0)
  expect_equal(fit$log_likelihood, as.numeric(logLik(reference)), tolerance = 1e-9)
  expect_equal(unname(fit$params$beta[, 1]), unname(stats::coef(reference)),
               tolerance = 1e-7)
  expect_equal(get_results(fit)$std_error, sqrt(unname(diag(stats::vcov(reference)))),
               tolerance = 1e-5)
})

test_that("count regression densities stay accurate for large counts and Poisson limits", {
  for (family in c("poisson", "negative_binomial")) {
    y <- if (family == "poisson") c(0, 1e6, 1e12, 1e15) else 0:15
    mu <- if (family == "poisson") pmax(y, 1) else rep(3, length(y))
    alpha <- latents:::.latents_min_dispersion
    spec <- list(y = y, x = matrix(1, length(y), 1), z = matrix(0, length(y), 0),
                 offset = log(mu), family = family, log_normalizer = -lgamma(y + 1))
    params <- list(beta = matrix(0, 1, 1), common = numeric(), sigma2 = alpha)
    expected <- if (family == "poisson") stats::dpois(y, mu, log = TRUE) else
      stats::dnbinom(y, mu = mu, size = 1 / alpha, log = TRUE)
    expect_equal(drop(latents:::.mixture_log_density(spec, params)), expected,
                 tolerance = 1e-12)
  }
})

test_that("robust regression covariance refuses too few independent clusters", {
  data <- review_regression_data()
  data$id <- 1L
  fit <- mixture_regression(y ~ time, data, n_classes = 1, id = "id",
                            n_starts = 0, vcov_type = "none")
  expect_error(vcov(fit, type = "robust"), class = "latents_too_few_groups")
  data$id <- rep(1:3, length.out = nrow(data))
  fit <- mixture_regression(y ~ time, data, n_classes = 1, id = "id",
                            n_starts = 0, vcov_type = "none")
  expect_error(vcov(fit, type = "robust"), class = "latents_too_few_groups")
})

test_that("multilevel growth generic BIC uses the independent clusters", {
  data <- review_regression_data()
  fit <- mixture_regression(y ~ time, data, 1, id = "id", class_level = "group",
                            random = "intercept", cluster = "school", n_starts = 0,
                            vcov_type = "none")
  expect_identical(nobs(fit), length(unique(data$school)))
  expect_identical(attr(logLik(fit), "nobs"), nobs(fit))
  expect_equal(stats::BIC(fit), get_results(fit, "fit")$bic)
})

test_that("enumerated growth and weighted models use the appropriate inference", {
  data <- review_regression_data()
  growth <- enumerate_regressions(y ~ time, data, n_classes = 1, id = "id",
                                  class_level = "group", random = "intercept",
                                  n_starts = 0)
  chosen <- get_results(growth, "model", n_classes = 1)
  expect_identical(names(chosen$inference$theta), names(coef(chosen)))
  expect_true(all(is.finite(get_results(chosen, "coefficients")$std_error)))
  weighted <- enumerate_regressions(y ~ time, data, n_classes = 1, id = "id",
                                    class_level = "group", weights = "weight",
                                    n_starts = 0)
  chosen <- get_results(weighted, "model", n_classes = 1)
  expect_identical(chosen$inference$vcov_type, "robust")
  expect_true(all(is.finite(get_results(chosen)$std_error)))
})

test_that("classification tables retain classes with no modal assignments", {
  posterior <- matrix(rep(c(0.7, 0.3), each = 10), 10, 2)
  fit <- list(spec = list(n_classes = 2L, nesting = "observation"),
              expectation = list(tau = posterior))
  table <- latents:::.mixture_classification_table(fit)
  expect_identical(nrow(table), 4L)
  expect_equal(table$mean_posterior[table$assigned == "class_1"], c(0.7, 0.3))
  expect_true(all(is.na(table$mean_posterior[table$assigned == "class_2"])))
  expect_identical(table$n_assigned, c(10L, 0L, 10L, 0L))
})

test_that("trajectory estimates include offsets and errors are on the response scale", {
  data <- review_regression_data()
  fit <- mixture_regression(y ~ time + offset(off), data, 1, id = "id",
                            class_level = "group", n_starts = 0)
  rows <- data[1:4, ]
  inference <- latents:::.mixture_resolve_inference(fit, NULL)
  band <- latents:::.trajectory_band(fit, rows, inference, 0.95)
  expect_equal(band$estimate, predict(fit, rows)$fitted, tolerance = 1e-12)
  logistic <- mixture_regression(binary ~ time, data, 1, family = "binomial",
                                 id = "id", class_level = "group", n_starts = 0)
  inference <- latents:::.mixture_resolve_inference(logistic, NULL)
  band <- latents:::.trajectory_band(logistic, rows, inference, 0.95)
  reference <- stats::glm(binary ~ time, data, family = stats::binomial())
  predicted <- stats::predict(reference, rows, type = "link", se.fit = TRUE)
  expected <- predicted$se.fit * stats::plogis(predicted$fit) *
    stats::plogis(predicted$fit, lower.tail = FALSE)
  expect_equal(band$std_error, unname(expected), tolerance = 1e-5)
})

test_that("growth residual sums survive a large response translation", {
  spec <- list(n = 8L, x = cbind(1, rep(0:3, 2)), z = matrix(0, 8, 0),
               y = 1e8 + c(0.1, 1.2, 1.9, 3.1, -0.2, 1, 2.2, 2.9),
               offset = rep(0, 8), random_design = matrix(1, 8, 1),
               group_index = rep(1:2, each = 4))
  beta <- c(1e8, 1)
  residual <- spec$y - drop(spec$x %*% beta)
  moments <- latents:::.growth_statistics(spec)
  computed <- latents:::.growth_residuals(moments, beta)
  expect_equal(computed$rr, unname(drop(rowsum(residual^2, spec$group_index))))
  expect_equal(drop(computed$zr), unname(drop(rowsum(residual, spec$group_index))))
  factor <- matrix(sqrt(0.5), 1, 1)
  terms <- latents:::.growth_class_terms(moments, beta, factor, 0.3)
  expected <- vapply(split(residual, spec$group_index), function(r) {
    covariance <- diag(0.3, length(r)) + 0.5
    -0.5 * (length(r) * log(2 * pi) +
              as.numeric(determinant(covariance, logarithm = TRUE)$modulus) +
              drop(crossprod(r, solve(covariance, r))))
  }, numeric(1))
  expect_equal(terms$log_density, unname(expected), tolerance = 1e-12)
})

test_that("a boundary coordinate affects only transforms that depend on it", {
  theta <- c(beta = 1, log_dispersion = log(1e-8))
  covariance <- matrix(c(4, NA, NA, NA), 2, 2)
  result <- latents:::.mixture_delta(theta, covariance, function(t) c(exp(t[1]), exp(t[2])))
  expect_equal(result$std_error[1], 2 * exp(1), tolerance = 1e-8)
  expect_true(is.na(result$std_error[2]))
})

test_that("failed and unconverged regression bootstrap fits withhold the p-value", {
  data <- data.frame(y = 1:5, x = 1:5)
  null <- list(log_likelihood = -10, converged = TRUE,
               spec = list(response_name = "y", model_data = data, formula = y ~ x))
  alternative <- null
  alternative$log_likelihood <- -8
  grid <- data.frame(n_classes = 1:2, n_group_classes = 1L)
  testthat::local_mocked_bindings(.mixture_draw = function(object) data$y,
                                .mixture_try_fit = function(...) errorCondition("refit failed"),
                                .package = "latents")
  expect_warning(result <- latents:::.mixture_blrt(
    list(null, alternative), grid, 2L, 3L, 0L, list()),
    class = "latents_failed_replicates")
  expect_true(is.na(result$blrt_p_value))
  expect_identical(result$blrt_replicates, 0L)
  expect_identical(result$blrt_flagged, 3L)
  testthat::local_mocked_bindings(.mixture_try_fit = function(...) {
    list(log_likelihood = -10, converged = FALSE)
  }, .package = "latents")
  expect_warning(result <- latents:::.mixture_blrt(
    list(null, alternative), grid, 2L, 3L, 0L, list()),
    class = "latents_failed_replicates")
  expect_true(is.na(result$blrt_p_value))
  testthat::local_mocked_bindings(.mixture_try_fit = function(formula, data, n, ...) {
    list(log_likelihood = -10 + 0.5 * (n - 1), converged = TRUE)
  }, .package = "latents")
  result <- latents:::.mixture_blrt(list(null, alternative), grid, 2L, 3L, 0L, list())
  expect_equal(result$blrt_p_value, 1 / 4)
  expect_identical(result$blrt_replicates, 3L)
  expect_identical(result$blrt_flagged, 0L)
})

test_that("alternative trajectory predictors print using their own column", {
  data <- review_regression_data()
  data$off <- data$id
  fit <- mixture_regression(y ~ time + off, data, 1, id = "id",
                            class_level = "group", random = "intercept", n_starts = 0)
  table <- get_results(fit, "trajectories", time = "off")
  expect_identical(names(table)[2L], "off")
  expect_true(any(grepl("off", capture.output(print(table)), fixed = TRUE)))
  individual <- get_results(fit, "individual", time = "off")
  expect_identical(names(individual)[2L], "off")
  expect_equal(individual$off, data$off)
  expect_s3_class(plot(fit, time = "off"), "ggplot")
  expect_s3_class(plot(fit, what = "individuals", time = "off"), "ggplot")
})

test_that("model comparisons reject different weights and likelihood scales", {
  data <- review_regression_data()
  fit <- mixture_regression(y ~ time, data, 1, n_starts = 0, vcov_type = "none")
  changed <- fit
  changed$spec$sampling_weights <- seq_len(nrow(data))
  expect_error(compare_models(a = fit, b = changed), class = "latents_incomparable_models")
  changed <- fit
  changed$spec$trials <- rep(10, nrow(data))
  expect_error(compare_models(a = fit, b = changed), class = "latents_incomparable_models")
})

test_that("nonfinite evaluated predictors and outcomes are refused explicitly", {
  data <- data.frame(y = 1:10, x = 0:9, id = rep(1:5, each = 2))
  expect_error(mixture_regression(y ~ log(x), data, 1), class = "latents_bad_data")
  for (family in c("gaussian", "binomial", "poisson", "ordinal")) {
    expect_error(latents:::.mixture_response(c(1, Inf), family), class = "latents_bad_data")
  }
  expect_error(mixture_regression(y ~ x, data, 1, membership = ~ log(x)),
               class = "latents_bad_data")
  data$y <- data$y + sin(data$x)
  fit <- mixture_regression(y ~ x, data, 1, n_starts = 0)
  data$x[1] <- Inf
  expect_error(predict(fit, data), class = "latents_bad_data")
})

test_that("profile prediction refuses matrix-valued indicator columns", {
  data <- data.frame(y = c(-2, -1, 0, 1, 2))
  fit <- lpa(data, "y", 1, n_starts = 1, seed = 1)
  data$y <- I(cbind(data$y, data$y + 10))
  expect_error(predict(fit, data), class = "latents_bad_data")
})

test_that("transition class counts reflect sampling weights", {
  fit <- list(n_profiles = 2L, group_probabilities = c(a = 0.6, b = 0.4),
              group_posteriors = matrix(c(0.9, 0.5, 0.1, 0.5), 2, 2),
              sampling_weights = c(0.5, 1.5))
  expect_equal(latents:::.lta_table(fit, "group_classes")$count, c(1.2, 0.8))
})

test_that("noise shares retain all cases and sequence plots preserve noise cells", {
  fit <- withr::with_seed(3, {
    data <- as.data.frame(rbind(matrix(stats::rnorm(120 * 2, -2), 120, 2),
                                matrix(stats::rnorm(120 * 2, 2), 120, 2),
                                matrix(stats::runif(30 * 2, -10, 10), 30, 2)))
    names(data) <- c("a", "b")
    data$id <- rep(seq_len(90), each = 3)
    data$time <- rep(1:3, 90)
    quietly(multilpa(data, c("a", "b"), "id", 2, 1, time = "time",
                     noise = TRUE, seed = 1, n_starts = 3))
  })
  expect_gt(sum(fit$subject_profiles == 0L), 0)
  key <- latents:::.gg_profile_key(fit)
  expect_equal(key$share, unname(colSums(fit$subject_posteriors)) / fit$n_observations)
  counts <- get_results(fit, "counts")
  counts <- counts[counts$level == "individuals", ]
  expect_equal(sum(counts$effective_count), fit$n_observations)
  expect_equal(counts$effective_count[counts$class == 0], sum(fit$noise_posteriors))
  expect_equal(latents:::.multilpa_count_frame(summary(fit)), get_results(fit, "counts"))
  sequences <- get_results(fit, "sequences", format = "wide")
  expect_equal(sum(as.matrix(sequences[-c(1, 2)]) == "0", na.rm = TRUE),
               sum(fit$subject_profiles == 0L))
  plot <- plot(fit, what = "sequences")
  expect_equal(sum(plot$data$profile_label == "Noise"), sum(fit$subject_profiles == 0L))
  expect_equal(nrow(plot$data), fit$n_observations)
  expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
})

test_that("observed growth means include person sampling weights", {
  spec <- list(model_data = data.frame(time = c(0, 1, 0, 1)),
               y = c(1, 2, 7, 8), family = "gaussian", n_groups = 2L,
               group_index = c(1L, 1L, 2L, 2L), sampling_weights = c(0.5, 1.5),
               n_classes = 1L, nesting = "group")
  fit <- structure(list(spec = spec, expectation = list(group_tau = matrix(1, 2, 1))),
                   class = "latents_mixture_regression")
  result <- latents:::.growth_observed_means(fit, "time", 0.95)
  expect_equal(result$observed_mean, c(5.5, 6.5))
})

test_that("NB dispersion derivatives stay accurate above the Poisson series switch", {
  # Independent 256-bit finite harmonic sums (digamma recurrence) at y = mu.
  reference <- c(-0.000499750333291542, -0.000499251165833542)
  result <- latents:::.latents_negative_binomial_dispersion_derivatives(
    500, 500, 2.001e-6)
  expect_equal(c(result$score, result$curvature), reference, tolerance = 1e-12)
  result <- latents:::.latents_negative_binomial_dispersion_derivatives(0, 1e20, 0.001)
  expect_equal(c(result$score, result$curvature),
               c(38143.9465808988, -37143.9465808988), tolerance = 1e-12)
})

test_that("random-effect spread varies only the chosen time predictor", {
  spec <- list(model_data = data.frame(time = c(0, 2), off = c(10, 20)),
               id = "id", random = ~ time + off, random_covariance = "equal",
               n_classes = 1L)
  fit <- list(spec = spec, params = list(random_covariance = list(diag(3))))
  bands <- data.frame(trajectory_time = c(10, 20), class = "class_1", estimate = 0)
  spread <- latents:::.growth_spread(fit, bands, "off", 0.8)
  expected_sd <- sqrt(1 + mean(spec$model_data$time)^2 + bands$trajectory_time^2)
  expect_equal(spread$spread_high, stats::qnorm(0.9) * expected_sd)
  expect_equal(spread$spread_low, -stats::qnorm(0.9) * expected_sd)
})
