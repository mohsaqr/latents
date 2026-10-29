# The comparisons with independent references are in tests/equivalence/test-review-0970-estimation.R.

test_that("joined categorical items estimate the full category grid", {
  data <- data.frame(a = rep(c(0, 0, 1), c(50, 25, 25)),
                     b = rep(c(0, 1, 0), c(50, 25, 25)))
  fit <- quietly(lca(data, c("a", "b"), 1, n_starts = 1,
                     local_dependence = "a:b"))
  expect_equal(fit$n_parameters, 3)
  expect_identical(fit$categorical_levels[["a:b"]], c("0:0", "0:1", "1:0", "1:1"))
  probability <- fit$response_probabilities[["a:b"]]
  expect_equal(sum(probability), 1, tolerance = 1e-12)
  expect_gt(unname(probability[1L, "1:1"]), 0)
  expect_equal(unname(probability[1L, "1:1"]), fit$min_probability, tolerance = 0)
  posterior <- predict(fit, newdata = data.frame(a = 1, b = 1), type = "posterior")
  expect_equal(posterior$probability_profile_1, 1)
  expect_equal(predict(fit, newdata = data, type = "posterior"),
               predict(fit, type = "posterior"))
  # Two classes of one joined four-category item are not identified (CORE-5).
  expect_warning(mixture <- quietly(lca(data, c("a", "b"), 2, n_starts = 1, seed = 1,
                                        local_dependence = "a:b")),
                 class = "latents_unidentified")
  expect_equal(mixture$n_parameters, 2 * 3 + 1)
  expect_equal(unname(mixture$response_probabilities[["a:b"]][, "1:1"]),
               rep(mixture$min_probability, 2), tolerance = 0)
  posterior <- predict(mixture, newdata = data.frame(a = 1, b = 1), type = "posterior")
  expect_equal(posterior$probability_profile_1 + posterior$probability_profile_2, 1)
})

test_that("extreme homogeneous transitions fall back to log-space messages", {
  layout <- .multilpa_sequence_layout(c(1L, 1L), 1:2, 1L)
  x <- matrix(c(0, 40), ncol = 1)
  parameters <- list(
    means = matrix(c(0, 40), 2), variances = matrix(1, 2),
    initial_probabilities = matrix(c(0.5, 0.5), 1),
    transition_probabilities = array(c(1, 1e-320, 1e-320, 1), c(2, 2, 1)),
    group_probabilities = 1)
  # The second matrix leaves the backward messages finite, but overflows
  # the backward-to-scale ratio used for transition counts.
  transitions <- list(matrix(c(1, 1e-320, 1e-320, 1), 2),
                       matrix(c(1, 1, 1e-320, 1e-320), 2))
  invisible(lapply(transitions, function(transition) {
    parameters$transition_probabilities <- array(transition, c(2, 2, 1))
    density <- .latents_measurement_log_density(x, parameters, offset = FALSE)
    emission <- .multilpa_sequence_emission(density$log_density, layout)
    reference <- .latents_forward_backward(
      emission$log_density, layout, log(parameters$initial_probabilities),
      list(array(log(transition), c(1, 2, 2))))
    moments <- .multilpa_sequence_moments(reference, emission$log_density, layout,
                                          transition, weights = 1, n_observations = 2)
    expectation <- .multilpa_transition_expectation(x, layout, parameters)
    expect_equal(expectation$log_likelihood,
                 sum(reference$log_scaled + emission$offset), tolerance = 1e-12)
    expect_equal(expectation$subject_posteriors, moments$posterior, tolerance = 1e-12)
    expect_equal(expectation$sequence[[1L]]$transition, moments$transition,
                 tolerance = 1e-12)
    expect_true(all(is.finite(unlist(expectation))))
    expect_equal(rowSums(expectation$subject_posteriors), c(1, 1), tolerance = 1e-12)
    expect_equal(sum(expectation$sequence[[1L]]$transition), 1, tolerance = 1e-12)
  }))
})

test_that("membership-covariate LMR comparisons retain a classed refusal", {
  set.seed(21)
  n <- 400
  z <- stats::rnorm(n)
  class <- stats::rbinom(n, 1, stats::plogis(z))
  data <- data.frame(a = stats::rnorm(n, 3 * class),
                     b = stats::rnorm(n, 3 * class), z = z)
  grid <- quietly(enumerate_lpa(
    data, c("a", "b"), n_profiles = 2:3, model = "EEI",
    profile_covariates = "z", n_starts = 2, seed = 1, tol = 1e-10))
  expect_error(.latents_vuong_parts(get_candidate(grid, 2)),
               class = "latents_invalid_lmr")
  result <- get_results(grid, "lmr")
  expect_equal(nrow(result), 1L)
  expect_match(result$note, "per-unit scores and observed information")
  expect_true(all(is.na(result[c("lr", "vlmr_p", "adjusted_lr", "almr_p",
                                 "reference_mean", "reference_sd")])))
})
