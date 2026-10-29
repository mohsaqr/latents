# The uniform noise component (`multilpa(noise = TRUE)`). The fourteen-structure
# comparison with mclust lives in equivalence/test-mclust-noise.R.
# The comparisons with independent references are in tests/equivalence/test-noise.R.

noise_frame <- function(seed = 3L) {
  set.seed(seed)
  x <- rbind(matrix(stats::rnorm(135L * 3L, 2), 135L, 3L),
             matrix(stats::rnorm(135L * 3L, -2), 135L, 3L),
             matrix(stats::runif(30L * 3L, -10, 10), 30L, 3L))
  colnames(x) <- c("a", "b", "c")
  as.data.frame(x)
}

noise_fit <- function(data, ...) {
  quietly(multilpa(data, c("a", "b", "c"), id = NULL, n_profiles = 2L,
                   noise = TRUE, n_starts = 3L, seed = 1L, tol = 1e-10, ...))
}

test_that("the noise component is visible through the tidy accessors", {
  fit <- noise_fit(noise_frame())
  assignments <- get_results(fit, "assignments")
  expect_true("posterior_noise" %in% names(assignments))
  expect_true(all(assignments$profile %in% c(0L, 1L, 2L)))
  expect_identical(sum(assignments$profile == 0L), fit$n_noise)
  posteriors <- get_results(fit, "posteriors")
  expect_setequal(unique(posteriors$profile), c(0L, 1L, 2L))
  expect_equal(as.vector(tapply(posteriors$posterior, posteriors$row, sum)),
               rep(1, fit$n_observations), tolerance = 1e-12)
  probabilities <- get_results(fit, "profile_probabilities")
  expect_equal(sum(probabilities$probability), 1, tolerance = 1e-12)
  expect_true(0L %in% probabilities$profile)
  classification <- get_results(fit, "classification")
  expect_true(0L %in% classification$class)
  expect_output(print(fit), "Noise component")
})

test_that("unsupported combinations are refused by class", {
  data <- noise_frame()
  data$unit <- rep(seq_len(60L), each = 5L)
  expect_error(quietly(multilpa(data, c("a", "b", "c"), "unit", n_profiles = 2L,
                                n_group_classes = 2L, noise = TRUE, n_starts = 1)),
               class = "latents_unsupported_noise")
  expect_error(quietly(multilpa(data, c("a", "b", "c"), id = NULL, n_profiles = 2L,
                                noise = "yes", n_starts = 1)),
               class = "latents_bad_argument")
  fit <- noise_fit(noise_frame())
  expect_error(parameter_inference(fit), class = "latents_unsupported_noise")
  expect_error(starting_values(fit), class = "latents_unsupported_noise")
  start <- list(means = fit$means, variances = fit$variances,
                profile_probabilities = fit$profile_probabilities,
                group_probabilities = 1)
  expect_error(noise_fit(noise_frame(), start = start),
               class = "latents_bad_start")
})

test_that("a noise fit is translation invariant", {
  data <- noise_frame()
  moved <- data
  moved[c("a", "b", "c")] <- sweep(as.matrix(data), 2L, c(100, -50, 7), "+")
  original <- noise_fit(data)
  shifted <- noise_fit(moved)
  expect_equal(shifted$log_likelihood, original$log_likelihood, tolerance = 1e-10)
  expect_equal(shifted$noise_probability, original$noise_probability,
               tolerance = 1e-8)
  expect_equal(shifted$hypervolume, original$hypervolume, tolerance = 1e-10)
})
