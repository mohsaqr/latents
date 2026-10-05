# The uniform noise component (`multilpa(noise = TRUE)`). The fourteen-structure
# comparison with mclust lives in equivalence/test-mclust-noise.R.

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

test_that("the likelihood is the mixture with a uniform component", {
  data <- noise_frame()
  fit <- noise_fit(data)
  x <- as.matrix(data)
  volume <- fit$hypervolume
  densities <- cbind(vapply(seq_len(2L), function(k) {
    fit$profile_probabilities[1L, k] * exp(rowSums(stats::dnorm(
      x, rep(fit$means[k, ], each = nrow(x)),
      rep(sqrt(fit$variances[k, ]), each = nrow(x)), log = TRUE)))
  }, numeric(nrow(x))), fit$noise_probability / volume)
  expect_equal(fit$log_likelihood, sum(log(rowSums(densities))), tolerance = 1e-10)
  # Posteriors over the Gaussian profiles and the noise sum to one, and so do
  # the proportions.
  expect_equal(unname(rowSums(fit$subject_posteriors)) + unname(fit$noise_posteriors),
               rep(1, nrow(x)), tolerance = 1e-12)
  expect_equal(sum(fit$profile_probabilities) + fit$noise_probability, 1,
               tolerance = 1e-12)
  # One proportion and the hypervolume: mclust's count.
  plain <- quietly(multilpa(data, c("a", "b", "c"), id = NULL, n_profiles = 2L,
                            n_starts = 1L, seed = 1L))
  expect_equal(fit$n_parameters, plain$n_parameters + 2)
  # The scattered points are the ones the component takes.
  expect_gte(mean(fit$subject_profiles[271:300] == 0L), 0.9)
})

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
                                n_group_classes = 2L, noise = TRUE)),
               class = "latents_unsupported_noise")
  expect_error(quietly(multilpa(data, c("a", "b", "c"), id = NULL, n_profiles = 2L,
                                noise = "yes")),
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
