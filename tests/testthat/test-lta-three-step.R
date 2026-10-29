# r3step(initial_covariates =, transition_covariates =): the three-step latent
# Markov model. The Latent GOLD comparison is in
# equivalence/test-latentgold-step3-markov.R; here, an independent likelihood
# that sums over every path, equivariance, and the refusals.
# The comparisons with independent references are in tests/equivalence/test-lta-three-step.R.

# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

markov_fixture <- function(n = 200L, occasions = 4L, seed = 12L) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(seed)
  w <- stats::rnorm(n)
  z <- matrix(stats::rnorm(n * occasions), n)
  state <- matrix(NA_integer_, n, occasions)
  state[, 1L] <- ifelse(stats::runif(n) < stats::plogis(0.3 + 0.8 * w), 1L, 2L)
  # Sequential by construction: each occasion's state depends on the last.
  for (t in seq_len(occasions)[-1L]) {  # simulation only
    logit <- ifelse(state[, t - 1L] == 1L, -1.2 + 0.8 * z[, t], -1.0 - 0.6 * z[, t])
    moves <- stats::runif(n) < stats::plogis(logit)
    state[, t] <- ifelse(moves, 3L - state[, t - 1L], state[, t - 1L])
  }
  data <- data.frame(id = rep(seq_len(n), each = occasions),
                     time = rep(seq_len(occasions), n), w = rep(w, each = occasions),
                     z = as.vector(t(z)))
  s <- as.vector(t(state))
  data$y1 <- stats::rnorm(length(s), c(-1.2, 1.2)[s], 1)
  data$y2 <- stats::rnorm(length(s), c(-1, 1)[s], 1)
  data
}

markov_fit <- function(data, ...) {
  lta(data, c("y1", "y2"), "id", n_profiles = 2, time = "time", n_starts = 2,
      seed = 1, ...)
}

test_that("step-three slopes are equivariant in the covariates' units and origin", {
  data <- markov_fixture()
  fit <- markov_fit(data)
  base <- r3step(fit, data, initial_covariates = "w", transition_covariates = "z")
  shifted <- r3step(fit, transform(data, z = 5 + 3 * z, w = -2 + 0.5 * w),
                    initial_covariates = "w", transition_covariates = "z")
  slope <- base$term != "(Intercept)"
  scale <- ifelse(base$term == "z", 3, 0.5)
  expect_equal(shifted$estimate[slope], base$estimate[slope] / scale[slope], tolerance = 1e-6)
  expect_equal(shifted$standard_error[slope], base$standard_error[slope] / scale[slope],
               tolerance = 1e-5)
  expect_equal(attr(shifted, "log_likelihood"), attr(base, "log_likelihood"),
               tolerance = 1e-9)
  observed <- r3step(fit, data, transition_covariates = "z", vcov_type = "observed")
  expect_identical(attr(observed, "vcov_type"), "observed")
  expect_identical(unique(observed$block), c("initial", "transition"))
  expect_true(all(is.finite(observed$standard_error)))
})

test_that("the three-step transition model refuses what it does not define", {
  data <- markov_fixture(n = 120L)
  fit <- markov_fit(data)
  expect_error(r3step(fit, data, covariates = "w", transition_covariates = "z"),
               class = "latents_bad_argument")
  expect_error(r3step(fit, data, transition_covariates = "z", level = "groups"),
               class = "latents_bad_argument")
  expect_error(r3step(fit, data, transition_covariates = "y1"),
               class = "latents_bad_covariate")
  expect_error(r3step(fit, data, initial_covariates = "z"), class = "latents_bad_covariate")
  expect_error(r3step(fit, transform(data, k = 1), transition_covariates = "k"),
               class = "latents_bad_covariate")
  expect_error(r3step(fit, data, transition_covariates = "missing_column"))
  expect_error(r3step(fit, subset(data, id > 1), transition_covariates = "z"),
               class = "latents_bad_inference_data")
  with_covariates <- markov_fit(data, transition_covariates = "z")
  expect_error(r3step(with_covariates, data, transition_covariates = "z"),
               class = "latents_unsupported_three_step")
  profiles <- multilpa(data, c("y1", "y2"), "id", n_profiles = 2, n_group_classes = 1,
                       n_starts = 1, seed = 1)
  expect_error(r3step(profiles, data, transition_covariates = "z"),
               class = "latents_bad_argument")
})

test_that("two group classes are refused for the three-step chain", {
  skip_on_cran()
  data <- markov_fixture(n = 160L)
  fit <- markov_fit(data, n_group_classes = 2)
  expect_error(r3step(fit, data, transition_covariates = "z"),
               class = "latents_unsupported_three_step")
})
