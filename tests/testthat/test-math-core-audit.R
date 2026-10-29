# Regression tests from past audits; CI runs them on every platform.
# The comparisons with independent references are in tests/equivalence/test-math-core-audit.R.
skip_on_cran()

test_that("extreme common Gaussian densities retain profile and group priors", {
  parameters <- list(means = matrix(0, 2L, 1L), variances = matrix(1, 2L, 1L),
    profile_probabilities = rbind(c(.2, .8), c(.7, .3)),
    group_probabilities = c(.3, .7))
  x <- matrix(c(1e8, -1e8, 2e8), 3L, 1L)
  groups <- c(1L, 2L, 1L)
  expectation <- .multilpa_expectation(x, groups, parameters)
  # Identical measurement distributions supply no class information, regardless
  # of how surprising the observations are under those distributions.
  expect_equal(expectation$group_posteriors,
               matrix(rep(c(.3, .7), each = 2L), 2L), tolerance = 1e-14)
  expect_equal(expectation$subject_posteriors,
               matrix(rep(c(.55, .45), each = 3L), 3L), tolerance = 1e-14)
  expect_equal(rowSums(expectation$subject_posteriors), rep(1, 3L))
  expect_equal(rowSums(expectation$group_posteriors), rep(1, 2L))
  expect_equal(expectation$log_likelihood, sum(dnorm(x, log = TRUE)))

  # A categorical observation must still provide information when the Gaussian
  # block is extremely unlikely. Its reference is the categorical-only model.
  codes <- matrix(c(1L, 2L, 1L), 3L, 1L)
  parameters$response_probabilities <- list(rbind(c(.8, .2), c(.1, .9)))
  mixed <- .multilpa_expectation(x, groups, parameters, codes)
  categorical <- parameters
  categorical$means <- categorical$variances <- matrix(numeric(0), 2L, 0L)
  reference <- .multilpa_expectation(matrix(numeric(0), 3L, 0L), groups,
                                    categorical, codes)
  expect_equal(mixed$subject_posteriors, reference$subject_posteriors,
               tolerance = 1e-14)
  expect_equal(mixed$group_posteriors, reference$group_posteriors,
               tolerance = 1e-14)
})

test_that("all-categorical starts have no Gaussian covariance requirement", {
  data <- data.frame(group = rep(1:4, each = 3), item = rep(c(1, 2, 2), 4))
  start <- list(response_probabilities = list(matrix(c(1 / 3, 2 / 3), 1L)),
                profile_probabilities = matrix(1), group_probabilities = 1)
  fits <- lapply(c("diagonal", "full"), function(covariance) {
    multilpa(data, "item", "group", 1, 1, categorical = "item",
             covariance_model = covariance, start = start, n_starts = 1)
  })
  expect_equal(fits[[1]]$log_likelihood, 4 * log(1 / 3) + 8 * log(2 / 3))
  expect_equal(fits[[2]]$log_likelihood, fits[[1]]$log_likelihood)
  expect_equal(fits[[2]]$response_probabilities, fits[[1]]$response_probabilities)
  expect_null(fits[[2]]$covariances)
  expect_equal(fits[[2]]$n_parameters, 1)
})

test_that("time uniqueness respects distinct native numeric group identifiers", {
  data <- data.frame(group = rep(c(1, 1 + 1e-15), each = 3),
                     time = rep(1:3, 2), score = c(-2, -1, -3, 2, 1, 3))
  expect_equal(length(unique(data$group)), 2L)
  fit <- multilpa(data, "score", "group", 1, 1, n_starts = 1, time = "time")
  expect_equal(fit$n_groups, 2L)
  expect_equal(fit$group_index, rep(1:2, each = 3))
  expect_equal(fit$time_values, data$time)
  data$time[2] <- data$time[1]
  expect_error(multilpa(data, "score", "group", 1, 1, n_starts = 1, time = "time"),
               class = "latents_bad_time")
})
