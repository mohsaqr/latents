# Comparisons with independent references moved from tests/testthat/test-multilpa.R,
# so that they run in CI's equivalence job rather than under R CMD check.

multilpa_small_data <- data.frame(
  group = c("z", "z", "a", "a", "b", "b"),
  x = c(-2, -1, 0, 1, 2, 3),
  y = c(0, 2, 1, 4, 3, 5)
)

test_that("explicit starts give the independently calculated initial likelihood", {
  start <- list(
    means = matrix(c(-1, 1, 2, 4), nrow = 2L, byrow = TRUE),
    variances = matrix(c(1.5, 2, 2.5, 1), nrow = 2L, byrow = TRUE),
    profile_probabilities = matrix(c(0.8, 0.2, 0.3, 0.7),
                                   nrow = 2L, byrow = TRUE),
    group_probabilities = c(0.4, 0.6)
  )
  # Small, moderate densities allow direct probability arithmetic independent
  # of the fitter's log-sum-exp implementation.
  density <- cbind(
    dnorm(multilpa_small_data$x, -1, sqrt(1.5)) *
      dnorm(multilpa_small_data$y, 1, sqrt(2)),
    dnorm(multilpa_small_data$x, 2, sqrt(2.5)) *
      dnorm(multilpa_small_data$y, 4, 1)
  )
  conditional_density <- density %*% t(start$profile_probabilities)
  likelihood_by_group_type <- t(vapply(unique(multilpa_small_data$group),
    function(group) {
      apply(conditional_density[multilpa_small_data$group == group, , drop = FALSE],
            2L, prod)
    }, numeric(2)))
  expected_log_likelihood <- sum(log(as.numeric(
    likelihood_by_group_type %*% start$group_probabilities
  )))

  expect_warning(
    expect_warning(fit <- multilpa(multilpa_small_data, c("x", "y"), "group", 2L, 2L,
                     n_starts = 1L, max_iter = 1L, start = start, seed = 92),
                   class = "latents_thin_profiles"),
    "converg|iteration"
  )
  expect_equal(fit$log_likelihood_history[1L], expected_log_likelihood,
               tolerance = 1e-10)
  expect_false(fit$converged)
  expect_equal(fit$iterations, 1L)
  expect_length(fit$log_likelihood_history, 2L)
  expect_gte(fit$log_likelihood, expected_log_likelihood - 1e-10)
})
