# Review 2026-10-26: the shared EM kernel and the Gaussian measurement model.
# Calibrations against independent arithmetic (brute-force enumeration,
# mvtnorm, numDeriv, mclust, constrained optimization) and invariants of the
# offset stabilization added in 8ac3576.
# The comparisons with independent references are in tests/equivalence/test-review-1026-kernel.R.

kernel_review_mixed <- function() {
  set.seed(3)
  n <- 40L
  k <- 3L
  x <- matrix(stats::rnorm(n * 2L), n, 2L)
  x[c(3L, 7L), 1L] <- NA
  codes <- matrix(sample.int(3L, n * 2L, TRUE), n, 2L)
  codes[5L, 2L] <- NA
  counts <- matrix(stats::rpois(n * 2L, 3), n, 2L)
  counts[9L, 1L] <- NA
  parameters <- list(
    means = matrix(c(-1, 0, 1, 0.5, 0, -0.5), k, 2L),
    variances = matrix(c(1, 2, 0.5, 1, 1, 1.5), k, 2L),
    response_probabilities = list(
      matrix(c(.2, .5, .3, .3, .3, .4, .5, .2, .3), k, 3L, byrow = TRUE),
      matrix(c(.6, .2, .2, .1, .8, .1, .3, .3, .4), k, 3L, byrow = TRUE)),
    ordinal_intercepts = list(c(1, 0.2, -0.8)),
    ordinal_locations = matrix(c(-1, 0, 1.2), k, 1L),
    count_means = matrix(c(1, 3, 6, 2, 2.5, 4), k, 2L),
    count_dispersion = matrix(c(.3, .5, .2, .4, .1, .6), k, 2L),
    profile_probabilities = matrix(c(.2, .3, .5, .6, .3, .1), 2L, k, byrow = TRUE),
    group_probabilities = c(.4, .6))
  list(x = x, codes = codes, counts = counts, ordinal = matrix(sample.int(4L, n, TRUE), n, 1L),
       parameters = parameters, group_index = rep(seq_len(10L), each = 4L))
}

test_that("measurement offsets leave the nested likelihood and posteriors unchanged", {
  fixture <- kernel_review_mixed()
  invisible(lapply(c("poisson", "negative_binomial"), function(model) {
    extra <- list(ordinal = fixture$ordinal, count = fixture$counts, count_model = model)
    on <- .latents_measurement_log_density(fixture$x, fixture$parameters,
                                           fixture$codes, extra, offset = TRUE)
    off <- .latents_measurement_log_density(fixture$x, fixture$parameters,
                                            fixture$codes, extra, offset = FALSE)
    expect_equal(on$log_density + on$offset, off$log_density, tolerance = 1e-13)
    stabilized <- .multilpa_nested_structure(on$log_density, on$offset,
                                             fixture$group_index, fixture$parameters)
    plain <- .multilpa_nested_structure(off$log_density, off$offset,
                                        fixture$group_index, fixture$parameters)
    # Brute force: L_j = sum_h omega_h prod_{i in j} sum_k pi_hk f_k(y_i).
    density <- exp(off$log_density)
    group_likelihood <- vapply(seq_len(10L), function(group) {
      rows <- which(fixture$group_index == group)
      sum(vapply(seq_len(2L), function(h) {
        fixture$parameters$group_probabilities[h] *
          prod(density[rows, , drop = FALSE] %*% fixture$parameters$profile_probabilities[h, ])
      }, numeric(1)))
    }, numeric(1))
    expect_equal(stabilized$log_likelihood, sum(log(group_likelihood)), tolerance = 1e-12)
    expect_equal(stabilized$log_likelihood, plain$log_likelihood, tolerance = 1e-13)
    expect_equal(stabilized$subject_posteriors, plain$subject_posteriors, tolerance = 1e-12)
    expect_equal(stabilized$group_posteriors, plain$group_posteriors, tolerance = 1e-12)
  }))
})

test_that("an extreme Gaussian density does not round away priors and categorical evidence", {
  x <- matrix(1e8, 6L, 1L)
  parameters <- list(means = matrix(0, 2L, 1L), variances = matrix(1, 2L, 1L),
    response_probabilities = list(matrix(c(.9, .1, .2, .8), 2L, 2L, byrow = TRUE)),
    profile_probabilities = matrix(c(.3, .7), 1L), group_probabilities = 1)
  codes <- matrix(c(1L, 1L, 2L, 2L, 1L, 2L), 6L, 1L)
  measurement <- .latents_measurement_log_density(x, parameters, codes, offset = TRUE)
  result <- .multilpa_nested_structure(measurement$log_density, measurement$offset,
                                       seq_len(6L), parameters)
  # Both profiles have log f = -5e15; the posterior is prior times P(code).
  truth <- t(vapply(codes[, 1L], function(code) {
    weight <- c(.3, .7) * parameters$response_probabilities[[1L]][, code]
    weight / sum(weight)
  }, numeric(2)))
  expect_equal(result$subject_posteriors, truth, tolerance = 1e-14)
})

test_that("SQUAREM reaches plain EM's maximum with a monotone history", {
  skip_if_not_installed("MASS")
  set.seed(2)
  x <- rbind(MASS::mvrnorm(120L, c(0, 0), matrix(c(1, .6, .6, 1), 2L)),
             MASS::mvrnorm(80L, c(1.2, 1), diag(2) * 0.6),
             MASS::mvrnorm(60L, c(-1, 2), diag(c(.3, 1))))
  frame <- data.frame(a = x[, 1L], b = x[, 2L])
  fits <- lapply(c("squarem", "none"), function(acceleration) {
    suppressMessages(multilpa(frame, c("a", "b"), NULL, 3L, n_starts = 1L, seed = 1L,
      tol = 1e-10, max_iter = 5000L, covariance_model = "full",
      acceleration = acceleration))
  })
  expect_equal(fits[[1L]]$log_likelihood, fits[[2L]]$log_likelihood, tolerance = 1e-8)
  expect_lt(fits[[1L]]$iterations, fits[[2L]]$iterations)
  expect_gte(min(diff(fits[[1L]]$log_likelihood_history)), 0)
})

test_that("the profile likelihood is equivariant under rescaling and translation", {
  skip_if_not_installed("MASS")
  set.seed(2)
  x <- rbind(MASS::mvrnorm(120L, c(0, 0), matrix(c(1, .6, .6, 1), 2L)),
             MASS::mvrnorm(80L, c(2, 1), diag(2) * 0.6))
  fit <- function(values) {
    suppressMessages(multilpa(data.frame(a = values[, 1L], b = values[, 2L]),
      c("a", "b"), NULL, 2L, n_starts = 3L, seed = 4L, tol = 1e-12,
      volume = "varying", shape = "equal"))
  }
  base <- fit(x)
  scaled <- fit(x * 1000)
  moved <- fit(x + 1e6)
  # log L(c y) = log L(y) - n d log c.
  expect_equal(scaled$log_likelihood, base$log_likelihood - nrow(x) * 2 * log(1000),
               tolerance = 1e-12)
  expect_equal(moved$log_likelihood, base$log_likelihood, tolerance = 1e-10)
  # Relabelling the profiles permutes the posteriors and leaves the likelihood.
  parameters <- list(means = matrix(c(0, 2, 0, 1), 2L),
    variances = matrix(c(1, .6, 1, .6), 2L),
    covariances = array(c(1, .6, .6, 1, .6, 0, 0, .6), c(2L, 2L, 2L)),
    profile_probabilities = matrix(c(.6, .4), 1L), group_probabilities = 1)
  swapped <- parameters
  swapped$means <- parameters$means[2:1, ]
  swapped$variances <- parameters$variances[2:1, ]
  swapped$covariances <- parameters$covariances[, , 2:1]
  swapped$profile_probabilities <- parameters$profile_probabilities[, 2:1, drop = FALSE]
  original <- .multilpa_expectation(x, seq_len(nrow(x)), parameters)
  relabelled <- .multilpa_expectation(x, seq_len(nrow(x)), swapped)
  expect_equal(relabelled$log_likelihood, original$log_likelihood, tolerance = 1e-14)
  expect_equal(relabelled$subject_posteriors[, 2:1], original$subject_posteriors,
               tolerance = 1e-14)
})
