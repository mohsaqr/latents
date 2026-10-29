# Comparisons with independent references moved from tests/testthat/test-review-1026-multilevel-r2.R,
# so that they run in CI's equivalence job rather than under R CMD check.

r2_two_level <- function(n_groups = 24L, n_rows = 6L, seed = 3L) {
  withr::with_seed(seed, {
    group_class <- sample(1:2, n_groups, TRUE)
    shares <- rbind(c(0.8, 0.2), c(0.25, 0.75))
    do.call(rbind, lapply(seq_len(n_groups), function(group) {
      profile <- sample(1:2, n_rows, TRUE, prob = shares[group_class[group], ])
      centre <- c(-1.5, 1.5)[profile]
      data.frame(g = group, y1 = stats::rnorm(n_rows, centre, 1),
                 y2 = stats::rnorm(n_rows, centre * 0.8, 1.2),
                 z = stats::rnorm(n_rows), wg = stats::rnorm(1))
    }))
  })
}

# The data, and the same data with five extra groups of weight zero whose
# values are moved far away: a zero-weight unit must change nothing.
r2_padded <- function() {
  data <- r2_two_level()
  weights <- withr::with_seed(9, sample(1:3, 24L, TRUE))
  data$w <- weights[data$g]
  extra <- subset(data, g <= 5)
  extra$g <- extra$g + 100L
  extra$w <- 0
  extra$y1 <- extra$y1 * 4 + 6
  extra$y2 <- extra$y2 * 3 - 5
  extra$z <- extra$z * 5
  list(data = data, padded = rbind(data, extra))
}

r2_quiet <- function(expr) {
  withCallingHandlers(expr, warning = function(condition) {
    if (inherits(condition, c("latents_boundary", "latents_thin_profiles",
                              "latents_small_classes", "latents_unconverged",
                              "latents_extreme_coefficients"))) {
      invokeRestart("muffleWarning")
    }
  })
}

test_that("restricted cross-level pseudo likelihood ignores zero-weight groups (CL-1)", {
  sets <- r2_padded()
  fit <- function(data) r2_quiet(multilpa(
    data, c("y1", "y2"), "g", 2, 2, family = "restricted_cross_level",
    weights = "w", n_starts = 3, seed = 1, tol = 1e-12))
  plain <- fit(sets$data)
  padded <- fit(sets$padded)
  expect_equal(padded$log_likelihood, plain$log_likelihood, tolerance = 1e-10)
  # The reported likelihood is sum_j w_j log L_j at the reported estimates.
  data <- sets$data
  index <- match(data$g, unique(data$g))
  weights <- .latents_sampling_weights(data, "w", index, 24L)
  x <- as.matrix(data[c("y1", "y2")])
  row_density <- vapply(1:2, function(k) {
    plain$profile_prevalences[k] *
      stats::dnorm(x[, 1], plain$profile_means[k, 1], sqrt(plain$profile_variances[k, 1])) *
      stats::dnorm(x[, 2], plain$profile_means[k, 2], sqrt(plain$profile_variances[k, 2]))
  }, numeric(nrow(x)))
  group_means <- rowsum(x, index) / as.vector(table(index))
  group_density <- vapply(1:2, function(h) {
    plain$group_probabilities[h] *
      stats::dnorm(group_means[, 1], plain$group_means[h, 1],
                   sqrt(plain$group_mean_variances[h, 1])) *
      stats::dnorm(group_means[, 2], plain$group_means[h, 2],
                   sqrt(plain$group_mean_variances[h, 2]))
  }, numeric(24L))
  direct <- sum(weights * (log(rowSums(group_density)) +
                             as.vector(rowsum(log(rowSums(row_density)), index))))
  expect_equal(plain$log_likelihood, direct, tolerance = 1e-10)
})
