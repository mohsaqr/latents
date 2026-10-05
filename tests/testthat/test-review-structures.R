# Every Gaussian structure, with weights, checked directly in probability space.
# Regression tests from past audits; CI runs them on every platform.
skip_on_cran()

review_gaussian_density <- function(x, mean, covariance) {
  residual <- sweep(x, 2L, mean, "-")
  exp(-.5 * rowSums((residual %*% solve(covariance)) * residual)) /
    sqrt((2 * pi)^ncol(x) * det(covariance))
}

test_that("all fourteen structures preserve weighted multilevel mixture likelihood", {
  set.seed(921)
  x <- matrix(rnorm(36), 12L, 3L)
  group <- rep(1:4, each = 3)
  weights <- c(.5, 1.5, .75, 1.25)
  responsibilities <- cbind(seq(.2, .8, length.out = 12), seq(.8, .2, length.out = 12))
  expectation <- list(subject_posteriors = responsibilities,
                      group_posteriors = matrix(c(.7, .3), 4L, 2L, byrow = TRUE),
                      joint = list(responsibilities * .7, responsibilities * .3))
  lapply(.multilpa_structures(), function(structure) {
    variance_model <- if (structure %in% c("EEI", "EEE")) "equal" else "varying"
    covariance_model <- if (.multilpa_is_ellipsoidal(structure)) "full" else "diagonal"
    parameters <- .multilpa_maximization(x, expectation, variance_model, 1e-10,
                                         covariance_model, structure = structure)
    parameters$profile_probabilities <- rbind(c(.8, .2), c(.25, .75))
    parameters$group_probabilities <- c(.65, .35)
    covariance <- lapply(1:2, function(k) {
      if (is.null(parameters$covariances)) diag(parameters$variances[k, ]) else
        parameters$covariances[, , k]
    })
    density <- vapply(1:2, function(k) {
      review_gaussian_density(x, parameters$means[k, ], covariance[[k]])
    }, numeric(12))
    group_likelihood <- vapply(1:4, function(j) {
      rows <- which(group == j)
      sum(vapply(1:2, function(h) {
        parameters$group_probabilities[h] *
          prod(drop(density[rows, , drop = FALSE] %*% parameters$profile_probabilities[h, ]))
      }, numeric(1)))
    }, numeric(1))
    actual <- .multilpa_expectation(x, group, parameters, weights = weights)
    expect_equal(actual$log_likelihood, sum(weights * log(group_likelihood)),
                 tolerance = 1e-12, info = structure)
    expect_equal(unname(rowSums(actual$subject_posteriors)), weights[group],
                 tolerance = 1e-12, info = structure)
    expect_equal(unname(rowSums(actual$group_posteriors)), weights,
                 tolerance = 1e-12, info = structure)
    # Integer unit replication is an independent M-step weighting reference.
    copies <- c(1L, 3L, 2L, 4L)
    row_copies <- copies[group]
    weighted <- .multilpa_weigh(expectation, copies, group)
    repeated <- expectation
    repeated$subject_posteriors <- responsibilities[rep(1:12, row_copies), ]
    repeated$group_posteriors <- expectation$group_posteriors[rep(1:4, copies), ]
    repeated$joint <- lapply(expectation$joint, function(block) block[rep(1:12, row_copies), ])
    a <- .multilpa_maximization(x, weighted, variance_model, 1e-10,
                                 covariance_model, structure = structure)
    b <- .multilpa_maximization(x[rep(1:12, row_copies), ], repeated,
                               variance_model, 1e-10, covariance_model, structure = structure)
    expect_equal(a$means, b$means, tolerance = 1e-12, info = structure)
    expect_equal(a$variances, b$variances, tolerance = 1e-10, info = structure)
    expect_equal(a$covariances, b$covariances, tolerance = 1e-10, info = structure)
  })
})
