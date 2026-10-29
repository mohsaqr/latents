# Calibration and invariance checks of the two-level (multilevel) profile
# engine from the 2026-10-26 review: the nested E-step against brute-force
# enumeration, the closed-form M-step, the single-level limit, the sampling
# weights, the membership-logit M-step and the covariate scores and
# information against numerical derivatives (Vermunt 2003, 2008).
# The comparisons with independent references are in tests/equivalence/test-review-1026-multilevel.R.

# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

review_ml_data <- function(seed = 11, n_groups = 12, size = 4) {
  withr::local_seed(seed)
  g <- rep(seq_len(n_groups), each = size)
  h <- rep(sample(1:2, n_groups, TRUE), each = size)
  z <- ifelse(stats::runif(n_groups * size) < ifelse(h == 1, 0.8, 0.2), 1, 2)
  data.frame(
    g = g,
    y1 = stats::rnorm(n_groups * size, c(0, 2.5)[z]),
    y2 = stats::rnorm(n_groups * size, c(0, 1.5)[z], c(1, 0.6)[z]),
    wt = rep(sample(1:3, n_groups, TRUE), each = size)
  )
}

review_ml_covariate_data <- function(seed = 3, n_groups = 20, size = 6) {
  withr::local_seed(seed)
  g <- rep(seq_len(n_groups), each = size)
  gz <- rep(stats::rnorm(n_groups), each = size)
  hc <- rep(ifelse(stats::runif(n_groups) < stats::plogis(1.2 * gz[!duplicated(g)]),
                   2, 1), each = size)
  x <- stats::rnorm(n_groups * size)
  z <- ifelse(stats::runif(n_groups * size) < stats::plogis(-1 + 2 * (hc == 2) + 0.8 * x),
              2, 1)
  data.frame(
    g = g, x = x, gz = gz,
    wt = rep(stats::runif(n_groups, 0.5, 2), each = size),
    y1 = stats::rnorm(n_groups * size, c(-1.5, 1.5)[z], 0.8),
    y2 = stats::rnorm(n_groups * size, c(-1, 1)[z], c(0.7, 1)[z])
  )
}

review_ml_quiet <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (inherits(w, c("latents_boundary", "latents_small_classes",
                      "latents_thin_profiles"))) {
      invokeRestart("muffleWarning")
    }
  })
}

test_that("a converged fit is a fixed point of the closed-form M-step", {
  data <- review_ml_data()
  fit <- review_ml_quiet(multilpa(data, c("y1", "y2"), "g", n_profiles = 2,
                                  n_group_classes = 2, n_starts = 3, seed = 1,
                                  tol = 1e-12))
  y <- as.matrix(data[c("y1", "y2")])
  subject <- unname(fit$subject_posteriors)
  group <- unname(fit$group_posteriors)
  means <- crossprod(subject, y) / colSums(subject)
  variances <- t(vapply(1:2, function(k) {
    colSums(subject[, k] * sweep(y, 2L, means[k, ])^2) / sum(subject[, k])
  }, numeric(2)))
  expect_equal(unname(as.matrix(fit$means)), unname(means), tolerance = 1e-7)
  expect_equal(unname(as.matrix(fit$variances)), unname(variances), tolerance = 1e-7)
  expect_equal(unname(fit$group_probabilities), colMeans(group), tolerance = 1e-7)
  expect_true(fit$converged)
})

test_that("the E-step is invariant to relabelling profiles and group classes", {
  data <- review_ml_data()
  fit <- review_ml_quiet(multilpa(data, c("y1", "y2"), "g", n_profiles = 2,
                                  n_group_classes = 2, n_starts = 2, seed = 1))
  y <- as.matrix(data[c("y1", "y2")])
  parameters <- list(means = unname(as.matrix(fit$means)),
                     variances = unname(as.matrix(fit$variances)),
                     profile_probabilities = unname(as.matrix(fit$profile_probabilities)),
                     group_probabilities = unname(fit$group_probabilities))
  relabelled <- list(means = parameters$means[2:1, ],
                     variances = parameters$variances[2:1, ],
                     profile_probabilities = parameters$profile_probabilities[2:1, 2:1],
                     group_probabilities = parameters$group_probabilities[2:1])
  original <- .multilpa_expectation(y, fit$group_index, parameters)
  permuted <- .multilpa_expectation(y, fit$group_index, relabelled)
  expect_equal(permuted$log_likelihood, original$log_likelihood, tolerance = 1e-13)
  expect_equal(permuted$group_posteriors[, 2:1], original$group_posteriors,
               tolerance = 1e-13)
  expect_equal(permuted$subject_posteriors[, 2:1], original$subject_posteriors,
               tolerance = 1e-13)
})

test_that("integer group weights equal replicated groups and weights are scale free", {
  data <- review_ml_data(seed = 21, n_groups = 15, size = 5)
  weighted <- review_ml_quiet(multilpa(data, c("y1", "y2"), "g", n_profiles = 2,
                                       n_group_classes = 2, weights = "wt",
                                       n_starts = 3, seed = 1, tol = 1e-12))
  replicated <- do.call(rbind, lapply(seq_len(15), function(j) {
    block <- data[data$g == j, ]
    do.call(rbind, lapply(seq_len(block$wt[1]), function(r) {
      within(block, g <- paste(j, r))
    }))
  }))
  copies <- review_ml_quiet(multilpa(replicated, c("y1", "y2"), "g", n_profiles = 2,
                                     n_group_classes = 2, n_starts = 3, seed = 1,
                                     tol = 1e-12))
  ordered <- function(fit, field) {
    unname(as.matrix(fit[[field]]))[order(fit$means[, "y1"]), ]
  }
  expect_equal(ordered(weighted, "means"), ordered(copies, "means"), tolerance = 1e-6)
  expect_equal(ordered(weighted, "variances"), ordered(copies, "variances"),
               tolerance = 1e-6)
  # The weights are normalized to sum to the number of groups (as Mplus does).
  scale <- sum(data$wt[!duplicated(data$g)]) / 15
  expect_equal(copies$log_likelihood, scale * weighted$log_likelihood, tolerance = 1e-8)
  rescaled <- within(data, wt <- wt * 1e6)
  big <- review_ml_quiet(multilpa(rescaled, c("y1", "y2"), "g", n_profiles = 2,
                                  n_group_classes = 2, weights = "wt",
                                  n_starts = 3, seed = 1, tol = 1e-12))
  expect_equal(big$log_likelihood, weighted$log_likelihood, tolerance = 1e-12)
  expect_equal(big$means, weighted$means, tolerance = 1e-12)
})

test_that("the parameter count and BIC use the group sample size", {
  data <- review_ml_data()
  fit <- review_ml_quiet(multilpa(data, c("y1", "y2"), "g", n_profiles = 2,
                                  n_group_classes = 2, n_starts = 2, seed = 1))
  # (H - 1) + H (K - 1) mixing, K means and K variances per indicator.
  expect_identical(as.integer(fit$n_parameters), 1L + 2L + 4L + 4L)
  expect_identical(nobs(fit), 12L)
  expect_equal(BIC(fit), fit$bic, tolerance = 1e-12)
  expect_equal(fit$bic, -2 * fit$log_likelihood + log(12) * fit$n_parameters,
               tolerance = 1e-12)
  expect_equal(fit$bic_individual,
               -2 * fit$log_likelihood + log(48) * fit$n_parameters, tolerance = 1e-12)
})

test_that("an all-missing group changes neither the fit nor the group sample size", {
  data <- review_ml_covariate_data()
  empty <- data.frame(g = 99, x = c(0.3, -0.2), gz = 0.1, wt = 1,
                      y1 = NA_real_, y2 = NA_real_)
  padded <- rbind(data, empty)
  fit_padded <- review_ml_quiet(multilpa(padded, c("y1", "y2"), "g", 2, 2,
                                         missing = "fiml", n_starts = 2, seed = 1,
                                         tol = 1e-12))
  fit <- review_ml_quiet(multilpa(data, c("y1", "y2"), "g", 2, 2, missing = "fiml",
                                  n_starts = 2, seed = 1, tol = 1e-12))
  expect_equal(fit_padded$log_likelihood, fit$log_likelihood, tolerance = 1e-9)
  expect_identical(nobs(fit_padded), nobs(fit))
  covariate_padded <- review_ml_quiet(multilpa(
    padded, c("y1", "y2"), "g", 2, 2, profile_covariates = "x",
    group_covariates = "gz", missing = "fiml", n_starts = 2, seed = 1, tol = 1e-10))
  expect_identical(nobs(covariate_padded), 20L)
  expect_equal(BIC(covariate_padded), covariate_padded$bic, tolerance = 1e-12)
  kept <- review_ml_quiet(multilpa(
    padded, c("y1", "y2"), "g", 2, 2, profile_covariates = "x",
    group_covariates = "gz", missing = "fiml", n_starts = 2, seed = 1, tol = 1e-10,
    drop_empty = FALSE))
  expect_identical(nobs(kept), 21L)
  expect_equal(BIC(kept), kept$bic, tolerance = 1e-12)
})

test_that("rows without expected mass do not move the membership step", {
  withr::local_seed(4)
  n <- 200
  design <- cbind(1, stats::rnorm(n), stats::rnorm(n, 5, 0.3))
  counts <- matrix(stats::runif(n * 3), n, 3)
  solved <- .multilpa_weighted_logits(design, counts, matrix(0, 3, 2))
  # Zero-mass rows with extreme covariates used to dominate the coordinate
  # basis, so the step stopped far from the maximum (coefficients off by 2.5).
  padded <- .multilpa_weighted_logits(
    rbind(design, cbind(1, c(1e8, -1e8, 3e7), c(1e9, 2e9, -1e9))),
    rbind(counts, matrix(0, 3, 3)), matrix(0, 3, 2))
  expect_true(padded$converged)
  expect_equal(padded$coefficients, solved$coefficients, tolerance = 1e-12)
  # The same rows carry no weight in the separation verdict either.
  coefficients <- matrix(c(0.2, 0.5), 2, 1)
  profile_design <- list(cbind(1, c(stats::rnorm(20), 1e7)))
  group_design <- matrix(1, 21, 1)
  weights <- c(rep(1, 20), 0)
  expect_true(.multilpa_cov_extreme(coefficients, matrix(0, 1, 1),
                                    profile_design, group_design))
  expect_false(.multilpa_cov_extreme(coefficients, matrix(0, 1, 1),
                                     profile_design, group_design,
                                     profile_weights = weights,
                                     group_weights = weights))
})

test_that("the covariate structure is invariant to relabelling the group classes", {
  data <- review_ml_covariate_data()
  fit <- review_ml_quiet(multilpa(data, c("y1", "y2"), "g", 2, 2,
                                  profile_covariates = "x", group_covariates = "gz",
                                  n_starts = 2, seed = 1))
  x <- sweep(as.matrix(data[c("y1", "y2")]), 2L, fit$center, "-")
  parameters <- list(means = sweep(unname(fit$means), 2L, fit$center, "-"),
                     variances = unname(fit$variances))
  beta <- unname(fit$profile_coefficients)
  # Swapping the two group classes swaps their profile intercepts and design
  # blocks and flips the sign of the group logit (the reference changes).
  swapped_design <- lapply(fit$profile_design[2:1], function(design) {
    design[, c(2L, 1L, 3L), drop = FALSE]
  })
  original <- .multilpa_cov_expectation(x, fit$group_index, parameters,
                                        fit$profile_design, fit$group_design,
                                        beta, unname(fit$group_coefficients))
  relabelled <- .multilpa_cov_expectation(x, fit$group_index, parameters,
                                          swapped_design, fit$group_design,
                                          beta[c(2L, 1L, 3L), , drop = FALSE],
                                          -unname(fit$group_coefficients))
  expect_equal(relabelled$log_likelihood, original$log_likelihood, tolerance = 1e-12)
  expect_equal(relabelled$group_posteriors[, 2:1], original$group_posteriors,
               tolerance = 1e-12)
  expect_equal(relabelled$subject_posteriors, original$subject_posteriors,
               tolerance = 1e-12)
})

test_that("seeded two-level fits are reproducible and leave the caller's stream alone", {
  data <- review_ml_covariate_data()
  withr::local_seed(77)
  before <- .Random.seed
  first <- review_ml_quiet(multilpa(data, c("y1", "y2"), "g", 2, 2, n_starts = 4,
                                    seed = 5, search = "emEM"))
  expect_identical(.Random.seed, before)
  second <- review_ml_quiet(multilpa(data, c("y1", "y2"), "g", 2, 2, n_starts = 4,
                                     seed = 5, search = "emEM"))
  expect_identical(first$starts, second$starts)
  expect_identical(first$means, second$means)
  covariate <- review_ml_quiet(multilpa(data, c("y1", "y2"), "g", 2, 2,
                                        profile_covariates = "x", n_starts = 3,
                                        seed = 5))
  again <- review_ml_quiet(multilpa(data, c("y1", "y2"), "g", 2, 2,
                                    profile_covariates = "x", n_starts = 3, seed = 5))
  expect_identical(covariate$profile_coefficients, again$profile_coefficients)
})
