skip_on_cran()

# Audit 2026-11-01, slice `inference`: regression tests for the fixes made in
# R/weights.R, R/inference-service.R and R/fixed-measurement.R, and the
# independent numerical checks the audit ran (validation/audit-2026-11-01/
# inference.md).

audit_inference_data <- function(seed = 11L, n_groups = 40L, size = 6L) {
  set.seed(seed)
  school <- rep(seq_len(n_groups), each = size)
  high <- rep(stats::rbinom(n_groups, 1L, 0.5), each = size)
  profile <- ifelse(stats::runif(n_groups * size) < stats::plogis(-0.5 + 1.5 * high), 2L, 1L)
  data <- data.frame(
    school = school,
    a = stats::rnorm(n_groups * size, ifelse(profile == 2L, 2, -1), 1),
    b = stats::rnorm(n_groups * size, ifelse(profile == 2L, 1, -1), 0.8),
    c = stats::rnorm(n_groups * size, ifelse(profile == 2L, 0.5, 0), 1.2))
  data$b <- data$b + 0.4 * data$a
  data$w <- rep(stats::runif(n_groups, 0.5, 2), each = size)
  data$w[data$school %in% c(3L, 7L)] <- 0
  data
}

# The objective, analytic score and free coordinates of a profile fit, built
# as parameter_inference() builds them.
audit_inference_pieces <- function(fit, data) {
  prepared <- latents:::.multilpa_inference_matrix(fit, data)
  centred <- fit
  centred$means <- sweep(fit$means, 2L, prepared$centers, "-")
  full <- latents:::.multilpa_coefficients(centred, "unconstrained")
  free <- latents:::.multilpa_free_index(fit, "unconstrained")
  restore <- latents:::.multilpa_restore_held(full, free)
  list(
    theta = unname(full)[free],
    objective = function(theta) {
      -latents:::.multilpa_expectation(
        prepared$x, fit$group_index, latents:::.multilpa_decode(restore(theta), fit),
        prepared$codes, unname(fit$sampling_weights), centred$extra_data)$log_likelihood
    },
    score = function(theta) {
      latents:::.multilpa_score(restore(theta), prepared$x, centred, prepared$codes)[free]
    },
    group_scores = function(theta) {
      latents:::.multilpa_group_scores(restore(theta), prepared$x, centred,
                                       prepared$codes)[, free, drop = FALSE]
    })
}

audit_inference_fit <- function(data, structure, weights = NULL) {
  do.call(multilpa, c(list(data, c("a", "b", "c"), "school", n_profiles = 2,
                           n_group_classes = 2, n_starts = 3, seed = 1, tol = 1e-10,
                           weights = weights),
                      latents:::.multilpa_structure_arguments(structure)))
}

test_that("analytic scores equal numerical gradients for every structure family, with and without weights", {
  skip_if_not_installed("numDeriv")
  data <- audit_inference_data()
  cases <- list(c("VVI", ""), c("VVV", ""), c("VEV", ""), c("VVI", "w"), c("VEV", "w"))
  invisible(lapply(cases, function(case) {
    fit <- audit_inference_fit(data, case[[1L]], if (nzchar(case[[2L]])) case[[2L]])
    pieces <- audit_inference_pieces(fit, data)
    set.seed(5)
    away <- pieces$theta + stats::rnorm(length(pieces$theta), 0, 0.05)
    numerical <- numDeriv::grad(pieces$objective, away)
    expect_lt(max(abs(numerical - pieces$score(away))) / max(abs(numerical)), 1e-5,
              label = paste(case, collapse = " "))
  }))
})

test_that("the Wald covariance is the inverse numerical Hessian, and the robust one the CR1 sandwich over positive-weight groups", {
  skip_if_not_installed("numDeriv")
  data <- audit_inference_data()
  plain <- audit_inference_fit(data, "VVV")
  pieces <- audit_inference_pieces(plain, data)
  inverse <- solve(numDeriv::hessian(pieces$objective, pieces$theta))
  wald <- attr(parameter_inference(plain, data), "covariance_unconstrained")
  expect_lt(max(abs(wald - inverse)) / max(abs(inverse)), 1e-6)

  weighted <- audit_inference_fit(data, "VVI", "w")
  pieces <- audit_inference_pieces(weighted, data)
  bread <- solve(numDeriv::hessian(pieces$objective, pieces$theta))
  scores <- pieces$group_scores(pieces$theta)
  positive <- sum(weighted$sampling_weights > 0)
  expect_identical(positive, 38L)
  sandwich <- bread %*% crossprod(scores) %*% bread * positive / (positive - 1)
  robust <- attr(parameter_inference(weighted, data), "covariance_unconstrained")
  expect_lt(max(abs(robust - sandwich)) / max(abs(sandwich)), 1e-5)
})

test_that("relative entropy is 1 - E / (n log K) with n the weighted classified rows; zero-weight rows keep their posterior row", {
  data <- audit_inference_data()
  fit <- multilpa(data, c("a", "b"), "school", n_profiles = 3, n_group_classes = 1,
                  n_starts = 2, seed = 1, weights = "w")
  posteriors <- fit$subject_posteriors
  row_weights <- fit$sampling_weights[fit$group_index]
  cells <- ifelse(posteriors > 0, posteriors * log(posteriors), 0)
  entropy <- -sum(row_weights * rowSums(cells))
  expect_equal(latents:::.multilpa_level_relative_entropy(fit, "individuals"),
               1 - entropy / (sum(row_weights) * log(3)), tolerance = 1e-12)
  expect_equal(latents:::.multilpa_complete_entropy(fit), entropy, tolerance = 1e-12)
  expect_identical(nrow(posteriors), nrow(data))
  expect_identical(nrow(latents:::.multilpa_classified_posteriors(fit)),
                   sum(row_weights > 0))
})

test_that("a weighted fit's refusal names the bootstrap only where its inference offers one", {
  refusal <- tryCatch(latents:::.latents_weighted_vcov(TRUE, "observed", FALSE),
                      latents_unsupported_weights = function(condition) condition)
  expect_s3_class(refusal, "latents_unsupported_weights")
  expect_false(grepl("bootstrap", conditionMessage(refusal), fixed = TRUE))
  offered <- tryCatch(latents:::.latents_weighted_vcov(TRUE, "opg", FALSE, bootstrap = TRUE),
                      latents_unsupported_weights = function(condition) condition)
  expect_true(grepl("method = \"bootstrap\"", conditionMessage(offered), fixed = TRUE))
  expect_identical(latents:::.latents_weighted_vcov(TRUE, "observed", TRUE), "robust")
  expect_identical(latents:::.latents_weighted_vcov(FALSE, "opg", FALSE), "opg")

  # Public paths: a mixture regression has no bootstrap; a profile fit has one.
  set.seed(3)
  rows <- data.frame(x = stats::rnorm(120), w = stats::runif(120, 0.5, 2))
  rows$y <- 1 + rows$x + stats::rnorm(120)
  regression <- mixture_regression(y ~ x, rows, n_classes = 1, weights = "w", n_starts = 0)
  refusal <- tryCatch(vcov(regression, type = "observed"),
                      latents_unsupported_weights = function(condition) condition)
  expect_s3_class(refusal, "latents_unsupported_weights")
  expect_false(grepl("bootstrap", conditionMessage(refusal), fixed = TRUE))

  profile <- audit_inference_fit(audit_inference_data(), "VVI", "w")
  refusal <- tryCatch(parameter_inference(profile, vcov_type = "observed"),
                      latents_unsupported_weights = function(condition) condition)
  expect_true(grepl("method = \"bootstrap\"", conditionMessage(refusal), fixed = TRUE))
})

test_that("an observed information that cannot be formed or inverted is refused with its class", {
  failing <- function(displacement) stop("the likelihood cannot be evaluated here")
  expect_error(latents:::.multilpa_observed_hessian(failing, failing, c(1, 1), 1e-4),
               class = "latents_singular_information")
})

test_that("fit_staged() classes its own argument errors", {
  data <- audit_inference_data()
  expect_error(fit_staged(data, c("a", "b"), "school", n_profiles = 2, n_group_classes = 1),
               class = "latents_bad_argument")
  expect_error(fit_staged(data, c("a", "b"), "school", n_profiles = 2, measurement = list()),
               class = "latents_bad_argument")
  expect_error(fit_staged(data, c("a", "b"), "school", n_profiles = 2,
                          variance_model = "bogus"),
               class = "latents_bad_argument")
  expect_error(fit_staged(data, c("a", "b"), "school", n_profiles = 2, missing = "listwise"),
               class = "latents_bad_argument")
})
