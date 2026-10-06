# SQUAREM acceleration of the multilpa EM: same maxima as plain EM, a
# monotone path, honest iteration counts, and plain EM where it must be.

acceleration_data <- function() {
  engagement <- c("browse", "lectures", "forum_read")
  list(data = subset(course_engagement, student <= 40), vars = engagement)
}

test_that("SQUAREM reaches plain EM's maximum from the same start", {
  skip_on_cran()
  fixture <- acceleration_data()
  fit <- function(acceleration, ...) {
    multilpa(fixture$data, fixture$vars, id = "student", n_profiles = 2,
             n_group_classes = 2, n_starts = 3, seed = 1, tol = 1e-12,
             acceleration = acceleration, ...)
  }
  plain <- fit("none")
  fast <- fit("squarem")
  expect_equal(fast$log_likelihood, plain$log_likelihood, tolerance = 1e-8)
  expect_equal(fast$means, plain$means, tolerance = 1e-4)
  expect_lte(fast$iterations, plain$iterations)
  full_plain <- fit("none", covariance_model = "full")
  full_fast <- fit("squarem", covariance_model = "full")
  expect_equal(full_fast$log_likelihood, full_plain$log_likelihood,
               tolerance = 1e-8)
})

test_that("the accelerated path is monotone and counts every EM step", {
  fixture <- acceleration_data()
  fit <- multilpa(fixture$data, fixture$vars, id = "student", n_profiles = 3,
                  n_group_classes = 2, n_starts = 1, seed = 2,
                  covariance_model = "full")
  history <- fit$log_likelihood_history
  expect_length(history, fit$iterations + 1L)
  expect_true(all(diff(history) >= -1e-10 * (1 + abs(history[-1L]))))
  expect_identical(fit$acceleration, "squarem")
})

test_that("a budget below one cycle takes plain steps", {
  fixture <- acceleration_data()
  expect_warning(
    one <- multilpa(fixture$data, fixture$vars, id = "student",
                    n_profiles = 2, n_group_classes = 2, n_starts = 1,
                    seed = 1, max_iter = 2),
    class = "latents_unconverged")
  expect_identical(one$iterations, 2L)
  expect_length(one$log_likelihood_history, 3L)
})

test_that("prior fits run plain EM whatever is requested", {
  fixture <- acceleration_data()
  x <- stats::na.omit(fixture$data[fixture$vars])
  quiet <- function(expression) {
    withCallingHandlers(expression, latents_single_level = function(w) {
      invokeRestart("muffleMessage")
    }, latents_unconverged = function(w) invokeRestart("muffleWarning"))
  }
  fit <- function(acceleration) {
    quiet(multilpa(x, fixture$vars, id = NULL, n_profiles = 2, n_starts = 1,
                   seed = 1, max_iter = 40, prior = prior_control(),
                   acceleration = acceleration))
  }
  expect_identical(fit("squarem")$log_likelihood_history,
                   fit("none")$log_likelihood_history)
})

test_that("a capped EVE/VVE orientation step never lowers the M-step objective", {
  # The objective the covariance M-step maximizes, given scatter W_k and
  # effective sizes n_k: -sum_k [n_k log|Sigma_k| + tr(Sigma_k^-1 W_k)] / 2.
  objective <- function(covariances, scatter, weights) {
    -sum(vapply(seq_along(weights), function(k) {
      block <- covariances[, , k]
      weights[k] * as.numeric(determinant(block, logarithm = TRUE)$modulus) +
        sum(solve(block) * scatter[[k]])
    }, numeric(1))) / 2
  }
  scatter_set <- latents:::.mixture_with_seed(4, {
    lapply(1:3, function(k) {
      z <- matrix(stats::rnorm(60 * 3), 60) %*% matrix(stats::runif(9), 3)
      crossprod(z)
    })
  })
  weights <- c(60, 60, 60)
  invisible(lapply(c("EVE", "VVE"), function(code) {
    start <- latents:::.multilpa_structure_covariances(
      scatter_set, weights, code, 1e-8, inner_steps = 1L)
    previous <- objective(start, scatter_set, weights)
    # Each capped step from the last answer must not go downhill.
    values <- vapply(1:15, function(i) {
      start <<- latents:::.multilpa_structure_covariances(
        scatter_set, weights, code, 1e-8, start = start, inner_steps = 5L)
      objective(start, scatter_set, weights)
    }, numeric(1))
    expect_true(all(diff(c(previous, values)) >= -1e-8), info = code)
    full <- latents:::.multilpa_structure_covariances(scatter_set, weights,
                                                      code, 1e-8)
    expect_gte(objective(full, scatter_set, weights) + 1e-6,
               max(values))
  }))
})

test_that("a converged EVE/VVE fit is a fixed point of the full orientation solve", {
  fixture <- acceleration_data()
  x <- stats::na.omit(fixture$data[fixture$vars])
  invisible(lapply(c("EVE", "VVE"), function(code) {
    fit <- withCallingHandlers(
      do.call(multilpa, c(list(data = x, vars = fixture$vars, id = NULL,
                               n_profiles = 2, n_starts = 1, seed = 1,
                               tol = 1e-12, max_iter = 20000),
                          latents:::.multilpa_structure_arguments(code))),
      latents_single_level = function(w) invokeRestart("muffleMessage"))
    expect_true(fit$converged, info = code)
    posteriors <- fit$subject_posteriors
    data_matrix <- as.matrix(x)
    scatter <- lapply(seq_len(ncol(posteriors)), function(k) {
      residuals <- sweep(data_matrix, 2L, fit$means[k, ], "-")
      crossprod(residuals, residuals * posteriors[, k])
    })
    solved <- latents:::.multilpa_structure_covariances(
      scatter, colSums(posteriors), code, fit$min_variance,
      start = fit$covariances, tol = 1e-12, max_iter = 100000L)
    expect_equal(unname(solved), unname(fit$covariances), tolerance = 1e-4,
                 info = code)
  }))
})

test_that("flattening and refilling parameters is lossless", {
  fixture <- acceleration_data()
  fit <- multilpa(fixture$data, fixture$vars, id = "student", n_profiles = 2,
                  n_group_classes = 2, n_starts = 1, seed = 1,
                  covariance_model = "full")
  parameters <- fit[c("means", "variances", "covariances",
                      "profile_probabilities", "group_probabilities")]
  flat <- latents:::.multilpa_flatten_parameters(parameters)
  expect_identical(latents:::.multilpa_restore_parameters(parameters, flat),
                   parameters)
  expect_error(latents:::.multilpa_restore_parameters(parameters, flat[-1L]))
})

test_that("an extrapolated point outside the parameter space is recognised", {
  fixture <- acceleration_data()
  fit <- multilpa(fixture$data, fixture$vars, id = "student", n_profiles = 2,
                  n_group_classes = 2, n_starts = 1, seed = 1,
                  covariance_model = "full")
  parameters <- fit[c("means", "variances", "covariances",
                      "profile_probabilities", "group_probabilities")]
  expect_true(latents:::.multilpa_valid_parameters(parameters))
  negative <- parameters
  negative$profile_probabilities[1L, ] <- c(1.2, -0.2)
  expect_false(latents:::.multilpa_valid_parameters(negative))
  indefinite <- parameters
  indefinite$covariances[, , 1L] <- -diag(3)
  expect_false(latents:::.multilpa_valid_parameters(indefinite))
})

test_that("an extrapolation that lowers the likelihood falls back to plain EM", {
  # A toy one-parameter map, theta -> 1.1 theta + 1, whose likelihood rises
  # with theta but collapses beyond 10. From theta = 0 the two EM steps reach
  # 1 and 2.1; the S3 steplength extrapolates to 30, lands at 34, where the
  # likelihood is -1000, so the cycle must keep the plain point 2.1.
  log_likelihood <- function(theta) if (theta < 10) theta else -1000
  evaluate <- function(point) {
    list(log_likelihood = log_likelihood(drop(point$means)))
  }
  step <- function(point, point_expectation) {
    moved <- point
    moved$means[] <- 1.1 * point$means + 1
    list(state = moved, expectation = evaluate(moved))
  }
  start <- list(means = matrix(0))
  accelerator <- list(flatten = latents:::.multilpa_flatten_parameters,
                      restore = latents:::.multilpa_restore_parameters,
                      valid = latents:::.multilpa_valid_parameters)
  cycle <- latents:::.latents_squarem_cycle(step, start, evaluate(start),
                                            evaluate, accelerator)
  expect_false(cycle$accelerated)
  expect_equal(drop(cycle$state$means), 2.1)
  expect_equal(cycle$history, c(1, 2.1, 2.1))
  # Where the landing point is better, it is taken.
  generous <- function(theta) theta
  evaluate_up <- function(point) list(log_likelihood = generous(drop(point$means)))
  step_up <- function(point, point_expectation) {
    moved <- point
    moved$means[] <- 1.1 * point$means + 1
    list(state = moved, expectation = evaluate_up(moved))
  }
  taken <- latents:::.latents_squarem_cycle(step_up, start, evaluate_up(start),
                                            evaluate_up, accelerator)
  expect_true(taken$accelerated)
  expect_equal(drop(taken$state$means), 34)
})
