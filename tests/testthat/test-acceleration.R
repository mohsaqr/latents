# SQUAREM acceleration of the multilpa EM: same maxima as plain EM, a
# monotone path, honest iteration counts, and plain EM where it must be.

acceleration_data <- function() {
  engagement <- c("browse", "lectures", "forum_read")
  list(data = subset(course_engagement, student <= 40), vars = engagement)
}

test_that("SQUAREM reaches plain EM's maximum from the same start", {
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

test_that("EVE, VVE and prior fits run plain EM whatever is requested", {
  fixture <- acceleration_data()
  x <- stats::na.omit(fixture$data[fixture$vars])
  paired <- function(...) {
    arguments <- list(data = x, vars = fixture$vars, id = NULL,
                      n_profiles = 2, n_starts = 1, seed = 1, max_iter = 40,
                      ...)
    quiet <- function(expression) {
      withCallingHandlers(expression, latents_single_level = function(w) {
        invokeRestart("muffleWarning")
      }, latents_unconverged = function(w) invokeRestart("muffleWarning"))
    }
    list(squarem = quiet(do.call(multilpa, c(arguments,
                                             acceleration = "squarem"))),
         none = quiet(do.call(multilpa, c(arguments, acceleration = "none"))))
  }
  eve <- paired(volume = "equal", shape = "varying", orientation = "equal")
  expect_identical(eve$squarem$log_likelihood_history,
                   eve$none$log_likelihood_history)
  map <- paired(prior = prior_control())
  expect_identical(map$squarem$log_likelihood_history,
                   map$none$log_likelihood_history)
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
    list(parameters = moved, expectation = evaluate(moved))
  }
  start <- list(means = matrix(0))
  cycle <- latents:::.multilpa_squarem_cycle(step, start, evaluate(start),
                                             evaluate)
  expect_false(cycle$accelerated)
  expect_equal(drop(cycle$parameters$means), 2.1)
  expect_equal(cycle$history, c(1, 2.1, 2.1))
  # Where the landing point is better, it is taken.
  generous <- function(theta) theta
  evaluate_up <- function(point) list(log_likelihood = generous(drop(point$means)))
  step_up <- function(point, point_expectation) {
    moved <- point
    moved$means[] <- 1.1 * point$means + 1
    list(parameters = moved, expectation = evaluate_up(moved))
  }
  taken <- latents:::.multilpa_squarem_cycle(step_up, start, evaluate_up(start),
                                             evaluate_up)
  expect_true(taken$accelerated)
  expect_equal(drop(taken$parameters$means), 34)
})
