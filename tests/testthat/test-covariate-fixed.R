# `start` and `fixed` in membership-covariate models. With the measurement
# held at an unconditional fit's estimates, the covariate model is the two-step
# estimator (Bakk and Kuha 2018): only the membership logits are free, so its
# likelihood and standard errors can be checked against a direct
# maximisation over the logits with the class densities fixed.

skip_on_cran()

two_step_case <- function() {
  set.seed(3)
  n <- 600L
  z <- stats::rnorm(n)
  state <- vapply(z, function(value) {
    sample.int(3L, 1L, prob = exp(c(0.8, 0, -0.8) * value))
  }, integer(1))
  centres <- rbind(c(-2, 0), c(0, 2), c(2, -2))
  data <- data.frame(y1 = centres[state, 1L] + stats::rnorm(n),
                     y2 = centres[state, 2L] + stats::rnorm(n), z = z)
  plain <- lpa(data, c("y1", "y2"), n_profiles = 3, model = "VVI", n_starts = 3, seed = 1)
  two_step <- lpa(data, c("y1", "y2"), n_profiles = 3, model = "VVI", profile_covariates = "z",
                  start = starting_values(plain), fixed = "measurement", n_starts = 3, seed = 1)
  list(data = data, plain = plain, two_step = two_step)
}

# The logit-only negative log likelihood with the plain fit's class densities.
logit_objective <- function(case) {
  densities <- vapply(seq_len(3L), function(class) {
    with(case$plain, stats::dnorm(case$data$y1, means[class, 1L], sqrt(variances[class, 1L])) *
           stats::dnorm(case$data$y2, means[class, 2L], sqrt(variances[class, 2L])))
  }, numeric(nrow(case$data)))
  design <- cbind(1, case$data$z)
  function(coefficients) {
    eta <- design %*% cbind(matrix(coefficients, 2L, 2L), 0)
    shares <- exp(eta - apply(eta, 1L, max))
    -sum(log(rowSums(shares / rowSums(shares) * densities)))
  }
}

test_that("a held measurement model leaves only the membership logits free", {
  case <- two_step_case()
  fit <- case$two_step
  expect_identical(fit$fixed, c("means", "variances"))
  expect_equal(unname(fit$means), unname(case$plain$means), tolerance = 1e-12)
  expect_equal(unname(fit$variances), unname(case$plain$variances), tolerance = 1e-12)
  expect_identical(as.integer(fit$n_parameters), 4L)
  expect_identical(as.integer(fit$n_parameters_with_measurement), 16L)
  # Its likelihood is the maximum over the logits alone.
  objective <- logit_objective(case)
  optimum <- stats::optim(rep(0, 4L), objective, method = "BFGS",
                          control = list(reltol = 1e-14, maxit = 1000L))
  expect_equal(fit$log_likelihood, -optimum$value, tolerance = 1e-8)
  # The one-step model frees the measurement too, so it fits at least as well.
  one_step <- lpa(case$data, c("y1", "y2"), n_profiles = 3, model = "VVI",
                  profile_covariates = "z", n_starts = 3, seed = 1)
  expect_gte(one_step$log_likelihood, fit$log_likelihood - 1e-8)
})

test_that("standard errors are conditional on the held measurement model", {
  skip_if_not_installed("numDeriv")
  case <- two_step_case()
  table <- as.data.frame(parameter_inference(case$two_step))
  expect_identical(unique(table$parameter), "coefficient")
  coefficients <- as.vector(case$two_step$profile_coefficients)
  expected <- sqrt(diag(solve(numDeriv::hessian(logit_objective(case), coefficients))))
  expect_equal(table$standard_error, expected, tolerance = 1e-4)
})

test_that("`start` alone sets where the covariate model begins", {
  case <- two_step_case()
  default <- lpa(case$data, c("y1", "y2"), n_profiles = 3, model = "VVI",
                 profile_covariates = "z", n_starts = 3, seed = 1)
  started <- lpa(case$data, c("y1", "y2"), n_profiles = 3, model = "VVI",
                 profile_covariates = "z", start = starting_values(case$plain), n_starts = 1)
  expect_equal(started$log_likelihood, default$log_likelihood, tolerance = 1e-6)
  expect_identical(started$fixed, character())
  expect_identical(as.integer(started$n_parameters), as.integer(default$n_parameters))
})

test_that("held blocks follow the rules of a staged fit", {
  case <- two_step_case()
  expect_error(lpa(case$data, c("y1", "y2"), n_profiles = 3, profile_covariates = "z",
                   fixed = "means", n_starts = 1), class = "latents_bad_fixed")
  # A constrained structure ties the profiles' variances together.
  expect_error(lpa(case$data, c("y1", "y2"), n_profiles = 3, model = "VEI",
                   profile_covariates = "z", start = starting_values(case$plain),
                   fixed = "variances", n_starts = 1), class = "latents_bad_argument")
  means_only <- lpa(case$data, c("y1", "y2"), n_profiles = 3, model = "VEI",
                    profile_covariates = "z", start = starting_values(case$plain),
                    fixed = "means", n_starts = 1)
  expect_equal(unname(means_only$means), unname(case$plain$means), tolerance = 1e-12)
  expect_identical(means_only$covariance_structure, "VEI")
})
