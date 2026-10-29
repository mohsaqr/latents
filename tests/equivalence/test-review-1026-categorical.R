# Comparisons with independent references moved from tests/testthat/test-review-1026-categorical.R,
# so that they run in CI's equivalence job rather than under R CMD check.

test_that("adjacent-category ordinal probabilities and their Newton M-step", {
  skip_if_not_installed("numDeriv")
  log_probability <- .latents_ordinal_log_probabilities(c(.5, -.2, 1), c(.3, -.4, 0))
  expect_equal(rowSums(exp(log_probability)), rep(1, 3), tolerance = 1e-14)
  # log P(k) / P(k - 1) = a_k - a_(k-1) + eta_c
  expect_equal(t(diff(t(log_probability))),
               matrix(diff(c(0, .5, -.2, 1)), 3, 3, byrow = TRUE) + c(.3, -.4, 0),
               tolerance = 1e-14)
  withr::local_seed(2)
  counts <- matrix(stats::rexp(12) * 10, 3, 4)
  objective <- function(theta) {
    sum(counts * .latents_ordinal_log_probabilities(theta[1:3], c(theta[4:5], 0)))
  }
  solution <- .latents_ordinal_maximize(counts, numeric(3), numeric(3))
  theta <- c(solution$intercepts, solution$locations[1:2])
  expect_lt(max(abs(numDeriv::grad(objective, theta))), 1e-7)
  reference <- stats::optim(numeric(5), objective, method = "BFGS",
                            control = list(fnscale = -1, reltol = 1e-14, maxit = 1000))
  expect_gte(objective(theta), reference$value - 1e-10)
  expect_equal(theta, reference$par, tolerance = 1e-5)
})

test_that("negative-binomial Hessian and M-step match numDeriv and optim", {
  skip_if_not_installed("numDeriv")
  withr::local_seed(3)
  y <- stats::rnbinom(80, size = 2, mu = 4)
  weight <- matrix(stats::runif(160), 80, 2)
  objective <- function(theta) {
    sum(weight * vapply(1:2, function(k) stats::dnbinom(
      y, size = 1 / exp(theta[2 + k]), mu = exp(theta[k]), log = TRUE), numeric(80)))
  }
  theta <- log(c(3, 5, .4, .7))
  hessian <- .latents_negative_binomial_hessian(
    list(means = exp(theta[1:2]), dispersion = exp(theta[3:4])), y, weight, FALSE)
  expect_equal(hessian, numDeriv::hessian(objective, theta), tolerance = 1e-7)
  solution <- .latents_negative_binomial_maximize(y, weight, c(3, 5), c(.4, .7))
  reference <- stats::optim(theta, objective, method = "BFGS",
                            control = list(fnscale = -1, reltol = 1e-15, maxit = 2000))
  expect_gte(objective(log(c(solution$means, solution$dispersion))),
             reference$value - 1e-10)
  # Under-dispersed counts put the dispersion on its Poisson boundary.
  under <- .latents_negative_binomial_maximize(rep(c(3, 4, 5), 30), matrix(1, 90, 1), 4, .3)
  expect_equal(under$dispersion, .latents_min_dispersion)
  expect_equal(under$means, 4)
})

test_that("ordinal and count scores are the gradient of the observed log likelihood", {
  skip_if_not_installed("numDeriv")
  withr::local_seed(8)
  n <- 120
  data <- data.frame(o = sample(1:4, n, TRUE), y = stats::rpois(n, 3),
                     z = stats::rnbinom(n, size = 2, mu = 4))
  data$o[3] <- NA
  data$y[5] <- NA
  prior <- matrix(log(c(.4, .6)), n, 2, byrow = TRUE)
  joint <- function(theta, extra) {
    .latents_extra_log_density(extra, .latents_extra_decode(theta, extra, 2), n, 2) + prior
  }
  log_likelihood <- function(theta, extra) {
    value <- joint(theta, extra)
    top <- apply(value, 1L, max)
    sum(top + log(rowSums(exp(value - top))))
  }
  checks <- lapply(list(c("poisson", "varying"), c("negative_binomial", "varying"),
                        c("negative_binomial", "equal")), function(model) {
    extra <- .latents_set_count_model(
      .latents_prepare_extra(data, "o", c("y", "z"), missing = "fiml"), model[1], model[2])
    theta <- stats::rnorm(.latents_extra_n_parameters(extra, 2), 0, .3)
    value <- joint(theta, extra)
    posterior <- exp(value - apply(value, 1L, max))
    posterior <- posterior / rowSums(posterior)
    analytic <- drop(.latents_extra_scores(extra, posterior,
                                           .latents_extra_decode(theta, extra, 2)))
    expect_equal(analytic, numDeriv::grad(log_likelihood, theta, extra = extra),
                 tolerance = 1e-7)
    length(theta)
  })
  # (K - 1) + (C - 1) ordinal, C means per count, then C or 1 dispersions.
  expect_identical(unlist(checks), c(8L, 12L, 10L))
})

test_that("a mixed-indicator fit's posteriors and likelihood match hand computation", {
  withr::local_seed(10)
  n <- 150
  class <- sample(1:2, n, TRUE)
  data <- data.frame(
    g = stats::rnorm(n, c(0, 2)[class]), c = stats::rbinom(n, 1, c(.2, .8)[class]),
    o = vapply(class, function(k) sample(1:3, 1L, prob = if (k == 1) c(.6, .3, .1) else
      c(.1, .3, .6)), integer(1)),
    k = stats::rpois(n, c(2, 6)[class]))
  fit <- lpa(data, vars = c("g", "c", "o", "k"), categorical = "c", ordinal = "o",
             count = "k", n_profiles = 2, n_starts = 2, seed = 1)
  ordinal <- .latents_ordinal_log_probabilities(fit$ordinal_intercepts[[1]],
                                                fit$ordinal_locations[, 1])
  joint <- vapply(1:2, function(k) {
    log(fit$profile_probabilities[k]) +
      stats::dnorm(data$g, fit$means[k, "g"], sqrt(fit$variances[k, "g"]), log = TRUE) +
      log(fit$response_probabilities$c[k, data$c + 1]) + ordinal[k, data$o] +
      stats::dpois(data$k, fit$count_means[k, 1], log = TRUE)
  }, numeric(n))
  top <- apply(joint, 1L, max)
  posterior <- exp(joint - top) / rowSums(exp(joint - top))
  expect_equal(unname(fit$subject_posteriors), unname(posterior), tolerance = 1e-12)
  expect_equal(fit$log_likelihood, sum(top + log(rowSums(exp(joint - top)))),
               tolerance = 1e-12)
  # 1 class weight, 2 x 2 Gaussian, 2 x 1 binary, (3 - 1) + (2 - 1) ordinal, 2 means.
  expect_identical(as.integer(fit$n_parameters), 12L)
})
