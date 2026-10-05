# Negative-binomial regression mixtures, mixture_regression(family =
# "negative_binomial") (phase E7 of validation/ENGINE_DESIGN.md): the GLM
# block with NB2 counts under the regression structures. Here: a likelihood
# written with dnbinom() alone and recovery of a simulated truth. The
# comparisons with MASS::glm.nb() and numerical derivatives are in
# tests/equivalence/.

nb_data <- function(seed = 5, n = 1000) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(seed)
  x <- stats::runif(n, 0, 2)
  class <- 1L + (stats::runif(n) < 0.4)
  mu <- exp(ifelse(class == 1L, 0.5 + 0.8 * x, 3 - x))
  data.frame(y = stats::rnbinom(n, size = 1 / c(0.2, 0.25)[class], mu = mu),
             x = x, class = class)
}

nb_dense <- function(fit, data) {
  params <- fit$params
  shares <- exp(params$gamma[1L, ]) / sum(exp(params$gamma[1L, ]))
  sum(log(rowSums(vapply(seq_len(fit$spec$n_classes), function(k) {
    mu <- exp(params$beta[1L, k] + params$beta[2L, k] * data$x)
    shares[k] * stats::dnbinom(data$y, size = 1 / params$sigma2[k], mu = mu)
  }, numeric(nrow(data))))))
}

test_that("a simulated NB mixture is recovered", {
  data <- nb_data()
  fit <- mixture_regression(y ~ x, data, n_classes = 2,
                            family = "negative_binomial", n_starts = 3, seed = 1)
  coefficients <- get_results(fit, "coefficients")
  slopes <- sort(subset(coefficients, term == "x")$estimate)
  expect_equal(slopes, c(-1, 0.8), tolerance = 0.15)
  classes <- get_results(fit, "classes")
  expect_equal(sort(classes$share), c(0.4, 0.6), tolerance = 0.06)
  expect_equal(classes$dispersion, c(0.2, 0.25), tolerance = 0.3)
  expect_true(all(is.finite(classes$dispersion_std_error)))
  # Invariant: posteriors are probabilities.
  expect_equal(unname(rowSums(fit$expectation$tau)), rep(1, nrow(data)), tolerance = 1e-12)
})

test_that("a shared dispersion is one parameter, and counts equal estimates", {
  data <- nb_data()
  fit <- mixture_regression(y ~ x, data, n_classes = 2, variance = "equal",
                            family = "negative_binomial", n_starts = 2, seed = 1)
  expect_identical(length(unique(fit$params$sigma2)), 1L)
  expect_identical(fit$n_parameters, length(coef(fit)))
  expect_identical(fit$n_parameters, 6L)
  draws <- simulate(fit, nsim = 2, seed = 4)
  expect_identical(simulate(fit, nsim = 2, seed = 4), draws)
  expect_true(all(unlist(draws) >= 0) && all(unlist(draws) == round(unlist(draws))))
})

test_that("Poisson counts put the dispersion on its boundary, with a warning", {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(9)
  data <- data.frame(x = stats::runif(400, 0, 2))
  data$y <- stats::rpois(400, exp(1 + 0.5 * data$x))
  expect_warning(
    fit <- mixture_regression(y ~ x, data, n_classes = 1,
                              family = "negative_binomial", n_starts = 1, seed = 1),
    class = "latents_boundary")
  poisson <- mixture_regression(y ~ x, data, n_classes = 1, family = "poisson",
                                n_starts = 1, seed = 1)
  expect_equal(fit$log_likelihood, poisson$log_likelihood, tolerance = 1e-6)
  coefficients <- get_results(fit, "coefficients")
  expect_true(all(is.finite(coefficients$std_error)))
  expect_true(is.na(get_results(fit, "classes")$dispersion_std_error))
})

test_that("NB requests outside the family's support are refused", {
  data <- nb_data(n = 200)
  data$y[1L] <- 2.5
  expect_error(mixture_regression(y ~ x, data, n_classes = 2,
                                  family = "negative_binomial", n_starts = 1),
               class = "latents_bad_data")
  data$y[1L] <- -1
  expect_error(mixture_regression(y ~ x, data, n_classes = 2,
                                  family = "negative_binomial", n_starts = 1),
               class = "latents_bad_data")
  expect_error(mixture_regression(y ~ x, nb_data(n = 200), n_classes = 2,
                                  family = "poisson", variance = "equal", n_starts = 1),
               class = "latents_bad_argument")
})
