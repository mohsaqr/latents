# Negative-binomial regression mixtures, mixture_regression(family =
# "negative_binomial") (phase E7 of validation/ENGINE_DESIGN.md): the GLM
# block with NB2 counts under the regression structures. References: one
# class is MASS::glm.nb(); a likelihood written with dnbinom() alone;
# numerical derivatives; recovery of a simulated truth.

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

test_that("one class is MASS::glm.nb", {
  skip_if_not_installed("MASS")
  data <- nb_data()
  fit <- mixture_regression(y ~ x, data, n_classes = 1,
                            family = "negative_binomial", n_starts = 1, seed = 1)
  reference <- MASS::glm.nb(y ~ x, data = data)
  expect_equal(fit$log_likelihood, as.numeric(stats::logLik(reference)),
               tolerance = 1e-8)
  expect_equal(unname(fit$params$beta[, 1L]), unname(stats::coef(reference)),
               tolerance = 1e-6)
  expect_equal(fit$params$sigma2, 1 / reference$theta, tolerance = 1e-5)
  # glm.nb's errors are Fisher's expected information with the dispersion
  # held fixed; ours are the observed information with it estimated, so they
  # agree to a few percent.
  coefficients <- get_results(fit, "coefficients")
  expect_equal(coefficients$std_error, unname(sqrt(diag(stats::vcov(reference)))),
               tolerance = 0.05)
})

test_that("the likelihood is the dense NB mixture and the scores its gradient", {
  data <- nb_data()
  fit <- mixture_regression(y ~ x, data, n_classes = 2,
                            family = "negative_binomial", n_starts = 3, seed = 1)
  expect_equal(fit$log_likelihood, nb_dense(fit, data), tolerance = 1e-10)
  skip_if_not_installed("numDeriv")
  theta <- latents:::.mixture_pack(fit$spec, fit$params) + 0.03
  value <- function(v) {
    latents:::.mixture_expectation(
      fit$spec, latents:::.mixture_unpack(fit$spec, v, fit$params))$log_likelihood
  }
  params <- latents:::.mixture_unpack(fit$spec, theta, fit$params)
  scores <- latents:::.mixture_unit_scores(fit$spec, params,
                                           latents:::.mixture_expectation(fit$spec, params))
  expect_equal(unname(colSums(scores)), numDeriv::grad(value, theta), tolerance = 1e-6)
})
