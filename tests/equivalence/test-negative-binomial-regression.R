# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Negative-binomial regression mixtures, mixture_regression(family =
# "negative_binomial") (phase E7 of validation/ENGINE_DESIGN.md): the GLM
# block with NB2 counts under the regression structures. References: one
# class is MASS::glm.nb(); a likelihood written with dnbinom() alone;
# numerical derivatives of that likelihood.

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

# The NB2 mixture log likelihood written with dnbinom() alone, from the
# estimates coef() reports: class coefficients and log dispersions, class
# shares from the multinomial logits (class 1 the reference).
nb_dense <- function(estimates, data, n_classes) {
  shares <- exp(c(0, estimates[sprintf("membership.class_%d.(Intercept)", seq_len(n_classes)[-1L])]))
  shares <- shares / sum(shares)
  sum(log(rowSums(vapply(seq_len(n_classes), function(k) {
    mu <- exp(estimates[[sprintf("coefficient.class_%d.(Intercept)", k)]] +
                estimates[[sprintf("coefficient.class_%d.x", k)]] * data$x)
    size <- 1 / exp(estimates[[sprintf("log_dispersion.class_%d", k)]])
    shares[k] * stats::dnbinom(data$y, size = size, mu = mu)
  }, numeric(nrow(data))))))
}

test_that("one class is MASS::glm.nb", {
  skip_without_reference("MASS")
  data <- nb_data()
  fit <- mixture_regression(y ~ x, data, n_classes = 1,
                            family = "negative_binomial", n_starts = 1, seed = 1)
  reference <- MASS::glm.nb(y ~ x, data = data)
  # The two optimisers stop at their own convergence criteria; testthat's
  # relative tolerance of the original test, kept.
  expect_reference(as.numeric(logLik(fit)), as.numeric(stats::logLik(reference)),
                   reference = "MASS", covers = c("mixture_regression()", "logLik.latents_mixture_regression()"),
                   quantity = "log likelihood, one class", tolerance = 1e-8, measure = "relative")
  # Coefficients to 1e-6 relative: glm.nb alternates IRLS and a theta step
  # to its own tolerance (1e-8 on the deviance), not to the joint optimum.
  expect_reference(coef(fit)[c("coefficient.class_1.(Intercept)", "coefficient.class_1.x")],
                   stats::coef(reference), reference = "MASS",
                   covers = c("mixture_regression()", "coef.latents_mixture_regression()"),
                   quantity = "coefficients, one class", tolerance = 1e-6, measure = "relative",
                   key = c("(Intercept)", "x"))
  # latents reports log(dispersion); glm.nb reports theta = 1 / dispersion.
  # 1e-5 relative: theta is the slowest parameter of glm.nb's alternation.
  expect_reference(exp(coef(fit)[["log_dispersion.class_1"]]), 1 / reference$theta, reference = "MASS",
                   covers = c("mixture_regression()", "coef.latents_mixture_regression()"),
                   quantity = "dispersion 1 / theta, one class", tolerance = 1e-5, measure = "relative")
  # glm.nb's errors are Fisher's expected information with the dispersion
  # held fixed; latents' are the observed information with it estimated.
  expect_known_difference(sqrt(diag(vcov(fit)))[c("coefficient.class_1.(Intercept)", "coefficient.class_1.x")],
                          sqrt(diag(stats::vcov(reference))), reference = "MASS",
                          covers = c("mixture_regression()", "vcov.latents_mixture_regression()"),
                          quantity = "standard errors of the coefficients, one class", tolerance = 1e-6,
                          measure = "relative", key = c("(Intercept)", "x"),
                          reason = "glm.nb's standard errors use the expected information with the dispersion held fixed; latents uses the observed information with the dispersion estimated. They agree to a few percent.")
})

test_that("the likelihood is the dense NB mixture and the scores its gradient", {
  data <- nb_data()
  fit <- mixture_regression(y ~ x, data, n_classes = 2,
                            family = "negative_binomial", n_starts = 3, seed = 1)
  # The same sum of logs computed twice; only the order of the additions differs.
  expect_reference(as.numeric(logLik(fit)), nb_dense(coef(fit), data, 2L),
                   reference = "independent: NB2 mixture likelihood from stats::dnbinom()",
                   covers = c("mixture_regression()", "logLik.latents_mixture_regression()"),
                   quantity = "log likelihood, two classes", tolerance = 1e-10, measure = "relative")
  skip_without_reference("numDeriv")
  # The analytic scores are the package's internal routine behind the standard
  # errors of mixture_regression(); it has no exported accessor.
  spec <- fit$spec
  theta <- latents:::.mixture_pack(spec, fit$params) + 0.03
  value <- function(v) {
    latents:::.mixture_expectation(spec, latents:::.mixture_unpack(spec, v, fit$params))$log_likelihood
  }
  params <- latents:::.mixture_unpack(spec, theta, fit$params)
  scores <- latents:::.mixture_unit_scores(spec, params, latents:::.mixture_expectation(spec, params))
  # numDeriv's Richardson extrapolation is accurate to about 1e-8 relative on
  # a smooth sum of 1000 terms; 1e-6 relative, the original test's tolerance.
  expect_reference(unname(colSums(scores)), numDeriv::grad(value, theta), reference = "numDeriv",
                   covers = "mixture_regression()", quantity = "score at a point off the optimum",
                   tolerance = 1e-6, measure = "scaled",
                   note = "Package side: the internal analytic scores (.mixture_unit_scores) behind vcov(); reference: numDeriv::grad() of the package's likelihood.")
})
