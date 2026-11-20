# mixture_regression(random = ) with one class is a linear mixed model with a
# random intercept and slope fitted by maximum likelihood, so its accessors
# have exact references: nlme::lme(method = "ML") for the log likelihood, the
# fixed effects, the residual standard deviation, the random-effect covariance
# and the number of persons; and the observed information of the marginal
# likelihood, written densely below and differentiated by numDeriv, for every
# standard error, Wald interval and z test. The observed information of a
# mixed model is not lme's (X'V^-1X)^-1, which drops the cross terms between
# fixed effects and variance components; that difference is declared.

growth_one_class_data <- function() {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(7)
  n <- 80
  out <- data.frame(id = rep(seq_len(n), each = 5), time = rep(0:4, n))
  effects <- cbind(stats::rnorm(n, 0, 1), stats::rnorm(n, 0, 0.4))
  out$y <- 2 + 0.5 * out$time + effects[out$id, 1] + effects[out$id, 2] * out$time +
    stats::rnorm(nrow(out), 0, 0.7)
  out[stats::runif(nrow(out)) > 0.1, ]
}

growth_data <- growth_one_class_data()
growth_fit <- mixture_regression(y ~ time, growth_data, n_classes = 1, id = "id",
                                 class_level = "group", random = ~ 1 + time,
                                 n_starts = 1, seed = 1, tol = 1e-12)
growth_terms <- c("(Intercept)", "time")
growth_beta <- c("coefficient.class_1.(Intercept)", "coefficient.class_1.time")

# The marginal log likelihood from its definition: each person's y is
# N(X beta, Z G Z' + sigma^2 I), with X = Z = (1, time).
growth_dense <- function(beta, sigma, covariance) {
  sum(vapply(split(growth_data, growth_data$id), function(person) {
    design <- cbind(1, person$time)
    v <- design %*% covariance %*% t(design) + diag(sigma^2, nrow(person))
    factor <- chol(v)
    r <- backsolve(factor, person$y - drop(design %*% beta), transpose = TRUE)
    -0.5 * (nrow(person) * log(2 * pi) + 2 * sum(log(diag(factor))) + sum(r^2))
  }, numeric(1)))
}

# The same likelihood in three coordinates: the package's unconstrained ones
# (log sigma; G = L L' with L = [exp(a), 0; c, exp(b)]), standard deviations
# and a correlation, and variances and a covariance.
growth_unconstrained <- function(theta) {
  factor <- matrix(c(exp(theta[4]), theta[6], 0, exp(theta[5])), 2L)
  growth_dense(theta[1:2], exp(theta[3]), tcrossprod(factor))
}
growth_by_sd <- function(theta) {
  sds <- theta[4:5]
  growth_dense(theta[1:2], theta[3],
               diag(sds) %*% matrix(c(1, theta[6], theta[6], 1), 2L) %*% diag(sds))
}
growth_by_variance <- function(theta) {
  growth_dense(theta[1:2], theta[3], matrix(c(theta[4], theta[6], theta[6], theta[5]), 2L))
}

# Wald inference from estimates and standard errors with normal quantiles.
growth_wald <- function(estimate, error) {
  data.frame(estimate = estimate, std_error = error, statistic = estimate / error,
             p_value = 2 * stats::pnorm(-abs(estimate / error)),
             conf_low = estimate - stats::qnorm(0.975) * error,
             conf_high = estimate + stats::qnorm(0.975) * error)
}

test_that("one class: log likelihood, estimates and persons equal nlme::lme() by maximum likelihood", {
  skip_without_reference("nlme")
  reference <- nlme::lme(y ~ time, random = ~ 1 + time | id, data = growth_data, method = "ML",
                         control = nlme::lmeControl(opt = "optim", msTol = 1e-12))
  expect_reference(as.numeric(logLik(growth_fit)), as.numeric(stats::logLik(reference)),
                   reference = "nlme",
                   covers = c("mixture_regression()", "logLik.latents_growth_mixture()"),
                   quantity = "log likelihood, one class", tolerance = 1e-8, measure = "relative")
  # Both optimisers stop on their own criteria; fixed effects are the GLS
  # solution given the variance components, which converge to about 1e-7.
  # 1e-6 relative.
  expect_reference(stats::setNames(coef(growth_fit)[growth_beta], growth_terms),
                   nlme::fixef(reference), reference = "nlme",
                   covers = c("mixture_regression()", "coef.latents_growth_mixture()"),
                   quantity = "fixed effects, one class", tolerance = 1e-6, measure = "relative")
  expect_reference(exp(coef(growth_fit)[["log_sigma.class_1"]]), reference$sigma, reference = "nlme",
                   covers = c("mixture_regression()", "coef.latents_growth_mixture()"),
                   quantity = "residual standard deviation, one class", tolerance = 1e-6,
                   measure = "relative")
  # Variance components are the slowest part of either optimiser (their
  # likelihood is flat to second order): 1e-4 relative.
  covariance <- as.matrix(nlme::getVarCov(reference))
  correlation <- covariance[1, 2] / sqrt(covariance[1, 1] * covariance[2, 2])
  expect_reference(subset(parameter_inference(growth_fit), block %in% c("random_variance", "random_sd", "random_correlation", "residual_sd")),
                   data.frame(block = c("residual_sd", "random_variance", "random_variance",
                                        "random_sd", "random_sd", "random_correlation"),
                              term = c("residual", "(Intercept)", "time", "(Intercept)", "time",
                                       "(Intercept):time"),
                              estimate = c(reference$sigma, diag(covariance), sqrt(diag(covariance)),
                                           correlation)),
                   reference = "nlme",
                   covers = c("mixture_regression()", "parameter_inference.latents_growth_mixture()"),
                   quantity = "variance components, one class", tolerance = 1e-4,
                   measure = "relative", key = c("block", "term"))
  # nobs() counts persons, the independent units: lme's groups.
  expect_reference(nobs(growth_fit), nrow(nlme::ranef(reference)), reference = "nlme",
                   covers = c("mixture_regression()", "nobs.latents_growth_mixture()"),
                   quantity = "number of persons", tolerance = 0)
  # lme's fixed-effect covariance is (X'V^-1X)^-1, the expected-information
  # block; latents inverts the full observed information.
  expect_known_difference(subset(as.data.frame(growth_fit), select = c(term, std_error)),
                          data.frame(term = growth_terms,
                                     std_error = unname(sqrt(diag(stats::vcov(reference))))),
                          reference = "nlme",
                          covers = c("mixture_regression()", "as.data.frame.latents_growth_mixture()"),
                          quantity = "fixed-effect standard errors, one class", tolerance = 1e-6,
                          measure = "relative", key = "term",
                          reason = "lme's standard errors are (X'V^-1X)^-1 with the variance components held fixed; latents inverts the full observed information, whose fixed-effect/variance cross terms are not zero away from expectation. They agree to about 2e-4 here.")
})

test_that("one class: standard errors, intervals and tests are the inverse numerical Hessian's", {
  skip_without_reference("numDeriv")
  theta <- unname(coef(growth_fit))
  unconstrained <- solve(-numDeriv::hessian(growth_unconstrained, theta))
  errors <- sqrt(diag(unconstrained))
  wald <- cbind(term = growth_terms, growth_wald(theta[1:2], errors[1:2]))
  # numDeriv's Richardson Hessian of a smooth sum over 80 persons is accurate
  # to about 1e-8; 1e-6 relative on everything derived from it.
  expect_reference(as.data.frame(growth_fit), wald, reference = "numDeriv",
                   covers = c("mixture_regression()", "as.data.frame.latents_growth_mixture()"),
                   quantity = "coefficient table", tolerance = 1e-6, measure = "relative",
                   key = "term")
  summarised <- summary(growth_fit)
  expect_reference(as.data.frame(summarised), wald, reference = "numDeriv",
                   covers = c("summary.latents_growth_mixture()",
                              "as.data.frame.summary_latents_growth_mixture()"),
                   quantity = "coefficient table, summary()", tolerance = 1e-6,
                   measure = "relative", key = "term")
  expect_reference(get_results(summarised), wald, reference = "numDeriv",
                   covers = c("summary.latents_growth_mixture()",
                              "get_results.summary_latents_growth_mixture()"),
                   quantity = "coefficient table, get_results(summary())", tolerance = 1e-6,
                   measure = "relative", key = "term")
  # confint(): Wald intervals of every unconstrained parameter.
  expect_reference(confint(growth_fit),
                   data.frame(parameter = names(coef(growth_fit)), estimate = theta,
                              conf_low = theta - stats::qnorm(0.975) * errors,
                              conf_high = theta + stats::qnorm(0.975) * errors),
                   reference = "numDeriv",
                   covers = c("mixture_regression()", "confint.latents_growth_mixture()"),
                   quantity = "95% Wald intervals of the unconstrained parameters",
                   tolerance = 1e-6, measure = "relative", key = "parameter")
  # parameter_inference(): fixed effects with their tests, and the variance
  # components' standard errors, each from the Hessian in its own coordinates
  # (at the maximum the inverse Hessian transforms as the delta method).
  covariance <- growth_fit$params$random_covariance[[1L]]
  sigma <- exp(theta[3])
  sds <- sqrt(diag(covariance))
  by_sd <- sqrt(diag(solve(-numDeriv::hessian(
    growth_by_sd, c(theta[1:2], sigma, sds, covariance[1, 2] / prod(sds))))))
  by_variance <- sqrt(diag(solve(-numDeriv::hessian(
    growth_by_variance, c(theta[1:2], sigma, diag(covariance), covariance[1, 2])))))
  inference <- parameter_inference(growth_fit)
  coefficients <- growth_wald(theta[1:2], errors[1:2])
  expect_reference(subset(inference, block == "coefficient"),
                   data.frame(term = growth_terms, estimate = coefficients$estimate,
                              standard_error = coefficients$std_error,
                              statistic = coefficients$statistic, p_value = coefficients$p_value,
                              conf_low = coefficients$conf_low, conf_high = coefficients$conf_high),
                   reference = "numDeriv",
                   covers = c("mixture_regression()", "parameter_inference.latents_growth_mixture()"),
                   quantity = "fixed-effect inference", tolerance = 1e-6, measure = "relative",
                   key = "term")
  expect_reference(subset(inference, block %in% c("residual_sd", "random_variance", "random_sd", "random_correlation")),
                   data.frame(block = c("residual_sd", "random_sd", "random_sd", "random_correlation",
                                        "random_variance", "random_variance"),
                              term = c("residual", "(Intercept)", "time", "(Intercept):time",
                                       "(Intercept)", "time"),
                              standard_error = c(by_sd[3:6], by_variance[4:5])),
                   reference = "numDeriv",
                   covers = c("mixture_regression()", "parameter_inference.latents_growth_mixture()"),
                   quantity = "standard errors of the variance components", tolerance = 1e-6,
                   measure = "relative", key = c("block", "term"))
})
