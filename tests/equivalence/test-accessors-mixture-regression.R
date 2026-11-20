# mixture_regression() with one class is a normal linear model fitted by
# maximum likelihood, so every accessor of its result has an exact reference
# in stats::lm(): the coefficients, the ML residual standard deviation
# sqrt(RSS / n), the log likelihood, the fitted values and predictions, and the
# observed-information covariance sigma_ML^2 (X'X)^-1 = vcov(lm) (n - p) / n,
# whose log-sigma variance is 1 / (2n). Wald intervals and z tests follow from
# those standard errors with normal quantiles.

one_class_data <- function() {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(3)
  n <- 300
  x <- stats::runif(n)
  data.frame(x = x, y = 1 + 2 * x + stats::rnorm(n, sd = 0.5))
}

regression_data <- one_class_data()
regression_fit <- mixture_regression(y ~ x, regression_data, n_classes = 1, n_starts = 1, seed = 1)
regression_lm <- stats::lm(y ~ x, data = regression_data)
regression_n <- nrow(regression_data)
beta_names <- c("coefficient.class_1.(Intercept)", "coefficient.class_1.x")

# The coefficient table of the linear model with ML (normal-quantile) Wald inference.
lm_wald <- function() {
  estimates <- stats::coef(regression_lm)
  errors <- sqrt(diag(stats::vcov(regression_lm)) * (regression_n - 2) / regression_n)
  data.frame(term = names(estimates), estimate = unname(estimates), std_error = unname(errors),
             statistic = unname(estimates / errors), p_value = unname(2 * stats::pnorm(-abs(estimates / errors))),
             conf_low = unname(estimates - stats::qnorm(0.975) * errors),
             conf_high = unname(estimates + stats::qnorm(0.975) * errors))
}

test_that("one class: estimates and log likelihood equal lm() fitted by maximum likelihood", {
  # The one-class model has a closed-form optimum; the package reaches it by its general optimiser,
  # which stops at a gradient of about 1e-8: 1e-8 relative on the estimates.
  expect_reference(coef(regression_fit)[beta_names], stats::coef(regression_lm), reference = "stats",
                   covers = c("mixture_regression()", "coef.latents_mixture_regression()"),
                   quantity = "coefficients, one class", tolerance = 1e-8, measure = "relative",
                   key = c("(Intercept)", "x"))
  expect_reference(exp(coef(regression_fit)[["log_sigma.class_1"]]),
                   sqrt(sum(stats::residuals(regression_lm)^2) / regression_n), reference = "stats",
                   covers = c("mixture_regression()", "coef.latents_mixture_regression()"),
                   quantity = "residual standard deviation (ML), one class", tolerance = 1e-8, measure = "relative")
  expect_reference(as.numeric(logLik(regression_fit)), as.numeric(stats::logLik(regression_lm)), reference = "stats",
                   covers = c("mixture_regression()", "logLik.latents_mixture_regression()"),
                   quantity = "log likelihood, one class", tolerance = 1e-10, measure = "relative")
  expect_reference(nobs(regression_fit), stats::nobs(regression_lm), reference = "stats",
                   covers = "nobs.latents_mixture_regression()", quantity = "number of observations", tolerance = 0)
})

test_that("one class: the covariance is the ML-scaled vcov(lm), and inference follows from it", {
  # The observed information at the optimum is exact here; 1e-6 relative allows for the package's
  # second derivatives and for the 1e-8 stopping of the optimiser.
  expected <- stats::vcov(regression_lm) * (regression_n - 2) / regression_n
  dimnames(expected) <- list(beta_names, beta_names)
  expect_reference(vcov(regression_fit)[beta_names, beta_names], expected, reference = "stats",
                   covers = "vcov.latents_mixture_regression()", quantity = "covariance of the coefficients",
                   tolerance = 1e-6, measure = "relative")
  expect_reference(vcov(regression_fit)[["log_sigma.class_1", "log_sigma.class_1"]], 1 / (2 * regression_n),
                   reference = "independent: variance of log sigma at the ML estimate, 1 / (2n)",
                   covers = "vcov.latents_mixture_regression()", quantity = "variance of log sigma",
                   tolerance = 1e-6, measure = "relative")
  wald <- lm_wald()
  intervals <- cbind(`2.5 %` = wald$conf_low, `97.5 %` = wald$conf_high)
  rownames(intervals) <- beta_names
  expect_reference(confint(regression_fit)[beta_names, ], intervals, reference = "stats",
                   covers = "confint.latents_mixture_regression()", quantity = "95% Wald intervals of the coefficients",
                   tolerance = 1e-6, measure = "relative")
  # The p-values run to 1e-77: relative.
  expect_reference(as.data.frame(regression_fit), wald, reference = "stats",
                   covers = c("mixture_regression()", "as.data.frame.latents_mixture_regression()"),
                   quantity = "coefficient table", tolerance = 1e-6, measure = "relative", key = "term")
  expect_reference(get_results(regression_fit, "coefficients"), wald, reference = "stats",
                   covers = "get_results.latents_mixture_regression()", quantity = "coefficient table, get_results()",
                   tolerance = 1e-6, measure = "relative", key = "term")
  expect_reference(as.data.frame(summary(regression_fit)), wald, reference = "stats",
                   covers = c("summary.latents_mixture_regression()", "as.data.frame.summary_latents_mixture_regression()"),
                   quantity = "coefficient table, summary()", tolerance = 1e-6, measure = "relative", key = "term")
  expect_reference(get_results(summary(regression_fit)), wald, reference = "stats",
                   covers = "get_results.summary_latents_mixture_regression()",
                   quantity = "coefficient table, get_results(summary())", tolerance = 1e-6, measure = "relative",
                   key = "term")
})

test_that("one class: fitted values and predictions equal lm()'s", {
  expect_reference(predict(regression_fit), data.frame(row = seq_len(regression_n), fitted = unname(stats::fitted(regression_lm))),
                   reference = "stats", covers = "predict.latents_mixture_regression()", quantity = "fitted values",
                   tolerance = 1e-8, measure = "relative", key = "row")
  new_data <- data.frame(x = c(-0.5, 0, 0.25, 1.5))
  expect_reference(predict(regression_fit, newdata = new_data),
                   data.frame(row = seq_len(nrow(new_data)), fitted = unname(stats::predict(regression_lm, new_data))),
                   reference = "stats", covers = "predict.latents_mixture_regression()",
                   quantity = "predictions for new data", tolerance = 1e-8, measure = "relative", key = "row")
})
