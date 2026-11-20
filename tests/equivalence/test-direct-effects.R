# Comparisons with independent references moved from tests/testthat/test-direct-effects.R,
# so that they run in CI's equivalence job rather than under R CMD check.

direct_data <- function(seed = 3, n = 1500) {
  set.seed(seed)
  class <- sample(1:2, n, TRUE, prob = c(0.55, 0.45))
  covariate <- stats::rnorm(n, 1, 2)
  data <- data.frame(z = round(covariate, 4),
                     y1 = round(stats::rnorm(n, c(0, 2)[class]), 4),
                     y2 = round(stats::rnorm(n, c(0, 1.5)[class]), 4), y3 = 0,
                     y4 = round(stats::rnorm(n, c(0, 1)[class]), 4))
  data$y3 <- round(stats::rnorm(n, c(0, 1.8)[class]) + 0.8 * data$z, 4)
  data
}

direct_fit <- function(data = direct_data(), ...) {
  quietly(lpa(data, paste0("y", 1:4), 2, model = "VVI", n_starts = 10, seed = 1,
              tol = 1e-12, direct_effects = "y3 ~ z", ...))
}

test_that("a direct effect agrees with Latent GOLD", {
  skip_on_cran()
  # Latent GOLD 6.1 on these data (equivalence/latentgold-extra/direct-effects):
  # `y3 <- 1 + Cluster + z`, class-specific variances. Log likelihood
  # -9396.9944 with 18 parameters, BIC 18925.6267, the effect 0.7864 (s.e.
  # 0.0140), the first class's y3 intercept 0.0161 (s.e. 0.0421) and the y3
  # variances 0.9775 and 1.0451 (s.e. 0.0549 and 0.0630).
  fit <- direct_fit()
  listing <- "Typed from Latent GOLD 6.1's listing equivalence/latentgold-extra/direct-effects/direct.lst."
  # 5e-5 absolute: half a unit of the last digit Latent GOLD prints.
  expect_reference(as.numeric(logLik(fit)), -9396.9944, reference = "Latent GOLD",
                   covers = c("lpa()", "logLik.multilpa()"),
                   quantity = "log likelihood, direct effect of z on y3",
                   tolerance = 5e-5, measure = "absolute", note = listing)
  expect_reference(fit$n_parameters, 18, reference = "Latent GOLD", covers = "lpa()",
                   quantity = "number of parameters, direct effect of z on y3",
                   tolerance = 0, measure = "absolute", note = listing)
  expect_reference(fit$bic, 18925.6267, reference = "Latent GOLD", covers = "lpa()",
                   quantity = "BIC, direct effect of z on y3",
                   tolerance = 5e-5, measure = "absolute", note = listing)
  effects <- get_results(fit, "direct_effects")
  expect_named(effects, c("indicator", "covariate", "estimate", "standard_error",
                          "statistic", "p_value", "conf_low", "conf_high"))
  expect_reference(effects,
                   data.frame(indicator = "y3", covariate = "z", estimate = 0.7864,
                              standard_error = 0.0140),
                   reference = "Latent GOLD", covers = c("lpa()", "get_results.multilpa()"),
                   quantity = "direct effect of z on y3", tolerance = 5e-5,
                   measure = "absolute", key = c("indicator", "covariate"), note = listing)
  # Latent GOLD numbers the classes as this fit does here.
  y3 <- subset(parameter_inference(fit),
               term == "y3" & (parameter == "variance" |
                                 (parameter == "mean" & outcome == "profile_1")),
               select = c(outcome, parameter, estimate, standard_error))
  expect_reference(y3,
                   data.frame(outcome = c("profile_1", "profile_1", "profile_2"),
                              parameter = c("mean", "variance", "variance"),
                              estimate = c(0.0161, 0.9775, 1.0451),
                              standard_error = c(0.0421, 0.0549, 0.0630)),
                   reference = "Latent GOLD", covers = c("lpa()", "parameter_inference.multilpa()"),
                   quantity = "y3 intercept and variances with a direct effect", tolerance = 5e-5,
                   measure = "absolute", key = c("outcome", "parameter"), note = listing)
})
