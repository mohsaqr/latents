# Direct effects of covariates on continuous indicators: the covariate moves
# an indicator's mean by the same amount in every class (Latent GOLD's
# `y <- 1 + Cluster + z`).

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
  within <- function(actual, printed) {
    expect_true(all(abs(actual - printed) <= 5e-5),
                info = paste(format(actual, digits = 8), collapse = " "))
  }
  within(fit$log_likelihood, -9396.9944)
  expect_equal(fit$n_parameters, 18)
  within(fit$bic, 18925.6267)
  effects <- get_results(fit, "direct_effects")
  expect_named(effects, c("indicator", "covariate", "estimate", "standard_error",
                          "statistic", "p_value", "conf_low", "conf_high"))
  within(effects$estimate, 0.7864)
  within(effects$standard_error, 0.0140)
  inference <- parameter_inference(fit)
  y3 <- subset(inference, term == "y3")
  # Latent GOLD numbers the classes as this fit does here.
  within(subset(y3, parameter == "mean" & outcome == "profile_1")$estimate, 0.0161)
  within(subset(y3, parameter == "mean" & outcome == "profile_1")$standard_error, 0.0421)
  within(subset(y3, parameter == "variance")$estimate, c(0.9775, 1.0451))
  within(subset(y3, parameter == "variance")$standard_error, c(0.0549, 0.0630))
})

test_that("the direct-effect score is the gradient of the likelihood", {
  data <- direct_data(n = 400)
  fit <- quietly(lpa(data, paste0("y", 1:4), 2, model = "VVI", n_starts = 3, seed = 1,
                     direct_effects = c("y3 ~ z", "y1 ~ z")))
  x <- fit$indicator_data
  point <- fit[c("means", "variances", "profile_probabilities", "group_probabilities",
                 "direct_effects")]
  point$direct_effects <- point$direct_effects + c(0.05, -0.03)
  likelihood <- function(effects) {
    point$direct_effects <- effects
    .multilpa_expectation(x, fit$group_index, point, extra = fit$extra_data)$log_likelihood
  }
  numeric <- vapply(seq_along(point$direct_effects), function(j) {
    step <- 1e-6
    up <- point$direct_effects
    up[j] <- up[j] + step
    down <- point$direct_effects
    down[j] <- down[j] - step
    (likelihood(up) - likelihood(down)) / (2 * step)
  }, numeric(1))
  posteriors <- .multilpa_expectation(x, fit$group_index, point,
                                      extra = fit$extra_data)$subject_posteriors
  analytic <- as.vector(.latents_direct_scores(x, fit$extra_data, posteriors, point))
  expect_equal(analytic, numeric, tolerance = 1e-6)
  # At the estimate the effects' score is zero.
  expect_lt(max(abs(.latents_direct_scores(x, fit$extra_data, fit$subject_posteriors,
                                           fit))), 1e-3)
})

test_that("EM converges as fast on a covariate far from zero", {
  # Fitted on centred covariates, the iterations do not depend on where the
  # covariate is centred; on raw values far from zero they crawled.
  data <- direct_data(n = 600)
  iterations <- vapply(c(0, 40), function(shift) {
    quietly(lpa(transform(data, z = z + shift), paste0("y", 1:4), 2, model = "VVI",
                n_starts = 1, seed = 1, tol = 1e-12, acceleration = "none",
                direct_effects = "y3 ~ z"))$iterations
  }, numeric(1))
  expect_identical(iterations[[2L]], iterations[[1L]])
})

test_that("shifting a covariate moves the means; scaling it scales the effect", {
  skip_on_cran()
  data <- direct_data(n = 600)
  base <- quietly(lpa(data, paste0("y", 1:4), 2, model = "VVI", n_starts = 5, seed = 1,
                      tol = 1e-12, direct_effects = "y3 ~ z"))
  shifted <- quietly(lpa(transform(data, z = z + 10), paste0("y", 1:4), 2, model = "VVI",
                         n_starts = 5, seed = 1, tol = 1e-12, direct_effects = "y3 ~ z"))
  expect_equal(shifted$log_likelihood, base$log_likelihood, tolerance = 1e-8)
  expect_equal(unname(shifted$direct_effects), unname(base$direct_effects), tolerance = 1e-6)
  expect_equal(unname(shifted$means[, "y3"]),
               unname(base$means[, "y3"] - 10 * base$direct_effects), tolerance = 1e-5)
  # Standard errors follow the covariate's units over many orders of magnitude.
  invisible(lapply(c(10, 1e4, 1e-4), \(factor) {
    scaled <- quietly(lpa(transform(data, z = z * factor), paste0("y", 1:4), 2,
                          model = "VVI", n_starts = 5, seed = 1, tol = 1e-12,
                          direct_effects = "y3 ~ z"))
    expect_equal(unname(scaled$direct_effects) * factor, unname(base$direct_effects),
                 tolerance = 1e-6)
    expect_equal(get_results(scaled, "direct_effects")$standard_error * factor,
                 get_results(base, "direct_effects")$standard_error, tolerance = 1e-4)
  }))
})

test_that("a direct effect nests the model without it and removes the covariate's association", {
  skip_on_cran()
  data <- direct_data(n = 800)
  data$y2 <- data$y2 + 0.8 * data$z
  plain <- quietly(lpa(data, paste0("y", 1:4), 2, model = "VVI", n_starts = 10, seed = 1,
                       tol = 1e-10))
  direct <- quietly(lpa(data, paste0("y", 1:4), 2, model = "VVI", n_starts = 10, seed = 1,
                        tol = 1e-10, direct_effects = c("y2 ~ z", "y3 ~ z")))
  expect_equal(direct$n_parameters, plain$n_parameters + 2)
  expect_gte(direct$log_likelihood, plain$log_likelihood - 1e-6)
  # The shared covariate makes y2 and y3 associated within the classes of the
  # plain model; given the covariate they are not.
  pair <- function(fit) {
    residuals <- get_results(fit, "residuals", by = "overall")
    max(subset(residuals, indicator_1 %in% c("y2", "y3") &
                 indicator_2 %in% c("y2", "y3"))$statistic)
  }
  expect_gt(pair(plain), 10)
  expect_lt(pair(direct), 2)
})

test_that("a direct-effect fit predicts, refits and simulates with its covariates", {
  skip_on_cran()
  data <- direct_data(n = 500)
  fit <- quietly(lpa(data, paste0("y", 1:4), 2, model = "VVI", n_starts = 5, seed = 1,
                     tol = 1e-10, direct_effects = "y3 ~ z"))
  expect_equal(predict(fit, newdata = data, type = "posterior"),
               predict(fit, type = "posterior"), tolerance = 1e-10)
  expect_error(predict(fit, newdata = subset(data, select = -z)), class = "latents_bad_data")
  arguments <- .multilpa_refit_arguments(fit)
  expect_identical(arguments$direct_effects, "y3 ~ z")
  refit <- quietly(do.call(multilpa, c(list(data = .multilpa_resolve_data(fit, NULL),
                                            n_starts = 5, seed = 1, tol = 1e-10),
                                       arguments)))
  expect_equal(refit$log_likelihood, fit$log_likelihood, tolerance = 1e-6)
  expect_equal(refit$n_parameters, fit$n_parameters)
  # Draws keep the covariates and move y3 with them.
  set.seed(4)
  drawn <- .multilpa_simulate(fit)
  expect_identical(drawn$z, data$z)
  expect_gt(stats::cor(drawn$y3, drawn$z), 0.6)
  # Two-level fits take direct effects, with robust standard errors.
  grouped <- transform(data, g = rep(seq_len(50), each = 10))
  two <- quietly(multilpa(grouped, paste0("y", 1:4), "g", 2, 1, n_starts = 3, seed = 1,
                          tol = 1e-10, direct_effects = "y3 ~ z"))
  robust <- parameter_inference(two, vcov_type = "robust")
  expect_true(is.finite(subset(robust, parameter == "direct_effect")$standard_error))
})

test_that("direct effects the model cannot hold are refused by class", {
  data <- transform(direct_data(n = 200), g = rep(seq_len(20), each = 10),
                    item = factor(stats::rbinom(200, 1, 0.5)), constant = 1,
                    gappy = replace(stats::rnorm(200), 3, NA))
  vars <- paste0("y", 1:4)
  fit_with <- function(effects, ...) {
    lpa(data, vars, 2, n_starts = 1, seed = 1, direct_effects = effects, ...)
  }
  expect_error(fit_with("y3 z"), class = "latents_bad_argument")
  expect_error(fit_with("y3 ~ "), class = "latents_bad_argument")
  expect_error(fit_with(c("y3 ~ z", "y3 ~ z")), class = "latents_bad_argument")
  expect_error(fit_with("y3 ~ y1"), class = "latents_bad_argument")
  expect_error(fit_with("y9 ~ z"), class = "latents_bad_argument")
  expect_error(fit_with("y3 ~ absent"), class = "latents_bad_data")
  expect_error(fit_with("y3 ~ constant"), class = "latents_bad_data")
  expect_error(fit_with("y3 ~ gappy"), class = "latents_bad_data")
  expect_error(fit_with("y3 ~ item"), class = "latents_bad_data")
  expect_error(fit_with("y3 ~ z", model = "VVV"), class = "latents_bad_argument")
  expect_error(fit_with("y3 ~ z", model = "VEI"), class = "latents_bad_argument")
  expect_error(fit_with("y3 ~ z", local_dependence = "y1:y2"), class = "latents_bad_argument")
  expect_error(fit_with("y3 ~ z", noise = TRUE), class = "latents_bad_argument")
  expect_error(multilpa(data, c(vars, "item"), NULL, 2, 1, categorical = "item",
                        direct_effects = "item ~ z", n_starts = 1, seed = 1),
               class = "latents_bad_argument")
  expect_error(multilpa(data, vars, "g", 2, 2, profile_covariates = "z",
                        direct_effects = "y3 ~ z", n_starts = 1, seed = 1),
               class = "latents_bad_argument")
  expect_error(multilpa(data, vars, "g", n_profiles = 2, family = "additive",
                        direct_effects = "y3 ~ z"), class = "latents_bad_argument")
  plain <- quietly(lpa(data, vars, 2, n_starts = 1, seed = 1))
  expect_error(get_results(plain, "direct_effects"), class = "latents_no_direct_effects")
  expect_false("direct_effects" %in% names(get_results(plain, "all")))
})

# The parameter count, number of direct effects and covariance pattern of every
# multilpa() fit made while `code` runs.
refits_during <- function(code) {
  seen <- new.env()
  seen$fits <- list()
  suppressMessages(trace("multilpa", where = asNamespace("latents"), print = FALSE,
                         exit = bquote({
                           fit <- returnValue()
                           assign("fits", c(get("fits", .(seen)), list(data.frame(
                             parameters = fit$n_parameters,
                             effects = length(fit$direct_effects),
                             pattern = !is.null(fit$covariance_pattern)))), envir = .(seen))
                         })))
  on.exit(suppressMessages(untrace("multilpa", where = asNamespace("latents"))))
  force(code)
  do.call(rbind, seen$fits)
}

test_that("bootstrap replicates are refitted with the models' direct effects and patterns", {
  skip_on_cran()
  data <- direct_data(n = 400)
  null <- quietly(lpa(data, c("y1", "y2", "y3"), 1, n_starts = 2, seed = 1,
                      direct_effects = "y3 ~ z"))
  alternative <- quietly(lpa(data, c("y1", "y2", "y3"), 2, n_starts = 2, seed = 1,
                             direct_effects = "y3 ~ z"))
  refits <- refits_during(quietly(bootstrap_lrt(null, alternative, iter = 2, n_starts = 2,
                                                seed = 8)))
  expect_setequal(refits$parameters, c(null$n_parameters, alternative$n_parameters))
  expect_true(all(refits$effects == 1L))
  paired_null <- quietly(lpa(data, c("y1", "y2", "y4"), 1, n_starts = 2, seed = 1,
                             local_dependence = "y1:y2"))
  paired <- quietly(lpa(data, c("y1", "y2", "y4"), 2, n_starts = 2, seed = 1,
                        local_dependence = "y1:y2"))
  refits <- refits_during(quietly(bootstrap_lrt(paired_null, paired, iter = 2,
                                                n_starts = 2, seed = 8)))
  expect_setequal(refits$parameters, c(paired_null$n_parameters, paired$n_parameters))
  expect_true(all(refits$pattern))
  # Models that differ in their effects are not nested in the test's sense.
  plain <- quietly(lpa(data, c("y1", "y2", "y3"), 2, n_starts = 2, seed = 1))
  expect_error(bootstrap_lrt(null, plain, iter = 2, n_starts = 1, seed = 8),
               class = "latents_incomparable_models")
})

test_that("a staged fit refuses a measurement with direct effects or chosen covariances", {
  data <- transform(direct_data(n = 200), g = rep(seq_len(20), each = 10))
  vars <- c("y1", "y2", "y3")
  effects <- quietly(multilpa(data, vars, "g", 2, 1, n_starts = 2, seed = 1,
                              direct_effects = "y3 ~ z"))
  expect_error(fit_staged(data, vars, "g", 2, 2, measurement = effects, n_starts = 2,
                          seed = 1), class = "latents_bad_stage")
  paired <- quietly(multilpa(data, vars, "g", 2, 1, n_starts = 2, seed = 1,
                             local_dependence = "y1:y2"))
  expect_error(fit_staged(data, vars, "g", 2, 2, measurement = paired, n_starts = 2,
                          seed = 1, covariance_model = "full"), class = "latents_bad_stage")
})
