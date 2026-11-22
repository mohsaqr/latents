# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# multilpa() with a profile covariate (two profiles, two indicators, diagonal
# class-varying variances, membership logit on z, one group class) and the
# accessors that report it. References: a likelihood written with dnorm() and
# plogis() alone; Bayes' rule at latents' estimates for the posteriors; one
# Newton step of that likelihood with numDeriv's gradient and Hessian for the
# maximum; numDeriv's Hessian inverted for the standard errors; Wald intervals
# from their definition; the number of groups for nobs().

covariate_data <- local({
  set.seed(6)
  n <- 300L
  frame <- data.frame(g = rep(seq_len(30L), each = 10L), z = stats::rnorm(n))
  state <- 1L + (stats::runif(n) > stats::plogis(-0.3 + 1.2 * frame$z))
  frame$y1 <- stats::rnorm(n, c(-1.5, 1.5)[state], c(1, 0.8)[state])
  frame$y2 <- stats::rnorm(n, c(-1, 1)[state], c(0.9, 1.1)[state])
  frame
})
covariate_fit <- multilpa(covariate_data, c("y1", "y2"), "g", n_profiles = 2, n_group_classes = 1,
                          profile_covariates = "z", n_starts = 3, seed = 1, tol = 1e-14,
                          max_iter = 5000)
covariate_coef <- coef(covariate_fit)

# The components of the likelihood at theta, in coef()'s order: four means
# (profile major), four variances (profile major), the profile-1 logit's
# intercept and slope on z (profile 2 the reference).
covariate_terms <- function(theta) {
  means <- matrix(theta[1:4], 2L, 2L, byrow = TRUE)
  variances <- matrix(theta[5:8], 2L, 2L, byrow = TRUE)
  first <- stats::plogis(theta[[9L]] + theta[[10L]] * covariate_data$z)
  shares <- cbind(first, 1 - first)
  vapply(1:2, \(k) shares[, k] *
           stats::dnorm(covariate_data$y1, means[k, 1L], sqrt(variances[k, 1L])) *
           stats::dnorm(covariate_data$y2, means[k, 2L], sqrt(variances[k, 2L])),
         numeric(nrow(covariate_data)))
}
covariate_loglik <- function(theta) sum(log(rowSums(covariate_terms(theta))))

test_that("log likelihood, posteriors and nobs follow their definitions", {
  theta <- unname(covariate_coef)
  # The same sum of logs computed twice; only the order of the additions differs.
  expect_reference(as.numeric(logLik(covariate_fit)), covariate_loglik(theta),
                   reference = "independent: covariate mixture likelihood from dnorm() and plogis()",
                   covers = c("multilpa()", "logLik.multilpa_covariates()"),
                   quantity = "log likelihood, profile covariate", tolerance = 1e-10, measure = "relative")
  terms <- covariate_terms(theta)
  posterior <- terms / rowSums(terms)
  # Bayes' rule at latents' estimates. Posteriors are probabilities, many near
  # 0 or 1, so the difference is absolute.
  expect_reference(get_results(covariate_fit, "posteriors"),
                   data.frame(row = rep(seq_len(nrow(posterior)), 2L), profile = rep(1:2, each = nrow(posterior)),
                              posterior = as.vector(posterior)),
                   reference = "independent: Bayes' rule at latents' estimates, base R",
                   covers = c("multilpa()", "get_results.multilpa_covariates()"),
                   quantity = "posterior profile probabilities, profile covariate", tolerance = 1e-10,
                   measure = "absolute", key = c("row", "profile"))
  # nobs() of a two-level fit counts the groups, the independent units.
  expect_reference(nobs(covariate_fit), length(unique(covariate_data$g)),
                   reference = "independent: number of groups in the data",
                   covers = c("multilpa()", "nobs.multilpa_covariates()"),
                   quantity = "number of observations (groups)", tolerance = 0)
})

test_that("estimates, standard errors and intervals agree with numDeriv on the dense likelihood", {
  skip_without_reference("numDeriv")
  theta <- unname(covariate_coef)
  gradient <- numDeriv::grad(covariate_loglik, theta)
  hessian <- numDeriv::hessian(covariate_loglik, theta)
  # One Newton step from latents' estimates: the maximum of the dense
  # likelihood to second order.
  newton <- theta - solve(hessian, gradient)
  errors <- sqrt(diag(solve(-hessian)))
  names(newton) <- names(errors) <- names(covariate_coef)
  profiles <- data.frame(profile = rep(1:2, each = 2L), indicator = rep(c("y1", "y2"), 2L),
                         mean = newton[1:4], variance = newton[5:8],
                         mean_standard_error = errors[1:4], variance_standard_error = errors[5:8])
  note <- "Estimates: one Newton step of the dense likelihood (numDeriv gradient and Hessian) from latents' estimates; standard errors: numDeriv's Hessian inverted."
  # 1e-6 relative, as for lpa(): EM converges linearly, so the estimates stop
  # short of the maximum by more than the likelihood does. At tol = 1e-12 they
  # stopped up to 7e-6 short (seen in a check while writing this file), hence
  # tol = 1e-14 above. The standard errors rest on numDeriv's Hessian.
  expect_reference(as.data.frame(covariate_fit), profiles, reference = "numDeriv",
                   covers = c("multilpa()", "as.data.frame.multilpa_covariates()"),
                   quantity = "profile means, variances and their standard errors", tolerance = 1e-6,
                   measure = "relative", key = c("profile", "indicator"), note = note)
  summarised <- summary(covariate_fit)
  expect_reference(as.data.frame(summarised), profiles, reference = "numDeriv",
                   covers = c("summary.multilpa_covariates()", "as.data.frame.summary_multilpa_covariates()"),
                   quantity = "profile means, variances and their standard errors, summary()",
                   tolerance = 1e-6, measure = "relative", key = c("profile", "indicator"), note = note)
  expect_reference(get_results(summarised), profiles, reference = "numDeriv",
                   covers = c("summary.multilpa_covariates()", "get_results.summary_multilpa_covariates()"),
                   quantity = "profile means, variances and their standard errors, get_results(summary())",
                   tolerance = 1e-6, measure = "relative", key = c("profile", "indicator"), note = note)
  expect_reference(get_results(covariate_fit, "coefficients"),
                   data.frame(level = "profile", outcome = "profile_1", term = c("(Intercept)", "z"),
                              estimate = newton[9:10]),
                   reference = "numDeriv", covers = c("multilpa()", "get_results.multilpa_covariates()"),
                   quantity = "membership logit coefficients", tolerance = 1e-6, measure = "relative",
                   key = c("level", "outcome", "term"), note = note)
  # Wald intervals from numDeriv's standard errors: symmetric for the means and
  # the logits, on the log scale for the variances.
  critical <- stats::qnorm(0.975)
  estimates <- covariate_coef
  variance <- grepl("variance", names(estimates))
  half <- critical * errors
  low <- ifelse(variance, estimates * exp(-half / estimates), estimates - half)
  high <- ifelse(variance, estimates * exp(half / estimates), estimates + half)
  intervals <- cbind(`2.5%` = low, `97.5%` = high)
  rownames(intervals) <- names(estimates)
  # The intervals carry the standard errors' 1e-6 relative error. Scaled: an
  # endpoint near zero (the intercept's interval crosses it) would magnify that
  # error under a relative measure.
  expect_reference(confint(covariate_fit), intervals, reference = "numDeriv",
                   covers = c("multilpa()", "confint.multilpa_covariates()"),
                   quantity = "95% Wald intervals, profile covariate", tolerance = 1e-6, measure = "scaled",
                   note = "Wald intervals built in the test from numDeriv's standard errors at latents' estimates.")
})

test_that("predict() of a membership-covariate fit is not implemented (a declared API gap)", {
  # The reference a result would be compared with: Bayes' rule with the
  # membership logits on the new rows (covariate_terms() above).
  expect_defect_error(predict(covariate_fit, newdata = covariate_data),
                      defect = "predict() is not implemented for membership-covariate models: it always raises latents_unsupported_prediction; FINDINGS.md (latents, API gaps)",
                      covers = "predict.multilpa_covariates()",
                      quantity = "posterior profile probabilities of new rows, profile covariate")
})
