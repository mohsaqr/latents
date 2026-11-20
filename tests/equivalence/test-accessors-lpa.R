# lpa() on the Old Faithful eruptions (datasets::faithful, 272 eruptions, two
# indicators, two profiles, diagonal class-varying covariances = mclust's VVI)
# and the accessors that report it. References: mclust's E-step at latents'
# estimates (posteriors and log likelihood) and its EM continued from latents'
# posteriors (the maximum); the information criteria, entropy, average
# posteriors and odds of correct classification from their definitions;
# the standard errors from numDeriv's Hessian of a likelihood written with
# dnorm() alone; Wald intervals and z tests from their definitions.

faithful_x <- as.matrix(datasets::faithful)
faithful_fit <- lpa(datasets::faithful, c("eruptions", "waiting"), n_profiles = 2,
                    n_starts = 5, seed = 1)
faithful_coef <- coef(faithful_fit)
faithful_n <- nrow(faithful_x)

# The estimates as matrices (indicator x profile), read from coef() by name.
faithful_means <- vapply(1:2, \(k) faithful_coef[sprintf("measurement.mean.profile_%d.%s", k, colnames(faithful_x))],
                         numeric(2))
faithful_variances <- vapply(1:2, \(k) faithful_coef[sprintf("measurement.variance.profile_%d.%s", k, colnames(faithful_x))],
                             numeric(2))
faithful_shares <- faithful_coef[sprintf("profile.probability.profile_%d.group_class_1", 1:2)]

# mclust's E-step for these parameters, written as its VVV model with diagonal covariances.
faithful_estep <- function() {
  sigma <- array(unlist(lapply(1:2, \(k) diag(faithful_variances[, k]))), c(2L, 2L, 2L))
  parameters <- list(pro = unname(faithful_shares), mean = unname(faithful_means),
                     variance = list(modelName = "VVV", d = 2L, G = 2L, sigma = sigma,
                                     cholsigma = array(apply(sigma, 3L, chol), c(2L, 2L, 2L))))
  mclust::estepVVV(data = faithful_x, parameters = parameters)
}

# Posteriors as one row per unit and profile, the layout of get_results(, "posteriors").
long_posteriors <- function(z) {
  data.frame(row = rep(seq_len(nrow(z)), ncol(z)), profile = rep(seq_len(ncol(z)), each = nrow(z)),
             posterior = as.vector(z))
}

test_that("posteriors and log likelihood equal mclust's E-step at latents' estimates", {
  skip_without_reference("mclust")
  step <- faithful_estep()
  # One E-step evaluated by two programs on the same parameters: rounding only.
  expect_reference(as.numeric(logLik(faithful_fit)), step$loglik, reference = "mclust",
                   covers = c("lpa()", "logLik.multilpa()", "coef.multilpa()"),
                   quantity = "log likelihood at latents' estimates", tolerance = 1e-10, measure = "relative")
  # Posteriors are probabilities in [0, 1]; many are below 1e-10, so the difference is absolute.
  expect_reference(get_results(faithful_fit, "posteriors"), long_posteriors(step$z), reference = "mclust",
                   covers = c("lpa()", "get_results()", "get_results.multilpa()"),
                   quantity = "posterior profile probabilities", tolerance = 1e-10, key = c("row", "profile"))
})

test_that("the estimates are the maximum mclust's EM reaches from latents' posteriors", {
  skip_without_reference("mclust")
  continued <- mclust::meVVI(data = faithful_x, z = faithful_estep()$z,
                             control = mclust::emControl(tol = c(1e-12, 1e-12)))
  variances <- apply(continued$parameters$variance$sigma, 3L, diag)
  expected <- data.frame(profile = rep(1:2, each = 2), indicator = rep(colnames(faithful_x), 2),
                         mean = as.vector(continued$parameters$mean), variance = as.vector(variances))
  # Both stop at a relative change of 1e-12 or below in the log likelihood; the
  # estimates then agree to about the square root of that: 1e-6 relative.
  expect_reference(as.data.frame(faithful_fit), expected, reference = "mclust",
                   covers = c("lpa()", "as.data.frame.multilpa()"),
                   quantity = "profile means and variances", tolerance = 1e-6, measure = "relative",
                   key = c("profile", "indicator"))
  expect_reference(as.data.frame(summary(faithful_fit)), expected, reference = "mclust",
                   covers = c("summary.multilpa()", "as.data.frame.summary_multilpa()"),
                   quantity = "profile means and variances, summary()", tolerance = 1e-6, measure = "relative",
                   key = c("profile", "indicator"))
  expect_reference(get_results(summary(faithful_fit)), expected, reference = "mclust",
                   covers = c("summary.multilpa()", "get_results.summary_multilpa()"),
                   quantity = "profile means and variances, get_results(summary())", tolerance = 1e-6,
                   measure = "relative", key = c("profile", "indicator"))
  expect_reference(unname(faithful_shares), continued$parameters$pro, reference = "mclust",
                   covers = c("lpa()", "coef.multilpa()"), quantity = "profile shares",
                   tolerance = 1e-6, measure = "relative", key = c("profile 1", "profile 2"))
  expect_reference(as.numeric(logLik(faithful_fit)), continued$loglik, reference = "mclust",
                   covers = c("lpa()", "logLik.multilpa()"), quantity = "maximised log likelihood",
                   tolerance = 1e-10, measure = "relative")
})

test_that("information criteria, entropy and classification summaries follow their definitions", {
  z <- matrix(get_results(faithful_fit, "posteriors")$posterior, ncol = 2)
  log_likelihood <- as.numeric(logLik(faithful_fit))
  n_parameters <- 4 + 4 + 1          # means, variances, one free share
  n <- faithful_n
  entropy_sum <- -sum(ifelse(z > 0, z * log(z), 0))
  deviance <- -2 * log_likelihood
  # Akaike (1973); Cavanaugh (1999) KIC; Schwarz (1978); Sclove (1987); Bozdogan (1987) CAIC;
  # Banfield and Raftery (1993) AWE; Biernacki, Celeux and Govaert (2000) ICL-BIC; Biernacki and
  # Govaert (1997) CLC. Sums of a handful of doubles: rounding only.
  expected <- data.frame(
    log_likelihood = log_likelihood, n_parameters = n_parameters,
    aic = deviance + 2 * n_parameters, kic = deviance + 3 * (n_parameters + 1),
    bic_individual = deviance + n_parameters * log(n),
    sabic_individual = deviance + n_parameters * log((n + 2) / 24),
    caic_individual = deviance + n_parameters * (log(n) + 1),
    awe_individual = deviance + 2 * entropy_sum + 2 * n_parameters * (3 / 2 + log(n)),
    icl_individual = deviance + n_parameters * log(n) + 2 * entropy_sum,
    clc_individual = deviance + 2 * entropy_sum, key = "lpa")
  expect_reference(cbind(get_results(faithful_fit, "information_criteria"), key = "lpa"), expected,
                   reference = "independent: classification summaries, criteria and counts from their definitions, base R",
                   covers = c("lpa()", "get_results.multilpa()"), quantity = "information criteria",
                   tolerance = 1e-10, measure = "relative", key = "key")
  # Ramaswamy, DeSarbo, Reibstein and Robinson (1993): 1 - EN / (n log K).
  expect_reference(subset(get_entropy(faithful_fit), level == "individuals", select = c(level, entropy_sum, relative_entropy)),
                   data.frame(level = "individuals", entropy_sum = entropy_sum,
                              relative_entropy = 1 - entropy_sum / (n * log(2))),
                   reference = "independent: classification summaries, criteria and counts from their definitions, base R",
                   covers = "get_entropy()", quantity = "entropy", tolerance = 1e-10, measure = "relative",
                   key = "level")
  # Nagin (2005): the average posterior of each profile among the units assigned to each profile.
  modal <- max.col(z, ties.method = "first")
  avepp <- expand.grid(assigned_class = 1:2, class = 1:2)
  avepp$n_assigned <- tabulate(modal, 2L)[avepp$assigned_class]
  avepp$average_posterior <- mapply(\(j, k) mean(z[modal == j, k]), avepp$assigned_class, avepp$class)
  expect_reference(get_avepp(faithful_fit), avepp,
                   reference = "independent: classification summaries, criteria and counts from their definitions, base R",
                   covers = "get_avepp()", quantity = "average posterior probabilities", tolerance = 1e-12,
                   measure = "relative", key = c("assigned_class", "class"))
  # Nagin (2005): modal counts, posterior-weighted sizes, the diagonal of the AvePP matrix, and the odds of
  # correct classification OCC_k = (AvePP_k / (1 - AvePP_k)) / (pi_k / (1 - pi_k)).
  shares <- colMeans(z)
  diagonal <- vapply(1:2, \(k) mean(z[modal == k, k]), numeric(1))
  classification <- data.frame(class = 1:2, n_modal = tabulate(modal, 2L), proportion_modal = tabulate(modal, 2L) / n,
                               estimated_n = colSums(z), estimated_proportion = shares, average_posterior = diagonal,
                               odds_correct_classification = (diagonal / (1 - diagonal)) / (shares / (1 - shares)))
  expect_reference(get_classification(faithful_fit), classification,
                   reference = "independent: classification summaries, criteria and counts from their definitions, base R",
                   covers = "get_classification()", quantity = "classification summary", tolerance = 1e-10,
                   measure = "relative", key = "class")
  expect_reference(nobs(faithful_fit), n, reference = "independent: classification summaries, criteria and counts from their definitions, base R",
                   covers = "nobs.multilpa()", quantity = "number of observations", tolerance = 0)
})

test_that("standard errors are the inverse of numDeriv's Hessian of the dense likelihood", {
  skip_without_reference("numDeriv")
  free <- names(faithful_coef)[1:9]   # means, variances, the first share (the second is 1 minus it)
  log_likelihood <- function(theta) {
    means <- matrix(theta[1:4], 2L, 2L)
    variances <- matrix(theta[5:8], 2L, 2L)
    shares <- c(theta[[9L]], 1 - theta[[9L]])
    density <- vapply(1:2, \(k) shares[k] * stats::dnorm(faithful_x[, 1L], means[1L, k], sqrt(variances[1L, k])) *
                        stats::dnorm(faithful_x[, 2L], means[2L, k], sqrt(variances[2L, k])), numeric(faithful_n))
    sum(log(rowSums(density)))
  }
  expect_reference(log_likelihood(unname(faithful_coef[free])), as.numeric(logLik(faithful_fit)),
                   reference = "independent: two-profile normal mixture likelihood from stats::dnorm()",
                   covers = c("lpa()", "logLik.multilpa()"), quantity = "log likelihood, dense",
                   tolerance = 1e-10, measure = "relative")
  covariance <- solve(-numDeriv::hessian(log_likelihood, unname(faithful_coef[free])))
  dimnames(covariance) <- list(free, free)
  # numDeriv's Richardson extrapolation gives the diagonal of the Hessian to about 1e-8 relative; its
  # inverse carries that error. The correlations rest on near-zero off-diagonal second derivatives,
  # which numDeriv gives to about 1e-6 absolute.
  expect_reference(sqrt(diag(vcov(faithful_fit)))[free], sqrt(diag(covariance)), reference = "numDeriv",
                   covers = c("vcov.multilpa()", "lpa()"), quantity = "standard errors (observed information)",
                   tolerance = 1e-6, measure = "relative")
  expect_reference(stats::cov2cor(vcov(faithful_fit)[free, free]), stats::cov2cor(covariance), reference = "numDeriv",
                   covers = "vcov.multilpa()", quantity = "correlations of the estimates",
                   tolerance = 1e-5, measure = "absolute")
  # Wald intervals and z tests from the standard errors of numDeriv's covariance.
  inference <- parameter_inference(faithful_fit)
  mean_rows <- subset(inference, parameter == "mean")
  mean_names <- sprintf("measurement.mean.%s.%s", mean_rows$outcome, mean_rows$term)
  errors <- sqrt(diag(covariance))[mean_names]
  estimates <- faithful_coef[mean_names]
  wald <- data.frame(outcome = mean_rows$outcome, term = mean_rows$term, estimate = unname(estimates),
                     standard_error = unname(errors), statistic = unname(estimates / errors),
                     conf_low = unname(estimates - stats::qnorm(0.975) * errors),
                     conf_high = unname(estimates + stats::qnorm(0.975) * errors))
  expect_reference(mean_rows, wald, reference = "numDeriv",
                   covers = c("parameter_inference()", "parameter_inference.multilpa()"),
                   quantity = "Wald inference for the profile means", tolerance = 1e-6, measure = "relative",
                   key = c("outcome", "term"))
  expect_reference(confint(faithful_fit)[mean_names, ], cbind(`2.5%` = wald$conf_low, `97.5%` = wald$conf_high),
                   reference = "numDeriv", covers = "confint.multilpa()", quantity = "95% intervals of the profile means",
                   tolerance = 1e-6, measure = "relative", key = c(paste(mean_names, "low"), paste(mean_names, "high")))
})
