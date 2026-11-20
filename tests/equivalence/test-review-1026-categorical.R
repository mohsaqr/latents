# Comparisons with independent references moved from tests/testthat/test-review-1026-categorical.R,
# so that they run in CI's equivalence job rather than under R CMD check.

# The adjacent-category logit written from its definition in base R:
# log P(k) - log P(1) = a_(k-1) + (k - 1) eta, a the K - 1 category intercepts
# and eta the class location; one row per class.
categorical_ordinal_reference <- function(intercepts, locations) {
  t(vapply(locations, \(eta) {
    unnormalised <- c(0, intercepts + seq_along(intercepts) * eta)
    top <- max(unnormalised)
    unnormalised - top - log(sum(exp(unnormalised - top)))
  }, numeric(length(intercepts) + 1L)))
}

ordinal_note <- "Package side: the internal adjacent-category routines (.latents_ordinal_log_probabilities, .latents_ordinal_maximize) of lpa()/multilpa(ordinal =)."
count_note <- "Package side: the internal negative-binomial routines (.latents_negative_binomial_hessian, .latents_negative_binomial_maximize) of lpa()/multilpa(count_model = \"negative_binomial\")."

test_that("adjacent-category ordinal probabilities and their Newton M-step", {
  skip_without_reference("numDeriv")
  log_probability <- latents:::.latents_ordinal_log_probabilities(c(.5, -.2, 1), c(.3, -.4, 0))
  expect_equal(rowSums(exp(log_probability)), rep(1, 3), tolerance = 1e-14)
  # log P(k) / P(k - 1) = a_k - a_(k-1) + eta_c
  expect_reference(unname(t(diff(t(log_probability)))),
                   matrix(diff(c(0, .5, -.2, 1)), 3, 3, byrow = TRUE) + c(.3, -.4, 0),
                   reference = "independent: adjacent-category logit definition, base R",
                   covers = c("lpa()", "multilpa()"),
                   quantity = "adjacent-category log odds, three classes",
                   tolerance = 1e-14, measure = "relative", note = ordinal_note)
  withr::local_seed(2)
  counts <- matrix(stats::rexp(12) * 10, 3, 4)
  objective <- function(theta) {
    sum(counts * latents:::.latents_ordinal_log_probabilities(theta[1:3], c(theta[4:5], 0)))
  }
  solution <- latents:::.latents_ordinal_maximize(counts, numeric(3), numeric(3))
  theta <- c(solution$intercepts, solution$locations[1:2])
  labels <- c(sprintf("intercept %d", 1:3), sprintf("location %d", 1:2))
  # The package's solution claims stationarity: a zero gradient. numDeriv's
  # gradient of the objective there, to the original test's absolute bound
  # (numDeriv's error near zero is absolute).
  expect_reference(rep(0, 5), numDeriv::grad(objective, theta), reference = "numDeriv",
                   covers = c("lpa()", "multilpa()"),
                   quantity = "gradient of the ordinal objective at the M-step solution",
                   tolerance = 1e-7, measure = "absolute", key = labels,
                   note = paste(ordinal_note, "The package side is the zero gradient its M-step solution claims."))
  reference <- stats::optim(numeric(5), objective, method = "BFGS",
                            control = list(fnscale = -1, reltol = 1e-14, maxit = 1000))
  expect_gte(objective(theta), reference$value - 1e-10)
  # BFGS stops at its own criterion: 1e-5 relative, the original test's tolerance.
  expect_reference(theta, reference$par, reference = "stats", covers = c("lpa()", "multilpa()"),
                   quantity = "ordinal M-step solution against optim() BFGS",
                   tolerance = 1e-5, measure = "relative", key = labels,
                   note = paste(ordinal_note, "Reference: stats::optim() maximising the same objective."))
})

test_that("negative-binomial Hessian and M-step match numDeriv and optim", {
  skip_without_reference("numDeriv")
  withr::local_seed(3)
  y <- stats::rnbinom(80, size = 2, mu = 4)
  weight <- matrix(stats::runif(160), 80, 2)
  objective <- function(theta) {
    sum(weight * vapply(1:2, function(k) stats::dnbinom(
      y, size = 1 / exp(theta[2 + k]), mu = exp(theta[k]), log = TRUE), numeric(80)))
  }
  theta <- log(c(3, 5, .4, .7))
  hessian <- latents:::.latents_negative_binomial_hessian(
    list(means = exp(theta[1:2]), dispersion = exp(theta[3:4])), y, weight, FALSE)
  # Scaled: the cross-class entries are exactly zero, where numDeriv's error is
  # absolute; a relative measure would judge its rounding noise. 1e-7, the
  # original test's tolerance (numDeriv's Richardson extrapolation).
  expect_reference(unname(hessian), unname(numDeriv::hessian(objective, theta)), reference = "numDeriv",
                   covers = c("lpa()", "multilpa()"),
                   quantity = "negative-binomial Hessian in log mean and log dispersion",
                   tolerance = 1e-7, measure = "scaled",
                   note = paste(count_note, "Reference: numDeriv::hessian() of a dnbinom() likelihood written in the test."))
  solution <- latents:::.latents_negative_binomial_maximize(y, weight, c(3, 5), c(.4, .7))
  reference <- stats::optim(theta, objective, method = "BFGS",
                            control = list(fnscale = -1, reltol = 1e-15, maxit = 2000))
  expect_gte(objective(log(c(solution$means, solution$dispersion))),
             reference$value - 1e-10)
  # Under-dispersed counts put the dispersion on its Poisson boundary.
  under <- latents:::.latents_negative_binomial_maximize(rep(c(3, 4, 5), 30), matrix(1, 90, 1), 4, .3)
  expect_equal(under$dispersion, latents:::.latents_min_dispersion)
  expect_equal(under$means, 4)
})

test_that("ordinal and count scores are the gradient of the observed log likelihood", {
  skip_without_reference("numDeriv")
  withr::local_seed(8)
  n <- 120
  data <- data.frame(o = sample(1:4, n, TRUE), y = stats::rpois(n, 3),
                     z = stats::rnbinom(n, size = 2, mu = 4))
  data$o[3] <- NA
  data$y[5] <- NA
  prior <- matrix(log(c(.4, .6)), n, 2, byrow = TRUE)
  joint <- function(theta, extra) {
    latents:::.latents_extra_log_density(extra, latents:::.latents_extra_decode(theta, extra, 2), n, 2) + prior
  }
  log_likelihood <- function(theta, extra) {
    value <- joint(theta, extra)
    top <- apply(value, 1L, max)
    sum(top + log(rowSums(exp(value - top))))
  }
  checks <- lapply(list(c("poisson", "varying"), c("negative_binomial", "varying"),
                        c("negative_binomial", "equal")), function(model) {
    extra <- latents:::.latents_set_count_model(
      latents:::.latents_prepare_extra(data, "o", c("y", "z"), missing = "fiml"), model[1], model[2])
    theta <- stats::rnorm(latents:::.latents_extra_n_parameters(extra, 2), 0, .3)
    value <- joint(theta, extra)
    posterior <- exp(value - apply(value, 1L, max))
    posterior <- posterior / rowSums(posterior)
    analytic <- drop(latents:::.latents_extra_scores(extra, posterior,
                                                     latents:::.latents_extra_decode(theta, extra, 2)))
    data.frame(model = paste(model, collapse = ", "), coordinate = seq_along(theta),
               analytic = analytic,
               numerical = numDeriv::grad(log_likelihood, theta, extra = extra))
  })
  scores <- do.call(rbind, checks)
  # Scaled: score entries can sit near zero, where numDeriv's error is
  # absolute. 1e-7, the original test's tolerance.
  expect_reference(data.frame(model = scores$model, coordinate = scores$coordinate, score = scores$analytic),
                   data.frame(model = scores$model, coordinate = scores$coordinate, score = scores$numerical),
                   reference = "numDeriv", covers = c("lpa()", "multilpa()"),
                   quantity = "ordinal and count scores of the observed log likelihood",
                   tolerance = 1e-7, measure = "scaled", key = c("model", "coordinate"),
                   note = "Package side: the internal analytic scores (.latents_extra_scores) behind the standard errors of lpa()/multilpa(ordinal =, count =); reference: numDeriv::grad() of the package's observed log likelihood.")
  # (K - 1) + (C - 1) ordinal, C means per count, then C or 1 dispersions.
  expect_identical(vapply(checks, nrow, integer(1)), c(8L, 12L, 10L))
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
  # The joint density written in base R from the fit's estimates; the ordinal
  # category probabilities from their definition above.
  ordinal <- categorical_ordinal_reference(fit$ordinal_intercepts[[1]],
                                           fit$ordinal_locations[, 1])
  joint <- vapply(1:2, function(k) {
    log(fit$profile_probabilities[k]) +
      stats::dnorm(data$g, fit$means[k, "g"], sqrt(fit$variances[k, "g"]), log = TRUE) +
      log(fit$response_probabilities$c[k, data$c + 1]) + ordinal[k, data$o] +
      stats::dpois(data$k, fit$count_means[k, 1], log = TRUE)
  }, numeric(n))
  top <- apply(joint, 1L, max)
  posterior <- exp(joint - top) / rowSums(exp(joint - top))
  reference <- "independent: mixed-indicator mixture density, base R (dnorm, dpois, adjacent-category logit)"
  # The same sums of logs computed twice; the smallest posterior is about 4e-7.
  expect_reference(unname(fit$subject_posteriors), unname(posterior), reference = reference,
                   covers = "lpa()", quantity = "posterior class probabilities, mixed indicators",
                   tolerance = 1e-12, measure = "relative")
  expect_reference(as.numeric(logLik(fit)), sum(top + log(rowSums(exp(joint - top)))),
                   reference = reference, covers = c("lpa()", "logLik.multilpa()"),
                   quantity = "log likelihood, mixed indicators",
                   tolerance = 1e-12, measure = "relative")
  # 1 class weight, 2 x 2 Gaussian, 2 x 1 binary, (3 - 1) + (2 - 1) ordinal, 2 means.
  expect_identical(as.integer(fit$n_parameters), 12L)
})
