# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Growth mixture models: mixture_regression(random = ).

growth_fixture_data <- function(n = 160L, seed = 21L) {
  set.seed(seed)
  times <- 0:5
  class <- sample(1:2, n, TRUE, prob = c(0.6, 0.4))
  covariance <- matrix(c(1, 0.2, 0.2, 0.1), 2L)
  effects <- t(chol(covariance)) %*% matrix(stats::rnorm(2L * n), 2L)
  out <- data.frame(id = rep(seq_len(n), each = length(times)),
                    time = rep(times, n), x = stats::rnorm(n * length(times)),
                    age = rep(stats::rnorm(n), each = length(times)))
  person_class <- class[out$id]
  out$y <- c(2, 6)[person_class] + c(0.5, -0.4)[person_class] * out$time +
    0.3 * out$x + effects[1L, out$id] + effects[2L, out$id] * out$time +
    stats::rnorm(nrow(out), 0, c(0.8, 0.6)[person_class])
  out[stats::runif(nrow(out)) > 0.15, ]
}

growth_fit <- function(data, ...) {
  mixture_regression(y ~ time + x, data, n_classes = 2, id = "id",
                     class_level = "group", common = ~ x, membership = ~ age,
                     random = ~ 1 + time, n_starts = 2, seed = 1, ...)
}

# The random-effect covariance of each class from the estimates, as the
# parametrisation defines it: one matrix per class ("varying"), one shared
# ("equal"), or one shared and scaled by a class factor ("proportional").
growth_class_covariances <- function(fit, params) {
  lapply(seq_len(fit$spec$n_classes), function(k) {
    switch(fit$spec$random_covariance,
           varying = params$random_covariance[[k]],
           equal = params$random_covariance[[1L]],
           proportional = params$random_covariance[[1L]] * params$random_scale[k]^2)
  })
}

# An independent likelihood: each person's marginal covariance built densely,
# the class shares a softmax of (1, age) times the membership logits.
growth_dense_log_likelihood <- function(fit, data, params) {
  covariances <- growth_class_covariances(fit, params)
  persons <- split(data, data$id)
  sum(vapply(persons, function(rows) {
    logits <- drop(c(1, rows$age[1L]) %*% params$gamma)
    log_prior <- logits - max(logits) - log(sum(exp(logits - max(logits))))
    x <- cbind(1, rows$time, rows$x)
    z <- cbind(1, rows$time)
    densities <- vapply(seq_len(fit$spec$n_classes), function(k) {
      v <- z %*% covariances[[k]] %*% t(z) + params$sigma2[k] * diag(nrow(rows))
      mu <- x %*% c(params$beta[, k], params$common)
      upper <- chol(v)
      r <- backsolve(upper, rows$y - mu, transpose = TRUE)
      -0.5 * (nrow(rows) * log(2 * pi) + 2 * sum(log(diag(upper))) + sum(r^2))
    }, numeric(1))
    top <- max(densities + log_prior)
    top + log(sum(exp(densities + log_prior - top)))
  }, numeric(1)))
}

test_that("analytic growth scores equal numerical derivatives", {
  skip_on_cran()
  skip_without_reference("numDeriv")
  data <- growth_fixture_data()
  configurations <- list(list(covariance = "varying", variance = "varying", diagonal = FALSE),
                         list(covariance = "equal", variance = "equal", diagonal = FALSE),
                         list(covariance = "proportional", variance = "equal", diagonal = FALSE),
                         list(covariance = "varying", variance = "equal", diagonal = TRUE))
  invisible(lapply(configurations, function(cfg) {
    fit <- growth_fit(data, random_covariance = cfg$covariance, variance = cfg$variance,
                      random_diagonal = cfg$diagonal, vcov_type = "none")
    set.seed(5)
    theta <- latents:::.growth_pack(fit$spec, fit$params)
    theta <- theta + stats::rnorm(length(theta), 0, 0.05)
    log_likelihood <- function(v) {
      latents:::.growth_expectation(fit$spec, fit$stats,
                                    latents:::.growth_unpack(fit$spec, v, fit$params))$log_likelihood
    }
    analytic <- colSums(latents:::.growth_scores(
      fit$spec, fit$stats, latents:::.growth_unpack(fit$spec, theta, fit$params)))
    numerical <- numDeriv::grad(log_likelihood, theta)
    # The original bound, kept: an absolute 1e-4 on the scale of the largest
    # gradient entry, since numDeriv's error at a point off the optimum grows
    # with the size of the gradient, not with each entry.
    expect_reference(unname(analytic), numerical, reference = "numDeriv",
                     covers = "mixture_regression()",
                     quantity = sprintf("score at a perturbed point, covariance %s, variance %s, diagonal %s",
                                        cfg$covariance, cfg$variance, cfg$diagonal),
                     tolerance = 1e-4 * (1 + max(abs(numerical))), measure = "absolute",
                     note = "Package side: the internal analytic growth scores (.growth_scores) behind vcov() of mixture_regression(random = ); reference: numDeriv::grad() of the package's likelihood.")
  }))
})

test_that("observed standard errors equal the inverse numerical Hessian", {
  skip_without_reference("numDeriv")
  data <- growth_fixture_data(n = 120L)
  fit <- mixture_regression(y ~ time, data, n_classes = 2, id = "id",
                            class_level = "group", random = ~ 1,
                            random_covariance = "equal", variance = "equal",
                            n_starts = 2, seed = 1, tol = 1e-12)
  theta <- coef(fit)
  log_likelihood <- function(v) {
    latents:::.growth_expectation(fit$spec, fit$stats,
                                  latents:::.growth_unpack(fit$spec, v, fit$params))$log_likelihood
  }
  numerical <- solve(-numDeriv::hessian(log_likelihood, theta))
  # 1e-4, the original tolerance: second differences of a sum over persons.
  expect_reference(sqrt(diag(vcov(fit))), sqrt(diag(numerical)), reference = "numDeriv",
                   covers = c("mixture_regression()", "vcov.latents_growth_mixture()"),
                   quantity = "standard errors, random intercept, equal covariance", tolerance = 1e-4,
                   measure = "relative",
                   note = "Reference: the inverse of numDeriv::hessian() of the package's likelihood at coef().")
})


test_that("the growth mixture likelihood equals an independent dense likelihood", {
  skip_on_cran()
  data <- growth_fixture_data()
  invisible(lapply(c("varying", "equal", "proportional"), function(covariance) {
    fit <- growth_fit(data, random_covariance = covariance, vcov_type = "none")
    expect_s3_class(fit, "latents_growth_mixture")
    # The same likelihood summed in another order: 1e-10, the original tolerance.
    expect_reference(as.numeric(logLik(fit)), growth_dense_log_likelihood(fit, data, fit$params),
                     reference = "independent: dense marginal growth mixture likelihood, base R",
                     covers = c("mixture_regression()", "logLik.latents_growth_mixture()"),
                     quantity = sprintf("log likelihood, %s covariance", covariance),
                     tolerance = 1e-10, measure = "relative")
    # Away from the maximum too.
    set.seed(3)
    theta <- latents:::.growth_pack(fit$spec, fit$params)
    moved <- latents:::.growth_unpack(fit$spec, theta + stats::rnorm(length(theta), 0, 0.1),
                                      fit$params)
    expect_reference(latents:::.growth_expectation(fit$spec, fit$stats, moved)$log_likelihood,
                     growth_dense_log_likelihood(fit, data, moved),
                     reference = "independent: dense marginal growth mixture likelihood, base R",
                     covers = "mixture_regression()",
                     quantity = sprintf("log likelihood off the maximum, %s covariance", covariance),
                     tolerance = 1e-10, measure = "relative",
                     note = "Package side: the internal likelihood (.growth_expectation) that mixture_regression(random = ) maximises, at perturbed parameters.")
  }))
})
