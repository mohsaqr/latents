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

# An independent likelihood: each person's marginal covariance built densely.
growth_dense_log_likelihood <- function(fit, data, params) {
  spec <- fit$spec
  covariances <- latents:::.growth_covariances(spec, params)
  prior <- exp(latents:::.mixture_log_softmax(spec$w, params$gamma))
  persons <- split(seq_len(nrow(data)), match(data$id, spec$group_levels))
  sum(vapply(seq_along(persons), function(i) {
    rows <- data[persons[[i]], ]
    x <- cbind(1, rows$time, rows$x)
    z <- cbind(1, rows$time)
    densities <- vapply(seq_len(spec$n_classes), function(k) {
      v <- z %*% covariances[[k]] %*% t(z) + params$sigma2[k] * diag(nrow(rows))
      mu <- x %*% c(params$beta[, k], params$common)
      upper <- chol(v)
      r <- backsolve(upper, rows$y - mu, transpose = TRUE)
      -0.5 * (nrow(rows) * log(2 * pi) + 2 * sum(log(diag(upper))) + sum(r^2))
    }, numeric(1))
    top <- max(densities + log(prior[i, ]))
    top + log(sum(exp(densities + log(prior[i, ]) - top)))
  }, numeric(1)))
}

test_that("analytic growth scores equal numerical derivatives", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
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
    expect_lt(max(abs(analytic - numerical)), 1e-4 * (1 + max(abs(numerical))))
  }))
})

test_that("observed standard errors equal the inverse numerical Hessian", {
  skip_if_not_installed("numDeriv")
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
  expect_equal(unname(sqrt(diag(vcov(fit)))), sqrt(diag(numerical)), tolerance = 1e-4)
})
