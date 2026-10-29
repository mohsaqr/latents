# Comparisons with independent references moved from tests/testthat/test-growth-correlated.R,
# so that they run in CI's equivalence job rather than under R CMD check.

correlated_fixture <- function(n = 120L, seed = 8L, process_sd = 0.8) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(seed)
  times <- 0:5
  class <- sample(1:2, n, TRUE)
  out <- data.frame(id = rep(seq_len(n), each = length(times)), time = rep(times, n),
                    age = rep(stats::rnorm(n), each = length(times)))
  kernel <- t(chol(process_sd^2 * exp(-0.7 * abs(outer(times, times, "-"))) +
                     diag(1e-12, length(times))))
  process <- unlist(lapply(seq_len(n), function(i) {
    as.vector(kernel %*% stats::rnorm(length(times)))
  }))
  out$y <- c(1, 4)[class[out$id]] + c(0.5, -0.3)[class[out$id]] * out$time +
    rep(stats::rnorm(n), each = length(times)) + process +
    stats::rnorm(nrow(out), 0, 0.6)
  out[stats::runif(nrow(out)) > 0.1, ]
}

correlated_fit <- function(data, kind = "ar", ...) {
  quietly(mixture_regression(y ~ time, data, n_classes = 2, id = "id",
                             class_level = "group", random = ~ 1, membership = ~ age,
                             residual_correlation = kind, n_starts = 2, seed = 1,
                             vcov_type = "none", ...), "latents_random_boundary")
}

dense_log_likelihood <- function(fit, data) {
  params <- fit$params
  process <- function(t) {
    if (identical(fit$spec$residual_correlation, "ar")) {
      params$correlation[["sd"]]^2 * exp(-params$correlation[["rate"]] *
                                            abs(outer(t, t, "-")))
    } else params$correlation[["sd"]]^2 * outer(t, t, pmin)
  }
  sum(vapply(split(data, data$id), function(person) {
    eta <- c(0, sum(params$gamma[, 2L] * c(1, person$age[1L])))
    prior <- exp(eta - max(eta)) / sum(exp(eta - max(eta)))
    log_f <- vapply(1:2, function(k) {
      mean <- params$beta[1L, k] + params$beta[2L, k] * person$time
      v <- params$random_covariance[[k]][1L, 1L] + process(person$time) +
        diag(params$sigma2[k], nrow(person))
      r <- person$y - mean
      -0.5 * (nrow(person) * log(2 * pi) +
                as.numeric(determinant(v, logarithm = TRUE)$modulus) + sum(r * solve(v, r)))
    }, numeric(1))
    top <- max(log(prior) + log_f)
    top + log(sum(exp(log(prior) + log_f - top)))
  }, numeric(1)))
}

test_that("the correlated likelihood equals a dense computation", {
  data <- correlated_fixture()
  for (kind in c("ar", "brownian")) {  # two settings of one check
    fit <- correlated_fit(data, kind)
    expect_equal(fit$log_likelihood, dense_log_likelihood(fit, data), tolerance = 1e-9)
    expect_identical(fit$n_parameters,
                     correlated_fit(data, "ar")$n_parameters - (kind == "brownian"))
  }
  plain <- mixture_regression(y ~ time, data, n_classes = 2, id = "id",
                              class_level = "group", random = ~ 1, membership = ~ age,
                              n_starts = 2, seed = 1, vcov_type = "none")
  expect_identical(correlated_fit(data, "ar")$n_parameters, plain$n_parameters + 2L)
  expect_gte(correlated_fit(data, "ar")$log_likelihood, plain$log_likelihood - 1e-6)
})

test_that("conditional predictions add the empirical Bayes effects with the process", {
  data <- correlated_fixture()
  fit <- correlated_fit(data)
  new <- subset(data, id <= 4)
  ours <- predict(fit, new, type = "class_conditional")
  params <- fit$params
  dense <- unlist(lapply(1:2, function(k) {
    unlist(lapply(split(new, new$id), function(person) {
      t <- person$time
      g <- params$random_covariance[[k]][1L, 1L]
      v <- g + params$correlation[["sd"]]^2 * exp(-params$correlation[["rate"]] *
                                                     abs(outer(t, t, "-"))) +
        diag(params$sigma2[k], nrow(person))
      mean <- params$beta[1L, k] + params$beta[2L, k] * t
      mean + g * sum(solve(v, person$y - mean))
    }), use.names = FALSE)
  }))
  expect_equal(ours$fitted, dense, tolerance = 1e-8)
  expect_identical(nrow(simulate(fit, nsim = 2, seed = 3)), nrow(data))
  expect_equal(simulate(fit, nsim = 1, seed = 3), simulate(fit, nsim = 1, seed = 3))
  expect_error(predict(fit, subset(new, select = -time), type = "conditional"),
               class = "latents_bad_data")
})
