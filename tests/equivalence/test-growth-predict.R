# Comparisons with independent references moved from tests/testthat/test-growth-predict.R,
# so that they run in CI's equivalence job rather than under R CMD check.

predict_fixture <- function(n = 80L, seed = 4L) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(seed)
  times <- 0:4
  out <- data.frame(id = rep(seq_len(n), each = length(times)),
                    time = rep(times, n), age = rep(stats::rnorm(n), each = length(times)))
  class <- sample(1:2, n, TRUE)
  out$y <- c(1, 4)[class[out$id]] + c(0.5, -0.3)[class[out$id]] * out$time +
    rep(stats::rnorm(n), each = length(times)) +
    rep(stats::rnorm(n, 0, 0.3), each = length(times)) * out$time +
    stats::rnorm(nrow(out), 0, 0.6)
  out[stats::runif(nrow(out)) > 0.1, ]
}

predict_fit <- function(data) {
  mixture_regression(y ~ time, data, n_classes = 2, id = "id", class_level = "group",
                     membership = ~ age, random = ~ 1 + time, n_starts = 2, seed = 1,
                     vcov_type = "none")
}

test_that("conditional curves match a dense empirical Bayes computation", {
  data <- predict_fixture()
  fit <- predict_fit(data)
  new <- subset(data, id <= 6)
  ours <- predict(fit, new, type = "class_conditional")
  params <- fit$params
  # b_ik = G_k Z_i' V_ik^-1 (y_i - X_i beta_k), V_ik = Z_i G_k Z_i' + s_k^2 I.
  dense <- unlist(lapply(seq_len(2L), function(k) {
    covariance <- params$random_covariance[[k]]
    unlist(lapply(split(new, new$id), function(person) {
      z <- cbind(1, person$time)
      mean <- params$beta[1L, k] + params$beta[2L, k] * person$time
      v <- z %*% covariance %*% t(z) + diag(params$sigma2[k], nrow(person))
      effects <- covariance %*% t(z) %*% solve(v, person$y - mean)
      as.vector(mean + z %*% effects)
    }), use.names = FALSE)
  }))
  expect_equal(ours$fitted, dense, tolerance = 1e-8)
  # Each person's dense posterior under the membership prior.
  dense_posterior <- t(vapply(split(new, new$id), function(person) {
    z <- cbind(1, person$time)
    prior <- exp(c(0, sum(params$gamma[, 2L] * c(1, person$age[1L]))))
    log_joint <- log(prior / sum(prior)) + vapply(seq_len(2L), function(k) {
      covariance <- params$random_covariance[[k]]
      mean <- params$beta[1L, k] + params$beta[2L, k] * person$time
      v <- z %*% covariance %*% t(z) + diag(params$sigma2[k], nrow(person))
      r <- person$y - mean
      -0.5 * (nrow(person) * log(2 * pi) +
                as.numeric(determinant(v, logarithm = TRUE)$modulus) +
                sum(r * solve(v, r)))
    }, numeric(1))
    exp(log_joint - max(log_joint)) / sum(exp(log_joint - max(log_joint)))
  }, numeric(2)))
  posterior <- predict(fit, new, type = "posterior")
  expect_equal(unname(as.matrix(posterior[c("posterior_class_1", "posterior_class_2")])),
               unname(dense_posterior), tolerance = 1e-8)
})
