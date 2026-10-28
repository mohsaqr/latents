# predict() for growth mixture models. The lcmm comparison is in
# equivalence/test-growth-predict-lcmm.R; these tests use an independent dense
# implementation (full marginal covariance per person) and invariants.

# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

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

test_that("fitted and re-supplied persons give the fit's own posteriors", {
  data <- predict_fixture()
  fit <- predict_fit(data)
  fitted <- predict(fit, type = "posterior")
  expect_equal(fitted, as.data.frame(get_results(fit, "assignments")), ignore_attr = TRUE)
  # The same persons passed as new data are classified from scratch.
  expect_equal(predict(fit, data, type = "posterior"), fitted, tolerance = 1e-10)
  expect_equal(predict(fit, data, type = "class_conditional"),
               predict(fit, type = "class_conditional"), tolerance = 1e-10)
})

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

test_that("predictions satisfy their mixing identities", {
  data <- predict_fixture()
  fit <- predict_fit(data)
  new <- subset(data, id <= 10)
  by_class <- predict(fit, new, type = "class_response")
  expect_equal(as.vector(tapply(by_class$prior, by_class$row, sum)),
               rep(1, nrow(new)), tolerance = 1e-12)
  mixed <- as.vector(tapply(by_class$prior * by_class$fitted, by_class$row, sum))
  expect_equal(predict(fit, new, type = "response")$fitted, mixed, tolerance = 1e-12)
  conditional <- predict(fit, new, type = "class_conditional")
  expect_equal(predict(fit, new, type = "conditional")$fitted,
               as.vector(tapply(conditional$posterior * conditional$fitted,
                                conditional$row, sum)), tolerance = 1e-12)
  effects <- predict(fit, new, type = "random_effects")
  expect_named(effects, c("id", "class", "posterior", "(Intercept)", "time"))
  expect_identical(nrow(effects), 2L * length(unique(new$id)))
})

test_that("a person's predictions do not depend on row order or other persons", {
  data <- predict_fixture()
  fit <- predict_fit(data)
  new <- subset(data, id <= 8)
  shuffled <- new[rev(seq_len(nrow(new))), ]
  one <- predict(fit, new, type = "posterior")
  two <- predict(fit, shuffled, type = "posterior")
  expect_equal(two, one, tolerance = 1e-12)
  alone <- predict(fit, subset(data, id == 3), type = "posterior")
  expect_equal(alone$posterior_class_1, subset(one, id == 3)$posterior_class_1,
               tolerance = 1e-12)
})

test_that("trajectory predictions need no outcome; posteriors do", {
  data <- predict_fixture()
  fit <- predict_fit(data)
  grid <- data.frame(id = rep(1:2, each = 3), time = rep(c(0, 2, 4), 2),
                     age = rep(c(-1, 1), each = 3))
  expect_identical(nrow(predict(fit, grid, type = "response")), 6L)
  expect_error(predict(fit, grid, type = "posterior"), class = "latents_bad_data")
  expect_error(predict(fit, subset(grid, select = -id), type = "response"),
               class = "latents_bad_data")
  expect_error(predict(fit, transform(grid, time = NA), type = "response"),
               class = "latents_bad_data")
  expect_warning(predict(fit, grid, kind = "x"), class = "latents_unused_argument")
})

test_that("new clusters of a multilevel growth mixture are classified", {
  skip_on_cran()
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(3)
  schools <- do.call(rbind, lapply(seq_len(14), function(s) {
    improving <- if (s %% 2L == 1L) 0.85 else 0.2
    do.call(rbind, lapply(seq_len(8), function(i) {
      k <- if (stats::runif(1) < improving) 1L else 2L
      data.frame(school = s, student = (s - 1) * 8 + i, wave = 0:3,
                 score = (if (k == 1L) 50 + 3 * 0:3 else rep(55, 4)) +
                   stats::rnorm(1, 0, 2) + stats::rnorm(4, 0, 1.5))
    }))
  }))
  fit <- mixture_regression(score ~ wave, schools, n_classes = 2, id = "student",
                            class_level = "group", random = ~ 1,
                            n_group_classes = 2, cluster = "school", n_starts = 2,
                            seed = 1, vcov_type = "none")
  expect_equal(predict(fit, schools, type = "posterior"),
               predict(fit, type = "posterior"), tolerance = 1e-10)
  one_school <- predict(fit, subset(schools, school == 2), type = "posterior")
  expect_identical(nrow(one_school), 8L)
  expect_equal(rowSums(one_school[c("posterior_class_1", "posterior_class_2")]),
               rep(1, 8), tolerance = 1e-12)
  expect_error(predict(fit, subset(schools, select = -school), type = "posterior"),
               class = "latents_bad_data")
})
