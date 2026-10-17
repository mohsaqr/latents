# Growth mixtures with correlated residuals: mixture_regression(
# residual_correlation = "ar" | "brownian"). The lcmm::hlme comparison is in
# equivalence/test-growth-correlated-lcmm.R; here, a dense likelihood written
# out independently, scores against numerical gradients, the empirical Bayes
# formula with the process, a boundary case and the refusals.

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

test_that("analytic scores equal the likelihood's gradient", {
  data <- correlated_fixture()
  for (kind in c("ar", "brownian")) {  # two settings of one check
    fit <- correlated_fit(data, kind)
    theta <- latents:::.growth_pack(fit$spec, fit$params)
    shifted <- theta + 0.02 * cos(seq_along(theta))
    log_likelihood <- function(v) {
      latents:::.growth_expectation(fit$spec, fit$stats,
                                    latents:::.growth_unpack(fit$spec, v, fit$params))$log_likelihood
    }
    numeric_gradient <- vapply(seq_along(shifted), function(j) {
      h <- 1e-5 * max(1, abs(shifted[j]))
      (log_likelihood(replace(shifted, j, shifted[j] + h)) -
         log_likelihood(replace(shifted, j, shifted[j] - h))) / (2 * h)
    }, numeric(1))
    analytic <- colSums(latents:::.growth_scores(
      fit$spec, fit$stats, latents:::.growth_unpack(fit$spec, shifted, fit$params)))
    expect_equal(unname(analytic), numeric_gradient, tolerance = 1e-6)
  }
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

test_that("correlated growth inference reports the process", {
  data <- correlated_fixture()
  fit <- correlated_fit(data)
  table <- parameter_inference(fit, boundary = "fix")
  process <- subset(table, block == "residual_process")
  expect_identical(process$term, c("rate", "sd"))
  expect_equal(process$estimate, unname(fit$params$correlation), tolerance = 1e-12)
  random <- as.data.frame(get_results(fit, "random"))
  expect_identical(subset(random, startsWith(parameter, "process_"))$parameter,
                   c("process_rate", "process_sd"))
  expect_output(print(fit), "Correlated residuals: continuous-time autoregressive")
})

test_that("white noise residuals leave the process on its floor, reported", {
  skip_on_cran()
  data <- correlated_fixture(process_sd = 0)
  expect_warning(
    fit <- mixture_regression(y ~ time, data, n_classes = 2, id = "id",
                              class_level = "group", random = ~ 1, membership = ~ age,
                              residual_correlation = "ar", n_starts = 2, seed = 1,
                              vcov_type = "none"),
    class = "latents_random_boundary")
  expect_true(fit$converged)
  expect_true(length(fit$residual_boundary) > 0L)
  expect_error(parameter_inference(fit), class = "latents_boundary_fit")
  held <- parameter_inference(fit, boundary = "fix")
  expect_true(all(is.finite(subset(held, block == "coefficient")$standard_error)))
})

test_that("correlated residuals refuse what they do not define", {
  data <- correlated_fixture(n = 40L)
  expect_error(mixture_regression(y ~ time, data, n_classes = 2, id = "id",
                                  class_level = "group", residual_correlation = "ar"),
               class = "latents_bad_argument")
  expect_error(mixture_regression(y ~ time, data, n_classes = 2, id = "id",
                                  class_level = "group", random = ~ 1,
                                  residual_time = "time"),
               class = "latents_bad_argument")
  expect_error(mixture_regression(y ~ time, data, n_classes = 2, id = "id",
                                  class_level = "group", random = ~ 1,
                                  residual_correlation = "ar", residual_time = "absent"),
               class = "latents_bad_argument")
  expect_error(mixture_regression(y ~ time, transform(data, time = time - 2),
                                  n_classes = 2, id = "id", class_level = "group",
                                  random = ~ 1, residual_correlation = "brownian"),
               class = "latents_bad_data")
})
