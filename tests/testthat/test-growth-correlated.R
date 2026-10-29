# Growth mixtures with correlated residuals: mixture_regression(
# residual_correlation = "ar" | "brownian"). The lcmm::hlme comparison is in
# equivalence/test-growth-correlated-lcmm.R; here, a dense likelihood written
# out independently, scores against numerical gradients, the empirical Bayes
# formula with the process, a boundary case and the refusals.
# The comparisons with independent references are in tests/equivalence/test-growth-correlated.R.

# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

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
