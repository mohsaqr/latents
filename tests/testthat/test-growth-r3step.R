# r3step() for trajectory models (latent class growth and growth mixtures):
# the person-level membership regression of Vermunt (2010) with the step-one
# classification errors held fixed. References: the same likelihood written
# out independently here, maximized by optim; equivariance in the covariate.

r3step_fixture <- function(n = 240L, seed = 8L) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(seed)
  times <- 0:4
  x <- stats::rnorm(n)
  class <- 1L + (stats::runif(n) < stats::plogis(-0.3 + 1.0 * x))
  out <- data.frame(id = rep(seq_len(n), each = length(times)), time = rep(times, n),
                    x = rep(x, each = length(times)),
                    w = stats::rnorm(n * length(times)))
  out$y <- c(2, 3.5)[class[out$id]] + c(0.6, -0.4)[class[out$id]] * out$time +
    rep(stats::rnorm(n, 0, 0.8), each = length(times)) +
    stats::rnorm(nrow(out), 0, 0.9)
  out
}

r3step_fit <- function(data, ...) {
  mixture_regression(y ~ time, data, n_classes = 2, id = "id", class_level = "group",
                     random = ~ 1, n_starts = 2, seed = 1, vcov_type = "none", ...)
}

test_that("r3step on a growth mixture maximizes the fixed-error likelihood", {
  data <- r3step_fixture()
  fit <- r3step_fit(data)
  ours <- r3step(fit, data, "x", vcov_type = "observed")
  expect_identical(ours$outcome, c("class_1", "class_1"))
  expect_identical(ours$level, c("persons", "persons"))
  # Written out from Vermunt (2010): the modal class m_i, error matrix
  # D[k, j] = sum_i p_ik [m_i = j] / sum_i p_ik, and the log likelihood
  # sum_i log sum_k pi_k(x_i) D[k, m_i], logits against the last class.
  posterior <- predict(fit, type = "posterior")
  p <- as.matrix(posterior[c("posterior_class_1", "posterior_class_2")])
  modal <- max.col(p, ties.method = "first")
  errors <- t(vapply(1:2, function(k) {
    vapply(1:2, function(j) sum(p[modal == j, k]), numeric(1)) / sum(p[, k])
  }, numeric(2)))
  x <- data$x[!duplicated(data$id)]
  negative_log_likelihood <- function(theta) {
    eta <- theta[1L] + theta[2L] * x
    prior <- cbind(stats::plogis(eta), 1 - stats::plogis(eta))
    -sum(log(rowSums(prior * t(errors)[modal, ])))
  }
  reference <- stats::optim(c(0, 0), negative_log_likelihood, method = "BFGS",
                            control = list(reltol = 1e-14), hessian = TRUE)
  expect_equal(ours$estimate, reference$par, tolerance = 1e-5)
  expect_equal(ours$standard_error, sqrt(diag(solve(reference$hessian))),
               tolerance = 1e-3)
})

test_that("r3step slopes are equivariant in the covariate's units and origin", {
  data <- r3step_fixture()
  fit <- r3step_fit(data)
  base <- r3step(fit, data, "x")
  rescaled <- r3step(fit, transform(data, x = 10 + 4 * x), "x")
  expect_equal(rescaled$estimate[2L], base$estimate[2L] / 4, tolerance = 1e-6)
  expect_equal(rescaled$standard_error[2L], base$standard_error[2L] / 4, tolerance = 1e-5)
  expect_equal(rescaled$p_value[2L], base$p_value[2L], tolerance = 1e-5)
})

test_that("r3step works for latent class growth (no random effects)", {
  data <- r3step_fixture()
  fit <- mixture_regression(y ~ time, data, n_classes = 2, id = "id",
                            class_level = "group", n_starts = 2, seed = 1)
  result <- r3step(fit, data, "x")
  expect_identical(nrow(result), 2L)
  expect_true(all(is.finite(result$standard_error)))
})

test_that("r3step refuses what a trajectory model does not define", {
  data <- r3step_fixture(n = 120L)
  fit <- r3step_fit(data)
  expect_error(r3step(fit, data, "x", level = "groups"), class = "latents_bad_argument")
  expect_error(r3step(fit, data, "x", by_group_class = TRUE), class = "latents_bad_argument")
  expect_error(r3step(fit, data, "w"), class = "latents_bad_covariate")
  expect_error(r3step(fit, data, "time"), class = "latents_bad_covariate")
  expect_error(r3step(fit, data, c("x", "x")), class = "latents_bad_covariate")
  expect_error(r3step(fit, subset(data, id > 3), "x"),
               class = "latents_bad_inference_data")
  with_membership <- r3step_fit(data, membership = ~ x)
  expect_error(r3step(with_membership, data, "x"), class = "latents_unsupported_three_step")
})
