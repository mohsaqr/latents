# External lme4 comparison retained outside the shipped package tests.
# Same fixture and assertions as the original round-two growth regression.

calendar_reference_data <- function(n = 90L, seed = 3L) {
  withr::with_seed(seed, {
    class <- stats::rbinom(n, 1L, 0.45) + 1L
    covariance <- list(matrix(c(1, 0.3, 0.3, 0.25), 2L),
                       matrix(c(0.6, -0.1, -0.1, 0.15), 2L))
    effects <- t(vapply(seq_len(n), function(i) {
      as.vector(t(chol(covariance[[class[i]]])) %*% stats::rnorm(2L))
    }, numeric(2L)))
    out <- data.frame(id = rep(seq_len(n), each = 5L), time = rep(0:4, n))
    out$x <- stats::rnorm(n)[out$id]
    person_class <- class[out$id]
    out$y <- c(2, 4.5)[person_class] + c(0.5, -0.3)[person_class] * out$time +
      effects[out$id, 1L] + effects[out$id, 2L] * out$time +
      stats::rnorm(nrow(out), 0, c(0.7, 0.5)[person_class])
    out[stats::runif(nrow(out)) > 0.15, ]
  })
}

test_that("a one-class calendar-year random slope matches lme4", {
  skip_if_not_installed("lme4")
  data <- calendar_reference_data()
  data$year <- data$time + 2000
  fit <- suppressWarnings(mixture_regression(
    y ~ year, data, n_classes = 1L, id = "id", class_level = "group",
    random = ~ 1 + year, n_starts = 1L, seed = 1L))
  # lme4 on time from zero: the same model, which it fits without warnings.
  reference <- lme4::lmer(y ~ time + (1 + time | id), data = data, REML = FALSE,
                          control = lme4::lmerControl(
                            optCtrl = list(ftol_abs = 1e-15, xtol_abs = 1e-15)))
  expect_true(fit$converged)
  expect_length(fit$degenerate_classes, 0L)
  expect_lt(abs(fit$log_likelihood - as.numeric(stats::logLik(reference))), 1e-6)
  fixed <- lme4::fixef(reference)
  coefficients <- as.data.frame(get_results(fit, "coefficients"))
  expect_equal(coefficients$estimate,
               unname(c(fixed[1L] - 2000 * fixed[2L], fixed[2L])), tolerance = 1e-5)
})
