# Comparisons with independent references moved from tests/testthat/test-astra-core.R,
# so that they run in CI's equivalence job rather than under R CMD check.

test_that("constant raw outcomes can have identified regression residuals (A03)", {
  .latents_local_seed(1027)
  data <- data.frame(y = rep(2, 60), x = seq(-2, 3, length.out = 60),
                     o = stats::rnorm(60), w = 1)
  invisible(lapply(list(y ~ 0 + x, y ~ 1 + offset(o)), function(formula) {
    reference <- stats::lm(formula, data)
    invisible(lapply(c(FALSE, TRUE), function(weighted) {
      fit <- mixture_regression(formula, data, 1, n_starts = 0, vcov_type = "none",
                                 weights = if (weighted) "w" else NULL)
      expect_true(fit$converged)
      expect_equal(as.numeric(fit$params$beta), unname(stats::coef(reference)),
                    tolerance = 1e-10)
      expect_equal(unname(fit$params$sigma2), mean(stats::residuals(reference)^2),
                    tolerance = 1e-10)
      expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(reference)),
                    tolerance = 1e-10)
    }))
  }))
  expect_error(mixture_regression(y ~ 1, data, 1, n_starts = 0),
                class = "latents_bad_data")
  # The residual-scale check must follow the public empty-sample refusal.
  expect_error(mixture_regression(y ~ x, data[0, ], 1, n_starts = 0),
                class = "latents_bad_data")
  # Varying excluded values cannot hide a constant carrying response.
  extra <- data[1:2, ]; extra$y <- c(-3, 4); extra$w <- 0
  expect_error(mixture_regression(y ~ 1, rbind(data, extra), 1,
                                   weights = "w", n_starts = 0),
                class = "latents_bad_data")
})
