# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Comparisons with independent references moved from tests/testthat/test-astra-core.R,
# so that they run in CI's equivalence job rather than under R CMD check.

test_that("constant raw outcomes can have identified regression residuals (A03)", {
  latents:::.latents_local_seed(1027)
  data <- data.frame(y = rep(2, 60), x = seq(-2, 3, length.out = 60),
                     o = stats::rnorm(60), w = 1)
  formulas <- list(y ~ 0 + x, y ~ 1 + offset(o))
  cases <- expand.grid(formula = seq_along(formulas), weighted = c(FALSE, TRUE))
  tables <- lapply(seq_len(nrow(cases)), function(row) {
    formula <- formulas[[cases$formula[row]]]
    weighted <- cases$weighted[row]
    reference <- stats::lm(formula, data)
    fit <- mixture_regression(formula, data, 1, n_starts = 0, vcov_type = "none",
                              weights = if (weighted) "w" else NULL)
    expect_true(fit$converged)
    case <- data.frame(formula = deparse(formula), weighted = weighted)
    # Each formula has one coefficient.
    list(actual = data.frame(case, coefficient = as.numeric(fit$params$beta),
                             residual_variance = unname(fit$params$sigma2),
                             log_likelihood = as.numeric(logLik(fit))),
         expected = data.frame(case, coefficient = unname(stats::coef(reference)),
                               residual_variance = mean(stats::residuals(reference)^2),
                               log_likelihood = as.numeric(stats::logLik(reference))))
  })
  # stats::lm(): the maximum-likelihood residual variance is the mean squared residual.
  expect_reference(do.call(rbind, lapply(tables, `[[`, "actual")),
                   do.call(rbind, lapply(tables, `[[`, "expected")),
                   reference = "stats",
                   covers = c("mixture_regression()", "logLik.latents_mixture_regression()"),
                   quantity = "one-class Gaussian regression of a constant outcome",
                   tolerance = 1e-10, measure = "relative", key = c("formula", "weighted"))
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
