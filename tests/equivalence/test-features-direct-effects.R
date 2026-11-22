# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Comparisons with independent references moved from tests/testthat/test-features-direct-effects.R,
# so that they run in CI's equivalence job rather than under R CMD check.

de_unit_data <- function(seed = 7, n = 400, b = c(0.8, -0.6)) {
  withr::with_seed(seed, {
    class <- sample.int(2L, n, replace = TRUE, prob = c(0.55, 0.45))
    z <- round(stats::rnorm(n, 0.5, 1.5), 4)
    w <- round(stats::rnorm(n, -1, 1), 4)
    y <- vapply(1:4, \(j) stats::rnorm(n, c(0, c(1.8, 1.5, 1.6, 1.2)[[j]])[class],
                                       c(1, c(0.8, 1.2, 0.9, 1.1)[[j]])[class]), numeric(n))
    data <- data.frame(id = seq_len(n), z = z, w = w, round(y, 4))
    names(data)[4:7] <- paste0("y", 1:4)
    data$y3 <- data$y3 + b[[1L]] * data$z
    data$y1 <- data$y1 + b[[2L]] * data$z
    data
  })
}

test_that("one profile is a linear regression of each indicator", {
  data <- de_unit_data()
  fit <- quietly(lpa(data, paste0("y", 1:4), 1, model = "VVI", n_starts = 1, seed = 1,
                     tol = 1e-12, direct_effects = c("y3 ~ z + w", "y1 ~ z")))
  models <- list(y1 = stats::lm(y1 ~ z, data), y2 = stats::lm(y2 ~ 1, data),
                 y3 = stats::lm(y3 ~ z + w, data), y4 = stats::lm(y4 ~ 1, data))
  estimates <- coef(fit)
  indicators <- names(models)
  # A one-profile fit is these four regressions; both solve the same normal
  # equations, so only rounding separates them (the original 1e-9).
  expect_reference(estimates[sprintf("measurement.mean.profile_1.%s", indicators)],
                   vapply(models, \(m) stats::coef(m)[["(Intercept)"]], numeric(1)),
                   reference = "stats", covers = c("lpa()", "coef.multilpa()"),
                   quantity = "intercepts, one profile with direct effects",
                   tolerance = 1e-9, measure = "relative", key = indicators)
  expect_reference(estimates[sprintf("measurement.direct_effect.shared.%s", c("y3:z", "y3:w", "y1:z"))],
                   c(stats::coef(models$y3)[c("z", "w")], stats::coef(models$y1)[["z"]]),
                   reference = "stats", covers = c("lpa()", "coef.multilpa()"),
                   quantity = "direct effects, one profile", tolerance = 1e-9,
                   measure = "relative", key = c("y3:z", "y3:w", "y1:z"))
  # ML residual variance: the residual sum of squares over n.
  expect_reference(estimates[sprintf("measurement.variance.profile_1.%s", indicators)],
                   vapply(models, \(m) mean(stats::residuals(m)^2), numeric(1)),
                   reference = "stats", covers = c("lpa()", "coef.multilpa()"),
                   quantity = "residual variances, one profile with direct effects",
                   tolerance = 1e-9, measure = "relative", key = indicators)
  expect_reference(as.numeric(logLik(fit)),
                   sum(vapply(models, \(m) as.numeric(stats::logLik(m)), numeric(1))),
                   reference = "stats", covers = c("lpa()", "logLik.multilpa()"),
                   quantity = "log likelihood, one profile with direct effects",
                   tolerance = 1e-12, measure = "relative")
})
