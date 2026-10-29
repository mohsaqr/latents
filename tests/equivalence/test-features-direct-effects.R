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
  expect_equal(unname(fit$means[1, ]),
               unname(vapply(models, \(m) stats::coef(m)[["(Intercept)"]], numeric(1))),
               tolerance = 1e-9)
  expect_equal(unname(fit$direct_effects[c("y3:z", "y3:w", "y1:z")]),
               unname(c(stats::coef(models$y3)[c("z", "w")], stats::coef(models$y1)[["z"]])),
               tolerance = 1e-9)
  # ML residual variance: the residual sum of squares over n.
  expect_equal(unname(fit$variances[1, ]),
               unname(vapply(models, \(m) mean(stats::residuals(m)^2), numeric(1))),
               tolerance = 1e-9)
  expect_equal(fit$log_likelihood,
               sum(vapply(models, \(m) as.numeric(stats::logLik(m)), numeric(1))),
               tolerance = 1e-12)
})
