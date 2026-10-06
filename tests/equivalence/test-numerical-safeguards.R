# Exhaustive property checks: they run on CI and locally, and are skipped on
# CRAN to keep the check within its time budget.
skip_on_cran()

test_that("an unbalanced-group update agrees with exhaustive assignments", {
  data <- data.frame(group = rep(c("large", "singleton", "middle", "pair"),
                                 times = c(4, 1, 3, 2)),
                     x = c(-1.2, 0.8, -0.6, 1.2, 0.3, -0.9, 0.9, 1.1, -0.5, 0.7))
  parameters <- list(means = matrix(c(-0.7, 0.8), 2L, 1L),
                      variances = matrix(c(0.7, 0.8), 2L, 1L),
                      profile_probabilities = rbind(c(0.6, 0.4), c(0.4, 0.6)),
                      group_probabilities = c(0.5, 0.5))
  reference <- enumerated_em_update(as.matrix(data["x"]), data$group, parameters)
  expect_warning(
    fit <- multilpa(data, "x", "group", 2L, 2L, n_starts = 1L,
                      start = parameters, max_iter = 1L), "did not converge")
  expect_equal(unname(fit$group_probabilities), reference$group_probabilities,
               tolerance = 1e-10)
  expect_equal(unname(fit$profile_probabilities), reference$profile_probabilities,
               tolerance = 1e-10)
  expect_equal(unname(fit$means), unname(reference$means), tolerance = 1e-10)
  expect_equal(unname(fit$variances), unname(reference$variances), tolerance = 1e-10)
})
