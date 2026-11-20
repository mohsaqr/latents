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
  # One EM update whose posteriors enumerate every latent assignment of each group.
  enumerated <- "independent: one EM update from exhaustive latent assignments (helper-independent-likelihood.R)"
  expect_reference(unname(fit$group_probabilities), reference$group_probabilities,
                   reference = enumerated, covers = "multilpa()",
                   quantity = "group-class probabilities after one update, unbalanced groups",
                   tolerance = 1e-10, measure = "relative")
  expect_reference(unname(fit$profile_probabilities), unname(reference$profile_probabilities),
                   reference = enumerated, covers = "multilpa()",
                   quantity = "profile probabilities after one update, unbalanced groups",
                   tolerance = 1e-10, measure = "relative")
  expect_reference(unname(fit$means), unname(reference$means),
                   reference = enumerated, covers = "multilpa()",
                   quantity = "means after one update, unbalanced groups",
                   tolerance = 1e-10, measure = "relative")
  expect_reference(unname(fit$variances), unname(reference$variances),
                   reference = enumerated, covers = "multilpa()",
                   quantity = "variances after one update, unbalanced groups",
                   tolerance = 1e-10, measure = "relative")
})
