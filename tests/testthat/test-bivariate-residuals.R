.dependent_data <- function(seed = 4L, n_groups = 20L, per = 10L) {
  set.seed(seed)
  school <- rep(seq_len(n_groups), each = per)
  profile <- rep(rep(c(1L, 2L), length.out = n_groups), each = per)
  n <- length(school)
  # `a` and `b` share a term the profiles cannot explain; `c` does not.
  shared <- stats::rnorm(n)
  data.frame(school = school,
             a = stats::rnorm(n, ifelse(profile == 2L, 2, -2), 0.8) + shared,
             b = stats::rnorm(n, ifelse(profile == 2L, 2, -2), 0.8) + shared,
             c = stats::rnorm(n, ifelse(profile == 2L, 1, -1), 0.8))
}

test_that("a planted dependence is found and the innocent pairs are not", {
  data <- .dependent_data()
  fit <- multilpa(data, c("a", "b", "c"), "school", n_profiles = 2,
                  n_group_classes = 1, n_starts = 4, seed = 1)
  residuals <- bivariate_residuals(fit, data)
  worst <- residuals[1L, ]

  # ordered worst first, and the worst pair is the one that shares a term
  expect_equal(worst$indicator_1, "a")
  expect_equal(worst$indicator_2, "b")
  expect_gt(worst$residual, 0.5)
  expect_lt(worst$p_value, 1e-6)
  innocent <- subset(residuals, !(indicator_1 == "a" & indicator_2 == "b"))
  expect_lt(max(abs(innocent$residual)), 0.25)
  expect_true(all(abs(residuals$residual) == sort(abs(residuals$residual),
                                                  decreasing = TRUE)))
})

test_that("estimating the association makes the residual vanish", {
  data <- .dependent_data()
  diagonal <- multilpa(data, c("a", "b", "c"), "school", n_profiles = 2,
                       n_group_classes = 1, n_starts = 4, seed = 1)
  full <- multilpa(data, c("a", "b", "c"), "school", n_profiles = 2,
                   n_group_classes = 1, n_starts = 4, seed = 1,
                   covariance_model = "full")

  flagged <- bivariate_residuals(diagonal, data)
  resolved <- bivariate_residuals(full, data)
  expect_gt(max(abs(flagged$residual)), 0.5)
  # The full model estimates what the diagonal one assumed away.
  expect_lt(max(abs(resolved$residual)), 1e-3)
  expect_true(all(resolved$p_value > 0.99))
  # and its expectation is the fitted correlation, not zero
  expect_true(all(flagged$expected == 0))
  expect_gt(max(abs(resolved$expected)), 0.5)
})

test_that("the table is tidy, complete and correctly weighted", {
  data <- .dependent_data()
  fit <- multilpa(data, c("a", "b", "c"), "school", n_profiles = 2,
                  n_group_classes = 1, n_starts = 4, seed = 1)
  residuals <- bivariate_residuals(fit, data)

  expect_named(residuals, c("profile", "indicator_1", "indicator_2", "kind",
                            "observed", "expected", "residual", "effective_n",
                            "statistic", "df", "p_value"))
  # three pairs for three indicators, assessed in each of two profiles
  expect_equal(nrow(residuals), 3L * 2L)
  expect_setequal(unique(residuals$profile), c("profile_1", "profile_2"))
  expect_true(all(residuals$kind == "gaussian"))
  # each profile's rows carry that profile's own effective size
  sizes <- unique(residuals[c("profile", "effective_n")])
  expect_equal(sizes$effective_n[order(sizes$profile)],
               unname(fit$effective_profile_counts))
  expect_equal(residuals$residual, residuals$observed - residuals$expected)
})

test_that("pooling gives one row per pair", {
  data <- .dependent_data()
  fit <- multilpa(data, c("a", "b", "c"), "school", n_profiles = 2,
                  n_group_classes = 1, n_starts = 4, seed = 1)
  pooled <- bivariate_residuals(fit, data, by = "overall")

  expect_equal(nrow(pooled), 3L)
  expect_true(all(pooled$profile == "overall"))
  expect_equal(unique(pooled$effective_n), nrow(data))
})

test_that("categorical pairs are assessed with a chi-square", {
  set.seed(6)
  n_groups <- 40L
  school <- rep(seq_len(n_groups), each = 12L)
  profile <- rep(rep(c(1L, 2L), length.out = n_groups), each = 12L)
  n <- length(school)
  # Five indicators, so a two-class model is not saturated: the table has 31
  # degrees of freedom against 11 free parameters. u1 and u2 additionally share
  # a coin the profiles cannot explain; u3 to u5 do not.
  chance <- function() stats::runif(n) < ifelse(profile == 2L, 0.8, 0.2)
  shared <- stats::runif(n) < 0.5
  data <- data.frame(
    school = school,
    u1 = as.integer(ifelse(shared, TRUE, chance())),
    u2 = as.integer(ifelse(shared, TRUE, chance())),
    u3 = as.integer(chance()), u4 = as.integer(chance()),
    u5 = as.integer(chance()))
  indicators <- c("u1", "u2", "u3", "u4", "u5")
  fit <- multilpa(data, indicators, "school", n_profiles = 2,
                  n_group_classes = 1, categorical = indicators,
                  n_starts = 6, seed = 3)
  pooled <- bivariate_residuals(fit, data, by = "overall")

  expect_true(all(pooled$kind == "categorical"))
  expect_true(all(pooled$df == 1L))
  expect_equal(nrow(pooled), 10L)
  worst <- pooled[1L, ]
  expect_setequal(c(worst$indicator_1, worst$indicator_2), c("u1", "u2"))
  expect_lt(worst$p_value, 0.01)
  # the pairs with no planted dependence are not flagged
  innocent <- subset(pooled, !(indicator_1 == "u1" & indicator_2 == "u2"))
  expect_gt(min(innocent$p_value), 0.01)
})

test_that("a covariate fit is assessed too", {
  data <- .dependent_data()
  data$x <- stats::rnorm(nrow(data))
  fit <- fit_covariates(data, c("a", "b", "c"), "school", n_profiles = 2,
                        n_group_classes = 2, profile_covariates = "x",
                        n_starts = 4, seed = 1)
  residuals <- bivariate_residuals(fit, data)

  expect_equal(nrow(residuals), 3L * 2L)
  expect_gt(max(abs(residuals$residual)), 0.4)
})

test_that("a broken contract is refused", {
  data <- .dependent_data()
  fit <- multilpa(data, c("a", "b", "c"), "school", n_profiles = 2,
                  n_group_classes = 1, n_starts = 3, seed = 1)

  expect_error(bivariate_residuals(fit, data[1:10, ]),
               "one row per observation")
  expect_error(bivariate_residuals(fit, subset(data, select = c(school, a))),
               "must contain the fitted indicators")
  expect_error(bivariate_residuals("not a fit", data),
               "must be a fitted model")
})

test_that("a single indicator has no pair to assess", {
  data <- .dependent_data()
  fit <- multilpa(data, "a", "school", n_profiles = 2, n_group_classes = 1,
                  n_starts = 3, seed = 1)
  residuals <- bivariate_residuals(fit, data)

  expect_s3_class(residuals, "data.frame")
  expect_equal(nrow(residuals), 0L)
  expect_named(residuals, c("profile", "indicator_1", "indicator_2", "kind",
                            "observed", "expected", "residual", "effective_n",
                            "statistic", "df", "p_value"))
})
