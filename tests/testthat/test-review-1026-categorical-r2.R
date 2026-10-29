# Regression tests for the round-2 review of 2026-10-26, categorical and
# extra indicators (validation/review-2026-10-26/round-2/categorical.md).
# The comparisons with independent references are in tests/equivalence/test-review-1026-categorical-r2.R.

# More than 1.5 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

# --- C2-1: a shared negative-binomial dispersion is not forced to its floor --

# --- C2-2: half-observed local-dependence pairs ------------------------------

r2_partial_data <- function(seed = 31, n = 600) {
  if (exists(".Random.seed", globalenv())) {
    saved <- get(".Random.seed", globalenv())
    on.exit(assign(".Random.seed", saved, globalenv()), add = TRUE)
  }
  set.seed(seed)
  class <- sample(1:2, n, TRUE, prob = c(0.55, 0.45))
  u1 <- vapply(class, function(k) sample(1:3, 1L, prob = if (k == 1) c(.6, .3, .1) else
    c(.1, .3, .6)), integer(1))
  u2 <- ifelse(stats::runif(n) < 0.5, as.integer(u1 >= 2), stats::rbinom(n, 1, c(.3, .7)[class]))
  data <- data.frame(u1 = u1, u2 = u2, u3 = stats::rbinom(n, 1, c(.2, .8)[class]),
                     u4 = stats::rbinom(n, 1, c(.3, .75)[class]))
  data$u2[stats::runif(n) < ifelse(data$u1 == 3, 0.5, 0.1)] <- NA
  data$u1[sample(which(!is.na(data$u2)), 40)] <- NA
  data$u4[sample(n, 30)] <- NA
  data
}

r2_partial_fit <- function(data = r2_partial_data()) {
  lca(data, vars = names(data), n_classes = 2, local_dependence = "u1:u2",
      missing = "fiml", n_starts = 10, seed = 1, tol = 1e-13, max_iter = 20000)
}

test_that("half-observed pairs under MAR give consistent estimates (C2-2)", {
  skip_on_cran()
  if (exists(".Random.seed", globalenv())) {
    saved <- get(".Random.seed", globalenv())
    on.exit(assign(".Random.seed", saved, globalenv()), add = TRUE)
  }
  set.seed(44)
  n <- 40000
  class <- sample(1:2, n, TRUE, prob = c(0.6, 0.4))
  u1 <- stats::rbinom(n, 1, c(0.2, 0.8)[class])
  copy <- stats::runif(n) < 0.6
  u2 <- ifelse(copy, u1, stats::rbinom(n, 1, c(0.3, 0.7)[class]))
  data <- data.frame(u1 = u1, u2 = u2, u3 = stats::rbinom(n, 1, c(0.25, 0.8)[class]),
                     u4 = stats::rbinom(n, 1, c(0.3, 0.75)[class]),
                     u5 = stats::rbinom(n, 1, c(0.2, 0.6)[class]))
  # MAR: whether u2 is seen depends on the observed u1.
  data$u2[stats::runif(n) < ifelse(data$u1 == 1, 0.7, 0.05)] <- NA
  fit <- lca(data, vars = names(data), n_classes = 2, local_dependence = "u1:u2",
             missing = "fiml", n_starts = 5, seed = 1)
  responses <- get_results(fit, "responses")
  third <- subset(responses, indicator == "u3" & category == "1")
  low <- third$profile[which.min(third$probability)]
  joint <- subset(responses, indicator == "u1:u2")
  truth_u1 <- c(0.2, 0.8)
  truth_u2 <- c(0.3, 0.7)
  joint$class <- ifelse(joint$profile == low, 1L, 2L)
  first <- as.integer(sub(":.*$", "", joint$category))
  second <- as.integer(sub("^.*:", "", joint$category))
  p1 <- ifelse(first == 1, truth_u1[joint$class], 1 - truth_u1[joint$class])
  p2 <- ifelse(second == 1, truth_u2[joint$class], 1 - truth_u2[joint$class])
  truth <- p1 * (0.6 * (first == second) + 0.4 * p2)
  # The old code (the pair dropped from half-observed rows) was more than 50
  # standard errors off on the u1 = 1 cells.
  expect_true(all(abs(joint$probability - truth) < 3 * joint$probability_standard_error))
})

test_that("predict(), enumeration and refits use half-observed pair rows (C2-2)", {
  skip_on_cran()
  data <- r2_partial_data()
  fit <- r2_partial_fit(data)
  posterior <- predict(fit, newdata = data, type = "posterior")
  expect_equal(unname(as.matrix(posterior[c("probability_profile_1", "probability_profile_2")])),
               unname(fit$subject_posteriors), tolerance = 1e-12)
  grid <- quietly(enumerate_lca(data, names(data), 1:2, local_dependence = "u1:u2",
                                missing = "fiml", n_starts = 10, seed = 1, tol = 1e-13,
                                max_iter = 20000))
  cells <- as.data.frame(grid)
  expect_equal(cells$log_likelihood[cells$n_profiles == 2], as.numeric(logLik(fit)),
               tolerance = 1e-8)
  # A refit from the fit's own rebuilt data (what bootstrap and sensitivity
  # refit) keeps the half-observed rows: same likelihood at the fit's seed.
  frame <- .multilpa_resolve_data(fit, NULL)
  refit <- quietly(do.call(multilpa, c(list(data = frame, n_starts = 10, seed = 1,
                                            tol = 1e-13, max_iter = 20000),
                                       .multilpa_refit_arguments(fit))))
  expect_equal(as.numeric(logLik(refit)), as.numeric(logLik(fit)), tolerance = 1e-9)
  boot <- quietly(parameter_inference(fit, method = "bootstrap", iter = 3L, n_starts = 2L,
                                      seed = 1, boundary = "fix"))
  expect_s3_class(boot, "data.frame", exact = FALSE)
  # The original items come back with the observed one kept.
  labels <- .latents_categorical_labels(fit)
  expect_identical(unname(labels[, "u1"]), as.character(data$u1))
  expect_identical(unname(labels[, "u2"]), as.character(data$u2))
})

# --- C2-3: every pair of indicators has a bivariate residual -----------------

test_that("bivariate residuals cover every pair, with mixed and count kinds (C2-3)", {
  if (exists(".Random.seed", globalenv())) {
    saved <- get(".Random.seed", globalenv())
    on.exit(assign(".Random.seed", saved, globalenv()), add = TRUE)
  }
  set.seed(3)
  n <- 300
  class <- rep(1:2, c(160, 140))
  data <- data.frame(g1 = stats::rnorm(n, c(0, 2)[class]), g2 = stats::rnorm(n, c(0, 1)[class]),
                     c1 = stats::rbinom(n, 1, c(.2, .8)[class]), c2 = stats::rbinom(n, 1, c(.3, .7)[class]),
                     o1 = sample(1:4, n, TRUE), k1 = stats::rpois(n, c(1, 4)[class]))
  fit <- lpa(data, vars = names(data), categorical = c("c1", "c2"), ordinal = "o1",
             count = "k1", n_profiles = 2, n_starts = 3, seed = 1)
  for (by in c("overall", "profile")) {
    residuals <- get_results(fit, "residuals", by = by)
    expect_identical(nrow(residuals), 15L * if (by == "profile") 2L else 1L)
  }
  residuals <- get_results(fit, "residuals", by = "profile")
  kinds <- with(residuals, tapply(kind, paste(indicator_1, indicator_2), unique))
  expect_identical(as.vector(kinds[c("g1 g2", "c1 c2", "c1 o1", "g1 c1", "g2 o1", "o1 k1", "g1 k1")]),
                   c("gaussian", "categorical", "categorical", "mixed", "mixed", "count", "count"))
  untested <- subset(residuals, kind %in% c("mixed", "count"))
  expect_true(all(is.na(untested$statistic) & is.na(untested$df) & is.na(untested$p_value)))
  expect_true(all(untested$expected == 0))
  # The mixed association is the posterior-weighted correlation in the profile.
  row <- subset(residuals, profile == "profile_1" & indicator_1 == "g1" & indicator_2 == "k1")
  weight <- fit$subject_posteriors[, 1]
  values <- cbind(.multilpa_center_like(fit, as.matrix(data["g1"]))[, 1], data$k1)
  centred <- sweep(values, 2L, colSums(values * weight) / sum(weight))
  scatter <- crossprod(centred * weight, centred)
  expect_equal(row$observed, scatter[1, 2] / sqrt(scatter[1, 1] * scatter[2, 2]),
               tolerance = 1e-12)
  # The ordinal pair is a Pearson table on the ordinal's categories.
  table_row <- subset(residuals, profile == "profile_2" & indicator_1 == "c1" & indicator_2 == "o1")
  probability <- exp(.latents_ordinal_log_probabilities(fit$ordinal_intercepts[[1]],
                                                        fit$ordinal_locations[, 1]))
  weight <- fit$subject_posteriors[, 2]
  observed <- as.matrix(stats::xtabs(weight ~ factor(data$c1, 0:1) + factor(data$o1, 1:4)))
  expected <- sum(weight) * outer(fit$response_probabilities$c1[2, ], probability[2, ])
  expect_equal(table_row$statistic, sum((observed - expected)^2 / expected), tolerance = 1e-10)
  expect_identical(table_row$df, 3L)
})
