Sys.setenv(NOT_CRAN = "true")
pkgload::load_all(quiet = TRUE)
library(testthat)

test_that("bootstrap data mismatch is reproducible and Wald rejects the same input", {
  set.seed(719)
  d <- data.frame(g = rep(seq_len(30), each = 4), a = rnorm(120), b = rnorm(120))
  x <- multilpa(d, c("a", "b"), "g", 1, 1, n_starts = 1, seed = 1)
  wrong <- d
  wrong$a <- 10 * wrong$a
  wrong$b <- 10 * wrong$b
  expect_error(parameter_inference(x, data = wrong), class = "latents_bad_inference_data")
  ordinary <- parameter_inference(x, method = "bootstrap", iter = 12, n_starts = 1, seed = 1)
  mismatched <- parameter_inference(x, data = wrong, method = "bootstrap", iter = 12, n_starts = 1, seed = 1)
  expect_equal(mismatched$estimate, ordinary$estimate)
  expect_equal(mismatched$standard_error[1:2], 10 * ordinary$standard_error[1:2], tolerance = 1e-10)
  expect_equal(mismatched$standard_error[3:4], 100 * ordinary$standard_error[3:4], tolerance = 1e-10)
  expect_gt(mismatched$conf_low[3], 50)
  expect_lt(mismatched$estimate[3], 1)
})

test_that("weighted residual diagnostic disagrees in sign with weighted crossproduct", {
  set.seed(210)
  d <- data.frame(g = 1:120, a = rnorm(120), b = rnorm(120), w = rep(c(1, 10), each = 60))
  d$b[1:60] <- d$a[1:60] + rnorm(60, sd = 0.1)
  d$b[61:120] <- -d$a[61:120] + rnorm(60, sd = 0.1)
  x <- multilpa(d, c("a", "b"), "g", 1, 1, weights = "w", n_starts = 1, seed = 1)
  r <- sweep(as.matrix(d[c("a", "b")]), 2, x$means[1, ], "-")
  expected <- crossprod(r * d$w, r)
  weighted <- expected[1, 2] / sqrt(prod(diag(expected)))
  observed <- get_results(x, "residuals")$observed
  expect_gt(observed, 0)
  expect_lt(weighted, -0.7)
})

test_that("one-class identical imputations expose zero-variance missing-information bug", {
  set.seed(45)
  d <- data.frame(g = rep(1:30, each = 4), a = rnorm(120), b = rnorm(120))
  pooled <- pool_imputations(list(d, d), c("a", "b"), "g", 1,
                             n_group_classes = 1, n_starts = 1, seed = 1)
  estimates <- get_results(pooled)
  constant <- estimates$parameter == "probability"
  expect_true(all(estimates$within[constant] == 0))
  expect_true(all(estimates$between[constant] == 0))
  expect_true(all(is.nan(estimates$fmi[constant])))
  expect_true(is.nan(summary(pooled)$fmi_max))
})

test_that("matrix descriptives double-count observations and invent a fourth group", {
  d <- data.frame(g = rep(1:3, each = 2))
  d$a <- matrix(1:12, ncol = 2)
  actual <- descriptives(d, vars = "a", id = "g")
  expect_equal(actual$n, 12L)
  expect_equal(nrow(d), 6L)
  expect_equal(actual$n_groups, 4L)
  expect_equal(length(unique(d$g)), 3L)
})

test_that("R3STEP small predictor units alone cause the recorded numerical refusal", {
  set.seed(606)
  z <- rnorm(400)
  state <- rbinom(400, 1, plogis(1.5 * z))
  d <- data.frame(g = 1:400, a = rnorm(400, 6 * state, 0.8),
                  b = rnorm(400, 5 * state, 0.8), z = z)
  x <- multilpa(d, c("a", "b"), "g", 2, 1, n_starts = 1, seed = 1)
  ordinary <- r3step(x, d, "z", vcov_type = "observed")
  expect_true(all(is.finite(ordinary$standard_error)))
  d$z <- d$z * 1e-9
  expect_error(r3step(x, d, "z", vcov_type = "observed"),
                class = "latents_singular_information")
})

# These are characterization tests for audit evidence. Their passing status
# confirms the defects remain reproducible; it does not endorse the behavior.
