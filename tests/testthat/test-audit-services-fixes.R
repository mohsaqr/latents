test_that("bootstrap inference validates original observations before resampling", {
  set.seed(719)
  d <- data.frame(g = rep(1:30, each = 4), a = rnorm(120), b = rnorm(120))
  fit <- multilpa(d, c("a", "b"), "g", 1, 1, n_starts = 1, seed = 1)
  wrong <- d
  wrong$a <- 10 * wrong$a
  wrong$b <- 10 * wrong$b
  state <- .Random.seed
  expect_error(parameter_inference(fit, data = wrong, method = "bootstrap",
                                    iter = 2, n_starts = 1, seed = 1),
               class = "latents_bad_inference_data")
  expect_identical(.Random.seed, state)
  expect_error(parameter_inference(fit, data = d[-1, ], method = "bootstrap",
                                    iter = 2, n_starts = 1),
               class = "latents_bad_inference_data")
  expect_error(parameter_inference(fit, data = d[120:1, ], method = "bootstrap",
                                    iter = 2, n_starts = 1),
               class = "latents_bad_inference_data")
  original <- parameter_inference(fit, method = "bootstrap", iter = 3,
                                  n_starts = 1, seed = 1)
  supplied <- parameter_inference(fit, data = d, method = "bootstrap", iter = 3,
                                  n_starts = 1, seed = 1)
  expect_identical(original, supplied)
})

test_that("bootstrap inference checks ordering and normalized sampling weights", {
  set.seed(612)
  d <- data.frame(g = rep(1:30, each = 4), t = rep(1:4, 30),
                  a = rnorm(120), b = rnorm(120), w = rep(1:30, each = 4))
  fit <- multilpa(d, c("a", "b"), "g", 1, 1, time = "t", weights = "w",
                  n_starts = 1, seed = 1)
  wrong <- d
  wrong$t[1:4] <- 4:1
  expect_error(parameter_inference(fit, data = wrong, method = "bootstrap",
                                    iter = 2, n_starts = 1),
               class = "latents_bad_inference_data")
  wrong <- d
  wrong$w[1:4] <- wrong$w[1:4] * 2
  expect_error(parameter_inference(fit, data = wrong, method = "bootstrap",
                                    iter = 2, n_starts = 1),
               class = "latents_bad_inference_data")
  wrong$w[1] <- NA_real_
  expect_error(parameter_inference(fit, data = wrong, method = "bootstrap",
                                    iter = 2, n_starts = 1),
               class = "latents_bad_inference_data")
  rescaled <- d
  rescaled$w <- rescaled$w * 10
  expect_no_error(latents:::.multilpa_validate_bootstrap_data(fit, rescaled))
})

test_that("single-level and centered original data remain valid for bootstrap", {
  set.seed(910)
  d <- data.frame(g = rep(1:30, each = 4), a = rnorm(120) + 100,
                  b = rnorm(120) - 50)
  single <- lpa(d, c("a", "b"), 1, n_starts = 1, seed = 1)
  expect_no_error(parameter_inference(single, data = d, method = "bootstrap",
                                      iter = 3, n_starts = 1, seed = 1))
  invisible(lapply(c("person", "grand"), function(centering) {
    fit <- multilpa(d, c("a", "b"), "g", 1, 1, centering = centering,
                    n_starts = 1, seed = 1)
    stored <- parameter_inference(fit, method = "bootstrap", iter = 3,
                                  n_starts = 1, seed = 1)
    explicit <- parameter_inference(fit, data = d, method = "bootstrap", iter = 3,
                                    n_starts = 1, seed = 1)
    expect_equal(stored, explicit, tolerance = 1e-10)
  }))
})

test_that("identical one-class imputations report zero missing information", {
  set.seed(45)
  d <- data.frame(g = rep(1:30, each = 4), a = rnorm(120), b = rnorm(120))
  pooled <- pool_imputations(list(d, d), c("a", "b"), "g", 1,
                             n_group_classes = 1, n_starts = 1, seed = 1)
  estimates <- get_results(pooled)
  expect_true(all(estimates$between == 0))
  expect_true(all(estimates$riv == 0))
  expect_true(all(estimates$fmi == 0))
  expect_equal(summary(pooled)$fmi_max, 0)
  expect_equal(summary(pooled)$fmi_median, 0)
})

test_that("descriptive vectors retain rows and reject matrix identifiers or indicators", {
  d <- data.frame(g = rep(1:3, each = 2), a = 1:6)
  actual <- descriptives(d, vars = "a", id = "g")
  expect_equal(actual$n, 6L)
  expect_equal(actual$n_groups, 3L)
  d$a <- matrix(1:12, ncol = 2)
  expect_error(descriptives(d, vars = "a", id = "g"), class = "latents_bad_data")
  expect_error(descriptives(d, id = "g"), class = "latents_bad_data")
  d$a <- 1:6
  d$g <- cbind(1:6, 6:1)
  expect_error(descriptives(d, vars = "a", id = "g"), class = "latents_bad_data")
  expect_error(descriptives(d, vars = "a", by = "g"), class = "latents_bad_data")
})

test_that("weighted residual tables withhold unsupported inferential columns", {
  set.seed(210)
  d <- data.frame(g = 1:80, a = rnorm(80), b = rnorm(80), w = rep(c(1, 10), each = 40))
  fit <- multilpa(d, c("a", "b"), "g", 1, 1, weights = "w", n_starts = 1, seed = 1)
  rows <- get_results(fit, "residuals")
  expect_true(all(is.finite(rows$observed)))
  expect_true(all(is.na(rows$statistic)))
  expect_true(all(is.na(rows$p_value)))
  expect_true(all(is.na(rows$p_adjusted)))
  categorical <- data.frame(g = 1:120, a = rep(c(0, 0, 1, 1), 30),
                            b = rep(c(0, 1, 0, 1), 30), w = rep(1:30, each = 4))
  fit <- multilca(categorical, c("a", "b"), "g", 1, 1, weights = "w",
                  n_starts = 1, seed = 1)
  rows <- get_results(fit, "residuals")
  expect_true(all(is.finite(rows$observed)))
  expect_true(all(is.na(rows$statistic)))
  expect_true(all(is.na(rows$df)))
  expect_true(all(is.na(rows$p_value)))
})
