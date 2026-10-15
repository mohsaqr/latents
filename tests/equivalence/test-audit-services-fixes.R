test_that("weighted Gaussian residuals use the fitted distribution's weights", {
  set.seed(210)
  d <- data.frame(g = 1:120, a = rnorm(120), b = rnorm(120), w = rep(c(1, 10), each = 60))
  d$b[1:60] <- d$a[1:60] + rnorm(60, sd = 0.1)
  d$b[61:120] <- -d$a[61:120] + rnorm(60, sd = 0.1)
  fit <- multilpa(d, c("a", "b"), "g", 1, 1, weights = "w", n_starts = 1, seed = 1)
  residual <- sweep(as.matrix(d[c("a", "b")]), 2, fit$means[1, ], "-")
  scatter <- crossprod(residual * d$w, residual)
  expected <- scatter[1, 2] / sqrt(prod(diag(scatter)))
  actual <- get_results(fit, "residuals")
  expect_equal(actual$observed, expected, tolerance = 1e-12)
  expect_lt(actual$observed, -0.7)
  expect_equal(actual$effective_n, sum(fit$sampling_weights))
})

test_that("weighted categorical residuals use weighted cell counts", {
  d <- data.frame(g = 1:120, a = rep(c(0, 0, 1, 1), 30),
                  b = rep(c(0, 1, 0, 1), 30), w = rep(c(10, 1, 1, 10), 30))
  fit <- multilca(d, c("a", "b"), "g", 1, 1, weights = "w", n_starts = 1, seed = 1)
  counts <- xtabs(w ~ a + b, data = d)
  independent <- outer(rowSums(counts), colSums(counts)) / sum(counts)
  expected <- sqrt(sum((counts - independent)^2 / independent) / sum(counts))
  actual <- get_results(fit, "residuals")
  expect_equal(actual$observed, expected, tolerance = 1e-12)
  expect_equal(actual$expected, 0, tolerance = 1e-12)
})

test_that("Rubin pooling handles both zero-variance limits", {
  table <- data.frame(level = "group", outcome = "group_1", term = NA_character_,
                      parameter = "probability", estimate = 1, standard_error = 0,
                      statistic = NA_real_)
  constant <- latents:::.multilpa_rubin(list(table, table), 0.95)
  expect_equal(constant$riv, 0)
  expect_equal(constant$fmi, 0)
  expect_equal(constant$df, Inf)
  expect_equal(constant$standard_error, 0)
  second <- table
  second$estimate <- 0.5
  changing <- latents:::.multilpa_rubin(list(table, second), 0.95)
  expect_equal(changing$riv, Inf)
  expect_equal(changing$fmi, 1)
  expect_equal(changing$df, 1)
  expect_true(is.finite(changing$standard_error))
})

test_that("R3STEP estimates and errors are invariant to predictor units", {
  set.seed(606)
  z <- rnorm(400)
  state <- rbinom(400, 1, plogis(1.5 * z))
  d <- data.frame(g = 1:400, a = rnorm(400, 6 * state, 0.8),
                  b = rnorm(400, 5 * state, 0.8), z = z)
  fit <- multilpa(d, c("a", "b"), "g", 2, 1, n_starts = 1, seed = 1)
  invisible(lapply(c("observed", "robust"), function(type) {
    standard <- r3step(fit, d, "z", vcov_type = type)
    scaled <- d
    scaled$z <- scaled$z * 1e-9
    tiny <- r3step(fit, scaled, "z", vcov_type = type)
    factor <- ifelse(standard$term == "z", 1e-9, 1)
    expect_equal(tiny$estimate * factor, standard$estimate, tolerance = 1e-7)
    expect_equal(tiny$standard_error * factor, standard$standard_error, tolerance = 1e-7)
    expect_equal(tiny$p_value, standard$p_value, tolerance = 1e-7)
  }))
})
