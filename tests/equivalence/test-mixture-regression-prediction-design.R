
test_that("a single group slope is optimized as a slope rather than an intercept", {
  set.seed(86)
  data <- data.frame(g = rep(1:30, each = 4), x = rnorm(120), y = rnorm(120),
                     w = rep(seq(-2, 2, length.out = 30), each = 4))
  spec <- .mixture_spec(y ~ x, data, 2, "gaussian", "g", "observation", 2,
                        NULL, ~1, ~0 + w, "varying", 1e-6, "error")
  expectation <- .mixture_seed_expectation(spec, rep(1:2, 60), rep(1:2, 15))
  probability <- plogis(1.3 * spec$v[, 1])
  expectation$rho <- cbind(1 - probability, probability)
  updated <- .mixture_update_mixing(spec, .mixture_null_params(spec), expectation)
  reference <- stats::glm(probability ~ 0 + spec$v[, 1],
                          family = stats::quasibinomial())
  expect_equal(unname(updated$delta[1, 2]), unname(coef(reference)),
               tolerance = 1e-8)
})
