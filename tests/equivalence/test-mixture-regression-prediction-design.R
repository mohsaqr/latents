
test_that("a single group slope is optimized as a slope rather than an intercept", {
  set.seed(86)
  data <- data.frame(g = rep(1:30, each = 4), x = rnorm(120), y = rnorm(120),
                     w = rep(seq(-2, 2, length.out = 30), each = 4))
  spec <- latents:::.mixture_spec(y ~ x, data, 2, "gaussian", "g", "observation", 2,
                                  NULL, ~1, ~0 + w, "varying", 1e-6, "error")
  expectation <- latents:::.mixture_seed_expectation(spec, rep(1:2, 60), rep(1:2, 15))
  probability <- plogis(1.3 * spec$v[, 1])
  expectation$rho <- cbind(1 - probability, probability)
  updated <- latents:::.mixture_update_mixing(spec, latents:::.mixture_null_params(spec),
                                              expectation)
  # The group-class M step is a weighted logistic regression of the group
  # posteriors on the group covariate: glm() with a quasibinomial family.
  reference <- stats::glm(probability ~ 0 + spec$v[, 1],
                          family = stats::quasibinomial())
  expect_reference(unname(updated$delta[1, 2]), unname(stats::coef(reference)),
                   reference = "stats", covers = "mixture_regression()",
                   quantity = "group-class slope of the mixing M step", tolerance = 1e-8,
                   measure = "relative",
                   note = "Package side: the internal M step of the group-class logits (.mixture_update_mixing) inside mixture_regression(); reference: stats::glm(family = quasibinomial()).")
})
