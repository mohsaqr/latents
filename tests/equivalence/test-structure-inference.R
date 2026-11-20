# Wald inference for the ten covariance structures that constrain the volume,
# the shape or the orientation across profiles. The reference is
# helper-structure-reference.R: its own likelihood, its own chart of each
# structure and numerical derivatives, sharing nothing with R/inference.R. The
# full fourteen-structure comparison with numDeriv lives in
# equivalence/test-structure-inference.R.

skip_on_cran()

structure_data <- function(n_groups = 40L, size = 5L, seed = 11L) {
  set.seed(seed)
  unit <- rep(seq_len(n_groups), each = size)
  kind <- rep(stats::rbinom(n_groups, 1L, 0.5), each = size)
  high <- stats::rbinom(length(unit), 1L, ifelse(kind == 1L, 0.8, 0.2))
  a <- stats::rnorm(length(unit), 2.5 * high, 1 + 0.5 * high)
  b <- stats::rnorm(length(unit), -1.5 * high, 1.2 - 0.4 * high) +
    0.5 * a * (1 - high)
  data.frame(unit = unit, a = a, b = b)
}

structure_fit_at <- function(code, data, n_group_classes = 1L) {
  quietly(do.call(multilpa, c(
    list(data = data, vars = c("a", "b"), id = "unit", n_profiles = 2L,
         n_group_classes = n_group_classes, n_starts = 3L, seed = 2L,
         tol = 1e-14, max_iter = 20000L),
    latents:::.multilpa_structure_arguments(code))))
}

test_that("constrained-structure standard errors match an independent Hessian", {
  data <- structure_data()
  cases <- expand.grid(code = c("VII", "VEI", "EVI", "VEE", "EVE", "EEV", "EVV"),
                       classes = c(1L, 2L), stringsAsFactors = FALSE)
  tables <- lapply(seq_len(nrow(cases)), function(row) {
    code <- cases$code[row]
    fit <- structure_fit_at(code, data, cases$classes[row])
    expect_true(fit$converged, label = code)
    ours <- parameter_inference(fit)
    reference <- reference_structure_se(fit, data)
    keys <- data.frame(code = code, group_classes = cases$classes[row],
                       ours[c("level", "outcome", "term", "parameter")])
    list(actual = data.frame(keys, estimate = ours$estimate,
                             standard_error = ours$standard_error),
         expected = data.frame(keys, estimate = reference$estimate,
                               standard_error = reference$standard_error),
         likelihood = data.frame(code = code, group_classes = cases$classes[row],
                                 actual = fit$log_likelihood,
                                 expected = reference$log_likelihood))
  })
  actual <- do.call(rbind, lapply(tables, `[[`, "actual"))
  expected <- do.call(rbind, lapply(tables, `[[`, "expected"))
  likelihood <- do.call(rbind, lapply(tables, `[[`, "likelihood"))
  case_key <- c("code", "group_classes")
  row_key <- c(case_key, "level", "outcome", "term", "parameter")
  chart <- "independent: structure chart, likelihood and finite-difference Hessian (helper-structure-reference.R)"
  # The reference's own mixture likelihood at the fitted parameters.
  expect_reference(data.frame(likelihood[case_key], log_likelihood = likelihood$actual),
                   data.frame(likelihood[case_key], log_likelihood = likelihood$expected),
                   reference = chart, covers = "multilpa()",
                   quantity = "log likelihood, constrained structures", tolerance = 1e-10,
                   measure = "relative", key = case_key)
  # Natural parameters decoded from the reference chart: a round trip.
  expect_reference(actual[c(row_key, "estimate")], expected[c(row_key, "estimate")],
                   reference = chart, covers = c("multilpa()", "parameter_inference.multilpa()"),
                   quantity = "estimates, constrained structures", tolerance = 1e-10,
                   measure = "relative", key = row_key)
  # Standard errors from a central-difference Hessian (step 2e-4) of the
  # reference chart, mapped by a central-difference Jacobian: 1e-5 relative,
  # the original test's tolerance for the truncation error of the differences.
  informative <- expected$standard_error > 1e-8
  expect_reference(actual[informative, c(row_key, "standard_error")],
                   expected[informative, c(row_key, "standard_error")],
                   reference = chart, covers = c("multilpa()", "parameter_inference.multilpa()"),
                   quantity = "standard errors, constrained structures", tolerance = 1e-5,
                   measure = "relative", key = row_key)
  # With one group class its probability is fixed at one: no sampling variance.
  expect_not_comparable(actual[!informative, c(row_key, "standard_error")],
                        expected[!informative, c(row_key, "standard_error")],
                        reference = chart, covers = c("multilpa()", "parameter_inference.multilpa()"),
                        quantity = "standard errors of fixed parameters, constrained structures",
                        key = row_key,
                        reason = "A single group class has probability one by definition; its standard error is zero, and the reference's numerical value is rounding noise around zero.")
})
