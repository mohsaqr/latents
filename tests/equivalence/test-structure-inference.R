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
  invisible(lapply(seq_len(nrow(cases)), function(row) {
    code <- cases$code[row]
    fit <- structure_fit_at(code, data, cases$classes[row])
    expect_true(fit$converged, label = code)
    ours <- parameter_inference(fit)
    reference <- reference_structure_se(fit, data)
    expect_equal(reference$log_likelihood, fit$log_likelihood, tolerance = 1e-10,
                 info = code)
    expect_equal(ours$estimate, reference$estimate, tolerance = 1e-10, info = code)
    informative <- reference$standard_error > 1e-8
    expect_equal(ours$standard_error[informative],
                 reference$standard_error[informative], tolerance = 1e-5,
                 info = sprintf("%s with %d group classes", code, cases$classes[row]))
  }))
})
