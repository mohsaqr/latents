# A membership covariate that separates the profiles drives its logit
# coefficient toward infinity. The fit says so with a classed warning, and Wald
# inference is refused rather than reported with a meaningless standard error
# (Hauck and Donner, 1977). Before the refusal was explicit, it held only when
# the observed information happened to be numerically singular: seed 5 below
# reported z = -504 with standard error 669 and p = 0.45.

skip_on_cran()

.covariate_rows <- function(seed, separated = TRUE, unit = 1) {
  set.seed(seed)
  groups <- rep(seq_len(40), each = 15)
  z <- stats::rnorm(length(groups))
  profile <- if (separated) ifelse(z > 0, 2L, 1L) else
    1L + stats::rbinom(length(groups), 1L, stats::plogis(1.2 * z))
  data.frame(id = groups, z = z * unit,
             a = stats::rnorm(length(groups), c(0, 3)[profile]),
             b = stats::rnorm(length(groups), c(0, 3)[profile]))
}

.covariate_fit <- function(data) {
  multilpa(data, vars = c("a", "b"), id = "id", n_profiles = 2,
           n_group_classes = 1, profile_covariates = "z", n_starts = 2,
           seed = 1)
}

test_that("a separating covariate warns on the fit and refuses Wald inference", {
  expect_warning(fit <- .covariate_fit(.covariate_rows(5)),
                 class = "latents_extreme_coefficients")
  expect_error(parameter_inference(fit), class = "latents_boundary_fit")
  expect_error(vcov(fit), class = "latents_boundary_fit")
})

test_that("the refusal holds across simulated datasets", {
  refused <- vapply(3:8, \(seed) {
    fit <- quietly(.covariate_fit(.covariate_rows(seed)))
    inherits(tryCatch(parameter_inference(fit), latents_boundary_fit = identity),
             "latents_boundary_fit")
  }, logical(1))
  expect_true(all(refused))
})

test_that("the verdict does not depend on the covariate's units", {
  verdicts <- vapply(c(1e-3, 1, 1e3), \(unit) {
    quietly(.covariate_fit(.covariate_rows(5, unit = unit)))$extreme_logits
  }, logical(1))
  expect_true(all(verdicts))
  # A real, finite effect in tiny units has a large natural coefficient but a
  # small unit-free one: no warning, and inference is reported.
  data <- .covariate_rows(2, separated = FALSE, unit = 1e-3)
  expect_no_warning(fit <- .covariate_fit(data),
                    class = "latents_extreme_coefficients")
  expect_false(fit$extreme_logits)
  result <- as.data.frame(parameter_inference(fit))
  slope <- subset(result, parameter == "coefficient" & term == "z")
  expect_gt(abs(slope$estimate), 20)
  expect_true(is.finite(slope$standard_error))
})

test_that("the helper measures coefficients in root-mean-square units", {
  design <- list(cbind(1, c(-2, 2, -2, 2)))
  group_design <- matrix(1, 2L, 1L)
  none <- matrix(numeric(0), 1L, 0L)
  # Column scale 2: a natural coefficient of 11 is 22 unit-free.
  expect_true(.multilpa_cov_extreme(matrix(c(0, 11)), none, design, group_design))
  expect_false(.multilpa_cov_extreme(matrix(c(0, 9)), none, design, group_design))
  expect_error(.multilpa_cov_extreme(matrix(0), none, list(), group_design))
})
