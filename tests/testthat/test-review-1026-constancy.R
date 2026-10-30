skip_on_cran()

test_that("person scaling refuses constant units whose mean rounds", {
  # Three values of 0.1 average to 0.10000000000000002: the computed standard
  # deviation is about 1e-17, not zero, yet the unit is exactly constant.
  x <- matrix(c(0.1, 0.1, 0.1, 1, 2, 3, 5, 6, 8), ncol = 1L,
              dimnames = list(NULL, "y"))
  groups <- rep(1:3, each = 3L)
  expect_error(.multilpa_indicator_scales(x, groups, 3L, "person"),
               class = "latents_bad_data")
  sevenths <- matrix(c(rep(0.7, 3L), 1, 2, 3), ncol = 1L,
                     dimnames = list(NULL, "y"))
  expect_error(.multilpa_indicator_scales(sevenths, rep(1:2, each = 3L), 2L,
                                          "person"),
               class = "latents_bad_data")
})

test_that("person scaling accepts variation at a large origin", {
  x <- matrix(1e9 + c(0, 1, 2, 0, 2, 5), ncol = 1L, dimnames = list(NULL, "y"))
  scales <- .multilpa_indicator_scales(x, rep(1:2, each = 3L), 2L, "person")
  expect_equal(drop(scales), c(1, sd(c(0, 2, 5))), tolerance = 1e-12)
})

test_that("person scaling ignores units with no observed value", {
  x <- matrix(c(NA, NA, 1, 2, 4, 7), ncol = 1L, dimnames = list(NULL, "y"))
  scales <- .multilpa_indicator_scales(x, rep(1:3, each = 2L), 3L, "person")
  expect_identical(unname(scales[1L, 1L]), 1)
  expect_true(all(is.finite(scales)))
})

test_that("the step-three logit basis refuses a constant covariate exactly", {
  constant <- cbind("(Intercept)" = 1, z = rep(0.1, 9L))
  expect_error(.lta_three_step_basis(constant), class = "latents_bad_covariate")
  # The refusal must not depend on the constant's value.
  vapply(c(0.3, 1 / 3, 1e12 + 0.1), function(value) {
    design <- cbind("(Intercept)" = 1, z = rep(value, 10L))
    expect_error(.lta_three_step_basis(design), class = "latents_bad_covariate")
    TRUE
  }, logical(1))
})

test_that("the step-three logit basis standardizes a varying covariate at any origin", {
  z <- c(-1, 0, 1, 2)
  design <- cbind("(Intercept)" = 1, z = 1e9 + z)
  basis <- .lta_three_step_basis(design)
  # Invariant: the slope scale is the covariate's RMS spread, whatever the
  # origin, and the intercept row removes the centre in those units.
  spread <- sqrt(mean((z - mean(z))^2))
  expect_equal(basis[2L, 2L], 1 / spread, tolerance = 1e-10)
  expect_equal(basis[1L, 2L], -(1e9 + mean(z)) / spread, tolerance = 1e-10)
  expect_equal(.lta_three_step_basis(cbind("(Intercept)" = 1, z = z))[2L, 2L],
               basis[2L, 2L], tolerance = 1e-10)
})
