# Comparisons with independent references moved from tests/testthat/test-enumeration.R,
# so that they run in CI's equivalence job rather than under R CMD check.

test_that("Kneedle reproduces a hand-worked curve and a published value", {
  # Rescaled: x = 0, 1/3, 2/3, 1 and y = 1, 1/9, 1/18, 0. The line is 1 - x,
  # so the points lie 0, 5/9, 5/18, 0 below it, and perpendicular to it
  # those distances shrink by sqrt(2).
  knee <- .latents_kneedle(1:4, c(10, 2, 1.5, 1))
  expect_identical(knee$at, 2L)
  expect_equal(knee$distance, c(0, 5 / 9, 5 / 18, 0) / sqrt(2))
  # A two-level LCA grid's BIC (groups) at 4 group classes, 2 to 8 classes:
  # Kneedle was reported at 4 classes with a perpendicular distance of 0.1915.
  bic <- c(1143453, 1101283, 1083371, 1070010, 1060205, 1051111, 1044007)
  published <- .latents_kneedle(2:8, bic)
  expect_identical((2:8)[published$at], 4L)
  expect_equal(round(published$distance[published$at], 4), 0.1915)
  # A straight line has no knee; nor do fewer than three points.
  expect_identical(.latents_kneedle(1:5, 10:6)$at, NA_integer_)
  expect_identical(.latents_kneedle(1:2, c(5, 1))$at, NA_integer_)
  # The knee does not depend on the criterion's origin or units.
  expect_identical(.latents_kneedle(1:4, 1000 + 7 * c(10, 2, 1.5, 1))$at, 2L)
})
