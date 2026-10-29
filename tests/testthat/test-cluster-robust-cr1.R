# Cluster-robust variances carry the CR1 small-sample factor G / (G - 1), the
# default of sandwich::vcovCL(), at every robust path.
# The comparisons with independent references are in tests/equivalence/test-cluster-robust-cr1.R.

test_that("the factor scales the variance by G / (G - 1) and needs two units", {
  expect_equal(.inference_cluster_factor(20), 20 / 19)
  bread <- diag(2)
  meat <- crossprod(matrix(c(1, -1, 2, -2, 0.5, -0.5), 3L, 2L))
  expect_equal(.inference_sandwich(bread, meat, 3) / meat, matrix(1.5, 2L, 2L))
  expect_error(.inference_cluster_factor(1))
  expect_error(.inference_sandwich(bread, meat, 1))
})

test_that("the correction vanishes as the number of units grows", {
  factors <- vapply(c(5, 50, 500, 5000), .inference_cluster_factor, numeric(1))
  expect_true(all(diff(factors) < 0))
  expect_lt(factors[[4L]] - 1, 1e-3)
})
