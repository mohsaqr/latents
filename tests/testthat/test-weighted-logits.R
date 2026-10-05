
test_that("three categories converge, and the solution ignores covariate units", {
  set.seed(22)
  n <- 300
  x <- rnorm(n)
  design <- cbind(1, x = x)
  probabilities <- exp(cbind(0.3 + x, -0.2 - 0.5 * x, 0))
  probabilities <- probabilities / rowSums(probabilities)
  counts <- probabilities * runif(n, 0.2, 3)
  solved <- .multilpa_weighted_logits(design, counts, matrix(0, 2L, 2L))
  expect_true(solved$converged)
  scaled <- .multilpa_weighted_logits(cbind(1, x = x * 1e6), counts,
                                      matrix(0, 2L, 2L))
  expect_true(scaled$converged)
  expect_equal(scaled$coefficients[2L, ] * 1e6, solved$coefficients[2L, ],
               tolerance = 1e-7)
  expect_equal(scaled$coefficients[1L, ], solved$coefficients[1L, ],
               tolerance = 1e-7)
  expect_error(.multilpa_weighted_logits(design, counts, matrix(0, 2L, 2L),
                                         max_newton = 0))
})
