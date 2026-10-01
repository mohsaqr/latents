# The membership M-step: weighted multinomial logits on soft counts.

test_that("binary soft counts give the weighted logistic regression", {
  set.seed(21)
  n <- 200
  design <- cbind(1, x = rnorm(n))
  success <- plogis(-0.4 + 1.1 * design[, "x"])
  counts <- cbind(success * runif(n, 0.5, 2), (1 - success) * runif(n, 0.5, 2))
  solved <- .multilpa_weighted_logits(design, counts, matrix(0, 2L, 1L))
  # The first column is the modelled category against the last (reference).
  totals <- rowSums(counts)
  reference <- stats::glm(cbind(counts[, 1L], counts[, 2L]) ~ design[, "x"],
                          family = stats::quasibinomial(),
                          control = stats::glm.control(epsilon = 1e-14, maxit = 100))
  expect_true(solved$converged)
  # The step stops at a score of 1e-8 per unit of count, which moves a
  # coefficient by about that much relative to its information.
  expect_equal(unname(as.vector(solved$coefficients)),
               unname(stats::coef(reference)), tolerance = 1e-6)
  expect_lte(solved$scaled_score, 1e-8 * (1 + sum(totals)))
})

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
