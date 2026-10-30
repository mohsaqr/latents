skip_on_cran()

test_that("negative-binomial dispersion identifies equal-mean classes", {
  reference <- list(n_profiles = 2L, count_means = matrix(1.5, 2, 1),
                    count_dispersion = matrix(c(.1, 5), 2, 1))
  candidate <- reference
  candidate$count_dispersion <- reference$count_dispersion[2:1, , drop = FALSE]
  signature <- function(x) .latents_measurement_signature(
    list(x), list(reference), mixing = c(.8, .2))
  a <- signature(reference)
  b <- signature(candidate)
  variance <- 1.5 + c(.1, 5) * 1.5^2
  expect_equal(a[, 2], log1p(variance / mean(variance)), tolerance = 1e-14)
  expect_equal(b[, 2], a[2:1, 2], tolerance = 1e-14)
  expect_identical(.multilpa_match_order(a, b), 2:1)
  # These are distinguishable count distributions even though the means agree.
  expect_gt(abs(dnbinom(0, mu = 1.5, size = 10) -
                  dnbinom(0, mu = 1.5, size = .2)), .2)
})

test_that("shared and varying NB dispersion produce the same signature", {
  shared <- list(count_means = matrix(c(1.5, 4, 3, 8), 2, 2),
                 count_dispersion = matrix(c(.1, 2), 1, 2))
  varying <- shared
  varying$count_dispersion <- shared$count_dispersion[c(1L, 1L), , drop = FALSE]
  a <- .multilpa_extra_signature(shared, shared)
  b <- .multilpa_extra_signature(varying, varying)
  expect_identical(a, b)
  moved <- varying
  moved$count_means <- moved$count_means[2:1, , drop = FALSE]
  moved$count_dispersion <- moved$count_dispersion[2:1, , drop = FALSE]
  c <- .multilpa_extra_signature(moved, shared)
  expect_identical(c$count, a$count[2:1, , drop = FALSE])
  expect_identical(c$count_spread, a$count_spread[2:1, , drop = FALSE])
})

test_that("zero count means keep finite NB signatures and valid permutations", {
  reference <- list(count_means = matrix(c(0, 1.5), 2, 1),
                    count_dispersion = matrix(c(.1, 5), 2, 1))
  moved <- reference
  moved$count_means <- moved$count_means[2:1, , drop = FALSE]
  moved$count_dispersion <- moved$count_dispersion[2:1, , drop = FALSE]
  a <- .latents_measurement_signature(list(reference), list(reference))
  b <- .latents_measurement_signature(list(moved), list(reference))
  expect_true(all(is.finite(a)))
  expect_identical(.multilpa_match_order(a, a), 1:2)
  expect_identical(.multilpa_match_order(a, b), 2:1)
  reference$count_means[] <- 0
  zero <- .multilpa_extra_signature(reference, reference)
  expect_identical(zero$count_spread, matrix(0, 2, 1))
})

test_that("Poisson count signatures retain their existing shape and arithmetic", {
  block <- list(count_means = matrix(c(1.5, 4, 3, 8), 2, 2))
  result <- .multilpa_extra_signature(block, block)
  expect_named(result, 'count')
  expect_identical(result$count,
                   sweep(block$count_means, 2, sqrt(colMeans(block$count_means)), '/'))
})

test_that("correlation alone identifies full-covariance Gaussian classes", {
  reference <- list(means = matrix(0, 2, 2), variances = matrix(1, 2, 2),
    covariances = array(c(1, -.8, -.8, 1, 1, .8, .8, 1), c(2, 2, 2)))
  candidate <- reference
  candidate$covariances <- reference$covariances[, , 2:1, drop = FALSE]
  a <- .latents_measurement_signature(list(reference), list(reference), c(.8, .2))
  b <- .latents_measurement_signature(list(candidate), list(reference), c(.8, .2))
  expect_identical(a[, seq_len(4)], b[, seq_len(4)])
  expect_equal(a[, 5], c(-.8, .8))
  expect_identical(.multilpa_match_order(a, b), 2:1)
})

test_that("fitted NB classes align by dispersion across crossed prevalences", {
  set.seed(781)
  make <- function(probability) {
    state <- rbinom(600, 1, probability) + 1L
    alpha <- c(.1, 5)[state]
    data.frame(a = rnbinom(600, mu = 1.5, size = 1 / alpha),
               b = rnbinom(600, mu = 1.5, size = 1 / alpha))
  }
  data <- list(make(.2), make(.8))
  fits <- lapply(data, function(d) quietly(lca(d, c('a', 'b'), 2,
    count = c('a', 'b'), count_model = 'negative_binomial', n_starts = 1,
    seed = 19, max_iter = 1000, tol = 1e-10)))
  expect_true(all(vapply(fits, function(f) f$converged, logical(1))))
  expect_identical(order(rowMeans(fits[[1]]$count_dispersion)), 1:2)
  expect_identical(order(rowMeans(fits[[2]]$count_dispersion)), 2:1)
  expect_identical(.multilpa_pool_align(fits[[2]], fits[[1]])$profile_order, 2:1)
})
