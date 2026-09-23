# One seed rule for every verb that takes a seed: any whole number that
# set.seed() accepts, of either sign. Before this, multilpa(), the covariate
# fit, bootstrap_lrt() and sensitivity() refused negative seeds that set.seed()
# takes, while the bootstrap in parameter_inference() accepted a fraction that
# set.seed() truncates without notice, so seeds 1.2 and 1.7 were one stream.

seed_data <- function() {
  set.seed(3)
  data.frame(g = rep(seq_len(12), each = 8),
             a = c(stats::rnorm(48, -1.5), stats::rnorm(48, 1.5)),
             b = c(stats::rnorm(48, -1), stats::rnorm(48, 1)),
             x = stats::rnorm(96))
}

seed_fit <- function(seed, n_profiles = 2L, data = seed_data(), n_starts = 3L,
                     ...) {
  multilpa(data, c("a", "b"), "g", n_profiles = n_profiles,
           n_group_classes = 1L, n_starts = n_starts, seed = seed, ...)
}

# Pure noise and a short run: every start stops somewhere different, so the
# per-start log likelihoods are a fingerprint of the random stream. On the
# separated data above every start reaches one maximum and no seed shows.
seed_fingerprint <- function(seed) {
  set.seed(5)
  noise <- data.frame(g = rep(seq_len(12), each = 8),
                      a = stats::rnorm(96), b = stats::rnorm(96))
  fit <- quietly(seed_fit(seed, n_profiles = 3L, data = noise, n_starts = 4L,
                          max_iter = 3L))
  get_results(fit, "starts")
}

test_that("the seed rule is exactly the set of values set.seed() takes", {
  is_seed <- latents:::.multilpa_is_seed
  largest <- .Machine$integer.max
  expect_identical(is_seed(c(0, 1, -1, largest, -largest, 7L)), rep(TRUE, 6L))
  expect_identical(is_seed(c(1.5, -0.5, NA, Inf, -Inf, NaN, largest + 1,
                             -largest - 1)), rep(FALSE, 8L))
  expect_identical(is_seed(c("1", "2")), c(FALSE, FALSE))
  expect_identical(is_seed(TRUE), FALSE)
  expect_identical(is_seed(numeric()), logical())
  # The boundary is set.seed()'s own: one past it, set.seed() refuses too.
  expect_error(suppressWarnings(set.seed(-largest - 1)), "not a valid integer")
})

test_that("a negative seed fits, reproduces itself and is its own stream", {
  first <- seed_fingerprint(-7)
  expect_identical(seed_fingerprint(-7), first)
  # set.seed(-7) is not set.seed(7) with a sign dropped: it is another stream.
  expect_false(isTRUE(all.equal(seed_fingerprint(7), first)))
})

test_that("a supplied seed leaves the caller's random stream where it was", {
  data <- seed_data()
  set.seed(99)
  expected <- stats::runif(3)
  set.seed(99)
  seed_fit(-3, data = data)
  expect_identical(stats::runif(3), expected)
})

test_that("every verb refuses a seed set.seed() would not take exactly", {
  data <- seed_data()
  bad <- list(1.5, -0.5, .Machine$integer.max + 1, NA_real_, c(1, 2), "1")
  lapply(bad, \(seed) expect_error(seed_fit(seed),
                                   class = "latents_bad_argument"))
  expect_error(seed_fit(1.5, profile_covariates = "x"),
               class = "latents_bad_argument")
  fit <- seed_fit(1)
  expect_error(parameter_inference(fit, method = "bootstrap", iter = 2L,
                                   seed = 1.5),
               class = "latents_bad_argument")
  one_class <- multilpa(data, c("a", "b"), "g", n_profiles = 1L,
                        n_group_classes = 1L, n_starts = 1L, seed = 1)
  expect_error(bootstrap_lrt(one_class, fit, iter = 2L, seed = 0.5),
               class = "latents_bad_argument")
  expect_error(sensitivity(fit, seeds = c(1.2, 2.4)),
               class = "latents_bad_argument")
  expect_error(sensitivity(fit, seeds = c(1, .Machine$integer.max + 1)),
               class = "latents_bad_argument")
})

test_that("the covariate fit and the bootstraps accept a negative seed", {
  skip_on_cran()
  covariate <- quietly(seed_fit(-2, profile_covariates = "x"))
  expect_s3_class(covariate, "multilpa_covariates")
  expect_identical(quietly(seed_fit(-2, profile_covariates = "x"))$starts,
                   covariate$starts)

  fit <- seed_fit(1)
  first <- quietly(parameter_inference(fit, method = "bootstrap", iter = 3L,
                                       n_starts = 2L, seed = -11))
  again <- quietly(parameter_inference(fit, method = "bootstrap", iter = 3L,
                                       n_starts = 2L, seed = -11))
  expect_identical(first, again)

  one_class <- seed_fit(1, n_profiles = 1L)
  lrt <- quietly(bootstrap_lrt(one_class, fit, iter = 3L, n_starts = 2L,
                               seed = -5))
  expect_identical(quietly(bootstrap_lrt(one_class, fit, iter = 3L,
                                         n_starts = 2L, seed = -5)), lrt)

  result <- sensitivity(fit, seeds = c(-1, -2), n_starts = 2L)
  expect_identical(result$seed, c(-1L, -2L))
  expect_false(anyNA(result$log_likelihood))
})
