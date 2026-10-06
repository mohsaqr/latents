# Inference for a fit that holds part of its measurement model fixed.
#
# The references here are deliberately independent of the package: the
# conditional log likelihood is rebuilt by exhaustive enumeration of the latent
# assignments with direct dnorm arithmetic (helper-independent-likelihood.R),
# the held blocks are pasted back in by hand, and the information matrix comes
# from second central differences of that likelihood. Nothing in the reference
# path calls .multilpa_score(), .multilpa_decode() or optimHess().

#' Two-level data with two well-separated profiles and two group classes
#' @return A list with the data frame, the indicator matrix and the group vector.
fixed_measurement_data <- function() {
  set.seed(4021)
  n_groups <- 24L
  size <- 6L
  dat <- data.frame(g = rep(seq_len(n_groups), each = size))
  class_of_group <- rep(c(1, 2), length.out = n_groups)
  probability <- c(0.25, 0.75)[class_of_group][dat$g]
  profile <- rbinom(nrow(dat), 1, probability)
  dat$a <- rnorm(nrow(dat), 5 * profile, 1)
  dat$b <- rnorm(nrow(dat), 3 * profile, 1)
  list(data = dat, x = as.matrix(dat[c("a", "b")]), group = dat$g)
}

test_that("a held-measurement fit reports exactly the parameters it estimated", {
  prepared <- fixed_measurement_data()
  stage <- multilpa(prepared$data, c("a", "b"), "g", 2, 1, n_starts = 4,
                    seed = 11, tol = 1e-12)
  held <- multilpa(prepared$data, c("a", "b"), "g", 2, 2, n_starts = 3,
                   start = starting_values(stage, what = "measurement"),
                   fixed = "measurement", seed = 11, tol = 1e-12)
  expect_true(held$converged)
  expect_equal(held$fixed, c("means", "variances"))
  # The held blocks are the stage's own values, unmoved by the second stage.
  expect_equal(held$means, stage$means)
  expect_equal(held$variances, stage$variances)

  information <- parameter_inference(held, prepared$data)
  # Nothing measurement-level is reported, because nothing measurement-level was
  # estimated here.
  expect_false(any(information$level == "measurement"))
  expect_equal(nrow(information), held$n_group_classes * held$n_profiles +
                 held$n_group_classes)
  expect_equal(attr(information, "fixed"), c("means", "variances"))

  # The estimation-scale covariance is exactly n_parameters square: a held block
  # contributes no row and no column to the information matrix.
  estimation <- vcov(held, data = prepared$data, scale = "unconstrained")
  expect_equal(dim(estimation), rep(held$n_parameters, 2L))
  expect_equal(length(.multilpa_free_index(held, "unconstrained")),
               held$n_parameters)
  expect_equal(nrow(vcov(held, data = prepared$data)), nrow(information))
  expect_equal(nrow(confint(held, data = prepared$data)), nrow(information))
  expect_true(all(is.finite(information$standard_error)))
  expect_true(all(information$standard_error > 0))
})

test_that("an unheld fit takes the unconstrained path unchanged", {
  set.seed(4100)
  dat <- data.frame(g = rep(seq_len(30), each = 8), y = rnorm(240, 3, 2))
  plain <- multilpa(dat, "y", "g", 1, 1, n_starts = 1, seed = 4)
  expect_equal(plain$fixed, character(0))
  # Every coordinate is free, on both scales, and the free map is the identity.
  expect_equal(.multilpa_free_index(plain, "unconstrained"),
               seq_along(coef(plain, scale = "unconstrained")))
  expect_equal(.multilpa_free_index(plain, "natural"),
               seq_along(coef(plain)))

  information <- parameter_inference(plain, dat)
  expect_equal(nrow(information), length(coef(plain)))
  # Closed-form maximum-likelihood errors for a single Gaussian: the fixed-block
  # machinery must not have moved the ordinary answer by so much as an ulp.
  # The single profile and single group class are certain, hence the two trailing
  # probabilities of one with no error at all.
  variance <- mean((dat$y - mean(dat$y))^2)
  expect_equal(information$estimate, c(mean(dat$y), variance, 1, 1))
  expect_equal(information$standard_error,
               c(sqrt(variance / nrow(dat)), variance * sqrt(2 / nrow(dat)), 0, 0),
               tolerance = 1e-7)
  expect_equal(unname(confint(plain, data = dat)),
               unname(cbind(information$conf_low, information$conf_high)))
  expect_equal(rownames(vcov(plain, data = dat)), names(coef(plain)))

  # A fit made with an explicit empty `fixed` is the same fit, not a different
  # code path: identical estimates, identical errors, identical intervals.
  explicit <- multilpa(dat, "y", "g", 1, 1, n_starts = 1, seed = 4,
                       fixed = character())
  expect_equal(parameter_inference(explicit, dat), information)
})

test_that("inference refuses what a held fit cannot answer", {
  prepared <- fixed_measurement_data()
  stage <- multilpa(prepared$data, c("a", "b"), "g", 2, 1, n_starts = 4,
                    seed = 11, tol = 1e-12)
  held <- multilpa(prepared$data, c("a", "b"), "g", 2, 2, n_starts = 3,
                   start = starting_values(stage, what = "measurement"),
                   fixed = "measurement", seed = 11, tol = 1e-12)

  # A held coefficient is in coef(), because it is part of the model, but it was
  # not estimated here and has no interval.
  expect_error(confint(held, parm = "measurement.mean.profile_1.a",
                       data = prepared$data),
               class = "latents_held_parameter")

  # Holding every parameter the model has leaves nothing to report.
  set.seed(4055)
  small <- data.frame(g = rep(seq_len(12), each = 5), a = rnorm(60), b = rnorm(60))
  single <- multilpa(small, c("a", "b"), "g", 1, 1, n_starts = 1, seed = 2)
  everything <- multilpa(small, c("a", "b"), "g", 1, 1, n_starts = 1, seed = 2,
                         start = starting_values(single, what = "measurement"),
                         fixed = "measurement")
  expect_equal(everything$n_parameters, 0L)
  expect_error(parameter_inference(everything, small),
               class = "latents_no_free_parameters")
  expect_error(vcov(everything, data = small),
               class = "latents_no_free_parameters")

  # A `fixed` field naming something that is not a measurement block is refused
  # rather than silently treated as free.
  corrupted <- held
  corrupted$fixed <- c("means", "profile_probabilities")
  expect_error(.multilpa_free_index(corrupted, "unconstrained"),
               class = "latents_bad_fixed")
  expect_error(parameter_inference(corrupted, prepared$data),
               class = "latents_bad_fixed")
})

test_that("fit_staged makes good on its promise of conditional standard errors", {
  skip_on_cran()
  set.seed(4200)
  dat <- data.frame(g = rep(seq_len(40), each = 8))
  class_of_group <- rep(c(1, 2), length.out = 40)
  profile <- rbinom(nrow(dat), 1, c(0.2, 0.8)[class_of_group][dat$g])
  dat$a <- rnorm(nrow(dat), 4 * profile, 1)
  dat$b <- rnorm(nrow(dat), 3 * profile, 1)
  staged <- fit_staged(dat, c("a", "b"), "g", n_profiles = 2,
                       n_group_classes = 2, n_starts = 3, seed = 6)
  information <- parameter_inference(staged, dat)
  expect_false(any(information$level == "measurement"))
  expect_true(all(is.finite(information$standard_error)))
  expect_equal(dim(vcov(staged, data = dat, scale = "unconstrained")),
               rep(staged$n_parameters, 2L))
  expect_equal(nrow(confint(staged, data = dat)), nrow(information))
  expect_true(staged$n_parameters_with_measurement > staged$n_parameters)
})
