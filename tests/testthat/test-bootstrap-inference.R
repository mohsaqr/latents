# Cluster bootstrap inference. What is tested here is the statistics -- that
# relabelling is undone, that the resampling unit is the group, and that the
# interval agrees with the sandwich where both exist -- rather than the
# plumbing of the data frame it comes back in.

activity <- c("browse", "lectures", "forum_read")

clustered_data <- function(n_groups = 60L, per_group = 4L, seed = 11L) {
  set.seed(seed)
  n <- n_groups * per_group
  high <- rep(c(TRUE, FALSE), each = n / 2L)
  data.frame(unit = rep(seq_len(n_groups), each = per_group),
             x = stats::rnorm(n, ifelse(high, 2, -2), 1),
             y = stats::rnorm(n, ifelse(high, 1, -1), 1))
}

small_fit <- function(...) {
  multilpa(clustered_data(), c("x", "y"), "unit", n_profiles = 2L,
           n_group_classes = 1L, n_starts = 3L, seed = 1L, ...)
}

test_that("every permutation of the profiles appears exactly once", {
  permutations <- multilpa:::.multilpa_permutations
  expect_identical(dim(permutations(1L)), c(1L, 1L))
  for (k in 2:5) {
    orders <- permutations(k)
    expect_identical(nrow(orders), as.integer(factorial(k)), info = paste("k =", k))
    expect_identical(ncol(orders), k, info = paste("k =", k))
    # Each row is a permutation, and no row repeats.
    expect_true(all(apply(orders, 1L, function(row) identical(sort(row), seq_len(k)))),
                info = paste("k =", k))
    expect_identical(nrow(unique(orders)), as.integer(factorial(k)), info = paste("k =", k))
  }
})

test_that("matching recovers a relabelling that was applied on purpose", {
  fit <- small_fit()
  scrambled <- multilpa:::.multilpa_permute_profiles(fit, c(2L, 1L))
  order <- multilpa:::.multilpa_match_order(
    multilpa:::.multilpa_profile_signature(fit, fit),
    multilpa:::.multilpa_profile_signature(scrambled, fit))
  expect_identical(order, c(2L, 1L))
  # Undoing it returns the means the fit started with.
  restored <- multilpa:::.multilpa_permute_profiles(scrambled, order)
  expect_equal(restored$means, fit$means)
})

test_that("permuting profiles is invertible", {
  fit <- small_fit()
  order <- c(2L, 1L)
  twice <- multilpa:::.multilpa_permute_profiles(
    multilpa:::.multilpa_permute_profiles(fit, order), order(order))
  expect_equal(twice$means, fit$means)
  expect_equal(twice$variances, fit$variances)
  expect_equal(twice$profile_probabilities, fit$profile_probabilities)
})

test_that("a resample carries one group per draw, whatever repeats", {
  data <- clustered_data(n_groups = 5L, per_group = 3L)
  rows_by_group <- split(seq_len(nrow(data)), data$unit)
  drawn <- c(1L, 1L, 3L, 5L, 5L)
  resampled <- multilpa:::.multilpa_resample_groups(data, rows_by_group, "unit", drawn)
  expect_identical(nrow(resampled), 15L)
  # A group drawn twice becomes two groups, not one of double the size.
  expect_identical(length(unique(resampled$unit)), 5L)
  expect_identical(unname(table(resampled$unit)), unname(table(rep(seq_len(5L), each = 3L))))
})

test_that("a refit reproduces the structure it was fitted with", {
  # The regression this guards: refitting from `variance_model` and
  # `covariance_model` alone returned a VEI fit as VVI, a wider model.
  data <- clustered_data()
  fit <- multilpa(data, c("x", "y"), "unit", n_profiles = 2L, n_group_classes = 1L,
                  n_starts = 3L, seed = 1L, volume = "varying", shape = "equal",
                  orientation = "axis")
  expect_identical(fit$covariance_structure, "VEI")
  arguments <- multilpa:::.multilpa_refit_arguments(fit)
  expect_false(any(c("variance_model", "covariance_model") %in% names(arguments)))
  again <- do.call(multilpa, c(list(data = data), arguments,
                               list(n_starts = 3L, seed = 1L)))
  expect_identical(again$covariance_structure, "VEI")
  expect_identical(again$n_parameters, fit$n_parameters)
})

test_that("a constrained structure reports intervals the Wald path refuses", {
  skip_on_cran()
  fit <- small_fit(volume = "varying", shape = "equal", orientation = "axis")
  expect_error(parameter_inference(fit), class = "multilpa_unsupported_inference")
  boot <- parameter_inference(fit, method = "bootstrap", iter = 25L,
                              n_starts = 2L, seed = 5L)
  expect_s3_class(boot, "data.frame")
  expect_identical(attr(boot, "method"), "bootstrap")
  expect_true(all(is.finite(boot$standard_error)))
  expect_true(all(boot$standard_error > 0 | boot$estimate == 1))
  expect_true(all(boot$conf_low <= boot$estimate & boot$estimate <= boot$conf_high))
  # The interval is the inference; no normal approximation is put back on top.
  expect_true(all(is.na(boot$statistic)), label = "statistic is not reported")
  expect_true(all(is.na(boot$p_value)), label = "p_value is not reported")
})

test_that("the bootstrap agrees with the sandwich where both are available", {
  skip_on_cran()
  # Groups are the resampling unit, so the cluster-robust sandwich is the
  # comparison that matters; the naive observed information is not.
  fit <- small_fit()
  robust <- parameter_inference(fit, vcov_type = "robust")
  boot <- parameter_inference(fit, method = "bootstrap", iter = 150L,
                              n_starts = 2L, seed = 9L)
  expect_identical(robust$term, boot$term)
  estimated <- robust$standard_error > 0
  ratio <- boot$standard_error[estimated] / robust$standard_error[estimated]
  ## Every parameter, not the median of them. Unaligned replicates inflate the
  ## profile means and leave the variances and the mixing probabilities alone,
  ## because those are near enough symmetric under a relabelling to survive it;
  ## a median over the whole table therefore passes while the means are fifteen
  ## times too wide.
  expect_gt(min(ratio), 0.6)
  expect_lt(max(ratio), 1.8)
})

test_that("a replicate is relabelled onto the original before it is read", {
  # Label switching is the norm, not the exception: on this data 38 of 60
  # resamples come back with the profiles the other way round, and reading them
  # unaligned inflates the standard error of every profile mean about fifteen
  # fold. This is the guard for the wiring, which a calibration band alone did
  # not catch.
  fit <- small_fit()
  scrambled <- multilpa:::.multilpa_permute_profiles(fit, c(2L, 1L))
  expect_false(isTRUE(all.equal(scrambled$means, fit$means)))
  aligned <- multilpa:::.multilpa_align_labels(scrambled, fit)
  expect_equal(aligned$means, fit$means)
  expect_equal(aligned$variances, fit$variances)
  expect_equal(aligned$profile_probabilities, fit$profile_probabilities)
})

test_that("group classes are relabelled once their profiles agree", {
  fit <- multilpa(clustered_data(), c("x", "y"), "unit", n_profiles = 2L,
                  n_group_classes = 2L, n_starts = 3L, seed = 1L)
  scrambled <- multilpa:::.multilpa_permute_group_classes(fit, c(2L, 1L))
  aligned <- multilpa:::.multilpa_align_labels(scrambled, fit)
  expect_equal(aligned$profile_probabilities, fit$profile_probabilities)
  expect_equal(aligned$group_probabilities, fit$group_probabilities)
})

test_that("a held measurement block is refused rather than resampled", {
  fit <- small_fit()
  staged <- fit_staged(clustered_data(), c("x", "y"), "unit", n_profiles = 2L,
                       n_group_classes = 2L, n_starts = 2L, seed = 1L)
  expect_error(
    parameter_inference(staged, method = "bootstrap", iter = 5L, n_starts = 1L),
    class = "multilpa_unsupported_inference")
})

test_that("the resample count is validated before anything is fitted", {
  fit <- small_fit()
  expect_error(parameter_inference(fit, method = "bootstrap", iter = 1L),
               "at least two")
  expect_error(parameter_inference(fit, method = "bootstrap", iter = 10L, n_starts = 0L),
               "positive integer")
  expect_error(parameter_inference(fit, method = "bootstrap", iter = 10L, seed = c(1, 2)),
               class = "multilpa_bad_argument")
  # set.seed() would truncate 1.5 to 1 without notice.
  expect_error(parameter_inference(fit, method = "bootstrap", iter = 10L, seed = 1.5),
               class = "multilpa_bad_argument")
})

test_that("a seed makes the table reproducible and leaves the stream alone", {
  skip_on_cran()
  fit <- small_fit()
  set.seed(99L)
  before <- stats::runif(1L)
  state <- get(".Random.seed", envir = .GlobalEnv)
  first <- parameter_inference(fit, method = "bootstrap", iter = 20L,
                               n_starts = 2L, seed = 4L)
  expect_identical(get(".Random.seed", envir = .GlobalEnv), state)
  second <- parameter_inference(fit, method = "bootstrap", iter = 20L,
                                n_starts = 2L, seed = 4L)
  expect_equal(first$standard_error, second$standard_error)
  expect_equal(first$conf_low, second$conf_low)
  # A different seed moves the replicates, so the two are not the same draw.
  third <- parameter_inference(fit, method = "bootstrap", iter = 20L,
                               n_starts = 2L, seed = 5L)
  expect_false(isTRUE(all.equal(first$standard_error, third$standard_error)))
  expect_identical(before, {
    set.seed(99L); stats::runif(1L)
  })
})

test_that("confint reports the interval parameter_inference reports", {
  skip_on_cran()
  # Not `estimate +/- z * se` rebuilt from the bootstrap standard error: that
  # would be a second, different interval for the same fit.
  fit <- small_fit(volume = "varying", shape = "equal", orientation = "axis")
  table <- parameter_inference(fit, method = "bootstrap", iter = 25L,
                               n_starts = 2L, seed = 6L)
  interval <- confint(fit, method = "bootstrap", iter = 25L, n_starts = 2L, seed = 6L)
  expect_equal(unname(interval[, 1L]), table$conf_low)
  expect_equal(unname(interval[, 2L]), table$conf_high)
})

test_that("vcov and confint reach the bootstrap for a constrained structure", {
  skip_on_cran()
  fit <- small_fit(volume = "equal", shape = "spherical")
  expect_identical(fit$covariance_structure, "EII")
  covariance <- vcov(fit, method = "bootstrap", iter = 20L, n_starts = 2L, seed = 2L)
  expect_true(is.matrix(covariance))
  expect_identical(nrow(covariance), ncol(covariance))
  # A covariance matrix of replicates is symmetric with a non-negative diagonal.
  expect_equal(covariance, t(covariance))
  expect_true(all(diag(covariance) >= 0))
  interval <- confint(fit, method = "bootstrap", iter = 20L, n_starts = 2L, seed = 2L)
  expect_identical(nrow(interval), nrow(covariance))
  expect_true(all(interval[, 1L] <= interval[, 2L]))
})

test_that("a covariate fit refuses the bootstrap rather than ignoring it", {
  # Its structures are the four the Wald coordinates already express, so there
  # is nothing to resample for; the argument exists because every method of the
  # generic carries the generic's formals.
  data <- clustered_data()
  data$w <- stats::rnorm(nrow(data))
  # Perfectly separated profiles make a covariate's logit extreme; that is the
  # fixture's doing, and it is the one warning this call legitimately raises.
  fit <- quietly(multilpa(data, c("x", "y"), "unit", n_profiles = 2L,
                          n_group_classes = 2L, profile_covariates = "w",
                          n_starts = 2L, seed = 1L))
  expect_s3_class(fit, "multilpa_covariates")
  expect_error(parameter_inference(fit, method = "bootstrap"),
               class = "multilpa_unsupported_inference")
})
