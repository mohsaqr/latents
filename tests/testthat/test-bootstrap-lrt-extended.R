# bootstrap_lrt() on fits with missing indicators (`missing = "fiml"`) and on
# membership-covariate fits. What is tested is the simulation each replicate
# is refitted on: that it carries the observed missingness cell for cell, that
# it holds the covariates at their observed values and draws memberships from
# the fitted logits, and that the nesting rules for covariate models are kept.

skip_on_cran()

.lrt_fixture <- function(n_groups = 40L, size = 8L, seed = 3L, share = 0.15,
                         group_classes = FALSE) {
  set.seed(seed)
  n <- n_groups * size
  frame <- data.frame(g = rep(seq_len(n_groups), each = size), z = stats::rnorm(n))
  frame$w <- stats::ave(stats::rnorm(n), frame$g)
  # Optionally two kinds of group, more likely the higher their `w`, which
  # shift the profile intercept: the structure a second group class describes.
  shift <- if (group_classes) {
    kind <- stats::rbinom(n_groups, 1L, stats::plogis(1.5 * frame$w[
      match(seq_len(n_groups), frame$g)]))
    c(-1.5, 1.5)[kind[frame$g] + 1L]
  } else 0
  state <- 2L - stats::rbinom(n, 1L, stats::plogis(shift + 1.2 * frame$z))
  frame$y1 <- stats::rnorm(n, c(0, 2.5)[state])
  frame$y2 <- stats::rnorm(n, c(0, 2)[state])
  if (share > 0) {
    frame$y1[stats::runif(n) < share] <- NA
    frame$y2[stats::runif(n) < share] <- NA
  }
  frame
}

.lrt_quietly <- function(expr) {
  quietly(expr, classes = c(.multilpa_expected_warnings, "latents_failed_replicates",
                            "latents_extreme_coefficients", "latents_boundary"))
}

test_that("simulated data carry the observed missing cells, and only those", {
  frame <- .lrt_fixture()
  fit <- .lrt_quietly(multilpa(frame, c("y1", "y2"), "g", 2L, 1L, missing = "fiml",
                               n_starts = 2, seed = 1))
  set.seed(4)
  simulated <- .multilpa_carry_missingness(.multilpa_simulate(fit), frame,
                                           c("y1", "y2"))
  expect_identical(is.na(simulated$y1), is.na(frame$y1))
  expect_identical(is.na(simulated$y2), is.na(frame$y2))
  expect_identical(simulated$g, frame$g)
})

test_that("a FIML comparison is bootstrapped with the same missingness", {
  frame <- .lrt_fixture()
  one <- .lrt_quietly(multilpa(frame, c("y1", "y2"), "g", 1L, 1L, missing = "fiml",
                               n_starts = 2, seed = 1))
  two <- .lrt_quietly(multilpa(frame, c("y1", "y2"), "g", 2L, 1L, missing = "fiml",
                               n_starts = 2, seed = 1))
  # Every replicate must go through the missingness mask; count the calls
  # made from inside the verb rather than trusting the helper's own test.
  carry <- .multilpa_carry_missingness
  calls <- 0L
  local_mocked_bindings(.multilpa_carry_missingness = function(simulated, observed, vars) {
    calls <<- calls + 1L
    carry(simulated, observed, vars)
  })
  test <- .lrt_quietly(bootstrap_lrt(one, two, iter = 9, n_starts = 2, seed = 1))
  expect_identical(calls, 9L)
  result <- as.data.frame(test)
  expect_identical(result$n_valid, 9L)
  # Two well separated profiles: the observed statistic beats every replicate.
  expect_equal(result$p_value, 0.1)
  # A FIML fit against a complete-data fit is not a comparison of one model.
  complete <- .lrt_quietly(multilpa(stats::na.omit(frame), c("y1", "y2"), "g",
                                    2L, 1L, n_starts = 2, seed = 1))
  expect_error(bootstrap_lrt(one, complete, iter = 2),
               class = "latents_incomparable_models")
})

test_that("covariate simulation holds the covariates and follows the logits", {
  frame <- .lrt_fixture(n_groups = 150L, share = 0)
  fit <- .lrt_quietly(multilpa(frame, c("y1", "y2"), "g", 2L, 1L,
                               profile_covariates = "z", n_starts = 2, seed = 1))
  set.seed(6)
  simulated <- .multilpa_cov_simulate(fit, frame)
  expect_identical(simulated$z, frame$z)
  expect_identical(simulated$g, frame$g)
  # Refitting the simulated data recovers the fitted slope, which a draw that
  # ignored the covariate would flatten to zero.
  refit <- .lrt_quietly(multilpa(simulated, c("y1", "y2"), "g", 2L, 1L,
                                 profile_covariates = "z", n_starts = 2, seed = 1))
  fitted_slope <- abs(fit$profile_coefficients["z", 1L])
  expect_lt(abs(abs(refit$profile_coefficients["z", 1L]) - fitted_slope), 0.35)
})

test_that("covariate models are bootstrapped, one group class against two", {
  frame <- .lrt_fixture(n_groups = 80L, share = 0.1, group_classes = TRUE)
  null <- .lrt_quietly(multilpa(frame, c("y1", "y2"), "g", 2L, 1L,
                                profile_covariates = "z", missing = "fiml",
                                n_starts = 2, seed = 1))
  alternative <- .lrt_quietly(multilpa(frame, c("y1", "y2"), "g", 2L, 2L,
                                       profile_covariates = "z",
                                       group_covariates = "w", missing = "fiml",
                                       n_starts = 2, seed = 1))
  # The null's own stored columns lack `w`; the verb reads the alternative's.
  test <- .lrt_quietly(bootstrap_lrt(null, alternative, iter = 5, n_starts = 2,
                                     seed = 1))
  expect_identical(nrow(get_results(test, "replicates")), 5L)
  expect_true(all(get_results(test, "replicates")$valid))
})

test_that("covariate models that are not nested are refused by class", {
  frame <- .lrt_fixture(share = 0)
  fit <- function(...) {
    .lrt_quietly(multilpa(frame, c("y1", "y2"), "g", n_starts = 2, seed = 1, ...))
  }
  with_z <- fit(n_profiles = 2L, n_group_classes = 1L, profile_covariates = "z")
  with_w <- fit(n_profiles = 3L, n_group_classes = 1L, profile_covariates = "w")
  expect_error(bootstrap_lrt(with_z, with_w, iter = 2),
               class = "latents_bad_nesting")
  two_shared <- fit(n_profiles = 2L, n_group_classes = 2L, profile_covariates = "z")
  three_varying <- fit(n_profiles = 2L, n_group_classes = 3L, profile_covariates = "z",
                       profile_slopes = "group_class")
  expect_error(bootstrap_lrt(two_shared, three_varying, iter = 2),
               class = "latents_bad_nesting")
  plain <- fit(n_profiles = 3L, n_group_classes = 1L)
  expect_error(bootstrap_lrt(with_z, plain, iter = 2))
})

test_that("a covariate replicate is valid once its likelihood has converged", {
  frame <- .lrt_fixture(n_groups = 50L, share = 0)
  null <- .lrt_quietly(multilpa(frame, c("y1", "y2"), "g", 2L, 1L,
                                profile_covariates = "z", n_starts = 3, seed = 1))
  alternative <- .lrt_quietly(multilpa(frame, c("y1", "y2"), "g", 3L, 1L,
                                       profile_covariates = "z", n_starts = 3,
                                       seed = 1))
  expect_true(isTRUE(null$likelihood_converged))
  test <- .lrt_quietly(bootstrap_lrt(null, alternative, iter = 12, n_starts = 3,
                                     max_iter = 3000, seed = 2))
  replicates <- get_results(test, "replicates")
  expect_true("logits_settled" %in% names(replicates))
  # A replicate is refused only if a likelihood failed to settle, not because
  # a logit of the over-fitted third profile kept drifting.
  expect_true(all(replicates$valid))
  # On this seed one replicate's logits did not settle (measured), so the test
  # exercises the rule rather than passing because no logit drifted.
  expect_gte(sum(!replicates$logits_settled), 1L)
  expect_true(is.finite(as.data.frame(test)$p_value))
})
