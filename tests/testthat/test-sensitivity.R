# Seed sensitivity. What is tested is the quantity reported -- that `agreement`
# survives label switching, and that `optimum` counts basins rather than
# distinct doubles -- not that the loop runs.
#
# The label-alignment assertion below is not hypothetical. When this verb was
# first written, `.multilpa_permute_profiles()` moved the parameter blocks and
# left `subject_profiles` behind, because the bootstrap it was written for reads
# only parameters. Seeds reaching an identical maximum reported `agreement` of
# exactly 0. That is the defect these tests reproduce.

separated_data <- function(n_groups = 14L, per_group = 10L, seed = 7L) {
  set.seed(seed)
  n <- n_groups * per_group
  data.frame(school = rep(seq_len(n_groups), each = per_group),
             a = c(stats::rnorm(n / 2L, -1.6), stats::rnorm(n / 2L, 1.6)),
             b = c(stats::rnorm(n / 2L, -1.2), stats::rnorm(n / 2L, 1.2)))
}

separated_fit <- function(...) {
  multilpa(separated_data(), c("a", "b"), "school", n_profiles = 2L,
           n_group_classes = 2L, n_starts = 4L, seed = 1L, ...)
}

test_that("permuting a fit moves its assignments, not only its parameters", {
  # The regression test for the defect above, stated on the helper itself so it
  # holds for every caller rather than only through sensitivity().
  fit <- separated_fit()
  swapped <- multilpa:::.multilpa_permute_profiles(fit, c(2L, 1L))
  expect_equal(unname(swapped$means), unname(fit$means[c(2L, 1L), ]))
  # A case on profile 1 must now be on profile 2, and the posteriors must move
  # with it, or the object describes one labelling and assigns another.
  expect_identical(swapped$subject_profiles,
                   match(fit$subject_profiles, c(2L, 1L)))
  expect_equal(unname(swapped$subject_posteriors),
               unname(fit$subject_posteriors[, c(2L, 1L)]))
  expect_equal(unname(swapped$effective_profile_counts),
               unname(fit$effective_profile_counts[c(2L, 1L)]))
  # Permuting twice is the identity, which no partial permutation satisfies.
  back <- multilpa:::.multilpa_permute_profiles(swapped, c(2L, 1L))
  expect_identical(back$subject_profiles, fit$subject_profiles)
  expect_equal(unname(back$means), unname(fit$means))
})

test_that("seeds reaching the same maximum agree, whatever they numbered it", {
  skip_on_cran()
  fit <- separated_fit()
  result <- sensitivity(fit, seeds = 1:6, n_starts = 3L)
  expect_identical(nrow(result), 6L)
  expect_identical(names(result),
                   c("seed", "log_likelihood", "converged", "iterations",
                     "optimum", "best", "agreement"))
  # Well-separated profiles: one basin, and every seed finds it.
  expect_identical(unique(result$optimum), 1L)
  expect_true(all(result$best))
  # The point of the verb: identical likelihoods must give agreement 1, which
  # they only do if the arbitrary labels were matched first.
  expect_equal(result$agreement, rep(1, 6L))
})

test_that("agreement is a proportion and the reference seed reproduces itself", {
  skip_on_cran()
  fit <- separated_fit()
  result <- sensitivity(fit, seeds = c(1L, 5L, 9L), n_starts = 3L)
  expect_true(all(result$agreement >= 0 & result$agreement <= 1))
  # Seed 1 is the reference fit's own seed with the same search effort, so it is
  # the arithmetic check that the refit rebuilds the same model.
  reference_row <- subset(result, seed == 1L)
  expect_equal(reference_row$log_likelihood, fit$log_likelihood)
  expect_equal(reference_row$agreement, 1)
})

test_that("optimum counts basins, not distinct doubles", {
  # Hand-computed against the definition: sorted downwards, a new optimum begins
  # where the drop exceeds the tolerance. Values within tolerance of each other
  # are one basin even when they are not equal, which is why `==` would be wrong.
  index <- multilpa:::.multilpa_optimum_index
  expect_identical(index(c(-100, -100, -100), tolerance = 1e-4),
                   c(1L, 1L, 1L))
  # A chain each within tolerance of the last is one basin, not three.
  expect_identical(index(c(-100, -100.00005, -100.0001), tolerance = 1e-4),
                   c(1L, 1L, 1L))
  # Ordering is by likelihood, not by position: the best value gets 1.
  expect_identical(index(c(-105, -100, -110), tolerance = 1e-4),
                   c(2L, 1L, 3L))
  # A failed fit is NA and does not consume an index.
  expect_identical(index(c(-100, NA, -105), tolerance = 1e-4),
                   c(1L, NA, 2L))
  expect_identical(index(c(NA_real_, NA_real_), tolerance = 1e-4),
                   c(NA_integer_, NA_integer_))
})

test_that("a family that cannot be refit is refused by class, not attempted", {
  skip_on_cran()
  activity <- c("browse", "lectures")
  transitions <- quietly(lta(course_engagement, vars = activity, id = "student",
                             time = "sequence", n_profiles = 2L,
                             n_starts = 2L, max_iter = 300L, seed = 1L))
  expect_error(sensitivity(transitions, seeds = 1:3),
               class = "multilpa_unsupported_sensitivity")
})

test_that("the argument contract is enforced before anything is refitted", {
  fit <- separated_fit()
  expect_error(sensitivity(fit, seeds = 1L), "at least two")
  expect_error(sensitivity(fit, seeds = c(1L, 1L)), "at least two")
  expect_error(sensitivity(fit, seeds = 1:3, n_starts = 0L), "positive whole")
  expect_error(sensitivity(fit, seeds = 1:3, tolerance = -1), "positive number")
  expect_error(sensitivity(data.frame(a = 1), seeds = 1:3), "fitted `multilpa`")
})
