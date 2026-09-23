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
  swapped <- latents:::.multilpa_permute_profiles(fit, c(2L, 1L))
  expect_equal(unname(swapped$means), unname(fit$means[c(2L, 1L), ]))
  expect_equal(unname(swapped$standard_deviations),
               unname(fit$standard_deviations[c(2L, 1L), ]))
  expect_identical(rownames(swapped$means), c("profile_1", "profile_2"))
  expect_identical(colnames(swapped$subject_posteriors),
                   c("profile_1", "profile_2"))
  # A case on profile 1 must now be on profile 2, and the posteriors must move
  # with it, or the object describes one labelling and assigns another.
  expect_identical(swapped$subject_profiles,
                   match(fit$subject_profiles, c(2L, 1L)))
  expect_equal(unname(swapped$subject_posteriors),
               unname(fit$subject_posteriors[, c(2L, 1L)]))
  expect_equal(unname(swapped$effective_profile_counts),
               unname(fit$effective_profile_counts[c(2L, 1L)]))
  # Permuting twice is the identity, which no partial permutation satisfies.
  back <- latents:::.multilpa_permute_profiles(swapped, c(2L, 1L))
  expect_identical(back$subject_profiles, fit$subject_profiles)
  expect_equal(unname(back$means), unname(fit$means))
})

test_that("permuting group classes keeps their posterior tables and labels aligned", {
  data <- separated_data(n_groups = 18L)
  fit <- quietly(multilpa(data, c("a", "b"), "school", n_profiles = 2L,
                           n_group_classes = 3L, n_starts = 2L,
                           max_iter = 80L, seed = 1L))
  permutation <- c(3L, 1L, 2L)
  moved <- latents:::.multilpa_permute_group_classes(fit, permutation)
  expect_equal(unname(moved$group_posteriors),
               unname(fit$group_posteriors[, permutation]))
  expect_identical(moved$group_classes,
                   match(fit$group_classes, permutation))
  expect_equal(unname(moved$effective_group_counts),
               unname(fit$effective_group_counts[permutation]))
  expect_identical(names(moved$group_probabilities),
                   paste0("group_class_", 1:3))
  expect_identical(colnames(moved$group_posteriors),
                   paste0("group_class_", 1:3))
  counts <- get_results(moved, "counts")
  expect_equal(counts$effective_count[counts$level == "groups"],
               unname(moved$effective_group_counts))
  expect_identical(get_results(moved, "assignments")$group_class,
                   moved$group_classes[moved$group_index])
  restored <- latents:::.multilpa_permute_group_classes(
    moved, order(permutation))
  expect_equal(restored$group_posteriors, fit$group_posteriors)
  expect_identical(restored$group_classes, fit$group_classes)
  expect_equal(restored$profile_probabilities, fit$profile_probabilities)
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

test_that("a short original run is replayed with its own iteration settings", {
  data <- separated_data()
  fit <- quietly(multilpa(data, c("a", "b"), "school", n_profiles = 2L,
                           n_group_classes = 2L, n_starts = 2L, seed = 1L,
                           max_iter = 2L, tol = 1e-6))
  expect_identical(fit$max_iter, 2L)
  expect_equal(fit$tol, 1e-6)
  result <- quietly(sensitivity(fit, seeds = 1:2))
  reference <- result[result$seed == 1L, ]
  expect_equal(reference$log_likelihood, fit$log_likelihood)
  expect_identical(reference$iterations, fit$iterations)
})

test_that("a categorical refit preserves declared factor coding", {
  data <- separated_data()
  data$category <- ordered(ifelse(data$a > 0, "high", "low"),
                           levels = c("high", "low"))
  fit <- quietly(multilpa(data, c("a", "category"), "school", 2L, 2L,
                           categorical = "category", n_starts = 2L,
                           max_iter = 60L, seed = 1L))
  result <- quietly(sensitivity(fit, seeds = 1:2))
  reference <- result[result$seed == 1L, ]
  expect_equal(reference$log_likelihood, fit$log_likelihood)
  expect_equal(reference$agreement, 1)
  expect_identical(levels(get_results(fit, "data")$category),
                   c("high", "low"))
})

test_that("optimum counts basins, not distinct doubles", {
  # Hand-computed against the definition: sorted downwards, a new optimum begins
  # where the drop exceeds the tolerance. Values within tolerance of each other
  # are one basin even when they are not equal, which is why `==` would be wrong.
  index <- latents:::.multilpa_optimum_index
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
  transitions <- quietly(lta(engagement_small, vars = activity, id = "student",
                             time = "sequence", n_profiles = 2L,
                             n_starts = 2L, max_iter = 300L, seed = 1L))
  expect_error(sensitivity(transitions, seeds = 1:3),
               class = "latents_unsupported_sensitivity")
})

test_that("the argument contract is enforced before anything is refitted", {
  fit <- separated_fit()
  expect_error(sensitivity(fit, seeds = 1L), class = "latents_bad_argument")
  expect_error(sensitivity(fit, seeds = c(1L, 1L)),
               class = "latents_bad_argument")
  expect_error(sensitivity(fit, seeds = c(1.2, 2.4)),
               class = "latents_bad_argument")
  expect_error(sensitivity(fit, seeds = c("1", "2")),
               class = "latents_bad_argument")
  expect_error(sensitivity(fit, seeds = 1:3, n_starts = 0L), "positive whole")
  expect_error(sensitivity(fit, seeds = 1:3, n_starts = 1.5), "positive whole")
  expect_error(sensitivity(fit, seeds = 1:3, tolerance = -1), "positive number")
  expect_error(sensitivity(data.frame(a = 1), seeds = 1:3), "fitted `multilpa`")
})

test_that("agreement refuses data whose rows no longer match the reference", {
  data <- separated_data()
  fit <- separated_fit()
  expect_error(sensitivity(fit, data = data[nrow(data):1L, ], seeds = 1:2),
               class = "latents_bad_inference_data")
  expect_error(sensitivity(fit, data = head(data, -1L), seeds = 1:2),
               class = "latents_bad_inference_data")
})

test_that("a staged sensitivity refits the conditional second stage", {
  data <- separated_data(n_groups = 18L)
  staged <- quietly(fit_staged(data, c("a", "b"), "school", n_profiles = 2L,
                               n_group_classes = 2L, n_starts = 2L,
                               max_iter = 80L, seed = 1L))
  result <- sensitivity(staged, seeds = 1:3, n_starts = 2L, max_iter = 80L)
  reference <- result[result$seed == 1L, ]
  expect_equal(reference$log_likelihood, staged$log_likelihood,
               tolerance = 1e-8)
  expect_equal(reference$agreement, 1)
  expect_true(all(is.finite(result$log_likelihood)))
})

test_that("a partially fixed fit is not silently refitted as a joint model", {
  data <- separated_data()
  initial <- separated_fit()
  held <- quietly(multilpa(data, c("a", "b"), "school", n_profiles = 2L,
                            n_group_classes = 2L, n_starts = 2L, seed = 1L,
                            start = starting_values(initial), fixed = "means"))
  expect_error(sensitivity(held, seeds = 1:2),
               class = "latents_unsupported_sensitivity")
})
