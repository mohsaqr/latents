# Labelling of mixture_regression() classes: relabelling is likelihood-neutral, and group
# classes are ordered by composition rather than by near-equal sizes.

test_that("relabelling classes and group classes leaves the likelihood unchanged", {
  fit <- mixture_regression(score ~ hours, study_hours, 2, id = "student",
                            n_group_classes = 2, membership = ~ sleep,
                            group_membership = ~ motivation, n_starts = 2, seed = 1,
                            vcov_type = "none")
  swapped <- latents:::.mixture_permute(fit$spec, fit$params, 2:1, 2:1)
  expect_equal(latents:::.mixture_expectation(fit$spec, swapped)$log_likelihood,
               fit$log_likelihood, tolerance = 1e-12)
  # Swapping twice restores the original parameters.
  back <- latents:::.mixture_permute(fit$spec, swapped, 2:1, 2:1)
  expect_equal(back$class_logits, fit$params$class_logits, tolerance = 1e-12)
  expect_equal(back$delta, fit$params$delta, tolerance = 1e-12)
})

test_that("group classes are ordered by their share of the first class", {
  fit <- mixture_regression(score ~ hours, study_hours, 2, id = "student",
                            n_group_classes = 2, n_starts = 2, seed = 1)
  composition <- subset(get_results(fit, "group_classes"), class == "class_1")
  expect_true(all(diff(composition$probability) < 0))
})

test_that("a three-class relabelling permutes every block consistently", {
  # Permutation invariance holds at any parameter point, converged or not;
  # three classes on two-class data converge slowly, so the fit is capped.
  expect_warning(
    fit <- mixture_regression(questions ~ hours, study_hours, 3, family = "poisson",
                              membership = ~ sleep, n_starts = 1, seed = 1,
                              max_iter = 30, vcov_type = "none"),
    class = "latents_unconverged")
  order <- c(3L, 1L, 2L)
  permuted <- latents:::.mixture_permute(fit$spec, fit$params, order)
  expect_equal(latents:::.mixture_expectation(fit$spec, permuted)$log_likelihood,
               fit$log_likelihood, tolerance = 1e-12)
  expect_equal(unname(permuted$beta), unname(fit$params$beta[, order]))
})
