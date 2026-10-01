# Comparing transition models: simulation, nesting, bootstrap LRT, enumeration.

compare_activity <- c("browse", "lectures")

test_that("simulation keeps groups, occasions and covariates and draws the model", {
  fit <- lta(course_engagement, compare_activity, "student", n_profiles = 2,
             time = "sequence", n_starts = 2, seed = 1)
  set.seed(3)
  simulated <- .lta_simulate(fit, course_engagement)
  expect_identical(dim(simulated), dim(course_engagement))
  expect_identical(simulated$student, course_engagement$student)
  expect_identical(simulated$sequence, course_engagement$sequence)
  expect_identical(simulated$previous_grade, course_engagement$previous_grade)
  expect_false(isTRUE(all.equal(simulated$browse, course_engagement$browse)))
  # A fit to data simulated from the fit recovers its transition matrix.
  set.seed(4)
  many <- .lta_simulate(fit, course_engagement)
  refit <- lta(many, compare_activity, "student", n_profiles = 2, time = "sequence",
               n_starts = 2, seed = 1)
  stay <- function(f) sort(diag(f$transition_probabilities[, , 1L]))
  expect_lt(max(abs(stay(refit) - stay(fit))), 0.1)
})

test_that("general fits simulate with covariates, occasion measurement and order 2", {
  # Fifteen occasions with few groups late on: some moves never occur, so the
  # fit is on a probability boundary by design.
  fit <- quietly(lta(course_engagement, compare_activity, "student", n_profiles = 2,
                     time = "sequence", transition_covariates = "previous_grade",
                     measurement = "occasion", order = 2, n_starts = 1, seed = 1,
                     max_iter = 50), "latents_boundary")
  set.seed(5)
  simulated <- .lta_simulate(fit, course_engagement)
  expect_identical(nrow(simulated), nrow(course_engagement))
  expect_true(all(is.finite(simulated$browse)))
})

test_that("the nesting rule relaxes options one way only", {
  base <- list(vars = "y", id = "id", time = "t", n_profiles = 2L, n_group_classes = 1L,
               variance_model = "varying", transitions = "homogeneous",
               transition_covariates = character(), initial_covariates = character(),
               measurement = "invariant", order = 1L)
  fake <- function(...) structure(list(arguments = modifyList(base, list(...))),
                                  class = "multilpa_lta")
  expect_true(.lta_nested(fake(), fake(transitions = "occasion")))
  expect_true(.lta_nested(fake(), fake(order = 2L, transition_covariates = "x")))
  expect_true(.lta_nested(fake(), fake(n_profiles = 3L, n_group_classes = 2L)))
  expect_false(.lta_nested(fake(transitions = "occasion"), fake()))
  expect_false(.lta_nested(fake(transition_covariates = "x"),
                           fake(transition_covariates = "w")))
  expect_false(.lta_nested(fake(variance_model = "varying"),
                           fake(variance_model = "equal")))
})

test_that("bootstrap LRT between nested transition fits", {
  null <- lta(course_engagement, compare_activity, "student", n_profiles = 2,
              time = "sequence", n_starts = 2, seed = 1)
  alternative <- lta(course_engagement, compare_activity, "student", n_profiles = 2,
                     time = "sequence", transition_covariates = "previous_grade",
                     n_starts = 2, seed = 1)
  set.seed(9)
  before <- .Random.seed
  result <- bootstrap_lrt(null, alternative, data = course_engagement, iter = 4,
                          n_starts = 1, seed = 2)
  expect_identical(.Random.seed, before)
  test <- get_results(result, "test")
  expect_equal(test$statistic, max(0, 2 * (alternative$log_likelihood -
                                             null$log_likelihood)))
  expect_identical(test$n_valid, 4L)
  expect_match(test$alternative_family, "transitions on previous_grade")
  expect_error(bootstrap_lrt(alternative, null, data = course_engagement, iter = 2),
               class = "latents_bad_argument")
  expect_error(bootstrap_lrt(null, alternative, iter = 2),
               class = "latents_bad_argument")
  profile <- suppressMessages(lpa(course_engagement, compare_activity, n_profiles = 2,
                                  n_starts = 1, seed = 1))
  expect_error(bootstrap_lrt(profile, null, data = course_engagement, iter = 2),
               class = "latents_bad_argument")
})

test_that("transition enumeration crosses profiles, classes and structures", {
  grid <- enumerate_classes(course_engagement, compare_activity, "student",
                            time = "sequence", n_profiles = 2:3, n_group_classes = 1:2,
                            n_starts = 2, seed = 1)
  expect_s3_class(grid, "latents_transition_enumeration")
  table <- get_results(grid)
  expect_identical(nrow(table), 4L)
  fit <- candidate_fit(grid, n_profiles = 3, n_group_classes = 1)
  expect_equal(subset(table, n_profiles == 3 & n_group_classes == 1)$bic, fit$bic)
  expect_equal(get_results(grid, "best")$bic, min(table$bic))
  expect_identical(as.data.frame(grid), table)
  structures <- enumerate_classes(course_engagement, compare_activity, "student",
                                  time = "sequence", n_profiles = 2,
                                  model = c("VVI", "VEI"), n_starts = 2, seed = 1)
  expect_identical(get_results(structures)$model, c("VVI", "VEI"))
  expect_error(candidate_fit(structures, n_profiles = 2, n_group_classes = 1),
               class = "latents_unknown_candidate")
  expect_error(enumerate_classes(course_engagement, compare_activity, "student",
                                 time = "sequence", n_profiles = 1:2),
               class = "latents_bad_transition")
  expect_error(enumerate_classes(course_engagement, compare_activity, "student",
                                 time = "sequence", no_such_argument = 1),
               class = "latents_bad_argument")
  skip_if_not_installed("ggplot2")
  expect_true(ggplot2::is_ggplot(plot(grid)))
})
