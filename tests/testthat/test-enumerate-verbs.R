# enumerate_lpa() and enumerate_lca() are enumerate_classes(id = NULL) under
# the names of the single-level analyses, as lpa() and lca() are for one fit.

srl_indicators <- c("cognitive_strategies", "intrinsic_value",
                    "self_efficacy", "self_regulation", "test_anxiety")
esm_activities <- c("time_with_friends", "on_social_media", "sports")

test_that("enumerate_lpa preserves an omitted model argument", {
  data <- iris[c("Sepal.Length", "Sepal.Width")]
  grid <- enumerate_lpa(data, names(data), n_profiles = 1,
                        variance_model = "equal", n_starts = 1, seed = 1)
  expect_identical(as.data.frame(grid)$model, "EEI")
  expect_equal(candidate_fit(grid, 1)$log_likelihood,
               lpa(data, names(data), 1, variance_model = "equal",
                   n_starts = 1, seed = 1)$log_likelihood)
  expect_error(enumerate_lpa(data, names(data), 1, model = "basic",
                            variance_model = "equal"),
               class = "latents_bad_argument")
  categorical <- data.frame(a = rep(c("yes", "no"), 10))
  grid <- enumerate_lpa(categorical, "a", 1, categorical = "a",
                        n_starts = 1, seed = 1)
  expect_true(is.na(as.data.frame(grid)$model))
  expect_true(is.na(as.data.frame(grid)$error))
})

test_that("enumerate_lpa() is the single-level grid, fitted without a notice", {
  expect_no_message(
    models <- enumerate_lpa(srl, srl_indicators, n_profiles = 1:2,
                            n_starts = 2, seed = 1))
  expect_s3_class(models, "multilpa_enumeration")
  table <- as.data.frame(models)
  # The basic structures crossed with the counts, one group class throughout.
  expect_identical(unique(table$model), c("EEI", "VVI", "EEE", "VVV"))
  expect_identical(nrow(table), 8L)
  expect_true(all(table$n_group_classes == 1L))
  # The same grid as the two-level verb with id = NULL.
  same <- suppressMessages(enumerate_classes(srl, srl_indicators, id = NULL,
    n_profiles = 1:2, n_starts = 2, seed = 1))
  expect_equal(table, as.data.frame(same))
  # Each candidate is the fit lpa() makes.
  alone <- lpa(srl, srl_indicators, n_profiles = 2, model = "EEE",
               n_starts = 2, seed = 1)
  expect_equal(candidate_fit(models, n_profiles = 2,
                             model = "EEE")$log_likelihood,
               alone$log_likelihood)
})

test_that("enumerate_lca() treats every indicator as categorical", {
  skip_on_cran()
  expect_no_message(
    models <- enumerate_lca(student_esm, esm_activities, n_classes = 1:3,
                            n_starts = 2, seed = 1))
  table <- as.data.frame(models)
  # No covariance structure is crossed, so one candidate per count.
  expect_identical(table$n_profiles, 1:3)
  expect_true(all(is.na(table$model)))
  same <- suppressMessages(enumerate_classes(student_esm, esm_activities,
    id = NULL, n_profiles = 1:3, categorical = esm_activities, n_starts = 2,
    seed = 1))
  expect_equal(table, as.data.frame(same))
  alone <- lca(student_esm, esm_activities, n_classes = 2, n_starts = 2,
               seed = 1)
  expect_equal(candidate_fit(models, n_profiles = 2)$log_likelihood,
               alone$log_likelihood)
})

test_that("the single-level verbs refuse two-level and misplaced arguments", {
  expect_error(enumerate_lpa(srl, srl_indicators, id = "x"),
               class = "latents_bad_argument")
  expect_error(enumerate_lpa(srl, srl_indicators, n_group_classes = 2),
               class = "latents_bad_argument")
  expect_error(enumerate_lpa(srl, srl_indicators, model = "most"),
               class = "latents_bad_argument")
  expect_error(enumerate_lca(student_esm, esm_activities, id = "student"),
               class = "latents_bad_argument")
  expect_error(enumerate_lca(student_esm, esm_activities,
                             categorical = esm_activities),
               class = "latents_bad_argument")
  expect_error(enumerate_lca(student_esm, esm_activities, model = "EEE"),
               class = "latents_bad_argument")
})

test_that("the two-level verbs name the analysis they fitted with id = NULL", {
  profile_notice <- tryCatch(
    multilpa(srl, srl_indicators, id = NULL, n_profiles = 1, n_starts = 1,
             seed = 1),
    latents_single_level = conditionMessage)
  expect_match(profile_notice, "single-level latent profile model",
               fixed = TRUE)
  class_notice <- tryCatch(
    multilca(student_esm, esm_activities, id = NULL, n_profiles = 1,
             n_starts = 1, seed = 1),
    latents_single_level = conditionMessage)
  expect_match(class_notice, "single-level latent class model", fixed = TRUE)
  grid_notice <- tryCatch(
    enumerate_classes(student_esm, esm_activities, id = NULL, n_profiles = 1,
                      categorical = esm_activities, n_starts = 1, seed = 1),
    latents_single_level = conditionMessage)
  expect_match(grid_notice, "single-level latent class model", fixed = TRUE)
})

test_that("a two-level verb without `id` points to the single-level verbs", {
  expect_error(enumerate_classes(srl, srl_indicators),
               "enumerate_lpa", class = "latents_bad_argument")
  expect_error(multilpa(srl, srl_indicators, n_profiles = 2),
               "lpa()", fixed = TRUE, class = "latents_bad_argument")
  expect_error(multilca(student_esm, esm_activities, n_profiles = 2),
               "lca()", fixed = TRUE, class = "latents_bad_argument")
})
