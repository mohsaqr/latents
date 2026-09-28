# multilca() is multilpa() with every indicator categorical. What is tested is
# that equivalence, the call it records, and the one argument it refuses.

.lca_items <- c("time_with_friends", "on_social_media", "tv_video_games")
.lca_data <- subset(student_esm, day <= 2)

test_that("multilca() is multilpa() with every indicator categorical", {
  lca <- quietly(multilca(.lca_data, vars = .lca_items, id = "student",
                          n_profiles = 2, n_group_classes = 2, n_starts = 1,
                          max_iter = 30, seed = 1))
  same <- quietly(multilpa(.lca_data, vars = .lca_items, id = "student",
                           n_profiles = 2, n_group_classes = 2,
                           categorical = .lca_items, n_starts = 1,
                           max_iter = 30, seed = 1))
  expect_s3_class(lca, "multilpa")
  expect_identical(lca$log_likelihood, same$log_likelihood)
  expect_identical(get_results(lca, "responses"), get_results(same, "responses"))
  expect_identical(lca$measurement_model, "categorical")
  # The object reports the call the caller wrote.
  expect_identical(as.character(lca$call[[1L]]), "multilca")
  # The fitted-model verbs work on it unchanged.
  expect_s3_class(get_results(lca, "profile_probabilities"), "data.frame")
  expect_identical(nrow(get_results(lca, "responses")),
                   2L * length(.lca_items) * 2L)
})

test_that("multilca() refuses `categorical` and checks `vars`", {
  expect_error(multilca(.lca_data, vars = .lca_items, id = "student",
                        n_profiles = 2, categorical = "sports"),
               class = "latents_bad_argument")
  expect_error(multilca(.lca_data, vars = character(), id = "student",
                        n_profiles = 2))
})

test_that("multilca() with id = NULL fits one group class without being told", {
  activities <- c("time_with_friends", "on_social_media", "tv_video_games")
  sample_data <- subset(student_esm, day <= 1)
  expect_message(
    implicit <- multilca(sample_data, vars = activities, id = NULL,
                         n_profiles = 2, n_starts = 1, seed = 1),
    class = "latents_single_level")
  explicit <- suppressMessages(
    multilca(sample_data, vars = activities, id = NULL, n_profiles = 2,
             n_group_classes = 1, n_starts = 1, seed = 1))
  expect_identical(implicit$n_group_classes, 1L)
  expect_equal(implicit$log_likelihood, explicit$log_likelihood)
  # Naming a second level on single-level data is still refused.
  expect_error(
    suppressWarnings(multilca(sample_data, vars = activities, id = NULL,
                              n_profiles = 2, n_group_classes = 2,
                              n_starts = 1, seed = 1)),
    class = "latents_bad_argument")
})
