# Correct-behavior regressions for independent audit probes of the public API.
audit_api_frame <- function() {
  withr::local_seed(101)
  data.frame(g = rep(1:20, each = 4), t = rep(1:4, 20),
             y = rnorm(80), z = rnorm(80), covariate = rnorm(80))
}

test_that("comparison classification follows remaining rows and their order", {
  d <- audit_api_frame()
  models <- quietly(enumerate_lpa(d, c("y", "z"), n_profiles = 2:3,
                                  model = "VVI", n_starts = 1, seed = 1))
  compared <- compare_models(models)
  one <- head(compared, 1)
  expect_identical(unique(get_results(one, "classification")$model), one$model)
  reversed <- compared[2:1, ]
  expect_identical(unique(get_results(reversed, "classification")$model), reversed$model)
  expect_identical(nrow(get_results(compared[FALSE, ], "classification")), 0L)
})

test_that("homogeneous and general transitions compare through grids and direct fits", {
  d <- audit_api_frame()
  homogeneous <- quietly(lta(d, c("y", "z"), "g", time = "t", n_profiles = 2,
                              n_starts = 1, seed = 1))
  models <- quietly(enumerate_classes(d, c("y", "z"), "g", time = "t",
                                      n_profiles = 2:3, n_group_classes = 1,
                                      transitions = "occasion", n_starts = 1,
                                      seed = 1, max_iter = 60))
  general <- get_candidate(models, 2)
  compared <- compare_models(homogeneous = homogeneous, general = general)
  expect_equal(compared$log_likelihood, c(homogeneous$log_likelihood, general$log_likelihood))
  expect_equal(compared$bic_groups, c(homogeneous$bic, general$bic))
  expect_equal(compared$bic_individual, c(homogeneous$bic_individual, general$bic_individual))
  expect_equal(compared$smallest_share[2], min(colMeans(general$subject_posteriors)))
  quality <- get_results(compared, "classification")
  expect_identical(unique(quality$model), compared$model)
  expect_equal(subset(quality, model == "general")$estimated_n,
               unname(colSums(general$subject_posteriors)))
  expect_identical(nrow(compare_models(models)), 2L)
  changed <- quietly(lta(transform(d, y = y + 1), c("y", "z"), "g", time = "t",
                        n_profiles = 2, transitions = "occasion", n_starts = 1,
                        seed = 1, max_iter = 60))
  expect_error(compare_models(general, changed), class = "latents_incomparable_models")
})

test_that("membership covariate enumeration carries replication counts", {
  d <- audit_api_frame()
  models <- quietly(enumerate_lpa(d, c("y", "z"), n_profiles = 2:3,
                                  profile_covariates = "covariate", model = "VVI",
                                  n_starts = 1, seed = 1, max_iter = 60))
  candidates <- get_results(models)
  expect_identical(nrow(candidates), 2L)
  # Covariate fits record replication since round 2 of 2026-10-26 (RP-1).
  expect_true(all(candidates$n_best_replicated >= 1L))
  expect_identical(candidates$model, rep("VVI", 2))
  compared <- compare_models(a = get_candidate(models, 2), b = get_candidate(models, 3))
  expect_equal(compared$bic_groups, candidates$bic_groups)
  expect_identical(unique(get_results(compared, "classification")$model), c("a", "b"))
})

test_that("all returns available tables for regression comparisons and other grids", {
  d <- audit_api_frame()
  one <- quietly(mixture_regression(y ~ z, d, n_classes = 1, seed = 1))
  two <- quietly(mixture_regression(y ~ z, d, n_classes = 2, seed = 1, max_iter = 60))
  compared <- compare_models(one = one, two = two)
  expect_named(get_results(compared, "all"), "comparison")
  expect_identical(get_results(compared, "all")$comparison, get_results(compared))
  families <- quietly(enumerate_classes(d, c("y", "z"), "g", family = "additive",
                                        n_group_classes = 1, n_starts = 1, seed = 1))
  transitions <- quietly(enumerate_classes(d, c("y", "z"), "g", time = "t",
                                           n_profiles = 2, n_group_classes = 1,
                                           n_starts = 1, seed = 1))
  for (models in list(families, transitions)) {
    tables <- get_results(models, "all")
    expect_named(tables, c("candidates", "best"))
    expect_identical(tables$candidates, get_results(models, "candidates"))
    expect_identical(tables$best, get_results(models, "best"))
  }
  failed <- transitions
  failed$table$log_likelihood <- failed$table$bic <- failed$table$aic <- NA_real_
  failed$fits <- rep(list(NULL), nrow(failed$table))
  expect_named(get_results(failed, "all"), "candidates")
})

test_that("comparison identity retains original centered values and category labels", {
  d <- audit_api_frame()
  original <- quietly(multilpa(d, c("y", "z"), "g", n_profiles = 2,
                               n_group_classes = 1, centering = "grand", n_starts = 1, seed = 1))
  shifted <- quietly(multilpa(transform(d, y = y + 100, z = z + 100),
                              c("y", "z"), "g", n_profiles = 2,
                              n_group_classes = 1, centering = "grand", n_starts = 1, seed = 1))
  expect_error(compare_models(original, shifted), class = "latents_incomparable_models")
  labels <- transform(d, y = ifelse(y > 0, "high", "low"),
                       z = ifelse(z > 0, "high", "low"))
  renamed <- transform(labels, y = ifelse(y == "high", "yes", "no"),
                        z = ifelse(z == "high", "yes", "no"))
  a <- quietly(lca(labels, c("y", "z"), n_classes = 2, n_starts = 1, seed = 1))
  b <- quietly(lca(renamed, c("y", "z"), n_classes = 2, n_starts = 1, seed = 1))
  expect_error(compare_models(a, b), class = "latents_incomparable_models")
})

test_that("noise entropy agrees between fits, candidates and comparison classification", {
  d <- audit_api_frame()
  models <- quietly(enumerate_lpa(d, c("y", "z"), n_profiles = 1:2, noise = TRUE,
                                  model = "VVI", n_starts = 1, seed = 1))
  compared <- compare_models(models)
  for (p in 1:2) {
    fit <- get_candidate(models, p)
    entropy <- subset(get_results(fit, "entropy"), level == "individuals")$relative_entropy
    expect_equal(compared$smallest_share[p],
                 min(fit$effective_profile_counts) / sum(.multilpa_class_posteriors(fit)))
    expect_equal(compared$profile_entropy[p], entropy)
    expect_equal(get_results(models)$profile_entropy[p], entropy)
    quality <- subset(get_results(compared, "classification"), model == compared$model[p])
    expect_equal(unique(quality$relative_entropy), entropy)
  }
})

test_that("smallest-class curves distinguish covariance structures", {
  skip_if_not_installed("ggplot2")
  d <- audit_api_frame()
  models <- quietly(enumerate_lpa(d, c("y", "z"), n_profiles = 2:3,
                                  model = c("EEI", "VVI"), n_starts = 1, seed = 1))
  drawn <- plot(models, what = "smallest")
  expect_setequal(unique(drawn$data$model), c("EEI", "VVI"))
  expect_identical(length(unique(drawn$data$line)), 2L)
  expect_true(all(vapply(split(drawn$data$model, drawn$data$line),
                         function(structure) length(unique(structure)) == 1L,
                         logical(1))))
  expect_identical(length(unique(ggplot2::ggplot_build(drawn)$data[[2]]$group)), 2L)
})


test_that("categorical transition comparisons ignore absent Gaussian storage types", {
  d <- audit_api_frame()
  d <- transform(d, y = ifelse(y > 0, "high", "low"),
                 z = ifelse(z > 0, "high", "low"))
  a <- quietly(lta(d, c("y", "z"), "g", time = "t", n_profiles = 2,
                   categorical = c("y", "z"), n_starts = 1, max_iter = 60, seed = 1))
  b <- quietly(lta(d, c("y", "z"), "g", time = "t", n_profiles = 2,
                   categorical = c("y", "z"), transitions = "occasion",
                   n_starts = 1, max_iter = 60, seed = 1))
  compared <- compare_models(a = a, b = b)
  expect_equal(compared$bic_groups, c(a$bic, b$bic))
  expect_true(all(is.na(compared$structure)))
  expect_identical(nrow(get_results(compared, "classification")), 4L)
})
