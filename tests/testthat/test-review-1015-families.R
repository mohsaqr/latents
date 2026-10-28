# Regressions from the 2026-10-15 review of the group-class families, the
# cross-level families and the ordinal/count indicators.

# More than 2 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

review_1015_weighted <- local({
  data <- additive_draw(
    list(means = rbind(c(0, 0), c(3, 2.5)), between = rbind(c(0.5, 0.5), c(0.5, 0.5)),
         within = c(1, 1), weights = c(0.6, 0.4)),
    sizes = rep(4L, 80L), seed = 2)
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(9)
  group_weight <- stats::setNames(stats::runif(80L, 0.3, 3), unique(data$group))
  data$w <- unname(group_weight[data$group])
  data
})

review_1015_fit <- function() {
  withCallingHandlers(
    multilpa(review_1015_weighted, c("y1", "y2"), "group", n_group_classes = 2,
             family = "additive", n_starts = 2, seed = 1, weights = "w"),
    latents_weak_class = function(w) invokeRestart("muffleWarning"))
}

test_that("weighted group-class counts are weighted, as the class weights are", {
  fit <- review_1015_fit()
  expect_false(any(fit$between_zero))
  classes <- get_results(fit, "group_classes")
  expect_equal(classes$count, unname(fit$group_probabilities) * fit$n_groups,
               tolerance = 1e-7)
  # The unweighted posterior sum is a different number on this design.
  expect_gt(max(abs(classes$count - colSums(fit$group_posteriors))), 0.1)
})

test_that("weighted effective groups use the robust covariance the warning uses", {
  fit <- review_1015_fit()
  robust <- .additive_inference(fit, "robust", warn_weak = FALSE)
  expected <- .additive_effective_groups(fit, robust)$effective_groups
  classes <- get_results(fit, "group_classes")
  expect_equal(classes$effective_groups, expected, tolerance = 1e-10)
  expect_equal(get_results(fit, "fit")$min_effective_groups, min(expected),
               tolerance = 1e-10)
})

test_that("a weighted fit's plot draws the robust intervals of its table", {
  skip_if_not_installed("ggplot2")
  fit <- review_1015_fit()
  table <- withCallingHandlers(get_results(fit, "parameters"),
                               latents_weak_class = function(w) invokeRestart("muffleWarning"))
  drawn <- withCallingHandlers(plot(fit, what = "means"),
                               latents_weak_class = function(w) invokeRestart("muffleWarning"))
  means <- subset(table, parameter == "mean")
  expect_equal(drawn$data$conf_low, means$conf_low, tolerance = 1e-10)
  expect_equal(drawn$data$conf_high, means$conf_high, tolerance = 1e-10)
})

test_that("an ordinal factor with one observed category is refused", {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(3)
  data <- data.frame(x = stats::rnorm(60),
                     o = factor(rep("a", 60), levels = c("a", "b")))
  expect_error(multilpa(data, c("x", "o"), NULL, n_profiles = 2, ordinal = "o",
                        n_starts = 1, seed = 1),
               class = "latents_bad_data")
  expect_error(.latents_prepare_extra(data, ordinal = "o"), class = "latents_bad_data")
  # Two observed categories of three levels are fitted, with a warning.
  data$o <- factor(rep(c("a", "b"), 30), levels = c("a", "b", "c"))
  expect_warning(prepared <- .latents_prepare_extra(data, ordinal = "o"),
                 class = "latents_unused_category")
  expect_identical(prepared$ordinal_categories, 3L)
})

test_that("a failed full cross-level start is surfaced, as in every other family", {
  data <- additive_draw(
    list(means = rbind(c(0, 0), c(2, 1.5)), between = rbind(c(0.5, 0.5), c(0.3, 0.3)),
         within = c(1, 1), weights = c(0.7, 0.3)),
    sizes = rep(4L, 30L), seed = 2)
  original <- .cross_level_initialize
  local_mocked_bindings(.cross_level_initialize = function(x, group_index, group_means,
                                                           n_profiles, n_group_classes,
                                                           variance_model, between_variance,
                                                           min_variance, start_index, ...) {
    if (start_index == 2L) stop("forced start failure")
    original(x, group_index, group_means, n_profiles, n_group_classes, variance_model,
             between_variance, min_variance, start_index, ...)
  })
  expect_warning(
    fit <- multilpa(data, c("y1", "y2"), "group", n_profiles = 2, n_group_classes = 2,
                    family = "full_cross_level", n_starts = 3, seed = 1),
    class = "latents_failed_starts")
  starts <- get_results(fit, "starts")
  expect_identical(starts$error[[2L]], "forced start failure")
  expect_true(is.na(starts$log_likelihood[[2L]]))
})

test_that("SQUAREM reaches the plain-EM fixed point in fewer steps", {
  data <- additive_draw(
    list(means = rbind(c(0, 0), c(2, 1.5)), between = rbind(c(0.5, 0.5), c(0.4, 0.4)),
         within = c(1, 1), weights = c(0.65, 0.35)),
    sizes = rep(5L, 60L), seed = 2)
  stats <- .additive_prepare(data, c("y1", "y2"), "group")
  structure <- .additive_structure("additive")
  start <- .additive_initialize(stats, 2L, structure, 1e-6, 1L)
  none <- matrix(FALSE, 2L, 2L)
  plain <- .additive_em(stats, start, structure, 1e-6, 5000L, 1e-12, none,
                        accelerate = FALSE)
  fast <- .additive_em(stats, start, structure, 1e-6, 5000L, 1e-12, none)
  expect_true(plain$converged && fast$converged)
  expect_lt(fast$iterations, plain$iterations)
  expect_equal(fast$expectation$log_likelihood, plain$expectation$log_likelihood,
               tolerance = 1e-10)
  expect_equal(fast$parameters$means, plain$parameters$means, tolerance = 1e-4)
  # Every SQUAREM cycle is monotone.
  expect_true(all(diff(fast$history) >= -1e-10 * (1 + abs(fast$history[-1L]))))
})

test_that("the additive bootstrap replicate that crawled now converges", {
  data <- additive_draw(
    list(means = rbind(c(0, 0), c(2, 1.5)), between = rbind(c(0.5, 0.5), c(0.3, 0.3)),
         within = c(1, 1), weights = c(0.7, 0.3)),
    sizes = rep(4L, 40L), seed = 2)
  null <- multilpa(data, c("y1", "y2"), "group", n_group_classes = 1,
                   family = "additive", n_starts = 2, seed = 1)
  alternative <- multilpa(data, c("y1", "y2"), "group", n_group_classes = 2,
                          family = "additive", n_starts = 2, seed = 1)
  test <- bootstrap_lrt(null, alternative, iter = 3, n_starts = 2, seed = 1)
  replicates <- get_results(test, "replicates")
  # Replicate 3 stopped unconverged at max_iter = 1000 under plain EM.
  expect_true(all(replicates$valid))
  expect_true(is.finite(test$p_value))
})

test_that("an ordinal factor with a never-observed level warns to drop it", {
  data <- data.frame(x = c(stats::qnorm(seq(0.01, 0.99, length.out = 60))),
                     o = factor(rep(c("a", "b"), 30), levels = c("a", "b", "c")))
  expect_warning(prepared <- .latents_prepare_extra(data, ordinal = "o"),
                 class = "latents_unused_category")
  expect_identical(prepared$ordinal_categories, 3L)
  data$o <- droplevels(data$o)
  expect_no_warning(.latents_prepare_extra(data, ordinal = "o"))
})

test_that("cross-level fits store and report the profile prevalences", {
  data <- additive_draw(
    list(means = rbind(c(0, 0), c(2, 1.5)), between = rbind(c(0.5, 0.5), c(0.3, 0.3)),
         within = c(1, 1), weights = c(0.7, 0.3)),
    sizes = rep(4L, 30L), seed = 2)
  restricted <- multilpa(data, c("y1", "y2"), "group", n_profiles = 2,
                         n_group_classes = 2, family = "restricted_cross_level",
                         n_starts = 2, seed = 1)
  prevalence <- restricted$profile_prevalences
  expect_named(prevalence, c("profile_1", "profile_2"))
  expect_equal(sum(prevalence), 1, tolerance = 1e-12)
  # Near the EM fixed point the prevalences are the mean profile posteriors
  # (the rating part stops at its relative-change rule, so not exactly).
  expect_equal(unname(prevalence), unname(colMeans(restricted$subject_posteriors)),
               tolerance = 1e-3)
  profiles <- get_results(restricted, "profiles")
  individual <- subset(profiles, level == "individual")
  group <- subset(profiles, level == "group")
  expect_equal(individual$prevalence, unname(prevalence[individual$profile]))
  expect_equal(group$prevalence, unname(restricted$group_probabilities[group$profile]))
  full <- multilpa(data, c("y1", "y2"), "group", n_profiles = 2, n_group_classes = 2,
                   family = "full_cross_level", n_starts = 2, seed = 1)
  expect_equal(unname(full$profile_prevalences),
               unname(drop(crossprod(full$group_probabilities, full$composition))),
               tolerance = 1e-12)
})

test_that("the variances plot keys classes by colour and shape, as the means plot", {
  skip_if_not_installed("ggplot2")
  data <- additive_draw(
    list(means = rbind(c(0, 0), c(2.5, 2)), between = rbind(c(0.5, 0.5), c(0.3, 0.3)),
         within = c(1, 1), weights = c(0.6, 0.4)),
    sizes = rep(5L, 60L), seed = 2)
  quiet <- function(expr) withCallingHandlers(expr, warning = function(w) {
    invokeRestart("muffleWarning")
  })
  fits <- lapply(c("additive_dispersion", "dispersion"), function(family) {
    quiet(multilpa(data, c("y1", "y2"), "group", n_group_classes = 2, family = family,
                   n_starts = 2, seed = 1))
  })
  key <- .additive_key(fits[[1L]])
  built <- ggplot2::ggplot_build(quiet(plot(fits[[1L]], what = "variances")))
  points <- built$data[[1L]]
  frame <- built$plot$data
  class_of <- as.character(frame$group_class)
  expect_identical(points$fill, key$colour[match(class_of, key$group_class)])
  expect_identical(points$shape, key$shape[match(class_of, key$group_class)])
  # Both levels are drawn, in their own panels, never told apart by colour.
  expect_setequal(levels(frame$panel), c("Within groups", "Between groups"))
  expect_identical(length(unique(built$layout$layout$PANEL)), 2L)
  # Shared blocks keep the reserved grey diamond, and no class takes it.
  shared <- ggplot2::ggplot_build(quiet(plot(fits[[2L]], what = "variances")))
  shared_frame <- shared$plot$data
  shared_points <- shared$data[[1L]][shared_frame$group_class == "shared", ]
  expect_true(all(shared_points$shape == 23L & shared_points$fill == "#999999"))
  many <- .additive_key(list(group_probabilities = stats::setNames(rep(0.125, 8L),
                                                                   paste0("group_class_", 1:8))))
  classes <- many$group_class != "shared"
  expect_false(any(many$shape[classes] == 23L))
  expect_false(any(many$colour[classes] == "#999999"))
})
