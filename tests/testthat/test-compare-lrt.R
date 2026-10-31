# compare_models(test = "lrt"): the chi-square test of a model against the
# next, and the nesting rule that decides when it applies. The comparisons
# with mclust, glca and Mplus are in equivalence/test-lrt.R.

skip_on_cran()

eruption_fit <- function(model, data = datasets::faithful, ...) {
  lpa(data, c("eruptions", "waiting"), n_profiles = 2, model = model, n_starts = 3,
      seed = 1, ...)
}

test_that("each row is tested against the one before, by the chi-square formula", {
  restricted <- eruption_fit("EEI")
  varying <- eruption_fit("VVI")
  full <- eruption_fit("VVV")
  tested <- compare_models(EEI = restricted, VVI = varying, VVV = full, test = "lrt")
  table <- as.data.frame(tested)
  expect_true(all(c("lr_statistic", "lr_df", "lr_p_value") %in% names(table)))
  expect_true(is.na(table$lr_statistic[1L]) && is.na(table$lr_p_value[1L]))
  statistic <- 2 * (c(varying$log_likelihood, full$log_likelihood) -
                      c(restricted$log_likelihood, varying$log_likelihood))
  df <- c(varying$n_parameters, full$n_parameters) -
    c(restricted$n_parameters, varying$n_parameters)
  expect_equal(table$lr_statistic[2:3], statistic, tolerance = 1e-12)
  expect_identical(table$lr_df[2:3], as.integer(df))
  expect_equal(table$lr_p_value[2:3], stats::pchisq(statistic, df, lower.tail = FALSE),
               tolerance = 1e-12)
  # The information criteria are the comparison's own, with or without a test.
  plain <- as.data.frame(compare_models(EEI = restricted, VVI = varying, VVV = full))
  expect_identical(table[names(plain)], plain, ignore_attr = "models")
  expect_output(print(tested), "LR")
})

test_that("the statistic does not depend on the units of the indicators", {
  rescaled <- transform(datasets::faithful, waiting = waiting / 60)
  original <- as.data.frame(compare_models(eruption_fit("EEI"), eruption_fit("VVI"),
                                           test = "lrt"))
  minutes <- as.data.frame(compare_models(eruption_fit("EEI", rescaled),
                                          eruption_fit("VVI", rescaled), test = "lrt"))
  expect_equal(minutes$lr_statistic, original$lr_statistic, tolerance = 1e-6)
  expect_identical(minutes$lr_df, original$lr_df)
})

test_that("pairs that are not nested, or differ in class counts, are refused", {
  shared <- eruption_fit("EEE")
  varying <- eruption_fit("VVI")
  expect_error(compare_models(varying, shared, test = "lrt"), class = "latents_bad_nesting")
  expect_error(compare_models(shared, varying, test = "lrt"), class = "latents_bad_nesting")
  # Most general first: the order is part of the question.
  expect_error(compare_models(varying, eruption_fit("EEI"), test = "lrt"),
               class = "latents_bad_nesting")
  three <- lpa(datasets::faithful, c("eruptions", "waiting"), n_profiles = 3,
               model = "VVI", n_starts = 3, seed = 1)
  expect_error(compare_models(varying, three, test = "lrt"), class = "latents_bad_nesting")
  # Without `test`, the same pairs compare on their criteria as before.
  expect_s3_class(compare_models(varying, three), "latents_comparison")
})

test_that("covariances for chosen pairs nest between the diagonal and the full matrix", {
  activity <- c("browse", "lectures", "forum_read")
  fit <- function(...) lpa(engagement_small, activity, n_profiles = 2, n_starts = 2,
                           seed = 1, ...)
  diagonal <- fit(model = "VVI")
  paired <- fit(local_dependence = "browse:lectures")
  full <- fit(model = "VVV")
  tested <- as.data.frame(compare_models(diagonal, paired, full, test = "lrt"))
  expect_identical(tested$lr_df[2:3], c(2L, 4L))
  expect_error(compare_models(paired, fit(model = "EEE"), test = "lrt"),
               class = "latents_bad_nesting")
  expect_error(compare_models(paired, diagonal, test = "lrt"), class = "latents_bad_nesting")
})

test_that("a pair of items modelled jointly nests the same items apart", {
  first_week <- subset(student_esm, day <= 3)
  items <- c("time_with_friends", "on_social_media", "sports", "walking")
  apart <- lca(first_week, items, n_classes = 2, n_starts = 2, seed = 1)
  joint <- lca(first_week, items, n_classes = 2, n_starts = 2, seed = 1,
               local_dependence = "sports:walking")
  tested <- as.data.frame(compare_models(apart, joint, test = "lrt"))
  # Two binary items: three free probabilities a class jointly, two apart.
  expect_identical(tested$lr_df[2L], 2L)
  expect_equal(tested$lr_statistic[2L],
               2 * (joint$log_likelihood - apart$log_likelihood), tolerance = 1e-12)
  expect_error(compare_models(joint, apart, test = "lrt"), class = "latents_bad_nesting")
})

test_that("direct effects and membership covariates may only be added", {
  activity <- c("browse", "lectures")
  fit <- function(...) lpa(engagement_small, activity, n_profiles = 2, model = "EEI",
                           n_starts = 2, seed = 1, ...)
  plain <- fit()
  direct <- fit(direct_effects = "browse ~ previous_grade")
  covariate <- fit(profile_covariates = "previous_grade")
  expect_identical(as.data.frame(compare_models(plain, direct, test = "lrt"))$lr_df[2L], 1L)
  expect_identical(as.data.frame(compare_models(plain, covariate, test = "lrt"))$lr_df[2L],
                   1L)
  expect_error(compare_models(direct, plain, test = "lrt"), class = "latents_bad_nesting")
  expect_error(compare_models(direct, covariate, test = "lrt"),
               class = "latents_bad_nesting")
})

test_that("measurement invariance across occasions is tested in lta()", {
  set.seed(11)
  persons <- 300L
  first <- sample.int(2L, persons, replace = TRUE)
  second <- ifelse(stats::runif(persons) < 0.8, first, 3L - first)
  answer <- function(state, shift) {
    factor(stats::rbinom(length(state), 1L, ifelse(state == 1L, 0.2, 0.8 - shift)))
  }
  panel <- rbind(
    data.frame(person = seq_len(persons), occasion = 1L, u1 = answer(first, 0),
               u2 = answer(first, 0), u3 = answer(first, 0)),
    data.frame(person = seq_len(persons), occasion = 2L, u1 = answer(second, 0.3),
               u2 = answer(second, 0), u3 = answer(second, 0)))
  fit <- function(...) lta(panel, c("u1", "u2", "u3"), "person", n_profiles = 2,
                           time = "occasion", categorical = c("u1", "u2", "u3"),
                           n_starts = 3, seed = 1, ...)
  invariant <- fit(measurement = "invariant")
  occasion <- fit(measurement = "occasion")
  tested <- as.data.frame(compare_models(invariant, occasion, test = "lrt"))
  expect_identical(tested$lr_df[2L], 6L)
  expect_equal(tested$lr_statistic[2L],
               2 * (occasion$log_likelihood - invariant$log_likelihood), tolerance = 1e-12)
  expect_error(compare_models(occasion, invariant, test = "lrt"),
               class = "latents_bad_nesting")
  # A stayer class has its null value on the boundary: the chi-square
  # reference does not hold, and the bootstrap is the test for it.
  expect_error(compare_models(invariant, fit(mover_stayer = TRUE), test = "lrt"),
               class = "latents_bad_nesting")
})

test_that("weighted, unconverged and regression fits are refused", {
  weighted_data <- transform(datasets::faithful, w = rep(c(1, 2), length.out = 272))
  weighted <- function(model) {
    lpa(weighted_data, c("eruptions", "waiting"), n_profiles = 2, model = model,
        n_starts = 2, seed = 1, weights = "w")
  }
  expect_error(compare_models(weighted("EEI"), weighted("VVI"), test = "lrt"),
               class = "latents_unsupported_weights")
  unfinished <- quietly(eruption_fit("VVI", max_iter = 2), "latents_unconverged")
  expect_error(compare_models(eruption_fit("EEI"), unfinished, test = "lrt"),
               class = "latents_no_converge")
  few <- subset(growth_scores, student <= 30)
  trajectories <- mixture_regression(score ~ wave, few, n_classes = 2, id = "student",
                                     class_level = "group", n_starts = 2, seed = 1)
  expect_error(compare_models(trajectories, trajectories, test = "lrt"),
               class = "latents_bad_argument")
  expect_error(compare_models(eruption_fit("EEI"), eruption_fit("VVI"), test = "wald"),
               class = "simpleError")
})

test_that("a model passed as `test` is still a model", {
  restricted <- eruption_fit("EEI")
  varying <- eruption_fit("VVI")
  table <- as.data.frame(compare_models(restricted, test = varying))
  expect_identical(table$model, c("restricted", "test"))
  expect_false("lr_statistic" %in% names(table))
})
