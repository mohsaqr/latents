informative_data <- function() {
  activity <- c("browse", "lectures", "forum_read")
  data <- subset(course_engagement, student <= 40,
                 c("student", "browse", "lectures", "forum_read"))
  set.seed(9)
  # A quarter of the rows observed on nothing, and some missing one indicator.
  empty <- sample(nrow(data), round(nrow(data) / 4))
  data[empty, activity] <- NA
  partial <- setdiff(sample(nrow(data), 40), empty)
  data$forum_read[partial] <- NA
  data$weight <- 1 + data$student %% 3
  list(data = data, activity = activity,
       observed = rowSums(!is.na(data[activity])) > 0)
}

test_that("classification quality is taken over rows with an observed indicator", {
  skip_on_cran()
  setup <- informative_data()
  fit <- multilpa(setup$data, setup$activity, "student", n_profiles = 3,
                  n_group_classes = 2, model = "VVI", missing = "fiml",
                  n_starts = 2, seed = 1)
  informative <- fit$subject_posteriors[setup$observed, , drop = FALSE]
  expect_identical(fit$informative, unname(setup$observed))
  expect_identical(fit$n_informative, sum(setup$observed))
  entropy <- subset(get_results(fit, "entropy"), level == "individuals")
  expect_identical(entropy$n_units, sum(setup$observed))
  expect_equal(entropy$relative_entropy,
               .multilpa_relative_entropy(informative), tolerance = 1e-12)
  # Every row with nothing observed would pull it down; they are left out.
  expect_gt(entropy$relative_entropy,
            .multilpa_relative_entropy(fit$subject_posteriors))
  classification <- subset(get_results(fit, "classification"),
                           level == "individuals")
  expect_identical(sum(classification$n_modal), sum(setup$observed))
  avepp <- diag(.multilpa_average_posterior_matrix(informative))
  expect_equal(classification$average_posterior, unname(avepp), tolerance = 1e-12)
  # The ICL charge is the entropy of the same rows its sample size counts.
  criteria <- get_results(fit, "information_criteria", format = "wide")
  expected_icl <- -2 * fit$log_likelihood +
    log(sum(setup$observed)) * fit$n_parameters +
    2 * .multilpa_entropy_sum(informative)
  expect_equal(criteria$icl_individual, expected_icl, tolerance = 1e-10)
  # Class counts are of the classified rows (drop_empty = TRUE, as Mplus).
  expect_equal(sum(fit$effective_profile_counts), sum(setup$observed),
               tolerance = 1e-8)
  expect_identical(nrow(plot_diagnostics(fit, what = "posteriors")$data),
                   sum(setup$observed))
})

test_that("rows with nothing observed change no estimate and no statistic", {
  skip_on_cran()
  setup <- informative_data()
  # Enough starts that both reach the global maximum: the two draw different
  # random starts, because their rows differ.
  fit <- function(data) multilpa(data, setup$activity, "student",
                                 n_profiles = 3, n_group_classes = 2,
                                 model = "VVI", missing = "fiml", n_starts = 20,
                                 max_iter = 5000, seed = 1, tol = 1e-10)
  with_empty <- fit(setup$data)
  dropped <- fit(setup$data[setup$observed, ])
  expect_equal(with_empty$log_likelihood, dropped$log_likelihood, tolerance = 1e-6)
  expect_equal(subset(get_results(with_empty, "entropy"), level == "individuals"),
               subset(get_results(dropped, "entropy"), level == "individuals"),
               tolerance = 1e-4, ignore_attr = TRUE)
})

test_that("the enumeration and comparison tables use the same rows", {
  skip_on_cran()
  setup <- informative_data()
  grid <- enumerate_classes(setup$data, setup$activity, "student",
                            n_profiles = 2:3, n_group_classes = 2,
                            model = "VVI", missing = "fiml", n_starts = 2,
                            seed = 1)
  fit <- get_candidate(grid, "profiles3_groups2")
  row <- subset(as.data.frame(grid), candidate == "profiles3_groups2")
  entropy <- subset(get_results(fit, "entropy"), level == "individuals")
  expect_equal(row$profile_entropy, entropy$relative_entropy, tolerance = 1e-12)
  compared <- get_results(compare_models(grid, n_group_classes = 2),
                          "comparison")
  expect_equal(subset(compared, model == "profiles3_groups2")$avepp_min,
               row$avepp_min, tolerance = 1e-12)
})

test_that("a weighted fit with empty rows weights the rows it classifies", {
  skip_on_cran()
  setup <- informative_data()
  fit <- multilpa(setup$data, setup$activity, "student", n_profiles = 2,
                  n_group_classes = 2, model = "VVI", missing = "fiml",
                  weights = "weight", n_starts = 2, seed = 1)
  classification <- subset(get_results(fit, "classification"),
                           level == "individuals")
  expect_identical(sum(classification$n_modal), sum(setup$observed))
  expect_true(all(is.finite(classification$average_posterior)))
})

test_that("complete data keep every row", {
  skip_on_cran()
  fit <- multilpa(subset(course_engagement, student <= 30),
                  c("browse", "lectures"), "student", n_profiles = 2,
                  n_group_classes = 2, n_starts = 1, seed = 1)
  expect_true(all(fit$informative))
  expect_equal(subset(get_results(fit, "entropy"), level == "individuals")$relative_entropy,
               .multilpa_relative_entropy(fit$subject_posteriors), tolerance = 1e-12)
})

test_that("several panels share a legend, one panel is labelled directly", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  grid <- enumerate_classes(subset(course_engagement, student <= 40),
                            c("browse", "lectures"), "student",
                            n_profiles = 2:3, n_group_classes = 2:3,
                            model = "VVI", n_starts = 1, seed = 1)
  text_layers <- function(p) sum(vapply(p$layers,
                                        \(l) inherits(l$geom, "GeomText"),
                                        logical(1)))
  lapply(list(plot(grid), plot(grid, what = "diagnostics")), \(p) {
    expect_identical(text_layers(p), 0L)
    expect_identical(p$theme$legend.position, "bottom")
  })
  single <- plot(grid, what = "bic_groups")
  expect_identical(text_layers(single), 1L)
  expect_identical(single$theme$legend.position, "none")
})

test_that("lta() and covariate fits count informative rows in the individual BIC", {
  skip_on_cran()
  activity <- c("browse", "lectures", "forum_read")
  set.seed(4)
  data <- subset(course_engagement, student <= 60,
                 c("student", "sequence", activity))
  data$z <- as.numeric(data$student %% 2)
  empty <- sample(nrow(data), round(nrow(data) / 5))
  data[empty, activity] <- NA
  observed <- sum(rowSums(!is.na(data[activity])) > 0)
  transition <- lta(data, activity, "student", n_profiles = 2, time = "sequence",
                    missing = "fiml", n_starts = 2, seed = 1)
  expect_identical(transition$n_informative, observed)
  expect_identical(sum(transition$informative), observed)
  expect_equal(transition$bic_individual,
               -2 * transition$log_likelihood +
                 log(observed) * transition$n_parameters, tolerance = 1e-12)
  covariate <- multilpa(data, activity, "student", n_profiles = 2,
                        n_group_classes = 2, group_covariates = "z",
                        missing = "fiml", n_starts = 2, seed = 1)
  expect_identical(covariate$n_informative, observed)
  # The stored criterion and the criteria table now agree.
  expect_equal(covariate$bic_individual,
               get_results(covariate, "information_criteria",
                           format = "wide")$bic_individual, tolerance = 1e-12)
})

test_that("a general lta() fit compares on its informative rows", {
  skip_on_cran()
  # Regression: the comparison's view of an lta() fit dropped `informative`, so
  # rows with nothing observed entered the entropy charge of ICL.
  set.seed(11)
  data <- data.frame(id = rep(1:20, each = 3), t = rep(1:3, 20), y = stats::rnorm(60))
  data$y[seq(2, 29, 3)] <- NA
  fit <- suppressWarnings(lta(data, "y", "id", n_profiles = 2, time = "t",
                              transitions = "occasion", missing = "fiml",
                              n_starts = 1, seed = 1))
  compared <- get_results(compare_models(a = fit, b = fit), "comparison")
  informative <- fit$subject_posteriors[fit$informative, , drop = FALSE]
  expect_equal(compared$icl_individual[1L],
               -2 * fit$log_likelihood + log(fit$n_informative) * fit$n_parameters +
                 2 * .multilpa_entropy_sum(informative), tolerance = 1e-10)
})

test_that("a group with nothing observed is left out of the group statistics", {
  skip_on_cran()
  data <- subset(course_engagement, student <= 30,
                 c("student", "browse", "lectures", "forum_read"))
  activity <- c("browse", "lectures", "forum_read")
  empty_students <- c(4L, 17L)
  data[data$student %in% empty_students, activity] <- NA
  data[c(1L, 50L), activity] <- NA  # single empty rows: their students stay
  fit <- multilpa(data, activity, "student", n_profiles = 2, n_group_classes = 2,
                  model = "VVI", missing = "fiml", n_starts = 2, seed = 1)
  kept <- !(fit$group_values %in% empty_students)
  expect_identical(fit$n_groups, 30L)
  expect_identical(fit$n_informative_groups, 28L)
  expect_identical(.multilpa_informative_groups(fit), kept)
  expect_equal(fit$bic_groups,
               -2 * fit$log_likelihood + log(28) * fit$n_parameters, tolerance = 1e-12)
  groups <- subset(get_results(fit, "entropy"), level == "groups")
  expect_identical(groups$n_units, 28L)
  expect_equal(groups$relative_entropy,
               .multilpa_relative_entropy(fit$group_posteriors[kept, , drop = FALSE]),
               tolerance = 1e-12)
  classification <- subset(get_results(fit, "classification"), level == "groups")
  expect_identical(sum(classification$n_modal), 28L)
  # Group class counts are of the classified groups (drop_empty = TRUE).
  expect_equal(sum(fit$effective_group_counts), 28, tolerance = 1e-8)
})

test_that("in a single-level fit an empty row is an empty group", {
  skip_on_cran()
  data <- srl[c("intrinsic_value", "test_anxiety")]
  data[1:20, ] <- NA
  fit <- lpa(data, names(data), n_profiles = 2, model = "VVI", missing = "fiml",
             n_starts = 3, seed = 1)
  expect_identical(fit$n_informative_groups, fit$n_informative)
  expect_equal(fit$bic_groups, fit$bic_individual, tolerance = 1e-12)
})

drop_empty_data <- function() {
  activity <- c("browse", "lectures", "forum_read")
  data <- subset(course_engagement, student <= 30, c("student", activity))
  set.seed(8)
  data[sample(nrow(data), 25), activity] <- NA
  data[data$student %in% c(4L, 17L), activity] <- NA
  list(data = data, activity = activity,
       observed = rowSums(!is.na(data[activity])) > 0)
}

test_that("drop_empty = TRUE counts only what was observed, FALSE counts everything", {
  skip_on_cran()
  setup <- drop_empty_data()
  fit <- function(drop) multilpa(setup$data, setup$activity, "student",
                                 n_profiles = 2, n_group_classes = 2, model = "VVI",
                                 missing = "fiml", n_starts = 2, seed = 1,
                                 drop_empty = drop)
  dropped <- fit(TRUE)
  kept <- fit(FALSE)
  # The estimates do not depend on the switch.
  expect_equal(dropped$log_likelihood, kept$log_likelihood, tolerance = 1e-12)
  expect_equal(dropped$means, kept$means, tolerance = 1e-12)
  expect_true(dropped$drop_empty)
  expect_false(kept$drop_empty)
  # Counts: the classified rows and groups, or every row and group.
  expect_equal(sum(dropped$effective_profile_counts), sum(setup$observed),
               tolerance = 1e-8)
  expect_equal(sum(dropped$effective_group_counts), 28, tolerance = 1e-8)
  expect_equal(sum(kept$effective_profile_counts), nrow(setup$data), tolerance = 1e-8)
  expect_equal(sum(kept$effective_group_counts), 30, tolerance = 1e-8)
  # Group sample size and classification statistics follow the switch.
  expect_equal(kept$bic_groups, -2 * kept$log_likelihood + log(30) * kept$n_parameters,
               tolerance = 1e-12)
  expect_equal(subset(get_results(kept, "entropy"), level == "individuals")$relative_entropy,
               .multilpa_relative_entropy(kept$subject_posteriors), tolerance = 1e-12)
  expect_identical(subset(get_results(kept, "entropy"), level == "groups")$n_units, 30L)
  # The individual sample size left empty rows out before 0.9.38 as well.
  expect_equal(kept$bic_individual, dropped$bic_individual, tolerance = 1e-12)
})

test_that("drop_empty reaches every candidate and refuses a wrong value by class", {
  skip_on_cran()
  setup <- drop_empty_data()
  grid <- enumerate_classes(setup$data, setup$activity, "student", n_profiles = 2,
                            n_group_classes = 2, model = "VVI", missing = "fiml",
                            n_starts = 2, seed = 1, drop_empty = FALSE)
  expect_false(get_candidate(grid, "profiles2_groups2")$drop_empty)
  expect_error(multilpa(setup$data, setup$activity, "student", 2, 2, missing = "fiml",
                        drop_empty = NA), class = "latents_bad_argument")
  expect_error(lta(subset(course_engagement, student <= 20), "browse", "student",
                   n_profiles = 2, time = "sequence", drop_empty = "no"),
               class = "latents_bad_argument")
})

test_that("fits with different drop_empty are not compared", {
  skip_on_cran()
  setup <- drop_empty_data()
  fit <- function(drop, k) multilpa(setup$data, setup$activity, "student",
                                    n_profiles = k, n_group_classes = 1, model = "VVI",
                                    missing = "fiml", n_starts = 2, seed = 1,
                                    drop_empty = drop)
  expect_error(compare_models(a = fit(TRUE, 2), b = fit(FALSE, 2)),
               class = "latents_incomparable_models")
  expect_error(bootstrap_lrt(fit(TRUE, 1), fit(FALSE, 2), data = setup$data, iter = 2),
               class = "latents_incomparable_models")
})
