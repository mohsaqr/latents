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
  # Class sizes are estimates of the population and keep every row.
  expect_equal(sum(fit$effective_profile_counts), nrow(setup$data),
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
