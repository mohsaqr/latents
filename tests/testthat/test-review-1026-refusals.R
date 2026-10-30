# Review 2026-10-26, refusal parity: an input that `multilpa()` refuses is
# refused by every entry point, with the same class, and for the right reason.

skip_on_cran()

refusal_data <- function() {
  withr::local_seed(20261026)
  data.frame(g = rep(seq_len(20L), each = 4L), y1 = stats::rnorm(80L),
             y2 = stats::rnorm(80L), x = stats::rnorm(80L),
             c1 = factor(sample(c("a", "b"), 80L, TRUE)),
             c2 = factor(sample(c("p", "q"), 80L, TRUE)))
}

test_that("an enumeration refuses data that every candidate refuses, as multilpa does", {
  data <- refusal_data()
  # A misspelled column used to come back as a grid of failed cells, with no
  # condition at all.
  expect_error(enumerate_lpa(data, c("y1", "nope"), n_profiles = 1:2, n_starts = 1,
                             seed = 1), class = "latents_bad_data")
  constant <- within(data, y1 <- 0.1)
  expect_error(multilpa(constant, vars = c("y1", "y2"), id = "g", n_profiles = 2,
                        n_group_classes = 1, n_starts = 1, seed = 1),
               class = "latents_bad_data")
  expect_error(enumerate_lpa(constant, c("y1", "y2"), n_profiles = 1:2, n_starts = 1,
                             seed = 1), class = "latents_bad_data")
  expect_error(enumerate_classes(constant, c("y1", "y2"), id = "g", n_profiles = 1:2,
                                 n_starts = 1, seed = 1), class = "latents_bad_data")
  one_level <- within(data, c1 <- factor(rep("a", nrow(data))))
  expect_error(enumerate_lca(one_level, c("c1", "c2"), n_classes = 1:2, n_starts = 1,
                             seed = 1), class = "latents_bad_categorical")
})

test_that("a candidate that fails on its own still leaves the rest of the grid", {
  data <- refusal_data()[seq_len(12L), ]
  # Too many profiles for twelve rows fails some candidates for reasons of the
  # model, not of the data, so the enumeration reports them and goes on.
  grid <- quietly(enumerate_lpa(data, c("y1", "y2"), n_profiles = c(1L, 6L),
                                n_starts = 1, seed = 1))
  candidates <- get_results(grid, "candidates")
  expect_true(any(!is.na(candidates$log_likelihood)))
})

test_that("the shared-error rule needs the same class and message in every cell", {
  bad <- errorCondition("the data are bad", class = "latents_bad_data")
  other <- errorCondition("something else", class = "latents_bad_data")
  model <- errorCondition("too many classes", class = "latents_no_valid_start")
  expect_error(.latents_refuse_shared_data_error(list(bad, bad)),
               class = "latents_bad_data")
  expect_null(.latents_refuse_shared_data_error(list(bad, other)))
  expect_null(.latents_refuse_shared_data_error(list(model, model)))
  expect_null(.latents_refuse_shared_data_error(list(bad, list())))
})

test_that("a constant Gaussian outcome is named as the reason, not the starts", {
  data <- within(refusal_data(), y1 <- 0.1)
  failure <- expect_error(mixture_regression(y1 ~ x, data, n_classes = 2, n_starts = 2,
                                             seed = 1), class = "latents_bad_data")
  expect_match(conditionMessage(failure), "same value")
  # Even one class could not fit it, so this is about the data.
  expect_error(mixture_regression(y1 ~ x, data, n_classes = 1, n_starts = 0, seed = 1),
               class = "latents_bad_data")
  # A constant outcome in new data for prediction is legitimate.
  fit <- mixture_regression(score ~ hours, study_hours, n_classes = 1, n_starts = 0,
                            seed = 1)
  expect_identical(nrow(predict(fit, newdata = transform(study_hours[1:5, ], score = 50))), 5L)
})

test_that("more cores than the machine has are capped with a classed warning", {
  available <- parallel::detectCores()
  skip_if(!is.finite(available) || available < 1L)
  expect_identical(.latents_check_cores(1L), 1L)
  expect_warning(capped <- .latents_check_cores(available + 14L),
                 class = "latents_too_many_cores")
  expect_identical(capped, as.integer(available))
})
