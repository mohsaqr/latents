# Review of the inference, enumeration and comparison verbs (2026-10-15).
# Each test reproduces a defect found in that review: a mistyped control that
# came back as a table of failures instead of an error, a NULL covariance, an
# error message that named the wrong count, unclassed refusals, an identifier
# check that compared storage types rather than values, and a summary that
# counted non-converged replicates as errors.

review_data <- function() {
  set.seed(7)
  data.frame(school = rep(seq_len(12), each = 10),
             score_a = c(stats::rnorm(60, -1.5), stats::rnorm(60, 1.5)),
             score_b = c(stats::rnorm(60, -1), stats::rnorm(60, 1)))
}

review_fit <- function(n_profiles = 2L, data = review_data()) {
  multilpa(data, c("score_a", "score_b"), "school", n_profiles = n_profiles,
           n_group_classes = 1L, n_starts = 2L, seed = 1)
}

test_that("sensitivity() refuses a bad control before refitting any seed", {
  fit <- review_fit()
  # Before the fix every refit failed, and the verb returned NA rows with a
  # latents_sensitivity_dropped warning instead of refusing the argument.
  expect_no_warning(expect_error(sensitivity(fit, seeds = 1:2, tol = -1),
                                 "`tol` must be"))
  expect_no_warning(expect_error(sensitivity(fit, seeds = 1:2, max_iter = "a"),
                                 "`max_iter` must be"))
  expect_no_warning(expect_error(sensitivity(fit, seeds = 1:2, max_iter = 2.5),
                                 "`max_iter` must be"))
  # Valid controls still reach every refit.
  table <- sensitivity(fit, seeds = 1:2, max_iter = 500L, tol = 1e-8)
  expect_identical(nrow(table), 2L)
  expect_false(anyNA(table$log_likelihood))
})

test_that("enumerate_classes() refuses a bad seed once, not per candidate", {
  data <- review_data()
  lapply(list(1.5, "1", c(1, 2), NA_real_), function(seed) {
    expect_error(enumerate_classes(data, c("score_a", "score_b"), "school",
                                   n_profiles = 1:2, n_group_classes = 1,
                                   model = "VVI", n_starts = 1, seed = seed),
                 class = "latents_bad_argument")
  })
})

test_that("parameter_inference() refuses a bad seed on the Wald path too", {
  fit <- review_fit()
  expect_error(parameter_inference(fit, seed = 1.5), class = "latents_bad_argument")
  expect_s3_class(parameter_inference(fit, seed = 1), "data.frame")
})

test_that("vcov() on the estimation scale refuses the bootstrap instead of NULL", {
  skip_on_cran()
  fit <- review_fit()
  expect_error(vcov(fit, scale = "unconstrained", method = "bootstrap",
                    iter = 3L, n_starts = 1L, seed = 1),
               class = "latents_unsupported_inference")
  natural <- vcov(fit, method = "bootstrap", iter = 3L, n_starts = 1L, seed = 1)
  expect_true(is.matrix(natural))
  expect_identical(rownames(natural), colnames(natural))
})

test_that("confint() refuses an unknown coefficient with a classed error", {
  fit <- review_fit()
  expect_error(confint(fit, parm = 999), class = "latents_bad_argument")
  expect_error(confint(fit, parm = 1.5), class = "latents_bad_argument")
  expect_error(confint(fit, parm = "not_a_coefficient"),
               class = "latents_bad_argument")
  expect_identical(nrow(confint(fit, parm = 1:2)), 2L)
})

test_that("bootstrap_lrt() matches groups by value, not by storage type", {
  skip_on_cran()
  data <- review_data()
  smaller <- review_fit(1L, data)
  larger <- review_fit(2L, data)
  as_double <- data
  as_double$school <- as.numeric(as_double$school)
  # Before the fix: latents_bad_inference_data, while parameter_inference()
  # accepted the same frame.
  supplied <- bootstrap_lrt(smaller, larger, data = as_double, iter = 2L,
                            n_starts = 1L, max_iter = 2000L, tol = 1e-6, seed = 1)
  stored <- bootstrap_lrt(smaller, larger, iter = 2L, n_starts = 1L,
                          max_iter = 2000L, tol = 1e-6, seed = 1)
  expect_equal(as.data.frame(supplied), as.data.frame(stored))
  # A real reordering of the groups is still refused.
  shuffled <- data
  shuffled$school <- rev(shuffled$school)
  expect_error(bootstrap_lrt(smaller, larger, data = shuffled, iter = 2L,
                             n_starts = 1L, seed = 1),
               class = "latents_bad_inference_data")
})

test_that("the bootstrap summary does not call a non-converged replicate an error", {
  skip_on_cran()
  data <- review_data()
  smaller <- review_fit(1L, data)
  larger <- review_fit(2L, data)
  # One EM step per refit: no replicate converges, none raises an error.
  result <- quietly(bootstrap_lrt(smaller, larger, iter = 2L, n_starts = 1L,
                                  max_iter = 1L, seed = 1))
  expect_true(is.na(as.data.frame(result)$p_value))
  replicates <- get_results(result, "replicates")
  expect_true(all(is.na(replicates$statistic)))
  printed <- capture.output(print(summary(result)))
  expect_false(any(grepl("raised an error", printed, fixed = TRUE)))
  expect_true(any(grepl("2 were invalid", printed, fixed = TRUE)))
})

test_that("get_candidate() on a transition grid checks counts before matching", {
  skip_on_cran()
  set.seed(2)
  sequences <- data.frame(student = rep(seq_len(20), each = 4),
                          wave = rep(seq_len(4), times = 20),
                          y = stats::rnorm(80, rep(c(-2, 2), 40)))
  grid <- quietly(enumerate_classes(sequences, "y", "student", time = "wave",
                                    n_profiles = 2:3, n_starts = 1, seed = 1,
                                    max_iter = 50),
                  c("latents_unconverged", "latents_boundary", "latents_small_classes"))
  expect_s3_class(get_candidate(grid, 2), "multilpa_transitions")
  # Before the fix: "0 candidates have 2 profiles", with the 2-profile
  # candidate present in the grid.
  expect_error(get_candidate(grid, 2.5), "`n_profiles` must be")
  expect_error(get_candidate(grid, 2, 1.5), "`n_group_classes` must be")
})

test_that("a transition or family grid's best is taken over converged candidates", {
  # The minimum of a criterion is the profile grid's rule: unconverged
  # candidates are reported, never returned as the minimum. Before the fix
  # the unconverged row below, with the lowest BIC, was returned as "best".
  table <- data.frame(candidate = c("a", "b", "c"), bic = c(10, 5, 12),
                      converged = c(TRUE, FALSE, TRUE))
  tables <- latents:::.latents_grid_result_tables
  expect_identical(tables(table, "best", "bic")$candidate, "a")
  expect_identical(tables(table, "all", "bic")$best$candidate, "a")
  none <- transform(table, converged = FALSE)
  expect_error(tables(none, "best", "bic"), class = "latents_failed_candidate")
  expect_named(tables(none, "all", "bic"), "candidates")
})

test_that("an all-unconverged transition grid names no best candidate", {
  skip_on_cran()
  grid <- quietly(enumerate_classes(subset(course_engagement, student <= 30),
                                    c("browse", "lectures"), "student",
                                    time = "sequence", n_profiles = 2:3,
                                    n_starts = 1, seed = 1, max_iter = 15))
  expect_false(any(get_results(grid, "candidates")$converged))
  expect_error(get_results(grid, "best"), class = "latents_failed_candidate")
  expect_named(get_results(grid, "all"), "candidates")
})
