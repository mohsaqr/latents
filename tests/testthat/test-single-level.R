# `id = NULL`. What is tested here is that the single-level fit is the same
# model the two-level one reduces to, not that a code path runs: an explicit
# one-row-per-unit fit and an `id = NULL` fit must agree to the last digit.
#
# `id` has no default, so a forgotten grouping still raises. `id = NULL` is a
# statement the caller makes, and it warns, because this is a package for the
# two-level model. `quietly()` muffles that expected warning at each call.

single_level_data <- function(n = 160L, seed = 1L) {
  set.seed(seed)
  data.frame(x = c(stats::rnorm(n / 2L, -2), stats::rnorm(n / 2L, 2)),
             y = stats::rnorm(n))
}

flat_fit <- function(..., data = single_level_data()) {
  quietly(multilpa(data, c("x", "y"), id = NULL, n_profiles = 2L, ...))
}

test_that("a forgotten id raises, and a stated one warns", {
  data <- single_level_data()
  # The hazard this guards: before `id = NULL` existed, omitting `id` raised.
  # It must keep raising, or a forgotten grouping silently becomes a different
  # model.
  expect_error(multilpa(data, c("x", "y"), n_profiles = 2L, n_starts = 2L),
               class = "multilpa_bad_argument")
  expect_error(multilpa(data, c("x", "y"), n_profiles = 2L, n_starts = 2L),
               "pass `id = NULL`")
  # Saying it explicitly works, and says so back.
  expect_warning(multilpa(data, c("x", "y"), id = NULL, n_profiles = 2L,
                          n_starts = 2L, seed = 1L),
                 class = "multilpa_single_level")
})

test_that("id = NULL is the fit you get by numbering the rows yourself", {
  data <- single_level_data()
  implicit <- flat_fit(n_starts = 5L, seed = 1L, data = data)
  explicit_data <- data
  explicit_data$unit <- seq_len(nrow(explicit_data))
  explicit <- multilpa(explicit_data, c("x", "y"), id = "unit", n_profiles = 2L,
                       n_group_classes = 1L, n_starts = 5L, seed = 1L)
  expect_equal(implicit$log_likelihood, explicit$log_likelihood)
  expect_equal(unname(implicit$means), unname(explicit$means))
  expect_equal(unname(implicit$variances), unname(explicit$variances))
  expect_identical(implicit$n_parameters, explicit$n_parameters)
  expect_equal(unname(implicit$subject_posteriors), unname(explicit$subject_posteriors))
})

test_that("a single-level fit records what it is", {
  fit <- flat_fit(n_starts = 3L, seed = 1L)
  expect_true(fit$single_level)
  expect_identical(fit$n_group_classes, 1L)
  expect_identical(fit$n_groups, fit$n_observations)
  # It says so, rather than reporting one group class and 160 groups of one.
  printed <- utils::capture.output(print(fit))
  expect_true(any(grepl("Latent profile analysis: 2 profiles", printed, fixed = TRUE)))
  expect_true(any(grepl("160 observations", printed, fixed = TRUE)))
  expect_false(any(grepl("group class", printed, fixed = TRUE)))
})

test_that("the fabricated unit column is not handed back", {
  fit <- flat_fit(n_starts = 3L, seed = 1L)
  # The caller never wrote `.observation`, so no table may return it.
  expect_identical(names(get_data(fit, "data")), c("x", "y"))
  expect_false(".observation" %in% names(get_data(fit, "assignments")))
  # A two-level fit still reports the identifier its caller did supply.
  data <- single_level_data()
  data$unit <- rep(seq_len(40L), each = 4L)
  grouped <- multilpa(data, c("x", "y"), id = "unit", n_profiles = 2L,
                      n_group_classes = 2L, n_starts = 3L, seed = 1L)
  expect_true("unit" %in% names(get_data(grouped, "data")))
})

test_that("group classes without an id are refused, not silently dropped", {
  expect_error(
    multilpa(single_level_data(), c("x", "y"), id = NULL, n_profiles = 2L,
             n_group_classes = 2L, n_starts = 2L),
    class = "multilpa_bad_argument")
  # Naming one explicitly is the same request as leaving it alone, so it passes.
  expect_s3_class(flat_fit(n_group_classes = 1L, n_starts = 2L, seed = 1L),
                  "multilpa")
})

test_that("a column already called .observation is refused", {
  data <- single_level_data()
  data$.observation <- seq_len(nrow(data))
  expect_error(
    multilpa(data, c("x", "y"), id = NULL, n_profiles = 2L, n_starts = 2L),
    class = "multilpa_bad_data")
})

test_that("the package's verbs work on a single-level fit", {
  skip_on_cran()
  fit <- flat_fit(n_starts = 5L, seed = 1L)
  # Inference, classification and the tables, on a fit with no second level.
  inference <- parameter_inference(fit)
  expect_true(all(is.finite(inference$standard_error)))
  expect_identical(nrow(get_data(fit, "profiles")), 4L)
  # Both levels are still reported: two profiles, and the single group class the
  # single-level fit collapses to, which carries the whole sample.
  counts <- get_data(fit, "counts")
  expect_identical(sum(counts$level == "individuals"), 2L)
  expect_equal(counts$effective_proportion[counts$level == "groups"], 1)
  expect_gt(get_data(fit, "entropy")$relative_entropy[1L], 0.5)
  # The cluster bootstrap degenerates to the ordinary one when every unit holds
  # one row, which is the right bootstrap for independent observations. On 160
  # rows with two starts a resample occasionally fails to converge, and the
  # warning saying how many were dropped is the one this fixture expects.
  boot <- quietly(
    parameter_inference(fit, method = "bootstrap", iter = 25L,
                        n_starts = 2L, seed = 2L),
    c("multilpa_bootstrap_dropped", "multilpa_unconverged"))
  expect_true(all(boot$conf_low <= boot$estimate & boot$estimate <= boot$conf_high))
})
