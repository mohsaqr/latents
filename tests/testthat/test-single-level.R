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
               class = "latents_bad_argument")
  expect_error(multilpa(data, c("x", "y"), n_profiles = 2L, n_starts = 2L),
               "pass `id = NULL`")
  # Saying it explicitly works, and says so back.
  expect_warning(multilpa(data, c("x", "y"), id = NULL, n_profiles = 2L,
                          n_starts = 2L, seed = 1L),
                 class = "latents_single_level")
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
  expect_identical(names(get_results(fit, "data")), c("x", "y"))
  expect_false(".observation" %in% names(get_results(fit, "assignments")))
  # A two-level fit still reports the identifier its caller did supply.
  data <- single_level_data()
  data$unit <- rep(seq_len(40L), each = 4L)
  grouped <- multilpa(data, c("x", "y"), id = "unit", n_profiles = 2L,
                      n_group_classes = 2L, n_starts = 3L, seed = 1L)
  expect_true("unit" %in% names(get_results(grouped, "data")))
})

test_that("group classes without an id are refused, not silently dropped", {
  expect_error(
    multilpa(single_level_data(), c("x", "y"), id = NULL, n_profiles = 2L,
             n_group_classes = 2L, n_starts = 2L),
    class = "latents_bad_argument")
  # Naming one explicitly is the same request as leaving it alone, so it passes.
  expect_s3_class(flat_fit(n_group_classes = 1L, n_starts = 2L, seed = 1L),
                  "multilpa")
})

test_that("a column already called .observation is refused", {
  data <- single_level_data()
  data$.observation <- seq_len(nrow(data))
  expect_error(
    multilpa(data, c("x", "y"), id = NULL, n_profiles = 2L, n_starts = 2L),
    class = "latents_bad_data")
})

test_that("the package's verbs work on a single-level fit", {
  skip_on_cran()
  fit <- flat_fit(n_starts = 2L, seed = 1L)
  # Inference, classification and the tables, on a fit with no second level.
  inference <- parameter_inference(fit)
  expect_true(all(is.finite(inference$standard_error)))
  expect_identical(nrow(get_results(fit, "profiles")), 4L)
  # Both levels are still reported: two profiles, and the single group class the
  # single-level fit collapses to, which carries the whole sample.
  counts <- get_results(fit, "counts")
  expect_identical(sum(counts$level == "individuals"), 2L)
  expect_equal(counts$effective_proportion[counts$level == "groups"], 1)
  expect_gt(get_results(fit, "entropy")$relative_entropy[1L], 0.5)
  # The cluster bootstrap degenerates to the ordinary one when every unit holds
  # one row, which is the right bootstrap for independent observations. On 160
  # rows with two starts a resample occasionally fails to converge, and the
  # warning saying how many were dropped is the one this fixture expects.
  boot <- quietly(
    parameter_inference(fit, method = "bootstrap", iter = 10L,
                        n_starts = 1L, seed = 2L),
    c("latents_bootstrap_dropped", "latents_unconverged"))
  expect_true(all(boot$conf_low <= boot$estimate & boot$estimate <= boot$conf_high))
})

test_that("a truth column cross-tabulates against profiles, not a single class", {
  # With one observation per group, "constant within every group" is vacuously
  # true of every column, so the level detection sent every truth column to the
  # group level and tabulated it against a group class that is the same for all
  # rows -- proportion 1.0 everywhere and nothing learned.
  data <- single_level_data()
  data$known <- rep(c("a", "b"), each = nrow(data) / 2L)
  fit <- quietly(multilpa(data, c("x", "y"), id = NULL, n_profiles = 2L,
                          n_starts = 3L, seed = 1L))
  recovery <- get_results(fit, "assignments", data = data, truth = "known")
  expect_identical(unique(recovery$assignment), "profile")
  expect_identical(nrow(recovery), 4L)
  # Within each truth value the proportions are a distribution over profiles.
  shares <- aggregate(proportion ~ value, data = recovery, FUN = sum)
  expect_equal(shares$proportion, rep(1, 2))
  # A two-level fit is unaffected: a column constant within groups still
  # describes the group.
  grouped_data <- single_level_data()
  grouped_data$unit <- rep(seq_len(40L), each = 4L)
  grouped_data$kind <- rep(c("a", "b"), each = nrow(grouped_data) / 2L)
  grouped <- multilpa(grouped_data, c("x", "y"), id = "unit", n_profiles = 2L,
                      n_group_classes = 2L, n_starts = 3L, seed = 1L)
  expect_identical(
    unique(get_results(grouped, "assignments", data = grouped_data,
                    truth = "kind")$assignment),
    "group_class")
})
