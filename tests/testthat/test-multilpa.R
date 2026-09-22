# Generate and inspect a deterministic, clearly separated two-level population.
multilpa_fixture <- local({
  set.seed(6123)
  n_groups <- 80L
  group_size <- 16L
  group_types <- rep(1:2, each = n_groups / 2)
  group_index <- rep(seq_len(n_groups), each = group_size)
  profiles <- 1L + as.integer(
    runif(length(group_index)) > c(0.9, 0.1)[group_types[group_index]]
  )
  data <- data.frame(
    group = sprintf("school_%03d", group_index),
    reading = rnorm(length(profiles), c(-2.3, 2.3)[profiles],
                    c(0.6, 0.8)[profiles]),
    maths = rnorm(length(profiles), c(-1.8, 1.8)[profiles],
                  c(0.7, 0.9)[profiles])
  )
  stopifnot(!anyNA(data), !anyDuplicated(data),
            all(is.finite(as.matrix(data[c("reading", "maths")]))))
  list(data = data, profiles = profiles, group_types = group_types)
})

multilpa_small_data <- data.frame(
  group = c("z", "z", "a", "a", "b", "b"),
  x = c(-2, -1, 0, 1, 2, 3),
  y = c(0, 2, 1, 4, 3, 5)
)

test_that("a multilevel Gaussian population is recovered with valid posteriors", {
  data <- multilpa_fixture$data
  original_data <- data
  fit <- multilpa(data, c("reading", "maths"), "group", 2L, 2L,
                   n_starts = 4L, max_iter = 600L, seed = 827)

  expect_s3_class(fit, "multilpa")
  expect_identical(data, original_data)
  expect_true(fit$converged)
  expect_type(fit$subject_profiles, "integer")
  expect_type(fit$group_classes, "integer")
  expect_equal(dim(fit$subject_posteriors), c(nrow(data), 2L))
  expect_equal(dim(fit$group_posteriors), c(80L, 2L))
  expect_equal(rowSums(fit$subject_posteriors), rep(1, nrow(data)),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(rowSums(fit$group_posteriors), rep(1, 80L),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_true(all(is.finite(fit$subject_posteriors)))
  expect_true(all(is.finite(fit$group_posteriors)))
  expect_true(all(fit$subject_posteriors >= 0 & fit$subject_posteriors <= 1 + 1e-12))
  expect_true(all(fit$group_posteriors >= 0 & fit$group_posteriors <= 1 + 1e-12))
  expect_equal(rowSums(fit$profile_probabilities), c(1, 1),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(sum(fit$group_probabilities), 1, tolerance = 1e-10)
  expect_true(all(fit$variances > 0))
  expect_false(fit$boundary)
  expect_identical(as.character(fit$group_ids), unique(data$group))
  expect_equal(fit$subject_profiles, max.col(fit$subject_posteriors),
               ignore_attr = TRUE)
  expect_equal(fit$group_classes, max.col(fit$group_posteriors),
               ignore_attr = TRUE)

  # Resolve arbitrary mixture labels before checking statistical recovery.
  profile_order <- order(fit$means[, "reading"])
  fitted_means <- fit$means[profile_order, , drop = FALSE]
  fitted_variances <- fit$variances[profile_order, , drop = FALSE]
  expect_equal(unname(fitted_means), matrix(c(-2.3, -1.8, 2.3, 1.8),
                                          nrow = 2L, byrow = TRUE),
               tolerance = 0.15)
  expect_equal(unname(fitted_variances), matrix(c(0.36, 0.49, 0.64, 0.81),
                                              nrow = 2L, byrow = TRUE),
               tolerance = 0.15)
  sorted_probabilities <- fit$profile_probabilities[, profile_order, drop = FALSE]
  group_order <- order(sorted_probabilities[, 1L], decreasing = TRUE)
  expect_equal(unname(sorted_probabilities[group_order, , drop = FALSE]),
               matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2L, byrow = TRUE),
               tolerance = 0.1)
  expect_gt(mean(match(fit$subject_profiles, profile_order) ==
                   multilpa_fixture$profiles), 0.97)
  expect_gt(mean(match(fit$group_classes, group_order) ==
                   multilpa_fixture$group_types), 0.9)

  expect_length(fit$log_likelihood_history, fit$iterations + 1L)
  expect_true(all(diff(fit$log_likelihood_history) >= -1e-7))
  expect_equal(tail(fit$log_likelihood_history, 1L), fit$log_likelihood,
               ignore_attr = TRUE)
  expect_equal(fit$n_parameters, 11L)
  expect_equal(fit$aic, -2 * fit$log_likelihood + 2 * fit$n_parameters)
  expect_equal(fit$bic, -2 * fit$log_likelihood + log(80) * fit$n_parameters)
  expect_equal(fit$bic_individual,
               -2 * fit$log_likelihood + log(nrow(data)) * fit$n_parameters)
  expect_s3_class(fit$starts, "data.frame")
  expect_equal(nrow(fit$starts), 4L)
  expect_true(all(c("start", "log_likelihood", "converged", "iterations",
                    "error", "boundary") %in% names(fit$starts)))
  expect_equal(fit$log_likelihood, max(fit$starts$log_likelihood, na.rm = TRUE))
})

test_that("one profile gives the analytic maximum likelihood solution", {
  fit <- multilpa(multilpa_small_data, c("x", "y"), "group", 1L, 1L,
                   n_starts = 1L, seed = 81)
  expected_means <- vapply(multilpa_small_data[c("x", "y")], mean, numeric(1))
  expected_variances <- vapply(multilpa_small_data[c("x", "y")], function(x) {
    mean((x - mean(x))^2)
  }, numeric(1))
  expected_log_likelihood <- sum(dnorm(multilpa_small_data$x, expected_means["x"],
                                       sqrt(expected_variances["x"]), log = TRUE)) +
    sum(dnorm(multilpa_small_data$y, expected_means["y"],
              sqrt(expected_variances["y"]), log = TRUE))

  expect_equal(as.numeric(fit$means), unname(expected_means), tolerance = 1e-10)
  expect_equal(as.numeric(fit$variances), unname(expected_variances),
               tolerance = 1e-10)
  expect_equal(fit$log_likelihood, expected_log_likelihood, tolerance = 1e-10)
  expect_equal(unname(fit$subject_posteriors), matrix(1, nrow = 6L, ncol = 1L))
  expect_equal(unname(fit$group_posteriors), matrix(1, nrow = 3L, ncol = 1L))
  expect_equal(fit$n_parameters, 4L)
  expect_equal(as.numeric(logLik(fit)), expected_log_likelihood)
  expect_equal(attr(logLik(fit), "df"), 4L)
  expect_equal(attr(logLik(fit), "nobs"), 3L)
  expect_equal(nobs(fit), 3L)
  expect_equal(AIC(fit), fit$aic)
  expect_equal(BIC(fit), fit$bic)
  expect_equal(AIC(fit, k = log(3)), fit$bic)
  expect_output(printed_fit <- withVisible(print(fit)), "6 individuals in 3 groups")
  expect_false(printed_fit$visible)
  expect_identical(printed_fit$value, fit)
  summary_fit <- summary(fit)
  expect_s3_class(summary_fit, "summary_multilpa")
  # A field the fit does not carry must take its documented default, not leave
  # an unnamed hole in the summary.
  expect_false(anyNA(names(summary_fit)))
  expect_equal(as.data.frame(summary_fit)$mean, as.vector(t(fit$means)))
  counts <- get_results(summary_fit, "counts")
  expect_equal(subset(counts, level == "individuals")$effective_count,
               colSums(fit$subject_posteriors), ignore_attr = TRUE)
  expect_equal(subset(counts, level == "groups")$effective_count,
               colSums(fit$group_posteriors), ignore_attr = TRUE)
  expect_output(printed_summary <- withVisible(print(summary_fit)),
                "-- counts", fixed = TRUE)
  expect_false(printed_summary$visible)
  expect_identical(printed_summary$value, summary_fit)
  expect_error(print(summary_fit, digits = 0))
})

test_that("equal variances are shared across profiles and counted correctly", {
  fit <- multilpa(multilpa_fixture$data, c("reading", "maths"), "group",
                   2L, 2L, variance_model = "equal", n_starts = 2L, seed = 712)
  expect_equal(fit$variances[1L, ], fit$variances[2L, ], ignore_attr = TRUE)
  expect_equal(fit$n_parameters, 9L)
  expect_true(all(diff(fit$log_likelihood_history) >= -1e-7))
})

test_that("explicit starts give the independently calculated initial likelihood", {
  start <- list(
    means = matrix(c(-1, 1, 2, 4), nrow = 2L, byrow = TRUE),
    variances = matrix(c(1.5, 2, 2.5, 1), nrow = 2L, byrow = TRUE),
    profile_probabilities = matrix(c(0.8, 0.2, 0.3, 0.7),
                                   nrow = 2L, byrow = TRUE),
    group_probabilities = c(0.4, 0.6)
  )
  # Small, moderate densities allow direct probability arithmetic independent
  # of the fitter's log-sum-exp implementation.
  density <- cbind(
    dnorm(multilpa_small_data$x, -1, sqrt(1.5)) *
      dnorm(multilpa_small_data$y, 1, sqrt(2)),
    dnorm(multilpa_small_data$x, 2, sqrt(2.5)) *
      dnorm(multilpa_small_data$y, 4, 1)
  )
  conditional_density <- density %*% t(start$profile_probabilities)
  likelihood_by_group_type <- t(vapply(unique(multilpa_small_data$group),
    function(group) {
      apply(conditional_density[multilpa_small_data$group == group, , drop = FALSE],
            2L, prod)
    }, numeric(2)))
  expected_log_likelihood <- sum(log(as.numeric(
    likelihood_by_group_type %*% start$group_probabilities
  )))

  expect_warning(
    fit <- multilpa(multilpa_small_data, c("x", "y"), "group", 2L, 2L,
                     n_starts = 1L, max_iter = 1L, start = start, seed = 92),
    "converg|iteration"
  )
  expect_equal(fit$log_likelihood_history[1L], expected_log_likelihood,
               tolerance = 1e-10)
  expect_false(fit$converged)
  expect_equal(fit$iterations, 1L)
  expect_length(fit$log_likelihood_history, 2L)
  expect_gte(fit$log_likelihood, expected_log_likelihood - 1e-10)
})

test_that("group identifiers preserve first occurrence without unused levels", {
  data <- multilpa_small_data[c(5L, 2L, 3L, 6L, 1L, 4L), , drop = FALSE]
  data$group <- factor(data$group, levels = c("unused", "a", "b", "z"))
  fit <- multilpa(data, c("x", "y"), "group", 1L, 1L,
                   n_starts = 1L, seed = 42)
  expect_identical(as.character(fit$group_ids), c("b", "z", "a"))
  expect_equal(nrow(fit$group_posteriors), 3L)
  expect_equal(length(fit$subject_profiles), nrow(data))
})

test_that("nearby numeric group identifiers are not merged during display conversion", {
  data <- multilpa_small_data
  data$group <- rep(c(1 + 1e-15, 1), each = 3L)
  stopifnot(length(unique(data$group)) == 2L)
  fit <- multilpa(data, c("x", "y"), "group", 1L, 1L,
                   n_starts = 1L, seed = 82)
  expect_equal(fit$n_groups, 2L)
  expect_identical(fit$group_index, rep(1:2, each = 3L))
  expect_equal(unname(fit$group_sizes), c(3L, 3L))
  expect_length(fit$group_ids, 2L)
  expect_false(anyDuplicated(fit$group_ids) > 0L)
  expect_equal(dim(fit$group_posteriors), c(2L, 1L))
})

test_that("a single group class permits singleton groups and one indicator", {
  data <- data.frame(group = seq_len(12L), x = c(seq(-3, -2, length.out = 6L),
                                                seq(2, 3, length.out = 6L)))
  fit <- multilpa(data, "x", "group", 2L, 1L,
                   n_starts = 2L, seed = 62)
  expect_equal(dim(fit$means), c(2L, 1L))
  expect_equal(dim(fit$variances), c(2L, 1L))
  expect_equal(dim(fit$group_posteriors), c(12L, 1L))
  expect_equal(fit$n_parameters, 5L)
  expect_equal(sort(as.numeric(fit$means)), c(-2.5, 2.5), tolerance = 1e-5)
})

test_that("a two-observation one-profile fit does not lose matrix dimensions", {
  data <- data.frame(group = c("a", "a"), score = c(2, 4))
  fit <- multilpa(data, "score", "group", 1L, 1L,
                   n_starts = 1L, seed = 6)
  expect_equal(dim(fit$means), c(1L, 1L))
  expect_equal(dim(fit$variances), c(1L, 1L))
  expect_equal(as.numeric(fit$means), 3)
  expect_equal(as.numeric(fit$variances), 1)
  expect_equal(fit$log_likelihood, sum(dnorm(c(2, 4), 3, 1, log = TRUE)))
})

test_that("large indicator offsets preserve variances and likelihood", {
  original_fit <- multilpa(multilpa_small_data, c("x", "y"), "group", 1L, 1L,
                            n_starts = 1L, seed = 62)
  shifted_data <- transform(multilpa_small_data, x = x + 1e10, y = y - 1e10)
  shifted_fit <- multilpa(shifted_data, c("x", "y"), "group", 1L, 1L,
                           n_starts = 1L, seed = 62)
  expect_equal(unname(shifted_fit$means - original_fit$means),
               matrix(c(1e10, -1e10), nrow = 1L))
  expect_equal(shifted_fit$variances, original_fit$variances, tolerance = 1e-10)
  expect_equal(shifted_fit$log_likelihood, original_fit$log_likelihood,
               tolerance = 1e-10)
})

test_that("seeded fitting is repeatable and preserves an existing RNG state", {
  set.seed(198)
  state_before <- .Random.seed
  fit_one <- multilpa(multilpa_fixture$data, c("reading", "maths"), "group", 2L, 2L,
                       n_starts = 2L, seed = 612)
  expect_identical(.Random.seed, state_before)
  fit_two <- multilpa(multilpa_fixture$data, c("reading", "maths"), "group", 2L, 2L,
                       n_starts = 2L, seed = 612)
  expect_identical(fit_one$means, fit_two$means)
  expect_identical(fit_one$variances, fit_two$variances)
  expect_identical(fit_one$log_likelihood_history, fit_two$log_likelihood_history)
  expect_identical(fit_one$starts, fit_two$starts)
  expect_identical(.Random.seed, state_before)
})

test_that("permuting input rows permutes posteriors without changing the fit", {
  data <- multilpa_fixture$data
  start <- list(
    means = matrix(c(-2.3, -1.8, 2.3, 1.8), nrow = 2L, byrow = TRUE),
    variances = matrix(c(0.36, 0.49, 0.64, 0.81), nrow = 2L, byrow = TRUE),
    profile_probabilities = matrix(c(0.9, 0.1, 0.1, 0.9),
                                   nrow = 2L, byrow = TRUE),
    group_probabilities = c(0.5, 0.5)
  )
  fit <- multilpa(data, c("reading", "maths"), "group", 2L, 2L,
                   n_starts = 1L, start = start, seed = 81)
  row_order <- rev(seq_len(nrow(data)))
  reordered_fit <- multilpa(data[row_order, , drop = FALSE],
                             c("reading", "maths"), "group", 2L, 2L,
                             n_starts = 1L, start = start, seed = 81)
  expect_equal(reordered_fit$log_likelihood, fit$log_likelihood, tolerance = 1e-9)
  expect_equal(reordered_fit$means, fit$means, tolerance = 1e-8)
  expect_equal(reordered_fit$subject_posteriors,
               fit$subject_posteriors[row_order, , drop = FALSE],
               tolerance = 1e-8, ignore_attr = TRUE)
  group_order <- match(reordered_fit$group_ids, fit$group_ids)
  expect_equal(reordered_fit$group_posteriors,
               fit$group_posteriors[group_order, , drop = FALSE],
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_identical(as.character(reordered_fit$group_ids),
                   unique(data$group[row_order]))
})

test_that("seeded fitting restores the absence of an RNG state", {
  previous_state <- if (exists(".Random.seed", envir = .GlobalEnv,
                               inherits = FALSE)) .GlobalEnv$.Random.seed else NULL
  on.exit({
    if (!is.null(previous_state)) {
      assign(".Random.seed", previous_state, envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  fit <- multilpa(multilpa_small_data, c("x", "y"), "group", 1L, 1L,
                   n_starts = 1L, seed = 27)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  expect_true(is.finite(fit$log_likelihood))
})

test_that("seed restoration also runs when fitting fails", {
  set.seed(481)
  state_before <- .Random.seed
  data <- data.frame(group = c("a", "a"), score = c(-1e200, 1e200))
  expect_error(multilpa(data, "score", "group", 1L, 1L,
                         n_starts = 1L, seed = 62), "overflow|rescale")
  expect_identical(.Random.seed, state_before)
})

test_that("invalid data and unidentifiable specifications fail explicitly", {
  args <- list(data = multilpa_small_data, vars = c("x", "y"),
               group = "group", n_profiles = 2L, n_group_classes = 2L,
               n_starts = 1L)
  invalid_data <- list(
    transform(multilpa_small_data, x = c(NA, -1, 0, 1, 2, 3)),
    transform(multilpa_small_data, x = c(Inf, -1, 0, 1, 2, 3)),
    transform(multilpa_small_data, x = c(NaN, -1, 0, 1, 2, 3)),
    transform(multilpa_small_data, x = as.character(x)),
    transform(multilpa_small_data, x = 1),
    transform(multilpa_small_data, group = c(NA, "z", "a", "a", "b", "b")),
    multilpa_small_data[FALSE, , drop = FALSE]
  )
  invisible(lapply(invalid_data, function(data) {
    invalid_args <- args
    invalid_args$data <- data
    expect_error(do.call(multilpa, invalid_args))
  }))
  invalid_specs <- list(
    list(vars = "missing"), list(vars = character()),
    list(vars = c("x", "x")), list(group = "missing"),
    list(group = c("group", "x")), list(vars = c("x", "group")),
    list(n_profiles = 7L), list(n_group_classes = 4L),
    list(n_profiles = 1L, n_group_classes = 2L),
    list(data = transform(multilpa_small_data, group = seq_len(nrow(multilpa_small_data))))
  )
  invisible(lapply(invalid_specs, function(spec) {
    expect_error(do.call(multilpa, utils::modifyList(args, spec)))
  }))
})

test_that("invalid numeric controls and variance models are rejected", {
  args <- list(data = multilpa_small_data, vars = c("x", "y"),
               group = "group", n_profiles = 2L, n_group_classes = 2L,
               n_starts = 1L)
  invalid_specs <- list(
    list(n_profiles = 0L), list(n_profiles = 1.5), list(n_profiles = NA),
    list(n_group_classes = 0L), list(n_group_classes = 1.5),
    list(n_starts = 0L), list(n_starts = 1.5), list(n_starts = Inf),
    list(max_iter = -1L), list(max_iter = 1.5), list(max_iter = NA),
    list(tol = 0), list(tol = -1), list(tol = Inf),
    list(min_variance = 0), list(min_variance = -1), list(min_variance = NA),
    list(seed = NA), list(seed = c(1, 2)), list(variance_model = "full")
  )
  invisible(lapply(invalid_specs, function(spec) {
    expect_error(do.call(multilpa, utils::modifyList(args, spec)))
  }))
})

test_that("malformed explicit starting values are rejected", {
  start <- list(
    means = matrix(c(-1, 1, 2, 4), nrow = 2L, byrow = TRUE),
    variances = matrix(1, nrow = 2L, ncol = 2L),
    profile_probabilities = matrix(c(0.8, 0.2, 0.3, 0.7),
                                   nrow = 2L, byrow = TRUE),
    group_probabilities = c(0.4, 0.6)
  )
  invalid_starts <- list(
    start[c("means", "variances")],
    utils::modifyList(start, list(means = matrix(1, nrow = 1L, ncol = 2L))),
    utils::modifyList(start, list(means = matrix(NA_real_, nrow = 2L, ncol = 2L))),
    utils::modifyList(start, list(variances = matrix(-1, nrow = 2L, ncol = 2L))),
    utils::modifyList(start, list(variances = matrix(1e-12, nrow = 2L, ncol = 2L))),
    utils::modifyList(start, list(profile_probabilities = matrix(1, 2L, 2L))),
    utils::modifyList(start, list(profile_probabilities = matrix(c(-0.1, 1.1, 0.3, 0.7),
                                                               2L, byrow = TRUE))),
    utils::modifyList(start, list(profile_probabilities = matrix(c(0, 1, 0.3, 0.7),
                                                               2L, byrow = TRUE))),
    utils::modifyList(start, list(group_probabilities = c(0.4, 0.5))),
    utils::modifyList(start, list(group_probabilities = c(-0.1, 1.1))),
    utils::modifyList(start, list(group_probabilities = c(0, 1))),
    utils::modifyList(start, list(group_probabilities = c(NA, 1)))
  )
  invisible(lapply(invalid_starts, function(invalid_start) {
    expect_error(multilpa(multilpa_small_data, c("x", "y"), "group", 2L, 2L,
                           n_starts = 1L, start = invalid_start))
  }))
  start$variances[1L, 1L] <- 2
  expect_error(multilpa(multilpa_small_data, c("x", "y"), "group", 2L, 2L,
                         variance_model = "equal", n_starts = 1L, start = start))
})

test_that("a binding variance floor is reported and uses constrained ML", {
  data <- data.frame(group = c("a", "a", "b", "b"),
                     x = c(-0.001, 0.001, -0.002, 0.002))
  expect_warning(
    fit <- multilpa(data, "x", "group", 1L, 1L, n_starts = 1L,
                     min_variance = 0.01, seed = 17),
    "variance|boundary|floor"
  )
  expect_true(fit$boundary)
  expect_equal(as.numeric(fit$variances), 0.01)
  expect_equal(fit$log_likelihood, sum(dnorm(data$x, 0, 0.1, log = TRUE)))
  expect_true(all(diff(fit$log_likelihood_history) >= -1e-7))
})
