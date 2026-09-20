# Holding measurement parameters at supplied values, and the staged workflow
# built on it. The constrained M-step is checked against the formula it is
# supposed to solve, not against itself.

.fixed_fixture <- function(n_groups = 40L, size = 8L, seed = 33L) {
  set.seed(seed)
  frame <- data.frame(school = rep(seq_len(n_groups), each = size))
  frame$a <- stats::rnorm(nrow(frame))
  frame$b <- stats::rnorm(nrow(frame))
  frame$q <- sample(c("no", "yes"), nrow(frame), replace = TRUE)
  frame$r <- sample(1:3, nrow(frame), replace = TRUE)
  frame
}

.fixed_structured <- function(n_groups = 60L, size = 10L, seed = 4L) {
  set.seed(seed)
  classes <- rep(1:2, length.out = n_groups)
  frame <- data.frame(school = rep(seq_len(n_groups), each = size))
  prevalence <- c(0.85, 0.2)[classes[frame$school]]
  state <- 1L + (stats::runif(nrow(frame)) > prevalence)
  frame$a <- stats::rnorm(nrow(frame), c(-2, 2)[state], 0.8)
  frame$b <- stats::rnorm(nrow(frame), c(-1.5, 1.5)[state], 0.8)
  frame
}

test_that("a held mean is the point the spread is measured around", {
  data <- .fixed_fixture(n_groups = 35L, size = 6L, seed = 21L)
  x <- as.matrix(data[, c("a", "b")])
  group_index <- match(data$school, unique(data$school))
  parameters <- list(
    means = matrix(c(-1, 1, 0.5, 2), 2L, 2L),
    variances = matrix(c(1, 2, 1.5, 0.8), 2L, 2L),
    profile_probabilities = matrix(c(0.6, 0.4, 0.3, 0.7), 2L, 2L),
    group_probabilities = c(0.5, 0.5))
  expectation <- .multilpa_expectation(x, group_index, parameters)
  held <- list(means = parameters$means)
  updated <- .multilpa_maximization(x, expectation, "varying", 1e-6, "diagonal",
                                    held = held)
  expect_identical(updated$means, held$means)
  # The M-step for a free variance under a held mean, written from the
  # definition. Centring on the free weighted mean instead would give a smaller
  # number, so this comparison discriminates between the two.
  reference <- t(vapply(seq_len(2L), function(profile) {
    responsibility <- expectation$subject_posteriors[, profile]
    colSums(sweep(x, 2L, held$means[profile, ], "-")^2 * responsibility) /
      sum(responsibility)
  }, numeric(2L)))
  expect_equal(updated$variances, reference, ignore_attr = TRUE)
  free_means <- sweep(crossprod(expectation$subject_posteriors, x), 1L,
                      colSums(expectation$subject_posteriors), "/")
  wrong <- t(vapply(seq_len(2L), function(profile) {
    responsibility <- expectation$subject_posteriors[, profile]
    colSums(sweep(x, 2L, free_means[profile, ], "-")^2 * responsibility) /
      sum(responsibility)
  }, numeric(2L)))
  expect_false(isTRUE(all.equal(updated$variances, wrong, check.attributes = FALSE)))
})

test_that("a held variance is returned as supplied and leaves the mean free", {
  data <- .fixed_fixture(n_groups = 35L, size = 6L, seed = 21L)
  x <- as.matrix(data[, c("a", "b")])
  group_index <- match(data$school, unique(data$school))
  parameters <- list(
    means = matrix(c(-1, 1, 0.5, 2), 2L, 2L),
    variances = matrix(c(1, 2, 1.5, 0.8), 2L, 2L),
    profile_probabilities = matrix(c(0.6, 0.4, 0.3, 0.7), 2L, 2L),
    group_probabilities = c(0.5, 0.5))
  expectation <- .multilpa_expectation(x, group_index, parameters)
  updated <- .multilpa_maximization(x, expectation, "varying", 1e-6, "diagonal",
                                    held = list(variances = parameters$variances))
  expect_identical(updated$variances, parameters$variances)
  expect_equal(updated$means,
               sweep(crossprod(expectation$subject_posteriors, x), 1L,
                     colSums(expectation$subject_posteriors), "/"))
})

test_that("holding nothing leaves the maximization step exactly as it was", {
  data <- .fixed_fixture(n_groups = 35L, size = 6L, seed = 21L)
  x <- as.matrix(data[, c("a", "b")])
  group_index <- match(data$school, unique(data$school))
  parameters <- list(
    means = matrix(c(-1, 1, 0.5, 2), 2L, 2L),
    variances = matrix(c(1, 2, 1.5, 0.8), 2L, 2L),
    profile_probabilities = matrix(c(0.6, 0.4, 0.3, 0.7), 2L, 2L),
    group_probabilities = c(0.5, 0.5))
  expectation <- .multilpa_expectation(x, group_index, parameters)
  expect_identical(
    .multilpa_maximization(x, expectation, "varying", 1e-6, "diagonal", held = NULL),
    .multilpa_maximization(x, expectation, "varying", 1e-6, "diagonal"))
  # And an ordinary fit is untouched by the new argument existing.
  plain <- multilpa(data, c("a", "b"), "school", n_profiles = 2L,
                    n_group_classes = 2L, n_starts = 2, seed = 1)
  explicit <- multilpa(data, c("a", "b"), "school", n_profiles = 2L,
                       n_group_classes = 2L, n_starts = 2, seed = 1,
                       fixed = character())
  expect_equal(plain$log_likelihood, explicit$log_likelihood)
  expect_equal(plain$means, explicit$means)
  expect_identical(plain$n_parameters, explicit$n_parameters)
})

test_that("every measurement model can be held, exactly and with the right count", {
  data <- .fixed_fixture()
  cases <- list(
    list(label = "gaussian diagonal", vars = c("a", "b"),
         categorical = character(), covariance_model = "diagonal",
         fixed = "measurement", blocks = c("means", "variances"), removed = 8L),
    list(label = "full covariance", vars = c("a", "b"),
         categorical = character(), covariance_model = "full",
         fixed = "measurement", blocks = c("means", "covariances"), removed = 10L),
    list(label = "mixed", vars = c("a", "b", "q", "r"),
         categorical = c("q", "r"), covariance_model = "diagonal",
         fixed = "measurement",
         blocks = c("means", "variances", "response_probabilities"), removed = 14L),
    list(label = "categorical only", vars = c("q", "r"),
         categorical = c("q", "r"), covariance_model = "diagonal",
         fixed = "response_probabilities",
         blocks = "response_probabilities", removed = 6L))
  invisible(lapply(cases, function(case) {
    fit_one <- function(classes, ...) {
      quietly(multilpa(data, case$vars, "school", n_profiles = 2L,
        n_group_classes = classes, categorical = case$categorical,
        covariance_model = case$covariance_model, n_starts = 3, seed = 1, ...))
    }
    stage_one <- fit_one(1L)
    held <- fit_one(2L, start = starting_values(stage_one, what = "measurement"),
                    fixed = case$fixed)
    free <- fit_one(2L)
    drift <- max(vapply(case$blocks, function(block) {
      max(abs(unlist(held[[block]]) - unlist(stage_one[[block]])))
    }, numeric(1)))
    expect_lt(drift, 1e-12)
    # A full-covariance count is a double, because its covariance term divides
    # by two; the comparison is on value, not on storage mode.
    expect_equal(free$n_parameters - held$n_parameters, case$removed,
                 ignore_attr = TRUE, label = case$label)
    expect_false(is.unsorted(held$log_likelihood_history))
    # A constrained maximum cannot exceed the unconstrained one.
    expect_lte(held$log_likelihood, free$log_likelihood + 1e-8)
  }))
})

test_that("incomplete indicators can be fitted with the measurement held", {
  data <- .fixed_fixture()
  data$a[c(3L, 9L, 40L)] <- NA
  data$b[c(5L, 9L)] <- NA
  stage_one <- quietly(multilpa(data, c("a", "b"), "school",
    n_profiles = 2L, n_group_classes = 1L, n_starts = 3, seed = 1, missing = "fiml"))
  held <- quietly(multilpa(data, c("a", "b"), "school", n_profiles = 2L,
    n_group_classes = 2L, n_starts = 3, seed = 1, missing = "fiml",
    start = starting_values(stage_one, what = "measurement"), fixed = "measurement"))
  expect_equal(held$means, stage_one$means)
  expect_equal(held$variances, stage_one$variances)
  expect_false(is.unsorted(held$log_likelihood_history))
})

test_that("a measurement-only start carries no mixing values", {
  data <- .fixed_fixture()
  fit <- quietly(multilpa(data, c("a", "b"), "school", n_profiles = 2L,
                                   n_group_classes = 1L, n_starts = 2, seed = 1))
  everything <- starting_values(fit)
  measurement <- starting_values(fit, what = "measurement")
  expect_true(all(c("profile_probabilities", "group_probabilities") %in%
                    names(everything)))
  expect_false(any(c("profile_probabilities", "group_probabilities") %in%
                     names(measurement)))
  expect_equal(measurement$means, everything$means)
})

test_that("a `fixed` request that cannot be met is refused by class", {
  data <- .fixed_fixture()
  fit <- quietly(multilpa(data, c("a", "b"), "school", n_profiles = 2L,
                                   n_group_classes = 1L, n_starts = 2, seed = 1))
  start <- starting_values(fit, what = "measurement")
  expect_error(multilpa(data, c("a", "b"), "school", n_profiles = 2L,
                        n_group_classes = 2L, n_starts = 1, seed = 1,
                        start = start, fixed = "profile_probabilities"),
               class = "multilpa_bad_fixed")
  expect_error(multilpa(data, c("a", "b"), "school", n_profiles = 2L,
                        n_group_classes = 2L, n_starts = 1, seed = 1,
                        fixed = "means"),
               class = "multilpa_bad_fixed")
  # The model has no categorical block to hold.
  expect_error(multilpa(data, c("a", "b"), "school", n_profiles = 2L,
                        n_group_classes = 2L, n_starts = 1, seed = 1,
                        start = start, fixed = "response_probabilities"),
               class = "multilpa_bad_fixed")
  # The start carries no variances to hold.
  bare <- start
  bare$variances <- NULL
  expect_error(multilpa(data, c("a", "b"), "school", n_profiles = 2L,
                        n_group_classes = 2L, n_starts = 1, seed = 1,
                        start = bare, fixed = "variances"),
               class = "multilpa_bad_fixed")
  expect_error(multilpa(data, c("a", "b"), "school", n_profiles = 2L,
                        n_group_classes = 2L, n_starts = 1, seed = 1,
                        start = start, fixed = TRUE),
               class = "multilpa_bad_fixed")
})

test_that("a staged fit holds its first stage and improves on it", {
  data <- .fixed_structured()
  staged <- fit_staged(data, c("a", "b"), "school", n_profiles = 2L,
                                  n_group_classes = 2L, n_starts = 5, seed = 2)
  joint <- multilpa(data, c("a", "b"), "school", n_profiles = 2L,
                    n_group_classes = 2L, n_starts = 5, seed = 2)
  expect_s3_class(staged, "multilpa")
  expect_identical(staged$means, staged$stage_one$means)
  expect_identical(staged$variances, staged$stage_one$variances)
  # Adding group classes must buy something, and cannot buy more than the
  # joint fit that was free to move the measurement as well.
  expect_gt(staged$log_likelihood, staged$stage_one$log_likelihood)
  expect_lte(staged$log_likelihood, joint$log_likelihood + 1e-8)
  expect_identical(staged$n_parameters_with_measurement, joint$n_parameters)
  expect_identical(staged$fixed, c("means", "variances"))
  expect_output(print(staged), "Held fixed")
})

test_that("a supplied first stage is used, and a mismatched one is refused", {
  data <- .fixed_structured()
  staged <- fit_staged(data, c("a", "b"), "school", n_profiles = 2L,
                                  n_group_classes = 2L, n_starts = 5, seed = 2)
  measurement <- multilpa(data, c("a", "b"), "school", n_profiles = 2L,
                          n_group_classes = 1L, n_starts = 5, seed = 2)
  again <- fit_staged(data, c("a", "b"), "school", n_profiles = 2L,
                                 n_group_classes = 2L, n_starts = 5, seed = 2,
                                 measurement = measurement)
  expect_equal(again$log_likelihood, staged$log_likelihood)
  expect_error(fit_staged(data, c("a", "b"), "school", n_profiles = 3L,
                                     n_group_classes = 2L, n_starts = 2, seed = 2,
                                     measurement = measurement),
               class = "multilpa_bad_stage")
  expect_error(fit_staged(data, "a", "school", n_profiles = 2L,
                                     n_group_classes = 2L, n_starts = 2, seed = 2,
                                     measurement = measurement),
               class = "multilpa_bad_stage")
  expect_error(fit_staged(data, c("a", "b"), "school", n_profiles = 2L,
                                     n_group_classes = 1L, n_starts = 2, seed = 2))
})

test_that("the stages table describes staged and ordinary fits alike", {
  data <- .fixed_structured()
  staged <- fit_staged(data, c("a", "b"), "school", n_profiles = 2L,
                                  n_group_classes = 2L, n_starts = 3, seed = 2)
  table_staged <- as.data.frame(staged, what = "stages")
  expect_named(table_staged, c("stage", "group_classes", "fixed",
                               "log_likelihood", "parameters",
                               "parameters_with_measurement", "converged"))
  expect_identical(table_staged$stage, c("measurement", "membership"))
  expect_identical(table_staged$group_classes, c(1L, 2L))
  expect_true(is.na(table_staged$fixed[1L]))
  expect_identical(table_staged$fixed[2L], "means, variances")

  joint <- multilpa(data, c("a", "b"), "school", n_profiles = 2L,
                    n_group_classes = 2L, n_starts = 3, seed = 2)
  table_joint <- as.data.frame(joint, what = "stages")
  expect_identical(nrow(table_joint), 1L)
  expect_identical(table_joint$stage, "joint")
  expect_true(is.na(table_joint$fixed))
  expect_identical(table_joint$parameters,
                   table_joint$parameters_with_measurement)
})

test_that("held values survive a fit that performs no update at all", {
  data <- .fixed_fixture()
  stage_one <- quietly(multilpa(data, c("a", "b"), "school",
    n_profiles = 2L, n_group_classes = 1L, n_starts = 2, seed = 1))
  start <- starting_values(stage_one, what = "measurement")
  # max_iter = 0 runs no maximization, so nothing can put the held values back
  # afterwards. They must already be in place, or this reports a measurement
  # solution nobody chose. `max_iter = 0` with a `start` now evaluates that
  # start alone, so `n_starts` is 1 here; the random-restart half of the old
  # contract moved to the test below, which needs a maximization to exercise it.
  evaluated <- quietly(multilpa(data, c("a", "b"), "school",
    n_profiles = 2L, n_group_classes = 2L, n_starts = 1, seed = 3,
    max_iter = 0L, start = start, fixed = "measurement"))
  expect_equal(evaluated$means, stage_one$means)
  expect_equal(evaluated$variances, stage_one$variances)
  expect_identical(evaluated$iterations, 0L)
  expect_identical(evaluated$best_start, 1L)
})

test_that("held values survive every random restart when EM does run", {
  data <- .fixed_fixture()
  stage_one <- quietly(multilpa(data, c("a", "b"), "school",
    n_profiles = 2L, n_group_classes = 1L, n_starts = 2, seed = 1))
  start <- starting_values(stage_one, what = "measurement")
  # With `max_iter > 0` every one of the `n_starts` restarts runs a real EM
  # sequence, so the held blocks have to be reinstated at each of them. This is
  # the half of the contract that `max_iter = 0` can no longer check, now that
  # it evaluates the supplied start alone.
  fitted <- quietly(multilpa(data, c("a", "b"), "school",
    n_profiles = 2L, n_group_classes = 2L, n_starts = 4, seed = 3,
    max_iter = 25L, start = start, fixed = "measurement"))
  expect_equal(fitted$means, stage_one$means)
  expect_equal(fitted$variances, stage_one$variances)
  expect_gt(fitted$iterations, 0L)
  # Every restart, not merely the winning one, must report the held solution.
  starts_table <- as.data.frame(fitted, what = "starts")
  expect_identical(nrow(starts_table), 4L)
})
