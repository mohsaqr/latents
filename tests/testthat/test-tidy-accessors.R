# The tidy surface of a fit, its summary, and a set of starting values:
# every verb returns a base data.frame whose columns are atomic, the summary
# never carries an unnamed field, and starting values still round-trip.

tidy_fixture <- local({
  set.seed(3141)
  group <- rep(seq_len(20L), each = 12L)
  # Two group classes that differ sharply in profile prevalence, so that a
  # two-class fit is identified and no class is effectively empty.
  group_type <- rep(1:2, each = 10L)
  profile <- 1L + as.integer(stats::runif(length(group)) >
                               c(0.9, 0.1)[group_type[group]])
  data <- data.frame(
    school = group,
    score_a = stats::rnorm(length(group), c(-2, 2)[profile], 0.7),
    score_b = stats::rnorm(length(group), c(-1.5, 1.5)[profile], 0.8)
  )
  stopifnot("the fixture must be complete" = !anyNA(data),
            "the fixture must be finite" =
              all(is.finite(as.matrix(data[c("score_a", "score_b")]))))
  data
})

tidy_categorical_fixture <- local({
  set.seed(2718)
  n <- 240L
  profile <- rep(seq_len(2L), length.out = n)
  draw <- function(high) {
    ifelse(stats::runif(n) < ifelse(profile == 1L, high, 1 - high), "yes", "no")
  }
  data.frame(school = rep(seq_len(12L), each = n / 12L),
             a = draw(0.85), b = draw(0.8), c = draw(0.75))
})

tidy_fit <- function(...) {
  multilpa(tidy_fixture, c("score_a", "score_b"), "school", n_profiles = 2L,
           n_group_classes = 2L, n_starts = 2L, seed = 11L, ...)
}

tidy_categorical_fit <- function(...) {
  multilpa(tidy_categorical_fixture, c("a", "b", "c"), "school",
           n_profiles = 2L, n_group_classes = 2L,
           categorical = c("a", "b", "c"), n_starts = 2L, seed = 11L, ...)
}

# A tidy table is a base data.frame with atomic columns: a matrix or a list
# column smuggled into a column is the defect this checks for.
expect_tidy_frame <- function(frame) {
  expect_s3_class(frame, "data.frame")
  expect_true(all(vapply(frame, is.atomic, logical(1))))
  expect_false(any(vapply(frame, is.matrix, logical(1))))
  expect_identical(anyDuplicated(names(frame)), 0L)
  expect_false(any(names(frame) == "" | is.na(names(frame))))
}

test_that("every table of a fit is tidy", {
  fit <- tidy_fit()
  tables <- c("profiles", "responses", "profile_probabilities", "posteriors",
              "group_posteriors", "starts", "stages", "information_criteria",
              "classification", "entropy")
  invisible(lapply(tables, function(what) {
    expect_tidy_frame(get_results(fit, what))
  }))
  expect_identical(nrow(get_results(fit, "profiles")), 4L)
  expect_identical(nrow(get_results(fit, "responses")), 0L)
  expect_identical(nrow(get_results(fit, "profile_probabilities")), 4L)
  expect_identical(nrow(get_results(fit, "posteriors")),
                   nrow(tidy_fixture) * fit$n_profiles)
  expect_identical(nrow(get_results(fit, "posteriors", format = "wide")),
                   nrow(tidy_fixture))
  expect_identical(nrow(get_results(fit, "group_posteriors")),
                   20L * fit$n_group_classes)
  expect_identical(nrow(get_results(fit, "starts")), 2L)
  expect_identical(nrow(get_results(fit, "stages")), 1L)
})

test_that("the standardized profile means the plot draws are available as data", {
  fit <- tidy_fit()
  raw <- get_results(fit, "profiles")
  standardized <- get_results(fit, "profiles", scale = "standardized")
  expect_tidy_frame(standardized)
  expect_identical(names(standardized), names(raw))
  expect_identical(standardized$profile, raw$profile)
  expect_identical(standardized$indicator, raw$indicator)
  # Computed independently from the input data, not from the fit, so that the
  # table is checked against the definition rather than against itself.
  centre <- vapply(tidy_fixture[c("score_a", "score_b")], mean, numeric(1))
  spread <- vapply(tidy_fixture[c("score_a", "score_b")], stats::sd, numeric(1))
  expect_equal(standardized$mean,
               (raw$mean - rep(centre, times = 2L)) / rep(spread, times = 2L),
               ignore_attr = TRUE)
  # The whole row moves to the new scale together: a variance divides by the
  # squared spread, a standard deviation by the spread.
  expect_equal(standardized$variance, raw$variance / rep(spread, times = 2L)^2,
               ignore_attr = TRUE)
  expect_equal(standardized$standard_deviation,
               sqrt(standardized$variance))
  # Standardizing is a change of units, so it cannot change the shape: the
  # profile separation on one indicator rescales by exactly that indicator.
  expect_equal(diff(subset(standardized, indicator == "score_a")$mean) /
                 diff(subset(raw, indicator == "score_a")$mean),
               1 / unname(spread["score_a"]))
})

test_that("standardizing errors where the plot would, and only where it applies", {
  fit <- tidy_fit()
  expect_error(get_results(fit, "posteriors", scale = "standardized"),
               class = "latents_bad_argument")
  expect_error(get_results(fit, "profiles", scale = "nonsense"))
  stripped <- fit
  stripped$indicator_data <- NULL
  expect_error(get_results(stripped, "profiles", scale = "standardized"),
               class = "latents_no_indicator_data")
  flat <- fit
  flat$indicator_data[, 1L] <- 1
  expect_error(get_results(flat, "profiles", scale = "standardized"),
               class = "latents_bad_scale")
})

test_that("an unestimated stage is missing, not a sentinel level", {
  fit <- tidy_fit()
  stages <- get_results(fit, "stages")
  expect_identical(stages$stage, "joint")
  expect_true(is.na(stages$fixed))
  expect_type(stages$fixed, "character")
  expect_false(any(stages$fixed %in% "none", na.rm = TRUE))
})

test_that("a summary carries no unnamed field, whatever the fit omits", {
  fits <- list(
    diagonal = tidy_fit(),
    equal = tidy_fit(variance_model = "equal"),
    full = tidy_fit(covariance_model = "full"),
    categorical = tidy_categorical_fit()
  )
  invisible(lapply(fits, function(fit) {
    summarized <- summary(fit)
    expect_s3_class(summarized, "summary_multilpa")
    expect_false(anyNA(names(summarized)))
    expect_false(any(names(summarized) == ""))
    expect_identical(anyDuplicated(names(summarized)), 0L)
  }))
  # The same names whatever the parameterization, so a caller need not branch.
  expect_identical(names(summary(fits$diagonal)), names(summary(fits$full)))
  expect_identical(names(summary(fits$diagonal)),
                   names(summary(fits$categorical)))
})

test_that("a summary defaults its absent fields rather than dropping them", {
  summarized <- summary(tidy_fit())
  expect_identical(summarized$fixed, character())
  expect_false(summarized$staged)
  expect_identical(summarized$n_parameters_with_measurement,
                   summarized$n_parameters)
  expect_null(summarized$covariances)
  expect_null(summarized$response_probabilities)
})

test_that("every table of a summary is tidy", {
  summarized <- summary(tidy_fit())
  tables <- c("profiles", "responses", "profile_probabilities", "counts",
              "covariances", "model", "starts")
  invisible(lapply(tables, function(what) {
    expect_tidy_frame(get_results(summarized, what))
  }))
  expect_identical(names(as.data.frame(summarized)),
                   c("profile", "indicator", "mean", "variance",
                     "standard_deviation"))
  expect_identical(names(get_results(summarized, "counts")),
                   c("level", "class", "effective_count",
                     "effective_proportion"))
  expect_identical(nrow(get_results(summarized, "model")), 1L)
})

test_that("the summary tables report the quantities the summary was built from", {
  fit <- tidy_fit()
  summarized <- summary(fit)
  expect_equal(as.data.frame(summarized)$mean, as.vector(t(fit$means)))
  expect_equal(as.data.frame(summarized)$standard_deviation,
               sqrt(as.vector(t(fit$variances))))
  counts <- get_results(summarized, "counts")
  individuals <- subset(counts, level == "individuals")
  groups <- subset(counts, level == "groups")
  expect_equal(individuals$effective_count, colSums(fit$subject_posteriors),
               ignore_attr = TRUE)
  expect_equal(groups$effective_count, colSums(fit$group_posteriors),
               ignore_attr = TRUE)
  # The effective memberships partition the sample at each level.
  expect_equal(sum(individuals$effective_count), fit$n_observations)
  expect_equal(sum(groups$effective_count), fit$n_groups)
  expect_equal(sum(individuals$effective_proportion), 1)
  expect_equal(sum(groups$effective_proportion), 1)
  fit_row <- get_results(summarized, "model")
  expect_identical(fit_row$n_observations, fit$n_observations)
  expect_identical(fit_row$n_informative, fit$n_informative)
  expect_equal(fit_row$log_likelihood, fit$log_likelihood)
  expect_equal(fit_row$bic_individual, fit$bic_individual)
  expect_identical(get_results(summarized, "starts"), fit$starts)
})

test_that("the summary reports a covariance matrix under either parameterization", {
  diagonal <- get_results(summary(tidy_fit()), "covariances")
  expect_identical(nrow(diagonal), 8L)
  expect_identical(names(diagonal),
                   c("profile", "indicator", "indicator_2", "covariance"))
  off_diagonal <- subset(diagonal, indicator != indicator_2)
  # The diagonal parameterization states zero residual covariance; that is a
  # value, not a missing one.
  expect_equal(off_diagonal$covariance, rep(0, 4L))
  expect_false(anyNA(diagonal$covariance))
  on_diagonal <- subset(diagonal, indicator == indicator_2)
  expect_equal(on_diagonal$covariance,
               as.data.frame(summary(tidy_fit()))$variance)

  full_fit <- tidy_fit(covariance_model = "full")
  full <- get_results(summary(full_fit), "covariances")
  expect_identical(nrow(full), 8L)
  expect_equal(full$covariance, as.vector(full_fit$covariances))
  # A covariance matrix is symmetric whichever way it is read.
  swapped <- merge(full, full,
                   by.x = c("profile", "indicator", "indicator_2"),
                   by.y = c("profile", "indicator_2", "indicator"))
  expect_equal(swapped$covariance.x, swapped$covariance.y)
})

test_that("a categorical summary reports responses and no empty Gaussian block", {
  summarized <- summary(tidy_categorical_fit())
  responses <- get_results(summarized, "responses")
  expect_tidy_frame(responses)
  expect_identical(names(responses),
                   c("profile", "indicator", "category", "probability",
                     "threshold"))
  expect_identical(nrow(responses), 12L)
  expect_identical(nrow(get_results(summarized, "profiles")), 0L)
  expect_identical(nrow(get_results(summarized, "covariances")), 0L)
  expect_output(print(summarized), "-- responses", fixed = TRUE)
  # The Gaussian block of an all-categorical fit is empty, and the summary
  # says so under its own heading rather than promising means it has not.
  printed <- capture.output(print(summarized))
  gaussian <- grep("-- profiles ", printed, fixed = TRUE)
  expect_length(gaussian, 1L)
  expect_match(printed[gaussian + 1L], "no rows")
})

test_that("the summary print method is a report and returns its input", {
  summarized <- summary(tidy_fit())
  expect_output(printed <- withVisible(print(summarized)),
                "-- counts", fixed = TRUE)
  expect_false(printed$visible)
  expect_identical(printed$value, summarized)
  # No print method may teach the banned idiom.
  printed_text <- capture.output(print(summarized))
  expect_false(any(grepl("$", printed_text, fixed = TRUE)))
})

test_that("the tidy accessors refuse an object of the wrong class", {
  fit <- tidy_fit()
  expect_error(as.data.frame.summary_multilpa(fit),
               "must be an object of class `summary_multilpa`")
  expect_error(get_results(summary(fit), "nonsense"))
  expect_error(print.summary_multilpa(summary(fit), digits = 0),
               "between 1 and 22")
  expect_error(print.multilpa_start(fit), "must be a `multilpa_start` object")
  expect_error(as.data.frame.multilpa_start(fit),
               "must be an object of class `multilpa_start`")
})

test_that("starting values are a classed object that still round-trips", {
  fit <- tidy_fit()
  start <- starting_values(fit)
  expect_s3_class(start, "multilpa_start")
  expect_type(start, "list")
  expect_setequal(names(start), c("means", "variances",
                                  "profile_probabilities",
                                  "group_probabilities"))
  # The class must not change what the fitter reads: scoring the start without
  # a single EM update must reproduce the likelihood it came from.
  evaluated <- multilpa(tidy_fixture, c("score_a", "score_b"), "school",
                        n_profiles = 2L, n_group_classes = 2L, n_starts = 1L,
                        start = start, max_iter = 0L)
  expect_identical(evaluated$iterations, 0L)
  expect_equal(as.numeric(logLik(evaluated)), as.numeric(logLik(fit)),
               tolerance = 1e-12)
  refit <- multilpa(tidy_fixture, c("score_a", "score_b"), "school",
                    n_profiles = 2L, n_group_classes = 2L, n_starts = 1L,
                    start = start, max_iter = 5000L, tol = 1e-13)
  expect_equal(refit$log_likelihood, fit$log_likelihood, tolerance = 1e-9)
})

test_that("a measurement-only start round-trips into a fit that holds it", {
  fit <- tidy_fit()
  measurement <- starting_values(fit, what = "measurement")
  expect_s3_class(measurement, "multilpa_start")
  expect_setequal(names(measurement), c("means", "variances"))
  held <- multilpa(tidy_fixture, c("score_a", "score_b"), "school",
                   n_profiles = 2L, n_group_classes = 3L, n_starts = 1L,
                   seed = 11L, start = measurement, fixed = "measurement")
  expect_equal(unname(held$means), unname(fit$means), tolerance = 1e-10)
  expect_identical(nrow(get_results(measurement, "profiles")), 4L)
  # The mixing blocks are deliberately absent, so their table is empty rather
  # than an error.
  expect_identical(nrow(get_results(measurement, "profile_probabilities")), 0L)
})

test_that("every table of a starting-value set is tidy", {
  start <- starting_values(tidy_fit())
  tables <- c("profiles", "responses", "profile_probabilities", "covariances")
  invisible(lapply(tables, function(what) {
    expect_tidy_frame(get_results(start, what))
  }))
  profiles <- as.data.frame(start)
  expect_identical(names(profiles),
                   c("profile", "indicator", "mean", "variance",
                     "standard_deviation"))
  expect_identical(nrow(profiles), 4L)
  expect_equal(profiles$mean, as.vector(t(start$means)))
  probabilities <- get_results(start, "profile_probabilities")
  expect_identical(nrow(probabilities), 4L)
  totals <- tapply(probabilities$probability, probabilities$group_class, sum)
  expect_equal(unname(as.vector(totals)), rep(1, 2L))
  expect_identical(nrow(get_results(start, "responses")), 0L)
})

test_that("a categorical starting-value set tabulates its response blocks", {
  start <- starting_values(tidy_categorical_fit())
  responses <- get_results(start, "responses")
  expect_tidy_frame(responses)
  expect_identical(names(responses),
                   c("profile", "indicator", "category", "probability"))
  expect_identical(nrow(responses), 12L)
  totals <- tapply(responses$probability,
                   list(responses$profile, responses$indicator), sum)
  expect_equal(as.vector(totals), rep(1, 6L))
  expect_identical(nrow(get_results(start, "profiles")), 0L)
})

test_that("a full-covariance start tabulates the variances it implies", {
  start <- starting_values(tidy_fit(covariance_model = "full"))
  expect_true("covariances" %in% names(start))
  expect_null(start$variances)
  covariances <- get_results(start, "covariances")
  expect_tidy_frame(covariances)
  expect_equal(covariances$covariance, as.vector(start$covariances))
  # The Gaussian table must still report a variance, read off the diagonal.
  profiles <- get_results(start, "profiles")
  expect_identical(nrow(profiles), 4L)
  expect_true(all(profiles$variance > 0))
  expect_equal(profiles$standard_deviation, sqrt(profiles$variance))
})

test_that("printing a starting-value set describes it and returns it", {
  start <- starting_values(tidy_fit())
  expect_output(printed <- withVisible(print(start)),
                "2 profile\\(s\\), 2 continuous indicator\\(s\\)")
  expect_false(printed$visible)
  expect_identical(printed$value, start)
  expect_output(print(start), "Residual covariance: diagonal")
  expect_output(print(starting_values(tidy_fit(covariance_model = "full"))),
                "Residual covariance: full")
  expect_output(print(starting_values(tidy_categorical_fit())),
                "Categorical indicators: 3")
  expect_output(print(starting_values(tidy_fit(), what = "measurement")),
                "measurement only")
  printed_text <- capture.output(print(start))
  expect_false(any(grepl("$", printed_text, fixed = TRUE)))
})

test_that("starting_values still refuses an incomplete object by class", {
  expect_error(starting_values(list(means = matrix(0, 2L, 2L))),
               class = "latents_bad_start")
  expect_error(starting_values(list(means = matrix(0, 2L, 2L),
                                    variances = matrix(1, 2L, 2L),
                                    profile_probabilities = matrix(0.5, 1L, 2L),
                                    group_probabilities = 1),
                               covariance = "keep"),
               class = "latents_bad_start")
})
