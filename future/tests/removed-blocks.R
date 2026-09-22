# Test blocks removed from the shipped suite when `fit_random_intercept()`
# and `lmr_lrt()` moved to `future/`. Kept verbatim so that restoring a
# feature restores its tests with it. Each block is headed by the file it
# came from.

# --- from test-convenience.R ---
test_that("a model family without residuals says so instead of failing", {
  # Quadrature over 1422 observations is the slowest fit in this file, so it
  # runs on three indicators from a single start: the assertions below are
  # about what `diagnostics()` does with this model family, not about how well
  # the fit itself recovers anything.
  skip_on_cran()
  fit <- fit_random_intercept(course_engagement,
                              c("browse", "lectures", "forum_read"),
                              id = "student", n_profiles = 2, n_starts = 1, seed = 1)
  quality <- diagnostics(fit)
  expect_s3_class(get_data(quality, "entropy"), "data.frame")
  expect_error(get_data(quality, "residuals"),
               class = "multilpa_no_group_classes")
  expect_output(print(quality), "not available for this model family")
})


# --- from test-math-results-audit.R ---
test_that("random-intercept diagnostics distinguish profiles from group classes", {
  set.seed(174)
  d <- data.frame(g = rep(1:20, each = 5),
                  y = rep(rnorm(20), each = 5) + rnorm(100))
  fit <- fit_random_intercept(d, "y", "g", n_profiles = 1, n_starts = 1)
  indices <- get_data(fit, "information_criteria", format = "long")
  expect_equal(subset(indices, criterion == "bic" & convention == "groups")$value,
               fit$bic)
  expect_equal(subset(indices, criterion == "bic" & convention == "individuals")$value,
               fit$bic_individual)
  expect_true(all(is.na(subset(indices, convention == "groups" &
                                criterion %in% c("awe", "icl"))$value)))
  expect_equal(get_data(fit, "classification", level = "both")$level, "individuals")
  expect_error(get_data(fit, "classification", level = "groups"), "no discrete group")
  expect_equal(get_data(fit, "entropy")$level, "individuals")
  expect_equal(get_data(fit, "entropy")$entropy_sum, 0)
})


# --- from test-tidy-diagnostics.R ---
test_that("a model with no discrete group classes refuses by condition class", {
  set.seed(174)
  data <- data.frame(g = rep(seq_len(20L), each = 5L),
                     y = rep(stats::rnorm(20L), each = 5L) + stats::rnorm(100L))
  fit <- fit_random_intercept(data, "y", "g", n_profiles = 1, n_starts = 1)
  expect_error(get_data(fit, "classification", level = "groups"),
               class = "multilpa_no_group_classes")
  expect_error(get_data(fit, "average_posteriors", level = "groups"),
               class = "multilpa_no_group_classes")
  expect_error(get_data(fit, "residuals", data = data),
               class = "multilpa_no_group_classes")
  # "both" degrades to the level the model actually has.
  expect_identical(unique(get_data(fit, "classification", level = "both")$level),
                   "individuals")
  expect_identical(unique(get_data(fit, "average_posteriors", level = "both")$level),
                   "individuals")
})


# --- from test-get-data.R ---
test_that("a level a fit has not is refused rather than silently dropped", {
  intercepts <- fit_random_intercept(course_engagement, activity, "student",
                                     n_profiles = 2, n_starts = 1, seed = 1)
  expect_error(get_data(intercepts, "classification", level = "groups"),
               class = "multilpa_no_group_classes")
})


# --- from test-report-surface.R ---
test_that("plot(diagnostics()) draws for a random-intercept fit", {
  skip_on_cran()
  fit <- fit_random_intercept(.surface_data(), c("a", "b"), "school",
                              n_profiles = 2L, n_starts = 2L, seed = 1)
  draw(expect_invisible(plot(diagnostics(fit))))
})


# --- from test-diagnostics.R ---
test_that("lmr_lrt reports a statistic and withholds a p-value", {
  dat <- make_two_level()
  smaller <- multilpa(dat, c("a", "b"), "g", 1, 1, n_starts = 2, seed = 5)
  larger <- multilpa(dat, c("a", "b"), "g", 2, 1, n_starts = 5, seed = 5)
  result <- lmr_lrt(smaller, larger)
  expect_s3_class(result, "data.frame")
  expect_identical(nrow(result), 1L)
  expect_equal(result$statistic,
               2 * (larger$log_likelihood - smaller$log_likelihood))
  expect_identical(result$df, larger$n_parameters - smaller$n_parameters)
  expect_equal(result$n, smaller$n_observations)
  expect_equal(result$adjustment_factor, 1 + 1 / (result$df * log(result$n)))
  expect_equal(result$adjusted_statistic, result$statistic / result$adjustment_factor)
  # The adjustment always shrinks a positive statistic.
  expect_lt(result$adjusted_statistic, result$statistic)
  expect_true(is.na(result$p_value))
  grouped <- lmr_lrt(smaller, larger, n = "groups")
  expect_equal(grouped$n, smaller$n_groups)
  expect_equal(grouped$statistic, result$statistic)
  expect_gt(grouped$adjustment_factor, result$adjustment_factor)
})


# --- from test-diagnostics.R ---
test_that("lmr_lrt rejects invalid comparisons by condition class", {
  dat <- make_two_level()
  smaller <- multilpa(dat, c("a", "b"), "g", 1, 1, n_starts = 2, seed = 5)
  larger <- multilpa(dat, c("a", "b"), "g", 2, 1, n_starts = 5, seed = 5)
  expect_error(lmr_lrt(larger, smaller), class = "multilpa_bad_nesting")
  other <- multilpa(make_two_level(seed = 12L, n_groups = 20L),
                      c("a", "b"), "g", 2, 1, n_starts = 3, seed = 5)
  expect_error(lmr_lrt(smaller, other),
               class = "multilpa_incomparable_models")
  expect_error(lmr_lrt(smaller, larger, n = 1), "greater than one")
  reversed <- larger
  reversed$log_likelihood <- smaller$log_likelihood - 1
  expect_error(lmr_lrt(smaller, reversed),
               class = "multilpa_reversed_likelihood")
})


# --- from test-math-extensions-audit.R ---
test_that("univariate random-intercept mixtures retain profile-by-indicator shapes", {
  set.seed(4)
  data <- data.frame(g = rep(seq_len(12L), each = 6L),
                     y = rnorm(72L, rep(c(-2, 2), 36L), 0.5))
  fits <- lapply(c("varying", "equal"), function(variance_model)
    fit_random_intercept(data, "y", "g", 2L, n_starts = 1L,
                         variance_model = variance_model, seed = 1L))
  invisible(lapply(fits, function(fit) {
    expect_identical(dim(fit$means), c(2L, 1L))
    expect_identical(dim(fit$variances), c(2L, 1L))
    expect_equal(rowSums(fit$subject_posteriors), rep(1, 72L), tolerance = 1e-12)
    expect_true(is.finite(fit$log_likelihood))
    expect_true(fit$converged)
    expect_equal(unname(sort(drop(fit$means))), c(-2, 2), tolerance = 0.2)
  }))
})


# --- from test-math-extensions-audit.R ---
test_that("mixture random-intercept posteriors agree with adaptive integration", {
  x <- matrix(c(-0.5, 1, 0.2, 0.7, 0.9, -0.3), 3L, 2L)
  parameters <- list(means = matrix(c(-0.6, 0.8, -0.2, 0.6), 2L),
    variances = matrix(c(0.5, 0.8, 0.6, 0.9), 2L), random_sd = 0.35,
    profile_probabilities = c(0.3, 0.7))
  actual <- .ri_evaluate(x, c(1L, 1L, 2L), parameters, .ri_quadrature(61L), TRUE)
  reference <- lapply(list(1:2, 3L), function(rows) {
    density <- function(intercept) {
      vapply(seq_len(2L), function(k) parameters$profile_probabilities[k] *
        dnorm(x[rows, 1L], parameters$means[k, 1L] + intercept,
              sqrt(parameters$variances[k, 1L])) *
        dnorm(x[rows, 2L], parameters$means[k, 2L] + intercept,
              sqrt(parameters$variances[k, 2L])), numeric(length(rows))) |>
        matrix(nrow = length(rows), ncol = 2L)
    }
    integrate_weighted <- function(weight) {
      integrate(function(intercepts) vapply(intercepts, function(intercept) {
        components <- density(intercept)
        likelihood <- prod(rowSums(components)) * dnorm(intercept, sd = parameters$random_sd)
        if (likelihood == 0) 0 else likelihood * weight(intercept, components)
      }, numeric(1)), -Inf, Inf, rel.tol = 1e-11)$value
    }
    normalizer <- integrate_weighted(function(intercept, components) 1)
    mean <- integrate_weighted(function(intercept, components) intercept) / normalizer
    second <- integrate_weighted(function(intercept, components) intercept^2) / normalizer
    posteriors <- vapply(seq_along(rows), function(row)
      integrate_weighted(function(intercept, components)
        components[row, 1L] / sum(components[row, ])) / normalizer, numeric(1))
    list(mean = mean, sd = sqrt(second - mean^2), posterior = posteriors)
  })
  expect_equal(unname(actual$random_intercept_mean),
               vapply(reference, `[[`, numeric(1), "mean"), tolerance = 1e-10)
  expect_equal(unname(actual$random_intercept_sd),
               vapply(reference, `[[`, numeric(1), "sd"), tolerance = 1e-10)
  expect_equal(actual$subject_posteriors[, 1L],
               unlist(lapply(reference, `[[`, "posterior")), tolerance = 1e-10)
})


# --- from test-math-extensions-audit.R ---
test_that("LMR rejects equally sized fits on different data or group layouts", {
  object <- structure(list(n_observations = 4L, n_groups = 2L,
    vars = "y", group = "g", group_values = 1:2,
    group_index = c(1L, 1L, 2L, 2L), indicator_data = matrix(1:4, 4L),
    n_parameters = 2L, log_likelihood = -10), class = "multilpa")
  alternative <- object
  alternative$n_parameters <- 3L
  alternative$log_likelihood <- -9
  expect_equal(lmr_lrt(object, alternative)$statistic, 2)
  alternative$indicator_data[1L, 1L] <- 9L
  expect_error(lmr_lrt(object, alternative), class = "multilpa_incomparable_models")
  alternative$indicator_data <- object$indicator_data
  alternative$group_index <- c(1L, 2L, 1L, 2L)
  expect_error(lmr_lrt(object, alternative), class = "multilpa_incomparable_models")
  alternative$group_index <- object$group_index
  invisible(lapply(c("covariance_model", "variance_model", "min_variance", "min_probability"),
    function(field) {
      changed <- alternative
      changed[[field]] <- "different"
      expect_error(lmr_lrt(object, changed), class = "multilpa_incomparable_models")
    }))
  object$n_profiles <- 1L
  object$n_group_classes <- 2L
  alternative$n_profiles <- 3L
  alternative$n_group_classes <- 1L
  expect_error(lmr_lrt(object, alternative), class = "multilpa_bad_nesting")
  object$n_groups <- alternative$n_groups <- 1L
  expect_error(lmr_lrt(object, alternative, n = "groups"), "greater than one")
})


# --- from test-tidy-classes.R ---
test_that("a random-intercept summary is classed, tidy and fully named", {
  fit <- .tidy_random_intercept_fit()
  summary_object <- summary.multilpa_random_intercept(fit)
  expect_s3_class(summary_object, "summary_multilpa_random_intercept")
  expect_false(anyNA(names(summary_object)))
  expect_output(print.summary_multilpa_random_intercept(summary_object),
                "Random-intercept LPA")
  model <- get_data(summary_object, "model")
  expect_identical(nrow(model), 1L)
  expect_equal(model$random_intercept_sd, fit$random_sd)
  expect_equal(model$log_likelihood, fit$log_likelihood)
  expect_true(model$quadrature_check_passed)
  intercepts <- get_data(summary_object, "random_intercepts")
  expect_identical(nrow(intercepts), fit$n_groups)
  expect_identical(names(intercepts), c("group", "group_size", "mean", "sd"))
  expect_identical(nrow(get_data(summary_object, "starts")), 2L)
})


# --- from test-tidy-classes.R ---
test_that("random-intercept coefficients are named, natural and complete", {
  fit <- .tidy_random_intercept_fit()
  estimates <- coef.multilpa_random_intercept(fit)
  expect_type(estimates, "double")
  # A named vector is the base generic's contract; the tidy table is
  # as.data.frame(). Every free parameter must appear exactly once.
  expect_identical(length(estimates), as.integer(fit$n_parameters))
  expect_false(anyDuplicated(names(estimates)) > 0L)
  expect_true(all(is.finite(estimates)))
  # Names follow the package-wide level.parameter.outcome.term grammar, so a
  # name reads back into the columns parameter_inference() reports.
  expect_true(all(vapply(strsplit(names(estimates), ".", fixed = TRUE), length,
                         integer(1)) >= 3L))
  profiles <- get_data(fit, "profiles")
  expect_equal(unname(estimates[["measurement.mean.profile_1.score_a"]]),
               subset(profiles, profile == 1L & indicator == "score_a")$mean)
  expect_equal(unname(estimates[["measurement.variance.profile_1.score_a"]]),
               subset(profiles, profile == 1L & indicator == "score_a")$variance)
  expect_equal(unname(estimates[["group.standard_deviation.random_intercept"]]),
               fit$random_sd)
  # The logit is relative to the final profile, so exponentiating it returns
  # the ratio of the two profile probabilities.
  expect_equal(exp(unname(estimates[["profile.logit.profile_1"]])),
               fit$profile_probabilities[1L] / fit$profile_probabilities[2L])
  shared <- .tidy_random_intercept_fit(variance_model = "equal")
  shared_estimates <- coef.multilpa_random_intercept(shared)
  # An equal-variance fit estimates one variance per indicator, not one per
  # profile and indicator, so it must not be reported twice.
  expect_identical(length(shared_estimates), as.integer(shared$n_parameters))
  expect_true(all(c("measurement.variance.shared.score_a",
                    "measurement.variance.shared.score_b") %in%
                    names(shared_estimates)))
  # Every reported variance must be a variance, not the log coordinate it is
  # searched on: a negative one would be a scale mix-up, not an estimate.
  variances <- shared_estimates[grepl("^measurement\\.variance\\.",
                                      names(shared_estimates))]
  expect_true(all(variances > 0))
  # A single-profile fit estimates no logit and must still name every
  # parameter exactly once.
  single <- fit_random_intercept(
    data.frame(group = rep(seq_len(12), each = 5), score_a = rnorm(60)),
    "score_a", "group", 1, n_starts = 1, seed = 4)
  single_estimates <- coef.multilpa_random_intercept(single)
  expect_identical(length(single_estimates), as.integer(single$n_parameters))
  expect_false(any(grepl("logit", names(single_estimates))))
})


# --- from test-tidy-classes.R ---
test_that("a random-intercept fit refuses standard errors by condition class", {
  fit <- .tidy_random_intercept_fit()
  expect_error(vcov.multilpa_random_intercept(fit),
               class = "multilpa_no_inference")
  expect_error(parameter_inference.multilpa_random_intercept(fit),
               class = "multilpa_no_inference")
  expect_error(confint.multilpa_random_intercept(fit),
               class = "multilpa_no_inference")
})


# --- from test-tidy-classes.R ---
test_that("random-intercept plots draw both panels and return the fit", {
  fit <- .tidy_random_intercept_fit()
  file <- tempfile(fileext = ".pdf")
  on.exit(unlink(file), add = TRUE)
  grDevices::pdf(file)
  on.exit(grDevices::dev.off(), add = TRUE, after = FALSE)
  expect_s3_class(plot.multilpa_random_intercept(fit),
                  "multilpa_random_intercept")
  expect_s3_class(plot.multilpa_random_intercept(fit, scale = "standardized"),
                  "multilpa_random_intercept")
  expect_s3_class(plot.multilpa_random_intercept(fit, what = "random_intercepts"),
                  "multilpa_random_intercept")
  expect_error(plot.multilpa_random_intercept(fit, what = "nonsense"))
})

