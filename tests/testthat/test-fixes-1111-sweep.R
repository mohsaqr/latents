# Repairs of 2026-11-11, found by the argument sweep
# (equivalence/argument-sweep/): option values refused with the package's
# class, count arguments that changed nothing, tables and methods the general
# transition fit lacked, and non-finite information reaching eigen().

sweep_fix_panel <- function() {
  withr::with_seed(707L, {
    # A two-state chain that moves, so that no transition sits on its bound.
    first <- sample(1:2, 80L, TRUE)
    second <- ifelse(stats::runif(80L) < 0.3, 3L - first, first)
    third <- ifelse(stats::runif(80L) < 0.3, 3L - second, second)
    state <- as.vector(rbind(first, second, third))
    data.frame(id = rep(seq_len(80L), each = 3L), wave = rep(1:3, 80L),
               y1 = stats::rnorm(240L, c(0, 3)[state]),
               y2 = stats::rnorm(240L, c(0, -2)[state]) + 0.3 * stats::rnorm(240L))
  })
}

test_that("an unknown option value is refused with latents_bad_argument", {
  panel <- sweep_fix_panel()
  transition <- function(...) lta(panel, c("y1", "y2"), id = "id", time = "wave",
                                  n_profiles = 2L, n_starts = 1L, seed = 1L, ...)
  expect_error(transition(transitions = "free"), class = "latents_bad_argument")
  expect_error(transition(measurement = "free"), class = "latents_bad_argument")
  expect_error(transition(occasions = "all"), class = "latents_bad_argument")
  expect_error(transition(search = "random"), class = "latents_bad_argument")
  expect_error(transition(missing = "drop"), class = "latents_bad_argument")
  expect_error(transition(variance_model = "pooled"), class = "latents_bad_argument")
  regression <- function(...) mixture_regression(score ~ hours, study_hours, n_classes = 2L, ...)
  expect_error(regression(family = "gamma"), class = "latents_bad_argument")
  expect_error(regression(class_level = "unit"), class = "latents_bad_argument")
  expect_error(regression(variance = "pooled"), class = "latents_bad_argument")
  expect_error(regression(missing = "fiml"), class = "latents_bad_argument")
  expect_error(regression(vcov_type = "hc3"), class = "latents_bad_argument")
  expect_error(regression(select_start = "first"), class = "latents_bad_argument")
  # The message names the argument and the accepted values.
  expect_error(regression(family = "gamma"), "`family` must be one of")
})

test_that("a matched option keeps its meaning: defaults and partial names", {
  panel <- sweep_fix_panel()
  default <- lta(panel, c("y1", "y2"), id = "id", time = "wave", n_profiles = 2L,
                 n_starts = 2L, seed = 1L)
  spelled <- lta(panel, c("y1", "y2"), id = "id", time = "wave", n_profiles = 2L,
                 n_starts = 2L, seed = 1L, transitions = "homogeneous", occasions = "obs",
                 measurement = "invariant", variance_model = "varying")
  expect_identical(default$log_likelihood, spelled$log_likelihood)
})

test_that("count arguments that change nothing warn, and used ones do not", {
  expect_warning(lpa(srl, names(srl), 2L, n_starts = 1L, seed = 1L, count_dispersion = "equal"),
                 class = "latents_unused_argument")
  expect_warning(lpa(srl, names(srl), 2L, n_starts = 1L, seed = 1L,
                     count_model = "negative_binomial"),
                 class = "latents_unused_argument")
  counts <- withr::with_seed(505L, data.frame(visits = stats::rnbinom(200L, mu = 4, size = 2),
                                              posts = stats::rpois(200L, 3)))
  # Poisson counts have no dispersion to share.
  expect_warning(lpa(counts, names(counts), 2L, count = names(counts), n_starts = 1L, seed = 1L,
                     count_dispersion = "equal"), class = "latents_unused_argument")
  # One profile of overdispersed counts: an interior fit, so nothing else warns.
  expect_no_warning(lpa(counts, "visits", 1L, count = "visits", n_starts = 1L, seed = 1L,
                        count_model = "negative_binomial", count_dispersion = "equal"))
  # The warning changes nothing else: the fit is the one without the argument.
  quiet <- suppressWarnings(lpa(srl, names(srl), 2L, n_starts = 1L, seed = 1L,
                                count_dispersion = "equal"), classes = "latents_unused_argument")
  expect_identical(quiet$log_likelihood,
                   lpa(srl, names(srl), 2L, n_starts = 1L, seed = 1L)$log_likelihood)
})

test_that("a general transition fit reports its covariances", {
  panel <- sweep_fix_panel()
  two <- panel[panel$wave <= 2L, ]
  # With two occasions there is one transition, so the occasion-specific model
  # (the general engine) is the homogeneous one (the other engine), whose
  # covariance table is the reference.
  general <- lta(two, c("y1", "y2"), id = "id", time = "wave", n_profiles = 2L, n_starts = 4L,
                 seed = 1L, covariance_model = "full", transitions = "occasion", tol = 1e-12)
  homogeneous <- lta(two, c("y1", "y2"), id = "id", time = "wave", n_profiles = 2L,
                     n_starts = 4L, seed = 1L, covariance_model = "full", tol = 1e-12)
  expect_s3_class(general, "multilpa_lta")
  expect_equal(general$log_likelihood, homogeneous$log_likelihood, tolerance = 1e-8)
  table <- get_results(general, "covariances")
  expect_named(table, c("occasion", "profile", "indicator", "indicator_2", "covariance"))
  expect_identical(nrow(table), 8L)
  expect_true("covariances" %in% names(get_results(general, "all")))
  reference <- get_results(homogeneous, "covariances")
  expect_equal(sort(table$covariance), sort(reference$covariance), tolerance = 1e-5)
  # Symmetric, and its diagonal is the variance column of "profiles".
  mirrored <- merge(table, table, by.x = c("profile", "indicator", "indicator_2"),
                    by.y = c("profile", "indicator_2", "indicator"))
  expect_equal(mirrored$covariance.x, mirrored$covariance.y)
  diagonal <- subset(table, indicator == indicator_2)
  profiles <- get_results(general, "profiles")
  expect_equal(diagonal$covariance[order(diagonal$profile, diagonal$indicator)],
               profiles$variance[order(profiles$profile, profiles$indicator)])
  # A diagonal fit has no such table, as before.
  diagonal_fit <- lta(two, c("y1", "y2"), id = "id", time = "wave", n_profiles = 2L,
                      n_starts = 1L, seed = 1L, transitions = "occasion")
  expect_false("covariances" %in% names(get_results(diagonal_fit, "all")))
  expect_error(get_results(diagonal_fit, "covariances"), class = "latents_bad_argument")
})

test_that("predict() refuses a general transition fit with the transition fits' class", {
  panel <- sweep_fix_panel()
  general <- lta(panel, c("y1", "y2"), id = "id", time = "wave", n_profiles = 2L, n_starts = 1L,
                 seed = 1L, transitions = "occasion")
  expect_error(predict(general), class = "latents_unsupported_prediction")
})

test_that("non-finite information is not positive definite rather than an eigen() error", {
  broken <- matrix(c(1, NaN, NaN, 1), 2L)
  expect_false(.inference_positive_definite(broken, "min_le_eps_maxabs"))
  expect_false(.inference_positive_definite(matrix(c(1, Inf, Inf, 1), 2L), "min_le_eps_maxabs"))
  expect_error(.inference_positive_definite(broken, "min_le_eps_max"),
               class = "latents_singular_information")
  expect_error(.inference_positive_definite(broken, "ratio_le_eps"),
               class = "latents_singular_information")
  # Finite matrices are judged as before.
  expect_true(.inference_positive_definite(diag(2), "min_le_eps_maxabs"))
  expect_false(.inference_positive_definite(diag(c(1, 1e-12)), "min_le_eps_maxabs"))
  # The mixture engine's inverse reports NA standard errors with a warning.
  expect_warning(inverse <- .mixture_invert(broken), class = "latents_singular_information")
  expect_true(all(is.na(inverse)))
})

test_that("a general transition fit reports posteriors, entropy and criteria", {
  skip_on_cran()
  set.seed(3)
  n <- 150L
  occasions <- 4L
  state <- Reduce(function(previous, t) {
    ifelse(stats::runif(n) < 0.8, previous, 3L - previous)
  }, seq_len(occasions - 1L), accumulate = TRUE, sample(1:2, n, TRUE))
  state <- as.vector(t(do.call(cbind, state)))
  data <- data.frame(id = rep(seq_len(n), each = occasions),
                     t = rep(seq_len(occasions), n),
                     y1 = stats::rnorm(n * occasions, c(0, 3)[state]),
                     y2 = stats::rnorm(n * occasions, c(0, 2)[state]))
  fit <- quietly(lta(data, c("y1", "y2"), id = "id", time = "t", n_profiles = 2,
                     transitions = "occasion", n_starts = 2, seed = 1))
  expect_s3_class(fit, "multilpa_lta")
  posterior <- fit$subject_posteriors
  entropy_sum <- -sum(ifelse(posterior > 0, posterior * log(posterior), 0))
  entropy <- get_results(fit, "entropy")
  individuals <- subset(entropy, level == "individuals")
  expect_equal(individuals$entropy_sum, entropy_sum, tolerance = 1e-10)
  expect_equal(individuals$relative_entropy,
               1 - entropy_sum / (nrow(posterior) * log(2)), tolerance = 1e-10)
  criteria <- get_results(fit, "information_criteria")
  deviance <- -2 * fit$log_likelihood
  expect_equal(criteria$aic, deviance + 2 * fit$n_parameters, tolerance = 1e-10)
  expect_equal(criteria$bic_individual, fit$bic_individual, tolerance = 1e-10)
  expect_equal(criteria$icl_individual,
               deviance + fit$n_parameters * log(nrow(posterior)) + 2 * entropy_sum,
               tolerance = 1e-10)
  long <- get_results(fit, "posteriors")
  expect_identical(nrow(long), nrow(posterior) * 2L)
  # Invariant: each row's posteriors sum to one.
  expect_equal(unname(tapply(long$posterior, long$row, sum)), rep(1, nrow(posterior)),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_identical(nrow(get_results(fit, "posteriors", format = "wide")), nrow(posterior))
  expect_true(all(c("posteriors", "entropy", "information_criteria") %in%
                    names(get_results(fit, "all"))))
  expect_error(get_results(fit, "entropy", format = "wide"), class = "latents_bad_argument")
  expect_error(get_results(fit, "posteriors", format = "tall"), class = "latents_bad_argument")
  expect_error(starting_values(fit), class = "latents_bad_start")
  expect_error(starting_values(fit), "takes no")
})

test_that("an unknown plot view of a regression or a pooled result is refused by class", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  regression <- quietly(mixture_regression(hours ~ week, study_hours, n_classes = 2,
                                           n_starts = 2, seed = 1))
  expect_error(plot(regression, what = "radar"), class = "latents_bad_argument")
  expect_s3_class(plot(regression, what = "coefficients"), "ggplot")
  set.seed(4)
  # Two well-separated profiles, so neither completed data set fits on a bound.
  profile <- rep(1:2, each = 100)
  separated <- data.frame(a = stats::rnorm(200, c(0, 4)[profile]),
                          b = stats::rnorm(200, c(0, -3)[profile]))
  completed <- lapply(1:2, function(i) {
    transform(separated, a = a + stats::rnorm(200, sd = 0.1))
  })
  pooled <- quietly(suppressMessages(pool_imputations(
    completed, c("a", "b"), NULL, n_profiles = 2, n_starts = 2, seed = 1)))
  expect_error(plot(pooled, what = "radar"), class = "latents_bad_argument")
  # The default view and each named one still draw.
  expect_s3_class(plot(pooled), "ggplot")
  expect_s3_class(plot(pooled, what = "measurement"), "ggplot")
  expect_s3_class(plot(pooled, level = "measurement"), "ggplot")
})
