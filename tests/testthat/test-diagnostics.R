make_two_level <- function(seed = 11L, n_groups = 30L, per_group = 8L) {
  set.seed(seed)
  group <- rep(seq_len(n_groups), each = per_group)
  group_class <- rep(c(1L, 2L), length.out = n_groups)
  prevalence <- rbind(c(0.85, 0.15), c(0.2, 0.8))
  profile <- vapply(group, function(j) {
    sample.int(2L, 1L, prob = prevalence[group_class[j], ])
  }, integer(1))
  means <- rbind(c(-2, -2), c(2, 2))
  y <- t(vapply(profile, function(k) stats::rnorm(2L, means[k, ], 1), numeric(2)))
  data.frame(g = group, a = y[, 1L], b = y[, 2L])
}

test_that("information criteria match their documented formulas", {
  dat <- make_two_level()
  fit <- fit_multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 6, seed = 5)
  indices <- information_criteria(fit)
  expect_s3_class(indices, "data.frame")
  expect_identical(names(indices),
    c("criterion", "convention", "n", "value", "penalty", "definition"))
  expect_identical(nrow(indices), 12L)
  pick <- function(criterion_name, convention_name) {
    indices$value[indices$criterion == criterion_name &
                    indices$convention == convention_name]
  }
  q <- fit$n_parameters
  log_likelihood <- fit$log_likelihood
  n_individuals <- fit$n_observations
  n_groups <- fit$n_groups
  entropy_individuals <- .multilpa_entropy_sum(fit$subject_posteriors)
  entropy_groups <- .multilpa_entropy_sum(fit$group_posteriors)
  expect_equal(pick("log_likelihood", "none"), log_likelihood)
  expect_equal(pick("aic", "none"), fit$aic)
  expect_equal(pick("bic", "groups"), fit$bic)
  expect_equal(pick("bic", "individuals"), fit$bic_individual)
  expect_equal(pick("sabic", "individuals"),
               -2 * log_likelihood + q * log((n_individuals + 2) / 24))
  expect_equal(pick("sabic", "groups"),
               -2 * log_likelihood + q * log((n_groups + 2) / 24))
  expect_equal(pick("caic", "individuals"),
               -2 * log_likelihood + q * (log(n_individuals) + 1))
  expect_equal(pick("icl", "individuals"),
               -2 * log_likelihood + q * log(n_individuals) + 2 * entropy_individuals)
  expect_equal(pick("icl", "groups"),
               -2 * log_likelihood + q * log(n_groups) + 2 * entropy_groups)
  expect_equal(pick("awe", "individuals"),
               -2 * (log_likelihood - entropy_individuals) +
                 2 * q * (1.5 + log(n_individuals)))
  # Every penalty is nonnegative here and the criteria exceed -2L accordingly.
  penalties <- indices$penalty[!is.na(indices$penalty)]
  expect_true(all(penalties > 0))
})

test_that("information criteria order candidate models sensibly", {
  dat <- make_two_level()
  one <- fit_multilpa(dat, c("a", "b"), "g", 1, 1, n_starts = 3, seed = 5)
  two <- fit_multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 6, seed = 5)
  criteria <- c("bic", "sabic", "caic", "icl", "awe")
  better <- vapply(criteria, function(criterion) {
    one_value <- information_criteria(one)
    two_value <- information_criteria(two)
    select <- function(indices) {
      indices$value[indices$criterion == criterion &
                      indices$convention == "individuals"]
    }
    select(two_value) < select(one_value)
  }, logical(1))
  # The data are generated with two well-separated profiles, so every
  # criterion must prefer the two-profile model.
  expect_true(all(better))
})

test_that("classification diagnostics are internally consistent", {
  dat <- make_two_level()
  fit <- fit_multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 6, seed = 5)
  summary_table <- classification_table(fit, level = "both")
  expect_identical(nrow(summary_table), 4L)
  expect_true(all(summary_table$average_posterior >= 0 &
                    summary_table$average_posterior <= 1))
  individuals <- subset(summary_table, level == "individuals")
  expect_equal(sum(individuals$n_modal), fit$n_observations)
  expect_equal(sum(individuals$estimated_n), fit$n_observations)
  expect_equal(sum(individuals$estimated_proportion), 1)
  groups <- subset(summary_table, level == "groups")
  expect_equal(sum(groups$n_modal), fit$n_groups)
  expect_equal(sum(groups$estimated_n), fit$n_groups)
  # Well-separated profiles must classify far better than chance.
  expect_true(all(individuals$odds_correct_classification > 5))
  detail <- classification_table(fit, level = "individuals", detail = TRUE)
  expect_identical(nrow(detail), 4L)
  # Average posteriors within each modal group are a probability distribution.
  totals <- tapply(detail$average_posterior, detail$assigned_class, sum)
  expect_equal(unname(as.vector(totals)), rep(1, 2))
})

test_that("classification diagnostics handle a class with no modal members", {
  set.seed(3)
  dat <- data.frame(g = rep(seq_len(15L), each = 6L), y = stats::rnorm(90L))
  fit <- fit_multilpa(dat, "y", "g", 1, 1, n_starts = 1, seed = 2)
  table <- classification_table(fit)
  expect_identical(nrow(table), 1L)
  expect_equal(table$average_posterior, 1)
  # A single class has proportion one, so the odds are undefined, not infinite.
  expect_true(is.na(table$odds_correct_classification))
})

test_that("entropy is bounded and undefined for a single class", {
  dat <- make_two_level()
  fit <- fit_multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 6, seed = 5)
  entropy <- entropy_table(fit)
  expect_identical(nrow(entropy), 2L)
  expect_true(all(entropy$relative_entropy > 0 & entropy$relative_entropy <= 1))
  expect_true(all(entropy$entropy_sum >= 0))
  single <- fit_multilpa(dat, c("a", "b"), "g", 1, 1, n_starts = 1, seed = 5)
  expect_true(all(is.na(entropy_table(single)$relative_entropy)))
  # A degenerate, perfectly separated posterior has zero entropy.
  expect_equal(.multilpa_entropy_sum(matrix(c(1, 0, 0, 1), 2L, 2L)), 0)
  expect_equal(.multilpa_relative_entropy(matrix(c(1, 0, 0, 1), 2L, 2L)), 1)
  # A maximally uncertain posterior has zero relative entropy.
  expect_equal(.multilpa_relative_entropy(matrix(0.5, 4L, 2L)), 0)
})

test_that("tidy accessors return the documented shapes", {
  dat <- make_two_level()
  fit <- fit_multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 6, seed = 5)
  profiles <- as.data.frame(fit)
  expect_identical(names(profiles),
    c("profile", "indicator", "mean", "variance", "standard_deviation"))
  expect_identical(nrow(profiles), 4L)
  expect_equal(profiles$mean, as.vector(t(fit$means)))
  expect_equal(profiles$standard_deviation, sqrt(profiles$variance))
  probabilities <- as.data.frame(fit, what = "profile_probabilities")
  expect_identical(nrow(probabilities), 4L)
  totals <- tapply(probabilities$probability, probabilities$group_class, sum)
  expect_equal(unname(as.vector(totals)), rep(1, 2))
  posteriors <- as.data.frame(fit, what = "posteriors")
  expect_identical(nrow(posteriors), fit$n_observations)
  expect_identical(posteriors$row, seq_len(fit$n_observations))
  expect_equal(posteriors$group, dat$g)
  group_posteriors <- as.data.frame(fit, what = "group_posteriors")
  expect_identical(nrow(group_posteriors), fit$n_groups)
  expect_equal(sum(group_posteriors$group_size), fit$n_observations)
  expect_equal(sum(group_posteriors$log_likelihood), fit$log_likelihood)
  expect_identical(nrow(as.data.frame(fit, what = "starts")), 6L)
  expect_identical(as.data.frame(fit, what = "information_criteria"),
                   information_criteria(fit))
  expect_identical(as.data.frame(fit, what = "entropy"), entropy_table(fit))
  expect_identical(as.data.frame(fit, what = "classification", level = "groups"),
                   classification_table(fit, level = "groups"))
})

test_that("enumeration and inference tidy and print", {
  dat <- make_two_level()
  candidates <- enumerate_multilpa(dat, c("a", "b"), "g", profiles = 1:2,
                                 group_classes = 1:2, n_starts = 3, seed = 3)
  expect_s3_class(candidates, "multilpa_enumeration")
  grid <- as.data.frame(candidates)
  expect_identical(nrow(grid), 4L)
  expect_true(all(c("sabic_groups", "sabic_individual", "caic_groups",
                    "caic_individual", "awe_groups", "awe_individual",
                    "icl_groups", "icl_individual") %in% names(grid)))
  # A candidate that fails must be retained with NA criteria, not dropped.
  failed <- subset(grid, !converged)
  expect_true(all(is.na(failed$sabic_individual)))
  expect_output(print(candidates), "Class enumeration")
  fit <- fit_multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 6, seed = 5)
  information <- inference_multilpa(fit, dat)
  expect_s3_class(information, "multilpa_inference")
  estimates <- as.data.frame(information)
  expect_identical(nrow(estimates), length(coef(fit)))
  expect_equal(estimates$estimate, unname(coef(fit)))
  expect_equal(estimates$conf_high - estimates$conf_low,
               2 * stats::qnorm(0.975) * estimates$standard_error)
  # Mixing probabilities are boundary hypotheses and carry no Wald test.
  expect_true(all(is.na(subset(estimates,
    grepl("probability", estimates$parameter))$p_value)))
  expect_false(any(is.na(subset(estimates,
    grepl("^mean", estimates$parameter))$p_value)))
  expect_output(print(information), "Multilevel LPA inference")
})

test_that("multilpa_start round-trips a fitted solution", {
  dat <- make_two_level()
  fit <- fit_multilpa(dat, c("a", "b"), "g", 2, 2, n_starts = 6, seed = 5)
  start <- multilpa_start(fit)
  expect_setequal(names(start), c("means", "variances",
    "profile_probabilities", "group_probabilities"))
  expect_null(dimnames(start$means))
  refit <- fit_multilpa(dat, c("a", "b"), "g", 2, 2, start = start, n_starts = 1,
                      max_iter = 5000, tol = 1e-13)
  # The fitted solution must be a fixed point of the EM engine.
  expect_equal(refit$log_likelihood, fit$log_likelihood, tolerance = 1e-9)
  expect_equal(unname(refit$means), unname(fit$means), tolerance = 1e-6)
  full <- fit_multilpa(dat, c("a", "b"), "g", 2, 2, covariance_model = "full",
                     n_starts = 4, seed = 5)
  full_start <- multilpa_start(full)
  expect_true("covariances" %in% names(full_start))
  expect_null(full_start$variances)
  expect_false("covariances" %in% names(multilpa_start(full, covariance = "drop")))
  expect_equal(dim(multilpa_start(full, covariance = "drop")$variances), c(2L, 2L))
})

test_that("multilpa_start rejects incomplete input by condition class", {
  expect_error(multilpa_start(list(means = matrix(0, 2, 2))),
               class = "multilpa_bad_start")
  expect_error(multilpa_start(list(means = matrix(0, 2, 2),
    profile_probabilities = matrix(0.5, 1, 2), group_probabilities = 1)),
    class = "multilpa_bad_start")
  expect_error(multilpa_start(list(means = matrix(0, 2, 2),
    variances = matrix(1, 2, 2), profile_probabilities = matrix(0.5, 1, 2),
    group_probabilities = 1), covariance = "keep"),
    class = "multilpa_bad_start")
})

test_that("the Lo-Mendell-Rubin adjustment reproduces genuine Mplus TECH11", {
  reference <- readRDS(test_path("..", "fixtures", "mplus", "lmr-tech11.rds"))
  expect_identical(nrow(reference), 2L)
  adjusted <- reference$statistic / (1 + 1 / (reference$df * log(reference$n)))
  # Mplus prints TECH11 to three decimals; that is the comparison precision.
  expect_lt(max(abs(adjusted - reference$mplus_adjusted)), 5e-3)
})

test_that("lmr_lrt_multilpa reports a statistic and withholds a p-value", {
  dat <- make_two_level()
  smaller <- fit_multilpa(dat, c("a", "b"), "g", 1, 1, n_starts = 2, seed = 5)
  larger <- fit_multilpa(dat, c("a", "b"), "g", 2, 1, n_starts = 5, seed = 5)
  result <- lmr_lrt_multilpa(smaller, larger)
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
  grouped <- lmr_lrt_multilpa(smaller, larger, n = "groups")
  expect_equal(grouped$n, smaller$n_groups)
  expect_equal(grouped$statistic, result$statistic)
  expect_gt(grouped$adjustment_factor, result$adjustment_factor)
})

test_that("lmr_lrt_multilpa rejects invalid comparisons by condition class", {
  dat <- make_two_level()
  smaller <- fit_multilpa(dat, c("a", "b"), "g", 1, 1, n_starts = 2, seed = 5)
  larger <- fit_multilpa(dat, c("a", "b"), "g", 2, 1, n_starts = 5, seed = 5)
  expect_error(lmr_lrt_multilpa(larger, smaller), class = "multilpa_bad_nesting")
  other <- fit_multilpa(make_two_level(seed = 12L, n_groups = 20L),
                      c("a", "b"), "g", 2, 1, n_starts = 3, seed = 5)
  expect_error(lmr_lrt_multilpa(smaller, other),
               class = "multilpa_incomparable_models")
  expect_error(lmr_lrt_multilpa(smaller, larger, n = 1), "greater than one")
  reversed <- larger
  reversed$log_likelihood <- smaller$log_likelihood - 1
  expect_error(lmr_lrt_multilpa(smaller, reversed),
               class = "multilpa_reversed_likelihood")
})

test_that("every result class has a working tidy accessor", {
  # Rule 0: no result object should make the caller reach in with `$`.
  set.seed(4)
  group <- rep(seq_len(30L), each = 10L)
  n <- length(group)
  profile <- sample.int(2L, n, replace = TRUE)
  dat <- data.frame(school = group,
                    a = stats::rnorm(n, c(-1.5, 1.5)[profile]),
                    b = stats::rnorm(n, c(-1, 1)[profile]),
                    age = stats::rnorm(n),
                    resources = rep(stats::rnorm(30L), each = 10L))
  covariate_fit <- suppressWarnings(fit_multilpa_covariates(
    dat, c("a", "b"), "school", 2, 2, profile_covariates = "age",
    group_covariates = "resources", n_starts = 4, seed = 1))
  profiles <- as.data.frame(covariate_fit)
  expect_identical(names(profiles),
    c("profile", "indicator", "mean", "variance", "standard_deviation"))
  expect_identical(nrow(profiles), 4L)
  expect_equal(profiles$standard_deviation, sqrt(profiles$variance))
  coefficients <- as.data.frame(covariate_fit, what = "coefficients")
  expect_identical(names(coefficients),
    c("level", "outcome", "term", "estimate"))
  expect_setequal(unique(coefficients$level), c("profile", "group"))
  expect_true("age" %in% coefficients$term)
  expect_true("resources" %in% coefficients$term)
  expect_identical(nrow(as.data.frame(covariate_fit, what = "posteriors")), n)
  expect_identical(nrow(as.data.frame(covariate_fit, what = "group_posteriors")),
                   30L)
  intercept_fit <- fit_multilpa_random_intercept(dat, c("a", "b"), "school", 2,
                                               n_starts = 3, seed = 1)
  expect_identical(nrow(as.data.frame(intercept_fit)), 4L)
  intercepts <- as.data.frame(intercept_fit, what = "random_intercepts")
  expect_identical(names(intercepts), c("group", "group_size", "mean", "sd"))
  expect_identical(nrow(intercepts), 30L)
  expect_equal(sum(intercepts$group_size), n)
  expect_true(all(intercepts$sd > 0))
  expect_identical(nrow(as.data.frame(intercept_fit, what = "posteriors")), n)
})

test_that("preparation helpers enforce their own contracts", {
  # Extracted from fit_multilpa; they must still reject what it used to reject.
  expect_error(.multilpa_prepare_groups(c(1, NA, 2)), "nonmissing")
  expect_error(.multilpa_prepare_groups(c(1, Inf, 2)), "nonmissing")
  expect_error(.multilpa_prepare_groups(list(1, 2)), "nonmissing")
  groups <- .multilpa_prepare_groups(c("b", "a", "b", "c", "a"))
  # First-occurrence order, not sorted order.
  expect_identical(groups$values, c("b", "a", "c"))
  expect_identical(groups$index, c(1L, 2L, 1L, 3L, 2L))
  expect_identical(groups$n, 3L)
  expect_identical(groups$sizes, c(2L, 2L, 1L))
  dat <- data.frame(g = rep(1:4, each = 3), y = stats::rnorm(12),
                    u = rep(c(1, 2), 6))
  prepared <- .multilpa_prepare_indicators(dat, c("y", "u"), "u", "error", 1e-10)
  expect_identical(prepared$continuous, "y")
  expect_identical(dim(prepared$x), c(12L, 1L))
  expect_identical(unname(prepared$n_categories), 2L)
  # Every indicator categorical leaves a zero-column Gaussian block.
  all_categorical <- .multilpa_prepare_indicators(dat, "u", "u", "error", 1e-10)
  expect_identical(dim(all_categorical$x), c(12L, 0L))
  expect_error(.multilpa_check_arguments(dat, "y", "g", 0, 1, 1, 1, 1e-8, 1e-6,
                                       1e-10, NULL, character()),
               "positive integer")
  expect_error(.multilpa_check_arguments(dat, "y", "g", 2, 1, 1, 1, -1, 1e-6,
                                       1e-10, NULL, character()),
               "finite positive")
  expect_error(.multilpa_check_arguments(dat, "y", "g", 2, 1, 1, 1, 1e-8, 1e-6,
                                       1e-10, -1, character()),
               "nonnegative integer")
})

test_that("the null-default operator behaves like base %||%", {
  expect_identical(NULL %||% "fallback", "fallback")
  expect_identical("value" %||% "fallback", "value")
  expect_identical(NA %||% "fallback", NA)
  # A zero-length vector is not NULL and must pass through unchanged.
  expect_identical(character() %||% "fallback", character())
})

test_that("extracted covariate and inference helpers keep their contracts", {
  set.seed(9)
  group <- rep(seq_len(24L), each = 10L)
  n <- length(group)
  profile <- sample.int(2L, n, replace = TRUE)
  dat <- data.frame(school = group,
                    a = stats::rnorm(n, c(-1.5, 1.5)[profile]),
                    b = stats::rnorm(n, c(-1, 1)[profile]),
                    age = stats::rnorm(n),
                    resources = rep(stats::rnorm(24L), each = 10L))
  # Covariate column contracts, previously inline in fit_multilpa_covariates().
  expect_error(.multilpa_cov_check_covariates(dat, "a", character(),
                                            c("a", "b"), "school"),
               "distinct from indicators")
  bad <- dat
  bad$age[3L] <- NA
  expect_error(.multilpa_cov_check_covariates(bad, "age", character(),
                                            c("a", "b"), "school"),
               "finite numeric")
  expect_null(.multilpa_cov_check_covariates(dat, "age", "resources",
                                           c("a", "b"), "school"))
  first_rows <- match(seq_len(24L), group)
  designs <- .multilpa_cov_designs(dat, c("a", "b"), "age", "resources",
                                 first_rows, 2L)
  expect_identical(dim(designs$x), c(n, 2L))
  # Centring is a translation: the centred matrix has zero column means.
  expect_equal(unname(colMeans(designs$x)), c(0, 0))
  expect_identical(nrow(designs$w), 24L)
  expect_identical(colnames(designs$w), c("(Intercept)", "resources"))
  expect_identical(nrow(designs$stacked_design), n * 2L)
  collinear <- dat
  collinear$copy <- collinear$age
  expect_error(.multilpa_cov_designs(collinear, c("a", "b"), c("age", "copy"),
                                   character(), first_rows, 2L),
               "rank deficient")

  # Inference contracts, previously inline in inference_multilpa().
  fit <- fit_multilpa(dat, c("a", "b"), "school", 2, 2, n_starts = 5, seed = 1)
  prepared <- .multilpa_inference_matrix(fit, dat)
  expect_identical(dim(prepared$x), c(n, 2L))
  expect_equal(unname(colMeans(prepared$x)), c(0, 0))
  expect_equal(prepared$centers, colMeans(as.matrix(dat[, c("a", "b")])))
  expect_error(.multilpa_inference_matrix(fit, dat[-1L, ]), "original numeric")
  shuffled <- dat[c(2L, 1L, seq(3L, n)), ]
  expect_error(.multilpa_inference_matrix(fit, shuffled), "row order")
  expect_null(.multilpa_check_regularity(fit, "observed"))
  unconverged <- fit
  unconverged$converged <- FALSE
  expect_error(.multilpa_check_regularity(unconverged, "observed"), "converged fit")
  bounded <- fit
  bounded$boundary <- TRUE
  expect_error(.multilpa_check_regularity(bounded, "observed"), "bound-active")
  categorical <- fit
  categorical$response_probabilities <- list(matrix(0.5, 2L, 2L))
  expect_error(.multilpa_check_regularity(categorical, "observed"),
               class = "multilpa_unsupported_inference")

  # A singular information matrix must be refused, not inverted.
  flat <- function(displacement) 0
  flat_gradient <- function(displacement) rep(0, length(displacement))
  expect_error(.multilpa_observed_hessian(flat, flat_gradient, c(1, 1), 1e-4),
               "positive definite")
})
