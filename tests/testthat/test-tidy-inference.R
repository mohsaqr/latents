# One naming scheme, one column set, one scale convention, across every class
# that reports parameters. These are contract tests: they are here so that a
# caller can write code against `coef()` or `parameter_inference()` once and
# have it keep working when the fit changes family.

.tidy_inference_columns <- c("level", "outcome", "term", "parameter",
                             "estimate", "standard_error", "statistic",
                             "p_value", "p_adjusted", "conf_low", "conf_high")

.tidy_plain_fit <- function() {
  set.seed(1)
  data <- data.frame(school = rep(seq_len(12), each = 8),
                     reading = stats::rnorm(96), maths = stats::rnorm(96))
  list(data = data,
       fit = multilpa(data, c("reading", "maths"), "school", n_profiles = 2,
                      n_group_classes = 2, n_starts = 2, seed = 1))
}

.tidy_covariate_fit <- function() {
  set.seed(5)
  school <- rep(seq_len(16), each = 8)
  high_class <- rep(rep(c(FALSE, TRUE), length.out = 16), each = 8)
  x <- stats::rnorm(128)
  profile <- ifelse(stats::runif(128) <
                      stats::plogis(-1 + 2 * high_class + 0.8 * x), 2L, 1L)
  data <- data.frame(
    school = school, x = x,
    y1 = stats::rnorm(128, ifelse(profile == 2L, 2, -2), 0.7),
    y2 = stats::rnorm(128, ifelse(profile == 2L, 1.5, -1.5), 0.7))
  list(data = data,
       fit = fit_covariates(data, c("y1", "y2"), "school", n_profiles = 2,
                            n_group_classes = 2, profile_covariates = "x",
                            n_starts = 2, seed = 1))
}

test_that("a coefficient name is its tidy decomposition, in every class", {
  plain <- .tidy_plain_fit()
  covariate <- .tidy_covariate_fit()
  expected <- function(fit, data) {
    inference <- parameter_inference(fit, data)
    paste(inference$level, inference$parameter, inference$outcome,
          inference$term, sep = ".")
  }
  # level.parameter.outcome.term, with a missing term simply left off.
  expect_identical(names(coef(plain$fit)),
                   sub("[.]NA$", "", expected(plain$fit, plain$data)))
  expect_identical(names(coef(covariate$fit)),
                   expected(covariate$fit, covariate$data))
})

test_that("the four name parts read back into the tidy columns", {
  plain <- .tidy_plain_fit()
  inference <- parameter_inference(plain$fit, plain$data)
  parts <- strsplit(sub("^([^.]+)[.]([^.]+)[.]([^.]+)[.]?", "\\1\r\\2\r\\3\r",
                        names(coef(plain$fit))), "\r", fixed = TRUE)

  expect_identical(vapply(parts, `[[`, character(1), 1L), inference$level)
  expect_identical(vapply(parts, `[[`, character(1), 2L), inference$parameter)
  expect_identical(vapply(parts, `[[`, character(1), 3L), inference$outcome)
  term <- vapply(parts, function(p) if (length(p) >= 4L) p[[4L]] else NA_character_,
                 character(1))
  expect_identical(term, inference$term)
})

test_that("the unconstrained coefficients carry names of the same shape", {
  plain <- .tidy_plain_fit()
  estimation <- coef(plain$fit, scale = "unconstrained")

  expect_length(estimation, plain$fit$n_parameters)
  expect_true(all(grepl("^(measurement|profile|group)[.]", names(estimation))))
  # The kind names the scale it is on, so a log variance cannot be mistaken
  # for a variance by a caller reading the name.
  labels <- multilpa:::.multilpa_coefficient_labels(plain$fit, "unconstrained")
  expect_setequal(unique(labels$parameter), c("mean", "log_variance", "logit"))
  expect_identical(names(estimation), multilpa:::.multilpa_parameter_names(labels))
})

test_that("vcov and confint are named exactly as coef is", {
  plain <- .tidy_plain_fit()
  covariate <- .tidy_covariate_fit()

  expect_identical(rownames(vcov(plain$fit, data = plain$data)),
                   names(coef(plain$fit)))
  expect_identical(colnames(vcov(plain$fit, data = plain$data)),
                   names(coef(plain$fit)))
  expect_identical(rownames(confint(plain$fit, data = plain$data)),
                   names(coef(plain$fit)))
  expect_identical(rownames(vcov(covariate$fit, covariate$data)),
                   names(coef(covariate$fit)))
  expect_identical(rownames(confint(covariate$fit, data = covariate$data)),
                   names(coef(covariate$fit)))
  expect_identical(
    rownames(vcov(plain$fit, data = plain$data, scale = "unconstrained")),
    names(coef(plain$fit, scale = "unconstrained")))
})

test_that("one generic returns one column set, whichever class it dispatches on", {
  plain <- .tidy_plain_fit()
  covariate <- .tidy_covariate_fit()
  from_plain <- parameter_inference(plain$fit, plain$data)
  from_covariate <- parameter_inference(covariate$fit, covariate$data)

  expect_named(from_plain, .tidy_inference_columns)
  expect_named(from_covariate, .tidy_inference_columns)
  expect_identical(vapply(from_plain, class, character(1)),
                   vapply(from_covariate, class, character(1)))
  expect_s3_class(from_plain, "data.frame")
  expect_s3_class(from_covariate, "data.frame")
  expect_equal(nrow(from_covariate), covariate$fit$n_parameters)
  # Both carry the level vocabulary and record the same diagnostics.
  expect_true(all(from_plain$level %in% c("measurement", "profile", "group")))
  expect_true(all(from_covariate$level %in% c("measurement", "profile", "group")))
  expect_identical(attr(from_plain, "vcov_type"), "observed")
  expect_identical(attr(from_covariate, "vcov_type"), "observed")
  expect_identical(attr(from_plain, "adjust"), "none")
  expect_identical(attr(from_covariate, "adjust"), "none")
  expect_equal(attr(from_plain, "level"), 0.95)
  expect_equal(attr(from_covariate, "level"), 0.95)
})

test_that("a correction is named, never applied by stealth", {
  plain <- .tidy_plain_fit()
  none <- parameter_inference(plain$fit, plain$data)
  corrected <- parameter_inference(plain$fit, plain$data, adjust = "bonferroni")
  tested <- !is.na(none$p_value)

  expect_equal(none$p_adjusted, none$p_value)
  expect_identical(attr(corrected, "adjust"), "bonferroni")
  expect_equal(corrected$p_adjusted[tested],
               pmin(1, none$p_value[tested] * sum(tested)))
  expect_true(all(is.na(corrected$p_adjusted[!tested])))
  expect_error(parameter_inference(plain$fit, plain$data, adjust = "nonesuch"))
})

test_that("the natural-scale covariance is what the reported intervals are built from", {
  plain <- .tidy_plain_fit()
  inference <- parameter_inference(plain$fit, plain$data)
  natural <- vcov(plain$fit, data = plain$data)
  estimation <- vcov(plain$fit, data = plain$data, scale = "unconstrained")

  expect_equal(inference$standard_error, unname(sqrt(pmax(diag(natural), 0))))
  expect_equal(inference$conf_low,
               inference$estimate - stats::qnorm(0.975) * inference$standard_error)
  # The natural matrix is singular because every set of probabilities sums to
  # one; the estimation-scale matrix drops those references and is not.
  expect_equal(nrow(natural), length(coef(plain$fit)))
  expect_equal(nrow(estimation), plain$fit$n_parameters)
  expect_lt(min(abs(eigen(natural, symmetric = TRUE, only.values = TRUE)$values)),
            1e-8)
  expect_gt(min(eigen(estimation, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("the delta method is invertible: the two covariate scales agree through the Jacobian", {
  covariate <- .tidy_covariate_fit()
  estimation <- vcov(covariate$fit, covariate$data, scale = "unconstrained")
  natural <- vcov(covariate$fit, covariate$data)
  theta <- multilpa:::.multilpa_cov_encode(covariate$fit)
  jacobian <- multilpa:::.multilpa_cov_natural_jacobian(
    covariate$fit, theta, multilpa:::.multilpa_cov_labels(covariate$fit))

  expect_equal(natural, jacobian %*% estimation %*% t(jacobian),
               ignore_attr = TRUE)
  # A numerical derivative of the estimation-to-natural map reproduces it, so
  # the Jacobian is the formula and not just a convenient matrix.
  labels <- multilpa:::.multilpa_cov_labels(covariate$fit)
  natural_of <- function(values) {
    multilpa:::.multilpa_cov_natural_estimate(covariate$fit, values, labels)
  }
  numerical <- vapply(seq_along(theta), function(index) {
    step <- numeric(length(theta))
    step[index] <- 1e-6
    (natural_of(theta + step) - natural_of(theta - step)) / 2e-6
  }, numeric(length(theta)))
  expect_equal(jacobian, numerical, tolerance = 1e-6)
})

test_that("coef reports the same numbers the tidy table does, on both spread models", {
  # Regression: a full-covariance covariate fit used to return raw log-Cholesky
  # coordinates under names that said `covariance`, so `coef()` reported a
  # negative residual variance while parameter_inference() reported the real one.
  set.seed(5)
  school <- rep(seq_len(20), each = 8)
  x <- stats::rnorm(160)
  profile <- ifelse(stats::runif(160) < stats::plogis(-0.5 + 0.8 * x), 2L, 1L)
  shared <- stats::rnorm(160)
  data <- data.frame(
    school = school, x = x,
    y1 = stats::rnorm(160, ifelse(profile == 2L, 2, -2), 0.7) + shared,
    y2 = stats::rnorm(160, ifelse(profile == 2L, 1.5, -1.5), 0.7) + 0.6 * shared)
  invisible(lapply(c("diagonal", "full"), function(covariance_model) {
    fit <- fit_covariates(data, c("y1", "y2"), "school", n_profiles = 2,
                          n_group_classes = 1, profile_covariates = "x",
                          covariance_model = covariance_model,
                          n_starts = 3, seed = 1)
    inference <- parameter_inference(fit, data)
    expect_equal(unname(coef(fit)), inference$estimate, info = covariance_model)
    expect_identical(names(coef(fit)),
                     paste(inference$level, inference$parameter,
                           inference$outcome, inference$term, sep = "."),
                     info = covariance_model)
    # Every reported residual variance is positive, which is the property the
    # raw coordinates did not have.
    spread <- inference$parameter %in% c("variance", "covariance")
    diagonal_entries <- spread & inference$term %in%
      paste(fit$vars, fit$vars, sep = ":") | inference$parameter == "variance"
    expect_true(all(coef(fit)[diagonal_entries] > 0), info = covariance_model)
    # The estimation scale says so in the name and never claims to be natural.
    estimation <- coef(fit, scale = "unconstrained")
    expect_length(estimation, length(coef(fit)))
    expect_false(any(grepl("[.]variance[.]|[.]covariance[.]", names(estimation))),
                 info = covariance_model)
    expect_identical(names(estimation),
                     rownames(vcov(fit, data, scale = "unconstrained")),
                     info = covariance_model)
  }))
})

test_that("a fit with no standard errors says so by class rather than by number", {
  set.seed(7)
  data <- data.frame(person = rep(seq_len(30), each = 5),
                     wave = rep(seq_len(5), times = 30))
  data$score_a <- stats::rnorm(150)
  data$score_b <- stats::rnorm(150)
  fit <- fit_transitions(data, c("score_a", "score_b"), "person",
                         n_profiles = 2, time = "wave", n_starts = 2, seed = 1)

  expect_error(parameter_inference(fit, data), class = "multilpa_no_inference")
  expect_error(vcov(fit), class = "multilpa_no_inference")
})
