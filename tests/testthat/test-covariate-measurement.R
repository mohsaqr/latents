# Covariate fits with categorical, mixed and full-covariance measurement.
#
# The likelihood is checked against the definition written out independently,
# and the analytic score against numerical differentiation of that likelihood,
# so neither is compared with the machinery under test.

.cov_measurement_fixture <- function(n_groups = 20L, size = 8L, seed = 6L) {
  set.seed(seed)
  frame <- data.frame(g = rep(seq_len(n_groups), each = size),
                      z = stats::rnorm(n_groups * size))
  state <- 1L + (stats::runif(nrow(frame)) > stats::plogis(frame$z))
  frame$y1 <- stats::rnorm(nrow(frame), c(-2, 2)[state], 0.9)
  frame$y2 <- stats::rnorm(nrow(frame), c(-1.5, 1.5)[state], 0.9)
  frame$q <- ifelse(stats::runif(nrow(frame)) < c(0.8, 0.2)[state], "no", "yes")
  frame$r <- sample(1:3, nrow(frame), replace = TRUE)
  frame
}

# The two-level covariate likelihood, from the definition.
.cov_independent_likelihood <- function(fit, data) {
  continuous <- .multilpa_continuous_names(fit)
  x <- if (length(continuous)) as.matrix(data[continuous]) else
    matrix(0, nrow(data), 0L)
  density <- vapply(seq_len(fit$n_profiles), function(profile) {
    gaussian <- if (ncol(x) == 0L) numeric(nrow(data)) else if (is.null(fit$covariances)) {
      rowSums(vapply(seq_along(continuous), function(j) {
        stats::dnorm(x[, j], fit$means[profile, j],
                     sqrt(fit$variances[profile, j]), log = TRUE)
      }, numeric(nrow(data))))
    } else {
      covariance <- matrix(fit$covariances[, , profile], ncol(x), ncol(x))
      residual <- sweep(x, 2L, fit$means[profile, ], "-")
      factor <- chol(covariance)
      -0.5 * (ncol(x) * log(2 * pi) + 2 * sum(log(diag(factor))) +
                rowSums((residual %*% backsolve(factor, diag(ncol(x))))^2))
    }
    categorical <- if (is.null(fit$categorical_data)) numeric(nrow(data)) else
      rowSums(vapply(seq_along(fit$categorical), function(j) {
        log(fit$response_probabilities[[j]][profile, fit$categorical_data[, j]])
      }, numeric(nrow(data))))
    gaussian + categorical
  }, numeric(nrow(data)))
  group_prior <- exp(.multilpa_log_softmax(fit$group_design, fit$group_coefficients))
  sum(vapply(seq_len(fit$n_groups), function(group) {
    rows <- which(fit$group_index == group)
    log(sum(vapply(seq_len(fit$n_group_classes), function(class) {
      prior <- exp(.multilpa_log_softmax(fit$profile_design[[class]],
                                         fit$profile_coefficients))
      group_prior[group, class] *
        exp(sum(log(rowSums(exp(density[rows, , drop = FALSE]) *
                              prior[rows, , drop = FALSE]))))
    }, numeric(1))))
  }, numeric(1)))
}

test_that("every measurement model reproduces the likelihood written from its definition", {
  data <- .cov_measurement_fixture()
  cases <- list(
    list(label = "gaussian diagonal", vars = c("y1", "y2"), extra = list()),
    list(label = "gaussian full", vars = c("y1", "y2"),
         extra = list(covariance_model = "full")),
    list(label = "mixed", vars = c("y1", "y2", "q", "r"),
         extra = list(categorical = c("q", "r"))),
    list(label = "mixed, full covariance", vars = c("y1", "y2", "q", "r"),
         extra = list(categorical = c("q", "r"), covariance_model = "full")),
    list(label = "categorical only", vars = c("q", "r"),
         extra = list(categorical = c("q", "r"))))
  invisible(lapply(cases, function(case) {
    fit <- quietly(do.call(multilpa, c(list(
      data, case$vars, "g", 2L, 2L, profile_covariates = "z",
      n_starts = 2, seed = 3, max_iter = 50), case$extra)))
    expect_equal(fit$log_likelihood, .cov_independent_likelihood(fit, data),
                 tolerance = 1e-10, label = case$label)
    expect_false(is.unsorted(fit$log_likelihood_history))
  }))
})

test_that("the measurement blocks and parameter counts match the model fitted", {
  data <- .cov_measurement_fixture()
  gaussian <- quietly(multilpa(data, c("y1", "y2"), "g", 2L, 1L,
    profile_covariates = "z", n_starts = 2, seed = 3, max_iter = 50))
  full <- quietly(multilpa(data, c("y1", "y2"), "g", 2L, 1L,
    profile_covariates = "z", n_starts = 2, seed = 3, max_iter = 50, covariance_model = "full"))
  mixed <- quietly(multilpa(data, c("y1", "y2", "q", "r"), "g", 2L, 1L,
    profile_covariates = "z", n_starts = 2, seed = 3, max_iter = 50, categorical = c("q", "r")))
  categorical <- quietly(multilpa(data, c("q", "r"), "g", 2L, 1L,
    profile_covariates = "z", n_starts = 2, seed = 3, max_iter = 50, categorical = c("q", "r")))

  expect_identical(gaussian$measurement_model, "gaussian")
  expect_identical(mixed$measurement_model, "mixed")
  expect_identical(categorical$measurement_model, "categorical")
  expect_null(gaussian$covariances)
  expect_identical(dim(full$covariances), c(2L, 2L, 2L))
  expect_named(mixed$response_probabilities, c("q", "r"))
  expect_equal(unname(rowSums(mixed$response_probabilities[["r"]])), rep(1, 2))

  # 4 means + 4 variances + 2 logits; the full model swaps 4 variances for 6
  # covariance entries; the mixed model adds 2 * ((2-1) + (3-1)) responses.
  expect_equal(gaussian$n_parameters, 10)
  expect_equal(full$n_parameters, 12)
  expect_equal(mixed$n_parameters, 16)
  expect_equal(categorical$n_parameters, 8)
})

test_that("with no covariates the covariate model is the covariate-free one", {
  data <- .cov_measurement_fixture()
  # Intercept-only multinomial logits are a reparameterization of free mixing
  # weights, so the two likelihoods are the same function.
  cases <- list(list(vars = c("y1", "y2"), extra = list()),
                list(vars = c("y1", "y2"),
                     extra = list(covariance_model = "full")),
                list(vars = c("y1", "y2", "q"),
                     extra = list(categorical = "q")))
  invisible(lapply(cases, function(case) {
    plain <- quietly(do.call(multilpa, c(list(
      data, case$vars, "g", 2L, 2L, n_starts = 1, seed = 4,
      max_iter = 3000, tol = 1e-12), case$extra)))
    # The covariate estimator with no covariates: `multilpa()` cannot be
    # asked for it, because naming no covariate requests the covariate-free
    # model, which is exactly what this compares it against.
    covariate <- quietly(do.call(multilpa:::.multilpa_fit_covariates, c(list(
      data, case$vars, "g", 2L, 2L, n_starts = 1, seed = 4,
      max_iter = 3000, tol = 1e-12), case$extra)))
    expect_equal(plain$log_likelihood, covariate$log_likelihood, tolerance = 1e-6)
    expect_equal(plain$n_parameters, covariate$n_parameters)
  }))
})

test_that("the full-covariance score is the derivative of the likelihood", {
  data <- .cov_measurement_fixture(n_groups = 24L)
  fit <- quietly(multilpa(data, c("y1", "y2"), "g", 2L, 1L,
    profile_covariates = "z", n_starts = 3, seed = 1, covariance_model = "full",
    max_iter = 2000, tol = 1e-11))
  theta <- .multilpa_cov_encode(fit)
  expect_equal(length(theta), fit$n_parameters)
  x <- sweep(as.matrix(data[c("y1", "y2")]), 2L, fit$center, "-")
  log_likelihood <- function(values) {
    pieces <- .multilpa_cov_decode(values, fit)
    .multilpa_cov_expectation(x, fit$group_index, pieces$parameters,
                              fit$profile_design, fit$group_design,
                              pieces$beta, pieces$gamma)$log_likelihood
  }
  # At the maximum the score is zero, and comparing two near-zero vectors with
  # a relative tolerance fails on floating-point noise that differs between
  # platforms. A fixed perturbation gives a clearly non-zero gradient.
  probe <- theta + 0.05 * rep_len(c(1, -1), length(theta))
  analytic <- colSums(.multilpa_cov_group_scores(probe, x, fit))
  step <- 1e-5
  numerical <- vapply(seq_along(probe), function(index) {
    up <- probe; down <- probe
    up[index] <- up[index] + step
    down[index] <- down[index] - step
    (log_likelihood(up) - log_likelihood(down)) / (2 * step)
  }, numeric(1))
  expect_equal(analytic, numerical, tolerance = 1e-5, ignore_attr = TRUE)
  # The decoded parameters must round-trip through the log-Cholesky coordinates.
  expect_equal(.multilpa_cov_decode(theta, fit)$parameters$covariances,
               unname(fit$covariances), tolerance = 1e-10)
})

test_that("full-covariance standard errors match an independent information matrix", {
  data <- .cov_measurement_fixture(n_groups = 24L)
  fit <- quietly(multilpa(data, c("y1", "y2"), "g", 2L, 1L,
    profile_covariates = "z", n_starts = 3, seed = 1, covariance_model = "full",
    max_iter = 2000, tol = 1e-11))
  theta <- .multilpa_cov_encode(fit)
  x <- sweep(as.matrix(data[c("y1", "y2")]), 2L, fit$center, "-")
  negative <- function(values) {
    pieces <- .multilpa_cov_decode(values, fit)
    -.multilpa_cov_expectation(x, fit$group_index, pieces$parameters,
                               fit$profile_design, fit$group_design,
                               pieces$beta, pieces$gamma)$log_likelihood
  }
  independent <- solve(stats::optimHess(theta, negative))
  # `theta` is the estimation-scale encoding, so the comparison must ask for
  # that scale. `vcov()` now defaults to natural units, agreeing with
  # `confint()` and `parameter_inference()` on the same fit.
  package <- vcov(fit, data, scale = "unconstrained")
  expect_equal(sqrt(diag(package)), sqrt(diag(independent)),
               tolerance = 1e-5, ignore_attr = TRUE)
  expect_true(isTRUE(all.equal(package, t(package))))
  expect_gt(min(eigen(package, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("covariance estimates are reported in natural units with their own tests", {
  data <- .cov_measurement_fixture(n_groups = 24L)
  fit <- quietly(multilpa(data, c("y1", "y2"), "g", 2L, 1L,
    profile_covariates = "z", n_starts = 3, seed = 1, covariance_model = "full",
    max_iter = 2000, tol = 1e-11))
  inference <- parameter_inference(fit, data)
  covariances <- inference[inference$parameter == "covariance", ]
  expect_identical(nrow(covariances), 6L)
  # The delta method must return exactly the covariances the fit holds.
  fitted <- unlist(lapply(seq_len(2L), function(profile) {
    block <- matrix(fit$covariances[, , profile], 2L, 2L)
    block[lower.tri(block, diag = TRUE)]
  }), use.names = FALSE)
  expect_equal(covariances$estimate, fitted, tolerance = 1e-10)
  # A diagonal entry is a variance, so it gets no Wald test against zero; an
  # off-diagonal covariance is a real hypothesis and keeps one.
  diagonal <- covariances$term %in% c("y1:y1", "y2:y2")
  expect_true(all(is.na(covariances$p_value[diagonal])))
  expect_false(any(is.na(covariances$p_value[!diagonal])))
  expect_identical(nrow(parameter_inference(fit, data, vcov_type = "robust")),
                   nrow(inference))
  expect_identical(nrow(confint(fit, data = data)), nrow(inference))
})

test_that("inference refuses a categorical covariate fit rather than miscounting it", {
  data <- .cov_measurement_fixture()
  fit <- quietly(multilpa(data, c("y1", "y2", "q"), "g", 2L, 1L,
    profile_covariates = "z", n_starts = 2, seed = 3, categorical = "q"))
  expect_error(parameter_inference(fit, data),
               class = "multilpa_unsupported_inference")
  expect_error(vcov(fit, data), class = "multilpa_unsupported_inference")
  # The fit itself is still usable; only the standard errors are withheld.
  expect_s3_class(get_results(fit, "posteriors"), "data.frame")
  expect_true(is.finite(fit$log_likelihood))
})

test_that("natural covariance errors agree with the covariate-free model's own", {
  data <- .cov_measurement_fixture()
  # Without covariates the two fits are the same model, and each reports the
  # covariance block in natural units through its own delta method. Agreement
  # is therefore a cross-implementation check on the transformation, which
  # comparing a fit with itself could never provide.
  plain <- quietly(multilpa(data, c("y1", "y2"), "g", 2L, 1L, n_starts = 1, seed = 4, max_iter = 4000, tol = 1e-12, covariance_model = "full"))
  covariate <- quietly(multilpa:::.multilpa_fit_covariates(data, c("y1", "y2"), "g", 2L, 1L,
    n_starts = 1, seed = 4, max_iter = 4000, tol = 1e-12,
    covariance_model = "full"))
  from_plain <- parameter_inference(plain, data)
  from_covariate <- parameter_inference(covariate, data)
  plain_rows <- from_plain[from_plain$parameter == "covariance", ]
  covariate_rows <- from_covariate[from_covariate$parameter == "covariance", ]
  expect_identical(nrow(plain_rows), nrow(covariate_rows))
  expect_equal(plain_rows$estimate, covariate_rows$estimate, tolerance = 1e-5)
  expect_equal(plain_rows$standard_error, covariate_rows$standard_error,
               tolerance = 1e-4)
})
