# Comparisons with independent references moved from tests/testthat/test-covariate-measurement.R,
# so that they run in CI's equivalence job rather than under R CMD check.

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
      n_starts = 1, seed = 3, max_iter = 50), case$extra)))
    expect_equal(fit$log_likelihood, .cov_independent_likelihood(fit, data),
                 tolerance = 1e-10, label = case$label)
    expect_false(is.unsorted(fit$log_likelihood_history))
  }))
})

test_that("full-covariance standard errors match an independent information matrix", {
  data <- .cov_measurement_fixture(n_groups = 24L)
  fit <- quietly(multilpa(data, c("y1", "y2"), "g", 2L, 1L,
    profile_covariates = "z", n_starts = 1, seed = 1, covariance_model = "full",
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
