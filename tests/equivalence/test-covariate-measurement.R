# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

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

# Log multinomial probabilities, the last category the reference: base R.
.cov_log_softmax <- function(design, coefficients) {
  scores <- cbind(design %*% coefficients, 0)
  largest <- apply(scores, 1L, max)
  scores - largest - log(rowSums(exp(scores - largest)))
}

# The two-level covariate likelihood, from the definition.
.cov_independent_likelihood <- function(fit, data) {
  continuous <- fit$continuous %||% fit$vars
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
  group_prior <- exp(.cov_log_softmax(fit$group_design, fit$group_coefficients))
  sum(vapply(seq_len(fit$n_groups), function(group) {
    rows <- which(fit$group_index == group)
    log(sum(vapply(seq_len(fit$n_group_classes), function(class) {
      prior <- exp(.cov_log_softmax(fit$profile_design[[class]],
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
  fits <- lapply(cases, function(case) {
    quietly(do.call(multilpa, c(list(
      data, case$vars, "g", 2L, 2L, profile_covariates = "z",
      n_starts = 1, seed = 3, max_iter = 50), case$extra)))
  })
  labels <- vapply(cases, \(case) case$label, character(1))
  names(fits) <- labels
  expect_false(any(vapply(fits, \(fit) is.unsorted(fit$log_likelihood_history), logical(1))))
  # The same sum of logs computed twice; only the order of the additions differs.
  expect_reference(vapply(fits, \(fit) as.numeric(logLik(fit)), numeric(1)),
                   vapply(fits, \(fit) .cov_independent_likelihood(fit, data), numeric(1)),
                   reference = "independent: two-level covariate mixture likelihood, base R",
                   covers = c("multilpa()", "logLik.multilpa_covariates()"),
                   quantity = "log likelihood by measurement model", tolerance = 1e-10,
                   measure = "relative")
})

test_that("full-covariance standard errors match an independent information matrix", {
  data <- .cov_measurement_fixture(n_groups = 24L)
  fit <- quietly(multilpa(data, c("y1", "y2"), "g", 2L, 1L,
    profile_covariates = "z", n_starts = 1, seed = 1, covariance_model = "full",
    max_iter = 2000, tol = 1e-11))
  theta <- latents:::.multilpa_cov_encode(fit)
  x <- sweep(as.matrix(data[c("y1", "y2")]), 2L, fit$center, "-")
  negative <- function(values) {
    pieces <- latents:::.multilpa_cov_decode(values, fit)
    -latents:::.multilpa_cov_expectation(x, fit$group_index, pieces$parameters,
                                         fit$profile_design, fit$group_design,
                                         pieces$beta, pieces$gamma)$log_likelihood
  }
  independent <- solve(stats::optimHess(theta, negative))
  # `theta` is the estimation-scale encoding, so the comparison must ask for
  # that scale. `vcov()` now defaults to natural units, agreeing with
  # `confint()` and `parameter_inference()` on the same fit.
  package <- vcov(fit, data, scale = "unconstrained")
  # 1e-5 relative, the original tolerance: optimHess() without a gradient
  # takes second finite differences of the likelihood.
  expect_reference(sqrt(diag(package)), unname(sqrt(diag(independent))), reference = "stats",
                   covers = c("multilpa()", "vcov.multilpa_covariates()"),
                   quantity = "standard errors, full covariance, unconstrained scale",
                   tolerance = 1e-5, measure = "relative",
                   note = "Reference: the inverse of stats::optimHess() of the package's internal covariate likelihood (.multilpa_cov_expectation) at the estimate.")
  expect_true(isTRUE(all.equal(package, t(package))))
  expect_gt(min(eigen(package, symmetric = TRUE, only.values = TRUE)$values), 0)
})
