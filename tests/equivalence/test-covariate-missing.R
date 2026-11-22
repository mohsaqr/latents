# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Comparisons with independent references moved from tests/testthat/test-covariate-missing.R,
# so that they run in CI's equivalence job rather than under R CMD check.

.cov_missing_fixture <- function(n_groups = 40L, size = 8L, seed = 3L,
                                 share = 0.15) {
  set.seed(seed)
  n <- n_groups * size
  frame <- data.frame(g = rep(seq_len(n_groups), each = size),
                      z = stats::rnorm(n))
  group_class <- rep(stats::rbinom(n_groups, 1L, 0.5), each = size)
  state <- 1L + stats::rbinom(n, 1L, stats::plogis(1.2 * frame$z +
                                                     c(-1.5, 1.5)[group_class + 1L]))
  frame$y1 <- stats::rnorm(n, c(0, 2)[state])
  frame$y2 <- stats::rnorm(n, c(0, 1.5)[state]) + 0.5 * frame$y1
  frame$q <- vapply(state, function(profile) {
    sample(c("a", "b", "c"), 1L,
           prob = list(c(0.6, 0.3, 0.1), c(0.1, 0.3, 0.6))[[profile]])
  }, character(1))
  blank <- function(value) replace(value, stats::runif(n) < share, NA)
  frame$y1 <- blank(frame$y1)
  frame$y2 <- blank(frame$y2)
  frame$q <- blank(frame$q)
  frame
}

# Log multinomial probabilities, the last category the reference: base R.
.cov_missing_log_softmax <- function(design, coefficients) {
  scores <- cbind(design %*% coefficients, 0)
  largest <- apply(scores, 1L, max)
  scores - largest - log(rowSums(exp(scores - largest)))
}

# The observed-data likelihood from its definition: each row's Gaussian
# density over its observed coordinates only, times the categorical
# probability of its observed category, mixed over profiles and group classes.
.cov_missing_likelihood <- function(fit, data) {
  continuous <- c("y1", "y2")
  x <- as.matrix(data[continuous])
  density <- vapply(seq_len(fit$n_profiles), function(profile) {
    covariance <- if (is.null(fit$covariances)) diag(fit$variances[profile, ]) else
      fit$covariances[, , profile]
    vapply(seq_len(nrow(x)), function(row) {
      observed <- which(!is.na(x[row, ]))
      gaussian <- if (length(observed) == 0L) 0 else {
        sigma <- covariance[observed, observed, drop = FALSE]
        residual <- x[row, observed] - fit$means[profile, observed]
        -0.5 * (length(observed) * log(2 * pi) +
                  as.numeric(determinant(sigma)$modulus) +
                  sum(residual * solve(sigma, residual)))
      }
      category <- data$q[row]
      gaussian + if (is.na(category)) 0 else
        log(fit$response_probabilities$q[profile, category])
    }, numeric(1))
  }, numeric(nrow(x)))
  group_prior <- exp(.cov_missing_log_softmax(fit$group_design, fit$group_coefficients))
  sum(vapply(seq_len(fit$n_groups), function(group) {
    rows <- which(fit$group_index == group)
    log(sum(vapply(seq_len(fit$n_group_classes), function(class) {
      prior <- exp(.cov_missing_log_softmax(fit$profile_design[[class]],
                                            fit$profile_coefficients))
      group_prior[group, class] *
        exp(sum(log(rowSums(exp(density[rows, , drop = FALSE]) *
                              prior[rows, , drop = FALSE]))))
    }, numeric(1))))
  }, numeric(1)))
}

.cov_missing_fit <- function(data, covariance_model, n_group_classes = 2L,
                             missing = "fiml", tol = 1e-12) {
  quietly(multilpa(data, c("y1", "y2", "q"), "g", 2L, n_group_classes,
                   categorical = "q", missing = missing,
                   covariance_model = covariance_model,
                   profile_covariates = "z", n_starts = 2, seed = 1,
                   tol = tol, max_iter = 5000))
}

# Analytic scores and standard errors of one fit beside numerical ones, one
# row per unconstrained coordinate (named as vcov() names them).
.cov_missing_derivatives <- function(fit, data) {
  x <- sweep(as.matrix(data[c("y1", "y2")]), 2L, fit$center, "-")
  codes <- fit$categorical_data
  theta <- latents:::.multilpa_cov_encode(fit)
  negative <- function(value) {
    pieces <- latents:::.multilpa_cov_decode(value, fit)
    -latents:::.multilpa_cov_expectation(x, fit$group_index, pieces$parameters,
                                         fit$profile_design, fit$group_design,
                                         pieces$beta, pieces$gamma, codes)$log_likelihood
  }
  # Off the maximum, where the score is not zero and so can be wrong.
  set.seed(9)
  away <- theta + 0.05 * stats::rnorm(length(theta))
  analytic <- colSums(latents:::.multilpa_cov_group_scores(away, x, fit, codes))
  numerical <- vapply(seq_along(away), function(index) {
    step <- replace(numeric(length(away)), index, 1e-5)
    -(negative(away + step) - negative(away - step)) / 2e-5
  }, numeric(1))
  # The package's information is the Jacobian of the analytic score; the
  # reference is the Hessian of the likelihood alone.
  reported <- sqrt(diag(vcov(fit, scale = "unconstrained")))
  reference <- sqrt(diag(solve(stats::optimHess(theta, negative))))
  data.frame(coordinate = names(reported), analytic_score = unname(analytic),
             numerical_score = numerical, reported_error = unname(reported),
             reference_error = unname(reference))
}

test_that("a FIML covariate fit reproduces the observed-data likelihood", {
  data <- .cov_missing_fixture()
  models <- c(diagonal = "diagonal", full = "full")
  fits <- lapply(models, \(covariance_model) .cov_missing_fit(data, covariance_model))
  expect_identical(vapply(fits, \(fit) fit$missing, character(1)),
                   c(diagonal = "fiml", full = "fiml"))
  # The same likelihood summed twice; only the order of the additions differs.
  expect_reference(vapply(fits, \(fit) as.numeric(logLik(fit)), numeric(1)),
                   vapply(fits, \(fit) .cov_missing_likelihood(fit, data), numeric(1)),
                   reference = "independent: observed-data covariate mixture likelihood, base R",
                   covers = c("multilpa()", "logLik.multilpa_covariates()"),
                   quantity = "FIML log likelihood by covariance model", tolerance = 1e-10,
                   measure = "relative")
})

test_that("scores and standard errors with missing indicators match numerical ones", {
  data <- .cov_missing_fixture()
  models <- c(diagonal = "diagonal", full = "full")
  fits <- lapply(models, \(covariance_model) .cov_missing_fit(data, covariance_model))
  skip_if_not(all(vapply(fits, \(fit) isTRUE(fit$converged), logical(1))),
              "no start converged on this fixture")
  derivatives <- do.call(rbind, lapply(models, \(model) {
    cbind(covariance_model = model, .cov_missing_derivatives(fits[[model]], data))
  }))
  # 1e-5 absolute, the original bound: a central difference with step 1e-5
  # carries a truncation and rounding error of that order.
  expect_reference(data.frame(derivatives[c("covariance_model", "coordinate")],
                              score = derivatives$analytic_score),
                   data.frame(derivatives[c("covariance_model", "coordinate")],
                              score = derivatives$numerical_score),
                   reference = "independent: central differences (step 1e-5) of the covariate likelihood",
                   covers = "multilpa()", quantity = "score off the maximum, missing indicators",
                   tolerance = 1e-5, measure = "absolute", key = c("covariance_model", "coordinate"),
                   note = "Package side: the internal analytic group scores (.multilpa_cov_group_scores) behind vcov(); reference: central differences of the package's internal likelihood (.multilpa_cov_expectation).")
  # 1e-4 relative, the original bound: optimHess() without a gradient takes
  # second finite differences of the likelihood.
  expect_reference(data.frame(derivatives[c("covariance_model", "coordinate")],
                              standard_error = derivatives$reported_error),
                   data.frame(derivatives[c("covariance_model", "coordinate")],
                              standard_error = derivatives$reference_error),
                   reference = "stats", covers = c("multilpa()", "vcov.multilpa_covariates()"),
                   quantity = "standard errors, missing indicators, unconstrained scale",
                   tolerance = 1e-4, measure = "relative", key = c("covariance_model", "coordinate"),
                   note = "Reference: the inverse of stats::optimHess() of the package's internal covariate likelihood (.multilpa_cov_expectation) at the estimate.")
})
