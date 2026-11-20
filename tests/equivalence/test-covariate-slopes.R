# Comparisons with independent references moved from tests/testthat/test-covariate-slopes.R,
# so that they run in CI's equivalence job rather than under R CMD check.

.slopes_fixture <- function(n_groups = 60L, size = 8L, seed = 4L) {
  set.seed(seed)
  n <- n_groups * size
  frame <- data.frame(g = rep(seq_len(n_groups), each = size),
                      z = stats::rnorm(n))
  group_class <- rep(stats::rbinom(n_groups, 1L, 0.5), each = size) + 1L
  # Opposite slopes in the two kinds of group: the cross-level interaction the
  # shared-slope model cannot express.
  logit <- c(0.5, -0.5)[group_class] + c(1.5, -1.5)[group_class] * frame$z
  state <- 2L - stats::rbinom(n, 1L, stats::plogis(logit))
  frame$y1 <- stats::rnorm(n, c(0, 2.5)[state])
  frame$y2 <- stats::rnorm(n, c(0, 2)[state])
  frame
}

.slopes_fit <- function(data, n_group_classes = 2L, profile_slopes = "group_class",
                        covariates = "z", tol = 1e-10) {
  quietly(multilpa(data, c("y1", "y2"), "g", 2L, n_group_classes,
                   profile_covariates = covariates,
                   profile_slopes = profile_slopes, n_starts = 3, seed = 1,
                   tol = tol, max_iter = 3000))
}

# The two-level likelihood from its definition, with each group class's profile
# prior written from the coefficient table rather than from the stored design.
.slopes_likelihood <- function(fit, data) {
  coefficients <- get_results(fit, "coefficients")
  logit_of <- function(term, level = "profile") {
    coefficients$estimate[coefficients$level == level &
                            coefficients$term == term]
  }
  density <- vapply(seq_len(2L), function(profile) {
    stats::dnorm(data$y1, fit$means[profile, 1L], sqrt(fit$variances[profile, 1L])) *
      stats::dnorm(data$y2, fit$means[profile, 2L], sqrt(fit$variances[profile, 2L]))
  }, numeric(nrow(data)))
  group_share <- stats::plogis(logit_of("(Intercept)", "group"))
  group_share <- if (length(group_share)) c(group_share, 1 - group_share) else 1
  rows_by_group <- split(seq_len(nrow(data)), data$g)
  sum(vapply(rows_by_group, function(rows) {
    log(sum(vapply(seq_len(fit$n_group_classes), function(class) {
      slope_term <- if (identical(fit$profile_slopes, "group_class"))
        sprintf("z:group_class_%d", class) else "z"
      first <- stats::plogis(logit_of(sprintf("group_class_%d", class)) +
                               logit_of(slope_term) * data$z[rows])
      group_share[class] *
        prod(first * density[rows, 1L] + (1 - first) * density[rows, 2L])
    }, numeric(1))))
  }, numeric(1)))
}

test_that("slopes by group class reproduce the likelihood written from its definition", {
  data <- .slopes_fixture()
  fit <- .slopes_fit(data)
  expect_identical(fit$profile_slopes, "group_class")
  shared <- .slopes_fit(data, profile_slopes = "shared")
  fits <- list(group_class = fit, shared = shared)
  # The same likelihood summed twice; only the order of the additions differs.
  expect_reference(vapply(fits, \(model) as.numeric(logLik(model)), numeric(1)),
                   vapply(fits, \(model) .slopes_likelihood(model, data), numeric(1)),
                   reference = "independent: two-level likelihood with group-class slopes, base R",
                   covers = c("multilpa()", "logLik.multilpa_covariates()"),
                   quantity = "log likelihood by slope model", tolerance = 1e-10,
                   measure = "relative")
})

test_that("scores and standard errors with slopes by group class match numerical ones", {
  data <- .slopes_fixture()
  fit <- .slopes_fit(data)
  skip_if_not(fit$converged, "no start converged on this fixture")
  x <- sweep(as.matrix(data[c("y1", "y2")]), 2L, fit$center, "-")
  theta <- latents:::.multilpa_cov_encode(fit)
  negative <- function(value) {
    pieces <- latents:::.multilpa_cov_decode(value, fit)
    -latents:::.multilpa_cov_expectation(x, fit$group_index, pieces$parameters,
                                         fit$profile_design, fit$group_design,
                                         pieces$beta, pieces$gamma, NULL)$log_likelihood
  }
  set.seed(9)
  away <- theta + 0.05 * stats::rnorm(length(theta))
  analytic <- colSums(latents:::.multilpa_cov_group_scores(away, x, fit, NULL))
  numerical <- vapply(seq_along(away), function(index) {
    step <- replace(numeric(length(away)), index, 1e-5)
    -(negative(away + step) - negative(away - step)) / 2e-5
  }, numeric(1))
  reported <- sqrt(diag(vcov(fit, scale = "unconstrained")))
  coordinates <- names(reported)
  # 1e-5 absolute, the original bound: a central difference with step 1e-5
  # carries a truncation and rounding error of that order.
  expect_reference(unname(analytic), numerical,
                   reference = "independent: central differences (step 1e-5) of the covariate likelihood",
                   covers = "multilpa()", quantity = "score off the maximum, slopes by group class",
                   tolerance = 1e-5, measure = "absolute", key = coordinates,
                   note = "Package side: the internal analytic group scores (.multilpa_cov_group_scores) behind vcov(); reference: central differences of the package's internal likelihood (.multilpa_cov_expectation).")
  reference <- sqrt(diag(solve(stats::optimHess(theta, negative))))
  # 1e-4 relative, the original bound: optimHess() without a gradient takes
  # second finite differences of the likelihood.
  expect_reference(reported, unname(reference), reference = "stats",
                   covers = c("multilpa()", "vcov.multilpa_covariates()"),
                   quantity = "standard errors, slopes by group class, unconstrained scale",
                   tolerance = 1e-4, measure = "relative",
                   note = "Reference: the inverse of stats::optimHess() of the package's internal covariate likelihood (.multilpa_cov_expectation) at the estimate.")
  expect_equal(parameter_inference(fit, data)$standard_error,
               parameter_inference(fit)$standard_error)
})
