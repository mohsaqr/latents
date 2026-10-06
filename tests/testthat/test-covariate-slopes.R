# Profile-covariate slopes that differ across group classes
# (`multilpa(profile_slopes = "group_class")`).
#
# The likelihood is checked against its definition, the analytic score against
# numerical differentiation of that likelihood, and the standard errors against
# a Hessian of the likelihood alone.

skip_on_cran()

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
  expect_equal(fit$log_likelihood, .slopes_likelihood(fit, data), tolerance = 1e-10)
  shared <- .slopes_fit(data, profile_slopes = "shared")
  expect_equal(shared$log_likelihood, .slopes_likelihood(shared, data),
               tolerance = 1e-10)
})

test_that("the coefficients are named by group class and counted once each", {
  data <- .slopes_fixture()
  fit <- .slopes_fit(data)
  shared <- .slopes_fit(data, profile_slopes = "shared")
  profile_terms <- subset(get_results(fit, "coefficients"), level == "profile")$term
  expect_setequal(profile_terms, c("group_class_1", "group_class_2",
                                   "z:group_class_1", "z:group_class_2"))
  # One extra slope per extra group class, per covariate and profile logit.
  expect_identical(fit$n_parameters - shared$n_parameters, 1L)
  expect_identical(length(.multilpa_cov_encode(fit)), as.integer(fit$n_parameters))
})

test_that("the recovered slopes carry the simulated cross-level interaction", {
  data <- .slopes_fixture(n_groups = 100L, size = 10L)
  fit <- .slopes_fit(data)
  skip_if_not(fit$converged, "no start converged on this fixture")
  profile <- subset(parameter_inference(fit), level == "profile" &
                      grepl("^z:", term))
  # Labels are the fit's, so the test reads magnitudes and the sign contrast,
  # not which group class carries which sign.
  expect_true(all(abs(profile$estimate) > 1 & abs(profile$estimate) < 2))
  expect_lt(prod(sign(profile$estimate)), 0)
  expect_true(all(profile$p_value < 1e-6))
  shared <- .slopes_fit(data, profile_slopes = "shared")
  expect_gt(fit$log_likelihood, shared$log_likelihood + 10)
})

test_that("scores and standard errors with slopes by group class match numerical ones", {
  data <- .slopes_fixture()
  fit <- .slopes_fit(data)
  skip_if_not(fit$converged, "no start converged on this fixture")
  x <- sweep(as.matrix(data[c("y1", "y2")]), 2L, fit$center, "-")
  theta <- .multilpa_cov_encode(fit)
  negative <- function(value) {
    pieces <- .multilpa_cov_decode(value, fit)
    -.multilpa_cov_expectation(x, fit$group_index, pieces$parameters,
                               fit$profile_design, fit$group_design,
                               pieces$beta, pieces$gamma, NULL)$log_likelihood
  }
  set.seed(9)
  away <- theta + 0.05 * stats::rnorm(length(theta))
  analytic <- colSums(.multilpa_cov_group_scores(away, x, fit, NULL))
  numerical <- vapply(seq_along(away), function(index) {
    step <- replace(numeric(length(away)), index, 1e-5)
    -(negative(away + step) - negative(away - step)) / 2e-5
  }, numeric(1))
  expect_lt(max(abs(analytic - numerical)), 1e-5)
  reported <- sqrt(diag(vcov(fit, scale = "unconstrained")))
  reference <- sqrt(diag(solve(stats::optimHess(theta, negative))))
  expect_lt(max(abs(reported - reference) / reference), 1e-4)
  expect_equal(parameter_inference(fit, data)$standard_error,
               parameter_inference(fit)$standard_error)
})

test_that("with one group class the slopes by group class are the shared slopes", {
  data <- .slopes_fixture()
  varying <- .slopes_fit(data, n_group_classes = 1L)
  shared <- .slopes_fit(data, n_group_classes = 1L, profile_slopes = "shared")
  expect_equal(varying$log_likelihood, shared$log_likelihood, tolerance = 1e-10)
  expect_equal(unname(varying$profile_coefficients),
               unname(shared$profile_coefficients), tolerance = 1e-8)
})

test_that("rescaling a covariate rescales its group-class slopes and nothing else", {
  data <- .slopes_fixture()
  fit <- .slopes_fit(data)
  scaled_data <- data
  scaled_data$z <- 2 * scaled_data$z
  scaled <- .slopes_fit(scaled_data)
  expect_equal(scaled$log_likelihood, fit$log_likelihood, tolerance = 1e-8)
  slopes <- grepl("^z:", rownames(fit$profile_coefficients))
  expect_equal(scaled$profile_coefficients[slopes, ],
               fit$profile_coefficients[slopes, ] / 2, tolerance = 1e-5)
  expect_equal(scaled$profile_coefficients[!slopes, ],
               fit$profile_coefficients[!slopes, ], tolerance = 1e-5)
})

test_that("slopes by group class without a profile covariate are refused by class", {
  data <- .slopes_fixture()
  data$w <- stats::ave(data$z, data$g)
  expect_error(multilpa(data, c("y1", "y2"), "g", 2L, 2L,
                        profile_slopes = "group_class", n_starts = 1),
               class = "latents_bad_argument")
  expect_error(multilpa(data, c("y1", "y2"), "g", 2L, 2L, group_covariates = "w",
                        profile_slopes = "group_class", n_starts = 1),
               class = "latents_bad_argument")
})
