# Audit 2026-11-01, slice covariates-diagnostics: regression tests for the
# defects fixed in R/covariates.R, R/cross-level.R and R/report.R, and the
# independent checks of the membership-logit M-step, the covariate standard
# errors and the bivariate residual statistic.

skip_on_cran()

audit_warning_classes <- function(expr) {
  seen <- character()
  value <- withCallingHandlers(expr, warning = function(condition) {
    seen <<- c(seen, class(condition)[[1L]])
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = seen)
}

# 150 units around zero and 7 far away: the second profile has a Gaussian
# membership mass of 7.
audit_thin_data <- function() {
  withr::local_seed(2)
  data <- data.frame(rbind(matrix(stats::rnorm(600), 150L),
                           matrix(stats::rnorm(28, 6, 0.5), 7L)))
  names(data) <- paste0("y", 1:4)
  data$z <- stats::rnorm(157L)
  data
}

test_that("a covariate fit counts a constrained structure's own parameters as multilpa() does", {
  data <- audit_thin_data()
  vars <- paste0("y", 1:4)
  # VEI: 4 means and one volume per profile = 5 own parameters, below the
  # mass of 7, so neither fit may call the profile thin. The covariate fit
  # counted 8 (VVI's means and variances) and warned.
  plain <- audit_warning_classes(lpa(data, vars, n_profiles = 2, model = "VEI",
                                     n_starts = 2, seed = 1))
  covariate <- audit_warning_classes(lpa(data, vars, n_profiles = 2, model = "VEI",
                                         profile_covariates = "z", n_starts = 2,
                                         seed = 1))
  expect_equal(.latents_profile_parameters(4L, "varying", "diagonal", "VEI"), 5)
  expect_false("latents_thin_profiles" %in% plain$warnings)
  expect_false("latents_thin_profiles" %in% covariate$warnings)
  expect_equal(sum(covariate$value$effective_profile_counts < 8), 1L)
  # VVI has 8 own parameters, so both still warn: the check was not removed.
  expect_warning(lpa(data, vars, n_profiles = 2, model = "VVI", profile_covariates = "z",
                     n_starts = 2, seed = 1), class = "latents_thin_profiles")
})

test_that("a covariate fit holding its variances has no thin profile to report", {
  data <- audit_thin_data()
  vars <- paste0("y", 1:4)
  plain <- suppressWarnings(lpa(data, vars, n_profiles = 2, model = "VVI", n_starts = 2,
                                seed = 1))
  held <- audit_warning_classes(lpa(data, vars, n_profiles = 2, model = "VVI",
                                    start = starting_values(plain), fixed = "variances",
                                    profile_covariates = "z", n_starts = 1, seed = 1))
  expect_identical(held$value$fixed, "variances")
  expect_false("latents_thin_profiles" %in% held$warnings)
})

audit_cross_level_data <- function() {
  withr::local_seed(11)
  n <- 400L
  state <- sample(3L, n, TRUE)
  data.frame(id = rep(seq_len(40L), each = 10L),
             y1 = stats::rnorm(n, c(-2, 0, 2)[state]),
             y2 = stats::rnorm(n, c(1, -1, 0)[state]))
}

test_that("the restricted cross-level replication counts starts, not its two factors", {
  data <- audit_cross_level_data()
  fit <- suppressWarnings(multilpa(data, c("y1", "y2"), id = "id", n_profiles = 3,
                                   n_group_classes = 2, family = "restricted_cross_level",
                                   n_starts = 5, seed = 1))
  # Hand count: each factor ran five starts; the joint maximum is replicated
  # as often as its less replicated factor.
  expect_identical(nrow(get_results(fit, "starts")), 2L)
  replicated <- min(get_results(fit, "starts")$n_best_replicated)
  replication <- .multilpa_replication_table(fit)
  expect_identical(replication$n_starts, 5L)
  expect_identical(replication$n_best_replicated, as.integer(replicated))
  expect_equal(replication$share, replicated / 5)
  # Per-start convergence is not recorded for this family, so it is not invented.
  expect_true(is.na(replication$n_converged))
  checks <- get_results(diagnostics(fit, plots = FALSE), "checks")
  expect_identical(checks$display[checks$check == "best log likelihood reached by"],
                   sprintf("%d of 5 starts", replicated))
  expect_false("starts failed or unfinished" %in% checks$check)
})

test_that("cross-level tables refuse a misspelt table and name an ignored argument", {
  data <- audit_cross_level_data()
  fit <- suppressWarnings(multilpa(data, c("y1", "y2"), id = "id", n_profiles = 2,
                                   n_group_classes = 2, family = "full_cross_level",
                                   n_starts = 2, seed = 1))
  expect_error(get_results(fit, "bogus"), class = "latents_bad_argument")
  expect_error(as.data.frame(fit, what = "bogus"), class = "latents_bad_argument")
  expect_error(plot(fit, what = "bogus"), class = "latents_bad_argument")
  expect_warning(table <- get_results(fit, "groups", level = "groups"),
                 class = "latents_unused_argument")
  expect_identical(table, get_results(fit, "groups"))
  # A unique partial name keeps working, as under match.arg().
  expect_identical(get_results(fit, "comp"), get_results(fit, "composition"))
  # Without bivariate residuals the table is absent from the diagnostics.
  expect_error(get_results(diagnostics(fit, plots = FALSE), "residuals"),
               class = "latents_bad_argument")
})

test_that("an invalid cross-level variance mode is refused, not counted as failed starts", {
  data <- audit_cross_level_data()
  expect_error(multilpa(data, c("y1", "y2"), id = "id", n_profiles = 2, n_group_classes = 2,
                        family = "full_cross_level", between_variance = "bogus",
                        n_starts = 1, seed = 1),
               class = "latents_bad_argument")
  expect_error(multilpa(data, c("y1", "y2"), id = "id", n_profiles = 2, n_group_classes = 2,
                        family = "restricted_cross_level", variance_model = "bogus",
                        n_starts = 1, seed = 1),
               class = "latents_bad_argument")
})

test_that("the membership-logit M-step is the weighted multinomial maximum", {
  withr::local_seed(11)
  n <- 60L
  z <- stats::rnorm(n, 5, 2)
  design <- cbind(1, z)
  # Two categories: the weighted binomial GLM (stats IRLS) is the same model.
  share <- stats::plogis(0.5 + 0.4 * (z - 5))
  posterior <- cbind(share, 1 - share) * stats::runif(n, 0.5, 2)
  solved <- .multilpa_weighted_logits(design, posterior, matrix(0, 2L, 1L))
  glm_fit <- stats::glm(cbind(posterior[, 1L], posterior[, 2L]) ~ z, family = stats::quasibinomial(),
                        control = stats::glm.control(epsilon = 1e-14, maxit = 100L))
  expect_true(solved$converged)
  expect_equal(as.vector(solved$coefficients), unname(stats::coef(glm_fit)), tolerance = 1e-8)
  # Three categories against a direct maximization of the same log likelihood.
  eta <- cbind(0.5 + 0.4 * (z - 5), -0.3 - 0.6 * (z - 5), 0)
  counts <- exp(eta) / rowSums(exp(eta)) * stats::runif(n, 0.5, 2)
  solved <- .multilpa_weighted_logits(design, counts, matrix(0, 2L, 2L))
  negative_log_likelihood <- function(values) {
    scores <- cbind(design %*% matrix(values, 2L, 2L), 0)
    -sum(counts * (scores - log(rowSums(exp(scores)))))
  }
  direct <- stats::optim(rep(0, 4L), negative_log_likelihood, method = "BFGS",
                         control = list(reltol = 1e-15, maxit = 1000L))
  expect_equal(as.vector(solved$coefficients), direct$par, tolerance = 1e-5)
  expect_lte(negative_log_likelihood(as.vector(solved$coefficients)), direct$value + 1e-10)
})

test_that("covariate coefficient standard errors are the inverse numerical Hessian", {
  skip_if_not_installed("numDeriv")
  withr::local_seed(5)
  n_groups <- 60L
  group <- rep(seq_len(n_groups), each = 8L)
  group_class <- sample(1:2, n_groups, TRUE)[group]
  z <- stats::rnorm(length(group))
  w <- stats::rnorm(n_groups)[group]
  profile <- ifelse(stats::runif(length(group)) <
                      stats::plogis(ifelse(group_class == 1L, 1, -1) + 1.2 * z), 1L, 2L)
  data <- data.frame(id = group, a = c(-1.5, 1.5)[profile] + stats::rnorm(length(group)),
                     b = c(1, -1)[profile] + stats::rnorm(length(group)), z = z, w = w)
  check <- function(fit) {
    theta <- .multilpa_cov_encode(fit)
    x <- sweep(as.matrix(fit$indicator_data), 2L, fit$center, "-")
    log_likelihood <- function(values) {
      pieces <- .multilpa_cov_decode(values, fit)
      .multilpa_cov_expectation(x, fit$group_index, pieces$parameters, fit$profile_design,
                                fit$group_design, pieces$beta, pieces$gamma)$log_likelihood
    }
    expect_equal(log_likelihood(theta), fit$log_likelihood, tolerance = 1e-10)
    covariance <- solve(-numDeriv::hessian(log_likelihood, theta))
    coefficient <- .multilpa_cov_labels(fit)$parameter == "coefficient"
    table <- parameter_inference(fit)
    expect_equal(table$standard_error[table$parameter == "coefficient"],
                 sqrt(diag(covariance))[coefficient], tolerance = 1e-6)
  }
  check(suppressWarnings(multilpa(data, c("a", "b"), "id", n_profiles = 2, n_group_classes = 2,
                                  profile_covariates = "z", group_covariates = "w",
                                  n_starts = 3, seed = 1)))
  check(suppressWarnings(multilpa(data, c("a", "b"), "id", n_profiles = 2, n_group_classes = 2,
                                  profile_covariates = "z", profile_slopes = "group_class",
                                  n_starts = 3, seed = 1)))
})

test_that("the bivariate residual of a 2x2 table is the hand Pearson chi-square, BH-adjusted in diagnostics()", {
  withr::local_seed(3)
  n <- 400L
  class <- sample(1:2, n, TRUE, prob = c(0.6, 0.4))
  probabilities <- rbind(c(0.8, 0.7, 0.75, 0.2), c(0.2, 0.3, 0.25, 0.7))
  data <- as.data.frame(lapply(1:4, function(j) {
    ifelse(stats::runif(n) < probabilities[class, j], "yes", "no")
  }), col.names = paste0("u", 1:4))
  fit <- lca(data, paste0("u", 1:4), n_classes = 2, n_starts = 3, seed = 1)
  posterior <- fit$subject_posteriors[, 1L]
  levels_1 <- fit$categorical_levels$u1
  levels_2 <- fit$categorical_levels$u2
  observed <- vapply(levels_2, function(b) vapply(levels_1, function(a) {
    sum(posterior[data$u1 == a & data$u2 == b])
  }, numeric(1)), numeric(2))
  expected <- sum(posterior) * outer(fit$response_probabilities$u1[1L, ],
                                     fit$response_probabilities$u2[1L, ])
  chi_square <- sum((observed - expected)^2 / expected)
  residuals <- get_results(fit, "residuals")
  row <- subset(residuals, profile == "profile_1" & indicator_1 == "u1" & indicator_2 == "u2")
  expect_equal(row$statistic, chi_square, tolerance = 1e-10)
  expect_identical(row$df, 1L)
  expect_equal(row$p_value, stats::pchisq(chi_square, 1, lower.tail = FALSE), tolerance = 1e-10)
  # Both diagnostic entry points apply BH; explicit none retains raw values.
  expect_equal(residuals$p_adjusted, stats::p.adjust(residuals$p_value, "BH"))
  raw <- get_results(fit, "residuals", adjust = "none")
  expect_equal(raw$p_adjusted, raw$p_value)
  checked <- diagnostics(fit, plots = FALSE)
  adjusted <- get_results(checked, "residuals")
  expect_equal(adjusted$p_adjusted, stats::p.adjust(adjusted$p_value, "BH"))
  checks <- get_results(checked, "checks")
  local <- checks[checks$domain == "Local independence", , drop = FALSE]
  expect_identical(local$criterion, "BH-adjusted p < 0.05")
  expect_identical(local$value, as.numeric(sum(adjusted$p_adjusted < 0.05)))
})

test_that("a bad `step` on a covariate fit is refused with a classed error", {
  fit <- multilpa(engagement_small, c("browse", "lectures", "forum_read"), "student",
                  n_profiles = 2, n_group_classes = 1,
                  profile_covariates = "previous_grade", n_starts = 1, seed = 1)
  expect_error(parameter_inference(fit, step = -1), class = "latents_bad_argument")
  expect_error(parameter_inference(fit, step = 0), class = "latents_bad_argument")
  expect_error(parameter_inference(fit, step = c(1e-4, 1e-4)), class = "latents_bad_argument")
})
