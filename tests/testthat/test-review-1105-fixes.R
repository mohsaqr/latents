# Fixes from the review of 2026-11-04 (validation/review-2026-11-04/SUMMARY.md).

two_step_lca_data <- function() {
  withr::with_seed(1, {
    n <- 300
    z <- stats::rnorm(n)
    class <- stats::rbinom(n, 1, stats::plogis(0.2 + 0.8 * z)) + 1L
    data <- data.frame(g = seq_len(n), z = z)
    data[paste0("u", 1:5)] <- lapply(1:5, function(item) {
      factor(ifelse(stats::runif(n) < c(0.75, 0.3)[class], "yes", "no"))
    })
    data
  })
}

test_that("a two-step class model's standard errors hold the response probabilities", {
  skip_if_not_installed("numDeriv")
  data <- two_step_lca_data()
  items <- paste0("u", 1:5)
  base <- multilpa(data, items, "g", n_profiles = 2, n_group_classes = 1,
                   categorical = items, n_starts = 5, seed = 1, tol = 1e-10,
                   max_iter = 5000)
  fit <- multilpa(data, items, "g", n_profiles = 2, n_group_classes = 1,
                  categorical = items, profile_covariates = "z",
                  start = starting_values(base, what = "measurement"),
                  fixed = "measurement", n_starts = 1, seed = 1, tol = 1e-10,
                  max_iter = 5000)
  held <- fit$response_probabilities
  density <- vapply(1:2, function(class) {
    exp(Reduce(`+`, lapply(items, function(item) {
      log(held[[item]][class, as.character(data[[item]])])
    })))
  }, numeric(nrow(data)))
  log_likelihood <- function(beta) {
    share <- stats::plogis(beta[1] + beta[2] * data$z)
    sum(log(share * density[, 1] + (1 - share) * density[, 2]))
  }
  estimates <- as.vector(fit$profile_coefficients)
  by_hand <- sqrt(diag(solve(-numDeriv::hessian(log_likelihood, estimates))))
  expect_equal(parameter_inference(fit)$standard_error, by_hand, tolerance = 1e-3)
  expect_identical(dim(stats::vcov(fit, scale = "unconstrained")),
                   c(fit$n_parameters, fit$n_parameters))
})

staged_and_joint <- function() {
  vars <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
  list(
    vars = vars,
    joint = multilpa(course_engagement, vars, "student", n_profiles = 2,
                     n_group_classes = 2, n_starts = 3, seed = 1),
    staged = fit_staged(course_engagement, vars, "student", n_profiles = 2,
                        n_group_classes = 2, n_starts = 3, seed = 1))
}

test_that("a staged fit is charged for its measurement against a joint fit", {
  fits <- staged_and_joint()
  table <- compare_models(joint = fits$joint, staged = fits$staged)
  staged <- subset(table, model == "staged")
  expect_identical(staged$n_parameters,
                   as.integer(fits$staged$n_parameters_with_measurement))
  expect_identical(staged$n_parameters, subset(table, model == "joint")$n_parameters)
  n_groups <- length(unique(course_engagement$student))
  expect_equal(staged$bic_groups, -2 * fits$staged$log_likelihood +
                 staged$n_parameters * log(n_groups))
  expect_equal(staged$aic, -2 * fits$staged$log_likelihood + 2 * staged$n_parameters)
  # Same model and measurement, lower likelihood: the staged fit cannot win.
  expect_gte(staged$bic_groups, subset(table, model == "joint")$bic_groups)
  expect_match(attr(table, "title"), "held measurement")
  # The stored fit is untouched; only the comparison pays.
  expect_lt(fits$staged$n_parameters, staged$n_parameters)
  expect_equal(fits$staged$bic, -2 * fits$staged$log_likelihood +
                 fits$staged$n_parameters * log(n_groups))
})

test_that("fits holding the same measurement are compared on what they estimate", {
  fits <- staged_and_joint()
  measurement <- fits$staged$stage_one
  other <- fit_staged(course_engagement, fits$vars, "student", n_profiles = 2,
                      n_group_classes = 3, measurement = measurement,
                      n_starts = 3, seed = 1)
  table <- compare_models(two = fits$staged, three = other)
  expect_identical(table$n_parameters,
                   as.integer(c(fits$staged$n_parameters, other$n_parameters)))
  expect_identical(attr(table, "title"), "Model comparison")
  expect_equal(table$bic_groups, c(fits$staged$bic, other$bic))
})

test_that("comparisons without a held measurement do not change", {
  fits <- staged_and_joint()
  other <- multilpa(course_engagement, fits$vars, "student", n_profiles = 3,
                    n_group_classes = 1, n_starts = 3, seed = 1)
  table <- compare_models(a = fits$joint, b = other)
  expect_identical(table$n_parameters,
                   as.integer(c(fits$joint$n_parameters, other$n_parameters)))
  expect_equal(table$bic_groups, c(fits$joint$bic, other$bic))
  expect_identical(attr(table, "title"), "Model comparison")
})

test_that("p_adjusted is Benjamini-Hochberg by default in both residual verbs", {
  fit <- lca(student_esm, c("time_with_friends", "on_social_media", "sports", "walking"),
             n_classes = 2, n_starts = 2, seed = 1)
  residuals <- get_results(fit, "residuals")
  expect_equal(residuals$p_adjusted, stats::p.adjust(residuals$p_value, method = "BH"))
  expect_false(isTRUE(all.equal(residuals$p_adjusted, residuals$p_value)))
  expect_equal(as.data.frame(get_residuals(fit))$p_adjusted, residuals$p_adjusted)
  raw <- get_results(fit, "residuals", adjust = "none")
  expect_equal(raw$p_adjusted, raw$p_value)
  expect_error(get_results(fit, "residuals", adjust = "nope"),
               class = "latents_bad_argument")
})

test_that("a cross-level fit with group classes warns about its criteria", {
  vars <- c("browse", "lectures")
  expect_warning(
    multilpa(course_engagement, vars, "student", n_profiles = 2, n_group_classes = 2,
             family = "restricted_cross_level", n_starts = 2, seed = 1),
    class = "latents_cross_level_criteria")
  seen <- character()
  withCallingHandlers(
    multilpa(course_engagement, vars, "student", n_profiles = 2, n_group_classes = 1,
             family = "full_cross_level", n_starts = 2, seed = 1),
    warning = function(condition) {
      seen <<- c(seen, class(condition)[1L])
      invokeRestart("muffleWarning")
    })
  expect_false("latents_cross_level_criteria" %in% seen)
})
