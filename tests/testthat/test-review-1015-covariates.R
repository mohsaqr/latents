# Regressions from the 2026-10-15 review of the covariate, three-step, staged
# and pooling verbs. Each test failed before its repair.

review_r3step_data <- function() {
  set.seed(21)
  g <- rep(seq_len(50), each = 12)
  x <- stats::rnorm(600)
  truth <- 1L + as.integer(stats::runif(600) < stats::plogis(-0.3 + 1.2 * x))
  data.frame(
    g = g, x = x,
    a = stats::rnorm(600, ifelse(truth == 2L, 1.2, -1.2)),
    b = stats::rnorm(600, ifelse(truth == 2L, 1.2, -1.2))
  )
}

review_covariate_data <- function() {
  set.seed(5)
  school <- rep(seq_len(40), each = 8)
  x <- stats::rnorm(320)
  profile <- ifelse(stats::runif(320) < stats::plogis(-0.5 + 0.8 * x), 2L, 1L)
  data.frame(
    school = school, x = x, x2 = 2 * x, w = rep(seq_len(40) / 10, each = 8),
    y1 = stats::rnorm(320, ifelse(profile == 2L, 2, -2), 0.7),
    y2 = stats::rnorm(320, ifelse(profile == 2L, 1.5, -1.5), 0.7)
  )
}

test_that("r3step() does not depend on where a covariate's zero lies", {
  data <- review_r3step_data()
  fit <- multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                  n_starts = 4, seed = 1)
  base <- r3step(fit, data, "x")
  shift <- 1e4
  shifted_data <- within(data, x <- x + shift)
  # Before the repair this returned a slope of about zero with a finite SE.
  shifted <- r3step(fit, shifted_data, "x")
  expect_equal(shifted$estimate[shifted$term == "x"],
               base$estimate[base$term == "x"], tolerance = 1e-6)
  expect_equal(shifted$standard_error[shifted$term == "x"],
               base$standard_error[base$term == "x"], tolerance = 1e-5)
  # The intercept moves exactly by the shift times the slope.
  expect_equal(shifted$estimate[shifted$term == "(Intercept)"] +
                 shift * shifted$estimate[shifted$term == "x"],
               base$estimate[base$term == "(Intercept)"], tolerance = 1e-6)
})

test_that("r3step() reaches the maximum of its likelihood", {
  data <- review_r3step_data()
  fit <- multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                  n_starts = 4, seed = 1)
  result <- r3step(fit, data, "x")
  pieces <- .multilpa_level_assignments(fit, "individuals")
  by_assigned <- t(.multilpa_error_matrix(pieces))[pieces$modal, , drop = FALSE]
  design <- cbind(1, data$x)
  objective <- function(theta) {
    prior <- .multilpa_logit_prior(matrix(theta, 2L, 1L), design)
    -sum(log(rowSums(prior * by_assigned)))
  }
  tight <- stats::optim(c(0, 0), objective, method = "BFGS",
                        control = list(maxit = 5000L, reltol = 1e-15))
  tight <- stats::optim(tight$par, objective, method = "Nelder-Mead",
                        control = list(maxit = 20000L, reltol = 1e-15))
  # Before the repair the default BFGS tolerance left it 1e-4 away.
  expect_equal(result$estimate, tight$par, tolerance = 1e-6)
})

test_that("r3step(by_group_class = TRUE) does not depend on a covariate's origin", {
  data <- review_r3step_data()
  fit <- multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 2,
                  n_starts = 4, seed = 1)
  base <- r3step(fit, data, "x", by_group_class = TRUE)
  # Before the repair this shift stopped the search with latents_no_converge.
  shifted <- r3step(fit, within(data, x <- x + 100), "x", by_group_class = TRUE)
  expect_equal(shifted$estimate[shifted$term == "x"],
               base$estimate[base$term == "x"], tolerance = 1e-6)
  expect_equal(attr(shifted, "log_likelihood"), attr(base, "log_likelihood"),
               tolerance = 1e-8)
  intercepts <- base$term %in% c("group_class_1", "group_class_2")
  expect_equal(shifted$estimate[intercepts] +
                 100 * shifted$estimate[shifted$term == "x"],
               base$estimate[intercepts], tolerance = 1e-6)
})

test_that("a covariate far from zero is not mistaken for separation", {
  data <- review_covariate_data()
  base <- multilpa(data, c("y1", "y2"), "school", n_profiles = 2,
                   n_group_classes = 1, profile_covariates = "x",
                   n_starts = 2, seed = 1)
  base_inference <- parameter_inference(base)
  # Before the repair the fit warned latents_extreme_coefficients and the
  # inference was refused with latents_boundary_fit.
  expect_no_warning(shifted <- multilpa(
    within(data, x <- x + 100), c("y1", "y2"), "school", n_profiles = 2,
    n_group_classes = 1, profile_covariates = "x", n_starts = 2, seed = 1))
  expect_false(shifted$extreme_logits)
  inference <- parameter_inference(shifted)
  expect_equal(inference$estimate[inference$term == "x"],
               base_inference$estimate[base_inference$term == "x"],
               tolerance = 1e-6)
  expect_equal(inference$standard_error[inference$term == "x"],
               base_inference$standard_error[base_inference$term == "x"],
               tolerance = 1e-5)
})

test_that("the one-step membership regression finds the maximum far from zero", {
  data <- review_covariate_data()
  base <- multilpa(data, c("y1", "y2"), "school", n_profiles = 2,
                   n_group_classes = 1, profile_covariates = "x",
                   n_starts = 2, seed = 1)
  # Before the repair the Newton search stopped 0.0024 below the maximum,
  # reported as converged.
  shifted <- multilpa(within(data, x <- x + 1e6), c("y1", "y2"), "school",
                      n_profiles = 2, n_group_classes = 1,
                      profile_covariates = "x", n_starts = 2, seed = 1)
  expect_equal(shifted$log_likelihood, base$log_likelihood, tolerance = 1e-9)
  expect_equal(shifted$profile_coefficients["x", 1L],
               base$profile_coefficients["x", 1L], tolerance = 1e-6)
})

test_that("covariate-model input errors carry classes", {
  data <- review_covariate_data()
  fit_with <- function(...) {
    multilpa(data, c("y1", "y2"), "school", n_starts = 1, seed = 1, ...)
  }
  expect_error(fit_with(n_profiles = 1, n_group_classes = 1,
                        profile_covariates = "x"),
               class = "latents_bad_argument")
  expect_error(fit_with(n_profiles = 2, n_group_classes = 1,
                        group_covariates = "w"),
               class = "latents_bad_argument")
  expect_error(fit_with(n_profiles = 2, n_group_classes = 1,
                        profile_covariates = "y1"),
               class = "multilpa_indicator_reused")
  expect_error(fit_with(n_profiles = 2, n_group_classes = 1,
                        profile_covariates = "school"),
               class = "latents_bad_covariate")
  expect_error(fit_with(n_profiles = 2, n_group_classes = 1,
                        profile_covariates = c("x", "x2")),
               class = "latents_bad_covariate")
  expect_error(fit_with(n_profiles = 2, n_group_classes = 1,
                        profile_covariates = c("x", "x")),
               class = "latents_bad_argument")
  expect_error(fit_with(n_profiles = 2, n_group_classes = 1,
                        profile_covariates = "absent"),
               class = "latents_bad_data")
  expect_error(multilpa(data, c("y1", "y2"), "school", n_profiles = 2,
                        n_group_classes = 1, profile_covariates = "x",
                        n_starts = 0),
               class = "latents_bad_argument")
})

test_that("confint() on a covariate fit refuses an unknown parameter by class", {
  data <- review_covariate_data()
  fit <- multilpa(data, c("y1", "y2"), "school", n_profiles = 2,
                  n_group_classes = 1, profile_covariates = "x",
                  n_starts = 2, seed = 1)
  expect_error(confint(fit, parm = "absent"), class = "latents_bad_argument")
  expect_error(confint(fit, parm = 0), class = "latents_bad_argument")
})

test_that("fit_staged() refuses a measurement it cannot hold faithfully", {
  set.seed(7)
  data <- data.frame(
    school = rep(seq_len(30), each = 8),
    a = c(stats::rnorm(120, -1), stats::rnorm(120, 1)),
    b = c(stats::rnorm(120, -1), stats::rnorm(120, 1)),
    wt = rep(stats::runif(30, 0.5, 2), each = 8),
    o = sample(seq_len(3), 240, TRUE)
  )
  measure <- function(...) {
    multilpa(data, n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1,
             ...)
  }
  stage <- function(vars, measurement) {
    fit_staged(data, vars, "school", n_profiles = 2, measurement = measurement,
               n_starts = 2, seed = 1)
  }
  # Before the repair each of these returned a second stage that held the
  # measurement under a different likelihood (unweighted, or raw scores
  # against person-centred means) or failed with an unrelated start error.
  expect_error(stage(c("a", "b"), measure(vars = c("a", "b"), id = "school",
                                          weights = "wt")),
               class = "latents_bad_stage")
  expect_error(stage(c("a", "b"), measure(vars = c("a", "b"), id = "school",
                                          centering = "person")),
               class = "latents_bad_stage")
  expect_error(stage(c("a", "b", "o"), measure(vars = c("a", "b", "o"),
                                               id = "school", ordinal = "o")),
               class = "latents_bad_stage")
  # An ordinary measurement is still carried forward unchanged.
  plain <- measure(vars = c("a", "b"), id = "school")
  staged <- stage(c("a", "b"), plain)
  expect_equal(unname(staged$means), unname(plain$means))
})

test_that("three_step() refuses level = 'groups' for a trajectory model", {
  fit <- mixture_regression(score ~ wave, growth_scores, n_classes = 2,
                            id = "student", class_level = "group",
                            n_starts = 1, seed = 1)
  # Before the repair the argument was silently ignored.
  expect_error(three_step(fit, growth_scores, "motivation", level = "groups"),
               class = "latents_bad_argument")
  expect_identical(unique(three_step(fit, growth_scores, "motivation")$level),
                   "persons")
})

test_that("person centring in the covariate model does not depend on units", {
  data <- review_covariate_data()
  data$y1 <- data$y1 + rep(seq_len(40) / 4, each = 8)
  fit_scaled <- function(units) {
    multilpa(within(data, {
      y1 <- y1 * units
      y2 <- y2 * units
    }), c("y1", "y2"), "school", n_profiles = 2, n_group_classes = 1,
    profile_covariates = "x", centering = "person", min_variance = 1e-30,
    n_starts = 2, seed = 1)
  }
  ordinary <- fit_scaled(1)
  # Before the repair this was refused with latents_bad_data ("Centring left
  # `y1`, `y2` constant").
  tiny <- fit_scaled(1e-9)
  # Agreement is to the EM stopping tolerance, not to rounding. Both starts
  # reach this maximum, with the profiles in opposite orders, so which order is
  # reported turns on the last digits of two equal likelihoods; with two
  # profiles a relabelling flips the coefficient's sign, which is compared
  # without it (review 2026-10-26 round 2: after the unit-free EM stop the
  # two scales reported opposite orders, 0.7763 and -0.7763).
  expect_equal(abs(tiny$profile_coefficients["x", 1L]),
               abs(ordinary$profile_coefficients["x", 1L]), tolerance = 1e-3)
  expect_equal(tiny$log_likelihood,
               ordinary$log_likelihood + 2 * 320 * log(1e9), tolerance = 1e-8)
})

test_that("covariate-fit inference does not depend on a covariate's origin", {
  data <- review_covariate_data()
  fit_shifted <- function(shift) {
    multilpa(within(data, x <- x + shift), c("y1", "y2"), "school",
             n_profiles = 2, n_group_classes = 1, profile_covariates = "x",
             n_starts = 2, seed = 1)
  }
  near <- fit_shifted(0)
  base <- parameter_inference(near)
  far <- fit_shifted(1e6)
  # Before the repair the finite-difference information was singular here.
  shifted <- parameter_inference(far)
  slope <- shifted$term == "x"
  expect_equal(shifted$standard_error[slope],
               base$standard_error[base$term == "x"], tolerance = 1e-6)
  robust_far <- parameter_inference(far, vcov_type = "robust")
  robust_near <- parameter_inference(near, vcov_type = "robust")
  expect_equal(robust_far$standard_error[slope],
               robust_near$standard_error[robust_near$term == "x"],
               tolerance = 1e-6)
})
