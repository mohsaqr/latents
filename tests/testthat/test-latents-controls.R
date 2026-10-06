# Three controls added together: the start a fit reports, the two-level form
# of R3STEP, and standard errors for a one-step covariate model whose
# indicators are categorical. Each is checked against the formula it is meant
# to solve, not against itself.

# A two-level latent class design with a student-level covariate, far from
# any boundary: two classes, five binary items, a slope of 0.8 on the logit.
.lca_design <- function(n_groups = 60, per = 8, continuous = FALSE) {
  set.seed(2024)
  g <- rep(seq_len(n_groups), each = per)
  cluster <- rep(sample(1:2, n_groups, TRUE), each = per)
  z <- rep(stats::rnorm(n_groups), each = per)
  class <- ifelse(stats::runif(length(g)) <
                    stats::plogis(ifelse(cluster == 1L, 1, -1) + 0.8 * z), 1L, 2L)
  draw <- function(p1, p2) {
    ifelse(stats::runif(length(g)) < ifelse(class == 1L, p1, p2), 1L, 2L)
  }
  data <- data.frame(g = g, z = z, y1 = draw(0.85, 0.2), y2 = draw(0.8, 0.25),
                     y3 = draw(0.75, 0.3), y4 = draw(0.2, 0.7),
                     y5 = draw(0.9, 0.35))
  if (continuous) data$w <- stats::rnorm(length(g), ifelse(class == 1L, 1, -1))
  data
}
.lca_items <- paste0("y", 1:5)

.categorical_covariate_fit <- function(...) {
  data <- .lca_design()
  fit <- multilca(data, vars = .lca_items, id = "g", n_profiles = 2,
                  n_group_classes = 2, profile_covariates = "z", n_starts = 2,
                  seed = 1, select_start = "converged", ...)
  list(fit = fit, data = data)
}

.numerical_gradient <- function(f, at, step = 1e-5) {
  vapply(seq_along(at), function(index) {
    shift <- replace(numeric(length(at)), index, step)
    (f(at + shift) - f(at - shift)) / (2 * step)
  }, numeric(1))
}

# ---- select_start ---------------------------------------------------------

test_that("the start rule takes the best converged start only when asked", {
  likelihood <- c(-10, -9.9999, -9.99992, -Inf)
  converged <- c(TRUE, FALSE, TRUE, FALSE)
  expect_identical(.multilpa_select_start(likelihood, converged), 2L)
  expect_identical(.multilpa_select_start(likelihood, converged, "converged"), 3L)
  # With nothing converged, both rules fall back to the highest likelihood.
  expect_identical(.multilpa_select_start(c(-10, -9.9), c(FALSE, FALSE),
                                          "converged"), 2L)
  # Among starts tied to within the tolerance, a converged one wins, then the
  # earliest, under either rule.
  tied <- c(-5, -5 - 1e-13, -5)
  expect_identical(.multilpa_select_start(tied, c(FALSE, TRUE, TRUE)), 2L)
  expect_identical(.multilpa_select_start(tied, c(TRUE, TRUE, TRUE)), 1L)
})

test_that("the start rule is invariant to failed starts and refuses bad input", {
  base <- .multilpa_select_start(c(-3, -2, -4), c(TRUE, TRUE, TRUE))
  padded <- .multilpa_select_start(c(-Inf, -3, -2, -4), c(FALSE, TRUE, TRUE, TRUE))
  expect_identical(padded, base + 1L)
  expect_error(.multilpa_select_start(c(-Inf, -Inf), c(FALSE, FALSE)),
               "finite likelihood")
  expect_error(multilpa(subset(course_engagement, student <= 5),
                        vars = c("browse", "lectures"), id = "student",
                        n_profiles = 2, n_group_classes = 1, n_starts = 1,
                        select_start = "best"),
               "should be one of")
})

test_that("select_start reaches every fitting path", {
  small <- subset(course_engagement, student <= 20)
  plain <- multilpa(small, vars = c("browse", "lectures"), id = "student",
                    n_profiles = 2, n_group_classes = 1, n_starts = 2,
                    seed = 1, select_start = "converged")
  expect_true(plain$converged)
  moves <- lta(small, vars = c("browse", "lectures"), id = "student",
               time = "sequence", n_profiles = 2, n_starts = 2, seed = 1,
               select_start = "converged")
  expect_true(moves$converged)
  covariate <- .categorical_covariate_fit()$fit
  expect_true(covariate$converged)
})

# ---- r3step(by_group_class = TRUE) ------------------------------------------

test_that("group-class R3STEP reduces to the pooled regression with one group class", {
  data <- subset(course_engagement, student <= 40)
  fit <- multilpa(data, c("browse", "lectures", "forum_read"), "student", 2, 1,
                  n_starts = 1, seed = 1)
  for (type in c("observed", "robust")) {
    pooled <- r3step(fit, data, "previous_grade", vcov_type = type)
    nested <- r3step(fit, data, "previous_grade", vcov_type = type,
                     by_group_class = TRUE)
    expect_true(all(nested$level == "individuals"))
    # The two BFGS paths use different starts and stop on likelihood change;
    # their coefficients agree to optimization precision, not bit for bit.
    expect_equal(nested$estimate, pooled$estimate, tolerance = 1e-4)
    expect_equal(nested$standard_error, pooled$standard_error, tolerance = 1e-5)
  }
})

.two_level_membership <- function(seed = 11, n_groups = 80, per = 10, slope = 1) {
  set.seed(seed)
  g <- rep(seq_len(n_groups), each = per)
  group_class <- rep(sample(1:2, n_groups, TRUE), each = per)
  z <- rep(stats::rnorm(n_groups), each = per)
  logit <- ifelse(group_class == 1L, 1.2, -1.2) + slope * z
  profile <- ifelse(stats::runif(length(g)) < stats::plogis(logit), 1L, 2L)
  shift <- ifelse(profile == 1L, 1.5, -1.5)
  data.frame(g = g, z = z, a = stats::rnorm(length(g), shift),
             b = stats::rnorm(length(g), shift),
             c = stats::rnorm(length(g), shift / 1.5))
}

test_that("the two-level R3STEP gradient is the derivative of its likelihood", {
  data <- .two_level_membership(n_groups = 40)
  fit <- multilpa(data, c("a", "b", "c"), "g", n_profiles = 2,
                  n_group_classes = 2, n_starts = 1, seed = 1)
  pieces <- .multilpa_level_assignments(fit, "individuals")
  model <- .multilpa_r3step_group_model(fit, as.matrix(data["z"]), pieces,
                                        .multilpa_error_matrix(pieces))
  probe <- model$start + 0.1 * rep_len(c(1, -1), length(model$start))
  expect_equal(unname(model$gradient(probe)),
               .numerical_gradient(model$objective, probe, 1e-6),
               tolerance = 1e-5)
  # The per-group scores sum to the gradient, which the sandwich relies on.
  expect_equal(-colSums(model$group_scores(probe)), model$gradient(probe))
})

test_that("the two-level R3STEP recovers the slope the pooled one attenuates", {
  data <- .two_level_membership()
  fit <- multilpa(data, c("a", "b", "c"), "g", n_profiles = 2,
                  n_group_classes = 2, n_starts = 1, seed = 1)
  nested <- r3step(fit, data, "z", by_group_class = TRUE)
  pooled <- r3step(fit, data, "z")
  slope <- abs(subset(nested, term == "z")$estimate)
  expect_lt(abs(slope - 1), 0.25)
  expect_lt(abs(subset(pooled, term == "z")$estimate), slope)
  # The group-class intercepts separate in the two directions they were built in.
  intercepts <- sort(subset(nested, grepl("^group_class_", term))$estimate)
  expect_lt(intercepts[1L], -0.6)
  expect_gt(intercepts[2L], 0.6)
  expect_true(isTRUE(attr(nested, "by_group_class")))
  expect_identical(attr(nested, "vcov_type"), "robust")
  expect_true(all(nested$standard_error > 0))
  # Only the slopes are tested and corrected.
  expect_identical(!is.na(nested$p_adjusted), nested$term == "z")
})

test_that("the two-level R3STEP is refused where it does not apply", {
  data <- .two_level_membership(n_groups = 30)
  fit <- multilpa(data, c("a", "b", "c"), "g", n_profiles = 2,
                  n_group_classes = 2, n_starts = 1, seed = 1)
  expect_error(r3step(fit, data, "z", level = "groups", by_group_class = TRUE),
               class = "latents_bad_argument")
  expect_error(r3step(fit, data, "z", by_group_class = NA), "TRUE or FALSE")
  small <- subset(course_engagement, student <= 20)
  moves <- lta(small, vars = c("browse", "lectures"), id = "student",
               time = "sequence", n_profiles = 2, n_starts = 1, seed = 1)
  expect_error(r3step(moves, small, "previous_grade", by_group_class = TRUE),
               class = "latents_bad_argument")
})

test_that("r3step now defaults to standard errors clustered on groups", {
  data <- .two_level_membership(n_groups = 30)
  fit <- multilpa(data, c("a", "b", "c"), "g", n_profiles = 2,
                  n_group_classes = 1, n_starts = 1, seed = 1)
  expect_identical(attr(r3step(fit, data, "z"), "vcov_type"), "robust")
})

# ---- categorical covariate inference ----------------------------------------

test_that("the categorical covariate score is the derivative of the likelihood", {
  fit <- .categorical_covariate_fit()$fit
  x <- matrix(numeric(0), fit$n_observations, 0L)
  codes <- fit$categorical_data
  log_likelihood <- function(theta) {
    pieces <- .multilpa_cov_decode(theta, fit)
    .multilpa_cov_expectation(x, fit$group_index, pieces$parameters,
                              fit$profile_design, fit$group_design,
                              pieces$beta, pieces$gamma, codes)$log_likelihood
  }
  theta <- .multilpa_cov_encode(fit)
  expect_equal(log_likelihood(theta), fit$log_likelihood)
  # Away from the maximum, where a wrong score cannot hide behind a zero.
  probe <- theta + 0.05 * rep_len(c(1, -1), length(theta))
  expect_equal(unname(colSums(.multilpa_cov_group_scores(probe, x, fit, codes))),
               .numerical_gradient(log_likelihood, probe), tolerance = 1e-5)
})

test_that("categorical covariate inference reports every category coherently", {
  skip_on_cran()
  pair <- .categorical_covariate_fit()
  fit <- pair$fit
  inference <- parameter_inference(fit)
  responses <- subset(inference, parameter == "response")
  expect_equal(nrow(responses), 2L * 2L * length(.lca_items))
  expect_true(all(is.na(responses$p_value)))
  expect_true(all(responses$standard_error > 0))
  # Two categories of one binary item are complements: equal standard errors.
  by_row <- split(responses$standard_error,
                  paste(responses$outcome, sub(":.*", "", responses$term)))
  expect_true(all(vapply(by_row, function(se) isTRUE(all.equal(se[1], se[2])),
                         logical(1))))
  # The natural covariance of a simplex sums to zero across its categories.
  covariance <- vcov(fit)
  first <- grep("response.profile_1.y1:", rownames(covariance), fixed = TRUE)
  expect_equal(unname(rowSums(covariance[first, first])), c(0, 0),
               tolerance = 1e-10)
  # Every accessor agrees on the parameter set and its order.
  expect_identical(rownames(covariance), names(coef(fit)))
  expect_identical(nrow(confint(fit)), nrow(inference))
  expect_equal(unname(coef(fit)), inference$estimate)
  expect_equal(nrow(vcov(fit, scale = "unconstrained")),
               length(.multilpa_cov_encode(fit)))
  # The stored frame and the supplied one give the same answer.
  expect_equal(parameter_inference(fit, pair$data), inference)
  robust <- parameter_inference(fit, vcov_type = "robust")
  expect_equal(robust$estimate, inference$estimate)
  expect_true(all(robust$standard_error > 0))
})

test_that("a mixed covariate fit carries both measurement blocks", {
  data <- .lca_design(continuous = TRUE)
  fit <- multilpa(data, vars = c(.lca_items, "w"), id = "g", n_profiles = 2,
                  n_group_classes = 2, profile_covariates = "z",
                  categorical = .lca_items, n_starts = 1, seed = 1,
                  select_start = "converged")
  inference <- parameter_inference(fit)
  expect_setequal(unique(inference$parameter),
                  c("mean", "variance", "response", "coefficient"))
  expect_equal(sum(inference$parameter == "response"), 2L * 2L * 5L)
  expect_true(all(inference$standard_error > 0))
})

test_that("OPG standard errors reproduce glca on the same data", {
  skip_on_cran()
  # glca 1.4.2, item(y1, ..., y5) ~ z, group = g, nclass = 2, ncluster = 2,
  # on this seeded design (150 groups of 10): log likelihood -4486.5647 and
  # standard errors 0.16569, 0.16736 (the two intercepts) and 0.13636 (z).
  data <- .lca_design(n_groups = 150, per = 10)
  fit <- multilca(data, vars = .lca_items, id = "g", n_profiles = 2,
                  n_group_classes = 2, profile_covariates = "z", n_starts = 1,
                  seed = 1)
  expect_equal(fit$log_likelihood, -4486.5647, tolerance = 1e-7)
  opg <- subset(parameter_inference(fit, vcov_type = "opg"), level == "profile")
  expect_equal(sort(opg$standard_error), sort(c(0.16569, 0.16736, 0.13636)),
               tolerance = 2e-3)
  # The observed information is what the curvature says; glca does not use it.
  observed <- subset(parameter_inference(fit), level == "profile" & term == "z")
  expect_lt(observed$standard_error, subset(opg, term == "z")$standard_error)
})

test_that("bound-active response probabilities are refused or held on request", {
  data <- .lca_design()
  # Item y6 never takes its second category in profile-defining rows, so one
  # probability is driven to the floor.
  data$y6 <- ifelse(data$y1 == 1L, 1L, sample(1:2, nrow(data), TRUE))
  fit <- multilca(data, vars = c(.lca_items, "y6"), id = "g", n_profiles = 2,
                  n_group_classes = 2, profile_covariates = "z", n_starts = 1,
                  seed = 1, select_start = "converged")
  floor_hit <- any(unlist(fit$response_probabilities) <= 1e-10 * (1 + 1e-7))
  skip_if_not(floor_hit && fit$converged, "no probability reached the floor")
  expect_error(parameter_inference(fit), class = "latents_boundary_fit")
  held <- parameter_inference(fit, boundary = "fix")
  expect_gt(length(attr(held, "fixed_at_bound")), 0L)
  flagged <- .multilpa_cov_parameter_names(fit) %in% attr(held, "fixed_at_bound")
  expect_true(all(is.na(held$standard_error[flagged])))
  expect_true(all(held$standard_error[!flagged] > 0))
})

test_that("an unconverged covariate fit is still refused for inference", {
  fit <- suppressWarnings(multilpa(
    .lca_design(), vars = .lca_items, id = "g", n_profiles = 2,
    n_group_classes = 2, profile_covariates = "z", categorical = .lca_items,
    n_starts = 1, max_iter = 1, seed = 1))
  skip_if(fit$converged, "one iteration happened to converge")
  expect_error(parameter_inference(fit), "converged")
})
