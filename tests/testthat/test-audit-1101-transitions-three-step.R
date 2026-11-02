# Audit 2026-11-01, slice transitions-three-step: lta(), the transition fit's
# methods, the three-step verbs and the label alignment. Each test guards one
# fix of the audit log validation/audit-2026-11-01/transitions-three-step.md.

skip_on_cran()

audit_engagement <- course_engagement
audit_vars <- c("browse", "lectures")

audit_transition_fit <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- lta(audit_engagement, audit_vars, "student", n_profiles = 2,
                     time = "sequence", n_starts = 1, seed = 1)
    }
    cached
  }
})

test_that("lta() classes its argument mistakes as multilpa() does", {
  bad <- function(...) {
    expect_error(lta(audit_engagement, ..., n_starts = 1), class = "latents_bad_argument")
  }
  bad(audit_vars, NULL, 2, "sequence")
  bad(audit_vars, 1, 2, "sequence")
  bad(audit_vars, c("student", "sequence"), 2, "sequence")
  bad(audit_vars, n_profiles = 2, time = "sequence")
  bad(audit_vars, "student", 2)
  bad(audit_vars, "student", 2, NULL)
  bad(1:2, "student", 2, "sequence")
  bad(audit_vars, "student", 2, "sequence", categorical = NA_character_)
  bad(audit_vars, "student", 2, "sequence", min_probability = 2)
  bad(audit_vars, "student", 2, "sequence", mover_stayer = "yes")
  # The general engine is reached after the same checks.
  bad(audit_vars, NULL, 2, "sequence", transitions = "occasion")
})

test_that("a valid lta() call is unchanged by the classed checks", {
  fit <- audit_transition_fit()
  expect_s3_class(fit, "multilpa_transitions")
  transitions <- get_results(fit, "transitions")
  # Every origin's row of transition probabilities sums to one.
  expect_equal(as.vector(tapply(transitions$probability, transitions$from, sum)),
               rep(1, fit$n_profiles), tolerance = 1e-12)
})

test_that("the transition table's flags and the Wald arguments are classed", {
  fit <- audit_transition_fit()
  expect_error(get_results(fit, "transitions", estimated = NA),
               class = "latents_bad_argument")
  expect_error(get_results(fit, "transitions", stable = c(TRUE, FALSE)),
               class = "latents_bad_argument")
  expect_error(parameter_inference(fit, level = 2), class = "latents_bad_argument")
  expect_error(parameter_inference(fit, step = -1), class = "latents_bad_argument")
  # vcov() reached the Hessian with any step: 0 divided by zero.
  expect_error(vcov(fit, step = 0), class = "latents_bad_argument")
  expect_error(vcov(fit, step = -1e-4), class = "latents_bad_argument")
  expect_error(vcov(fit, step = "a"), class = "latents_bad_argument")
  expect_error(confint(fit, level = 1.5), class = "latents_bad_argument")
  # A valid step still gives the covariance parameter_inference() reports.
  expect_equal(vcov(fit), attr(parameter_inference(fit), "covariance"),
               ignore_attr = TRUE, tolerance = 1e-12)
})

test_that("the ML polish stops at a singular Hessian and surfaces a malformed one", {
  set.seed(3)
  truth <- rep(1:2, each = 100)
  pieces <- list(posteriors = cbind(ifelse(truth == 1L, 0.9, 0.1), ifelse(truth == 1L, 0.1, 0.9)),
                 modal = truth, n_classes = 2L)
  values <- stats::rnorm(200, c(0, 3)[truth])
  fit <- .multilpa_step_ml_fit(pieces, values, "varying")
  standardized <- (values - fit$outcome_coordinates$center) / fit$outcome_coordinates$scale
  em <- fit
  em$means <- fit$outcome_coordinates$means
  em$variances <- fit$variances / fit$outcome_coordinates$scale^2
  # Move off the stationary point so the polish has a Newton step to take.
  em$means <- em$means + c(0.05, -0.05)
  polished <- .multilpa_step_ml_polish(em, standardized, "varying")
  expect_lt(abs(polished$means[1] - fit$outcome_coordinates$means[1]), 1e-6)
  # A singular Hessian gives no direction: the estimate is returned as it came.
  testthat::local_mocked_bindings(
    .multilpa_step_ml_hessian = function(gradient, parameters, scale) {
      matrix(0, length(parameters), length(parameters))
    })
  stopped <- .multilpa_step_ml_polish(em, standardized, "varying")
  expect_equal(stopped$means, em$means)
  # A malformed one is a mistake in the call, which is no longer swallowed.
  testthat::local_mocked_bindings(
    .multilpa_step_ml_hessian = function(gradient, parameters, scale) {
      diag(length(parameters) + 1L)
    })
  expect_error(.multilpa_step_ml_polish(em, standardized, "varying"))
})
