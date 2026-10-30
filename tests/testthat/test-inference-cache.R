# The observed information is the slowest step of Wald inference: two score
# evaluations over all the data per parameter. It is computed once per
# session for identical inputs, and anew whenever anything it depends on
# differs.

skip_on_cran()

clear_inference_cache <- function() {
  rm(list = ls(.latents_cache), envir = .latents_cache)
}

# Counts the score evaluations of the profile-model path.
count_score_calls <- function(env = parent.frame()) {
  original <- .multilpa_score
  counter <- new.env()
  counter$calls <- 0L
  testthat::local_mocked_bindings(.multilpa_score = function(...) {
    counter$calls <- counter$calls + 1L
    original(...)
  }, .env = env)
  counter
}

cache_fit <- function() {
  quietly(multilpa(engagement_small, c("browse", "lectures"), "student",
                   n_profiles = 2, n_group_classes = 2, n_starts = 1, seed = 1,
                   tol = 1e-10))
}

test_that("a repeated inference reuses the information and gives the same result", {
  clear_inference_cache()
  fit <- cache_fit()
  counter <- count_score_calls()
  first <- parameter_inference(fit)
  first_calls <- counter$calls
  counter$calls <- 0L
  again <- parameter_inference(fit)
  expect_identical(again, first)
  # The second call evaluates the score once (the stationarity check) and
  # differentiates nothing.
  expect_identical(counter$calls, 1L)
  expect_gt(first_calls, 2L * fit$n_parameters)
  # A covariance built on the same information shares it.
  counter$calls <- 0L
  parameter_inference(fit, vcov_type = "robust")
  expect_identical(counter$calls, 1L)
})

test_that("anything the information depends on computes it anew", {
  clear_inference_cache()
  fit <- cache_fit()
  parameter_inference(fit)
  counter <- count_score_calls()
  parameter_inference(fit, step = 2e-4)
  expect_gt(counter$calls, fit$n_parameters)
  counter$calls <- 0L
  moved <- fit
  moved$means[1L, 1L] <- moved$means[1L, 1L] + 1e-9
  parameter_inference(moved)
  expect_gt(counter$calls, fit$n_parameters)
  # A value computed anew equals one computed without the cache.
  clear_inference_cache()
  fresh <- parameter_inference(fit)
  expect_identical(parameter_inference(fit), fresh)
})

test_that("the memo keeps the most recent entries and stores no failure", {
  clear_inference_cache()
  computed <- 0L
  value_of <- function(key) {
    .latents_memo("probe", list(key), function() {
      computed <<- computed + 1L
      key * 10
    }, size = 2L)
  }
  expect_identical(value_of(1), 10)
  expect_identical(value_of(2), 20)
  expect_identical(value_of(1), 10)
  expect_identical(computed, 2L)
  # A third key pushes out the least recently used one (2), not 1.
  value_of(3)
  value_of(1)
  expect_identical(computed, 3L)
  value_of(2)
  expect_identical(computed, 4L)
  failures <- 0L
  failing <- function() {
    .latents_memo("probe", list("fails"), function() {
      failures <<- failures + 1L
      stop(errorCondition("no", class = "latents_probe_failure"))
    })
  }
  expect_error(failing(), class = "latents_probe_failure")
  expect_error(failing(), class = "latents_probe_failure")
  expect_identical(failures, 2L)
})

test_that("covariate and transition inference reuse their information too", {
  clear_inference_cache()
  covariate <- quietly(multilpa(engagement_small, c("browse", "lectures"),
                                "student", n_profiles = 2, n_group_classes = 1,
                                profile_covariates = "forum_read",
                                n_starts = 1, seed = 1))
  first <- parameter_inference(covariate)
  again <- parameter_inference(covariate)
  expect_identical(again, first)
  expect_identical(length(ls(.latents_cache)), 1L)
  entries <- get("observed_information", envir = .latents_cache)
  expect_identical(entries[[1L]]$key$path, "covariates")
  transition <- quietly(lta(subset(course_engagement, student <= 30),
                            c("browse", "lectures"), "student",
                            time = "sequence", n_profiles = 2, n_starts = 1,
                            seed = 1))
  first_transition <- parameter_inference(transition)
  expect_identical(parameter_inference(transition), first_transition)
  paths <- vapply(get("observed_information", envir = .latents_cache),
                  \(entry) entry$key$path, character(1))
  expect_setequal(paths, c("covariates", "transitions"))
})
