.categorical_indicators <- function() c("u1", "u2", "u3", "u4", "u5")

.categorical_fit <- function(data, ...) {
  multilpa(data, .categorical_indicators(), "clus", n_profiles = 2,
           n_group_classes = 2, categorical = .categorical_indicators(),
           n_starts = 8, seed = 4, ...)
}

test_that("categorical coordinates round-trip through encode and decode", {
  fixture <- readRDS(test_path("..", "fixtures", "mplus", "twolevel-categorical.rds"))
  fit <- .categorical_fit(fixture$data)
  theta <- multilpa:::.multilpa_coefficients(fit, "unconstrained")
  rebuilt <- multilpa:::.multilpa_decode(theta, fit)

  expect_length(theta, fit$n_parameters)
  expect_equal(rebuilt$response_probabilities, fit$response_probabilities)
  expect_equal(unname(rebuilt$profile_probabilities),
               unname(fit$profile_probabilities))
  expect_equal(unname(rebuilt$group_probabilities),
               unname(fit$group_probabilities))
  # every response block is a set of simplexes
  expect_true(all(vapply(rebuilt$response_probabilities,
                         function(block) max(abs(rowSums(block) - 1)) < 1e-12,
                         logical(1))))
})

test_that("the categorical score agrees with a numerical gradient", {
  fixture <- readRDS(test_path("..", "fixtures", "mplus", "twolevel-categorical.rds"))
  fit <- .categorical_fit(fixture$data)
  codes <- multilpa:::.multilpa_encode_categorical(
    fixture$data[, .categorical_indicators(), drop = FALSE])$codes
  x <- matrix(numeric(0), nrow(fixture$data), 0L)
  theta <- multilpa:::.multilpa_coefficients(fit, "unconstrained")
  likelihood <- function(parameters) {
    multilpa:::.multilpa_expectation(x, fit$group_index,
      multilpa:::.multilpa_decode(parameters, fit), codes)$log_likelihood
  }
  analytic <- unname(-multilpa:::.multilpa_score(theta, x, fit, codes))
  numerical <- vapply(seq_along(theta), function(j) {
    step <- 1e-5 * max(1, abs(theta[j]))
    up <- theta; up[j] <- up[j] + step
    down <- theta; down[j] <- down[j] - step
    (likelihood(up) - likelihood(down)) / (2 * step)
  }, numeric(1))

  expect_equal(likelihood(theta), fit$log_likelihood)
  expect_equal(analytic, numerical, tolerance = 1e-4)
  expect_lt(max(abs(analytic)), 0.01)
})

test_that("the tidy table carries the response block at the measurement level", {
  fixture <- readRDS(test_path("..", "fixtures", "mplus", "twolevel-categorical.rds"))
  fit <- .categorical_fit(fixture$data)
  inference <- parameter_inference(fit, fixture$data)
  responses <- subset(inference, parameter == "response")

  expect_true("response" %in% inference$parameter)
  expect_true(all(responses$level == "measurement"))
  # one row per profile, indicator and category
  expect_equal(nrow(responses), 2L * 5L * 2L)
  expect_true(all(grepl("^u[0-9]:", responses$term)))
  expect_setequal(unique(responses$outcome), c("profile_1", "profile_2"))
  # probabilities lie inside the unit interval and their errors are usable
  expect_true(all(responses$estimate > 0 & responses$estimate < 1))
  expect_true(all(is.finite(responses$standard_error) & responses$standard_error > 0))
  # each profile-by-indicator simplex still sums to one
  sums <- tapply(responses$estimate, paste(responses$outcome, responses$term),
                 function(v) v)
  expect_equal(sum(responses$estimate), 2 * 5)
})

test_that("mixed continuous and categorical measurement is supported", {
  set.seed(9)
  n_groups <- 24L
  group <- rep(seq_len(n_groups), each = 10L)
  profile <- 1L + as.integer(runif(length(group)) > 0.5)
  data <- data.frame(
    clus = group,
    y = rnorm(length(group), ifelse(profile == 2L, 2, -2), 0.8),
    u = ifelse(runif(length(group)) < ifelse(profile == 2L, 0.8, 0.2), 1L, 0L))
  fit <- multilpa(data, c("y", "u"), "clus", n_profiles = 2, n_group_classes = 1,
                  categorical = "u", n_starts = 5, seed = 2)
  inference <- parameter_inference(fit, data)

  expect_identical(fit$measurement_model, "mixed")
  expect_setequal(unique(subset(inference, level == "measurement")$parameter),
                  c("mean", "variance", "response"))
  expect_true(all(subset(inference, parameter == "mean")$term == "y"))
  expect_true(all(grepl("^u:", subset(inference, parameter == "response")$term)))
  # With one group class the mixing weight is fixed at 1 and carries no
  # uncertainty; every estimated parameter does.
  fixed <- inference$level == "group" & inference$estimate == 1
  expect_equal(inference$standard_error[fixed], 0)
  expect_true(all(inference$standard_error[!fixed] > 0))
})

test_that("the wrong categorical data is refused", {
  fixture <- readRDS(test_path("..", "fixtures", "mplus", "twolevel-categorical.rds"))
  fit <- .categorical_fit(fixture$data)
  altered <- fixture$data
  altered$u1 <- 1L - altered$u1

  expect_error(parameter_inference(fit, altered),
               class = "multilpa_bad_inference_data")
})

test_that("the robust sandwich is available for categorical measurement", {
  fixture <- readRDS(test_path("..", "fixtures", "mplus", "twolevel-categorical.rds"))
  fit <- .categorical_fit(fixture$data)
  observed <- parameter_inference(fit, fixture$data)
  robust <- parameter_inference(fit, fixture$data, vcov_type = "robust")

  expect_equal(robust$estimate, observed$estimate)
  expect_true(all(robust$standard_error > 0))
  expect_gt(attr(robust, "scaling_correction"), 0)
  expect_identical(attr(robust, "vcov_type"), "robust")
})
