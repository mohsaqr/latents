# The two criteria tidyLPA reports that this package did not, and reading the
# measurement model complete from one call.

.release_fixture <- function(seed = 2L) {
  set.seed(seed)
  frame <- data.frame(g = rep(seq_len(20), each = 8))
  state <- 1L + (stats::runif(160) > 0.5)
  frame$a <- stats::rnorm(160, c(-2, 2)[state], 0.8)
  frame$b <- stats::rnorm(160, c(-1.5, 1.5)[state], 0.8)
  frame$q <- ifelse(stats::runif(160) < c(0.8, 0.2)[state], "no", "yes")
  frame
}

test_that("KIC and CLC are the criteria their definitions say they are", {
  data <- .release_fixture()
  fit <- quietly(multilpa(data, c("a", "b"), "g", 2L, 2L,
                                   n_starts = 3, seed = 1))
  criteria <- get_results(fit, "information_criteria", format = "long")
  q <- fit$n_parameters
  log_likelihood <- fit$log_likelihood

  kic <- criteria[criteria$criterion == "kic", ]
  expect_identical(nrow(kic), 1L)
  # "not applicable" is NA, not a sentinel level that would become a third value
  # of a two-level variable and silently survive a convention filter.
  expect_true(is.na(kic$convention))
  expect_true(is.na(kic$n))
  expect_equal(kic$value, -2 * log_likelihood + 3 * (q + 1))

  # CLC uses no sample size, but its convention picks the level whose
  # classification uncertainty it penalizes, so there is one row per level.
  clc <- criteria[criteria$criterion == "clc", ]
  expect_identical(nrow(clc), 2L)
  expect_identical(clc$convention, c("groups", "individuals"))
  entropies <- c(.multilpa_entropy_sum(fit$group_posteriors),
                 .multilpa_entropy_sum(fit$subject_posteriors))
  expect_equal(clc$value, -2 * log_likelihood + 2 * entropies)
  # It must not vary with n, unlike every other row that carries a convention.
  expect_equal(clc$value[1L] - clc$value[2L], 2 * (entropies[1L] - entropies[2L]))
  # The formula is documentation, so it is off by default and arrives on request.
  annotated <- subset(get_results(fit, "information_criteria", format = "long", definitions = TRUE),
                      criterion == "clc")
  expect_identical(annotated$definition, rep("-2L + 2 EN", 2L))
  expect_false("definition" %in% names(criteria))
})

test_that("the criteria table keeps its shape as criteria are added", {
  data <- .release_fixture()
  fit <- quietly(multilpa(data, c("a", "b"), "g", 2L, 2L,
                                   n_starts = 2, seed = 1))
  criteria <- get_results(fit, "information_criteria", format = "long")
  expect_named(criteria, c("criterion", "convention", "n", "value"))
  expect_setequal(unique(criteria$criterion),
                  c("deviance", "aic", "kic", "bic", "sabic", "caic",
                    "awe", "icl", "clc"))
  expect_false(anyNA(criteria$value))
  # Reporting the deviance rather than the log likelihood makes the `value`
  # column point one way on EVERY row, with no row carved out.
  deviance <- -2 * fit$log_likelihood
  expect_equal(subset(criteria, criterion == "deviance")$value, deviance)
  # The complexity charge is not always positive: sabic subtracts when the
  # sample size is below 22, because log((n + 2) / 24) is then negative. So a
  # criterion may fall BELOW the deviance, and the table must not assume
  # otherwise now that the charge is no longer a column of its own.
  expect_true(any(criteria$value < deviance))

  # The reported shape: one row, and the columns a comparison table publishes.
  wide <- get_results(fit, "information_criteria")
  expect_identical(nrow(wide), 1L)
  expect_identical(names(wide)[1:2], c("log_likelihood", "n_parameters"))
  expect_false("penalty" %in% names(wide))
  expect_equal(wide$log_likelihood, fit$log_likelihood)
})

test_that("the measurement tables carry their own standard errors when asked", {
  data <- .release_fixture()
  fit <- quietly(multilpa(data, c("a", "b", "q"), "g", 2L, 1L,
    categorical = "q", n_starts = 3, seed = 1))
  bare <- as.data.frame(fit)
  with_errors <- get_results(fit, "profiles", data = data)
  # Adding errors must not change the estimates or the shape.
  expect_identical(nrow(bare), nrow(with_errors))
  expect_equal(bare$mean, with_errors$mean)
  expect_named(with_errors, c(names(bare), "mean_standard_error",
                              "variance_standard_error"))
  expect_false(anyNA(with_errors$mean_standard_error))

  # The errors must be the ones parameter_inference() reports, put in the
  # shape of the measurement table rather than recomputed.
  inference <- parameter_inference(fit, data)
  means <- inference[inference$parameter == "mean", ]
  key <- paste(sub("^profile_", "", means$outcome), means$term)
  expect_equal(with_errors$mean_standard_error,
               means$standard_error[match(paste(with_errors$profile,
                                                with_errors$indicator), key)])

  responses <- get_results(fit, "responses", data = data)
  expect_true("probability_standard_error" %in% names(responses))
  expect_false(anyNA(responses$probability_standard_error))
  reported <- inference[inference$parameter == "response", ]
  response_key <- paste(sub("^profile_", "", reported$outcome), reported$term)
  expect_equal(responses$probability_standard_error,
               reported$standard_error[match(
                 paste(responses$profile,
                       paste(responses$indicator, responses$category, sep = ":")),
                 response_key)])
})

test_that("asking for errors a fit cannot supply raises rather than returning blanks", {
  data <- .release_fixture()
  transitions_data <- data
  transitions_data$t <- rep(seq_len(8), times = 20)
  transition_fit <- quietly(lta(transitions_data,
    c("a", "b"), "g", n_profiles = 2L, time = "t", n_starts = 2, seed = 1))
  # The transition family has no standard errors, so the measurement table
  # must refuse the request instead of filling the columns with NA.
  expect_error(get_results(transition_fit, "profiles", data = data),
               class = "multilpa_no_inference")
  expect_s3_class(get_results(transition_fit, "profiles"), "data.frame")
})
