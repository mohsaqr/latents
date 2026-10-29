empty_mask_pair_data <- function() {
  withr::local_seed(381)
  n <- 180L
  class <- rep(1:2, each = n / 2L)
  data <- data.frame(
    a = stats::rbinom(n, 1, c(.15, .85)[class]),
    b = stats::rbinom(n, 1, c(.2, .8)[class]),
    c = stats::rbinom(n, 1, c(.1, .9)[class]),
    e = stats::rbinom(n, 1, c(.2, .8)[class]))
  data[1, c("b", "c", "e")] <- NA
  data[2, c("a", "c", "e")] <- NA
  data[3, ] <- NA
  data
}

test_that("partially observed joined items remain classified", {
  data <- empty_mask_pair_data()
  fit <- quietly(lca(data, names(data), n_classes = 2,
    local_dependence = "a:b", missing = "fiml", n_starts = 2,
    seed = 19, max_iter = 300))
  expect_identical(unname(fit$informative[1:4]), c(TRUE, TRUE, FALSE, TRUE))
  expect_equal(nobs(fit), nrow(data) - 1L)
  expect_equal(sum(get_results(fit, "classification")$n_modal), nrow(data) - 1L)
  assignments <- get_results(fit, "assignments")
  expect_false(anyNA(assignments$profile[1:2]))
  expect_true(is.na(assignments$profile[3]))
  expect_equal(unname(as.matrix(assignments[1:2, c("posterior_profile_1",
                                                 "posterior_profile_2")])),
               unname(fit$subject_posteriors[1:2, ]))
  expect_equal(.multilpa_categorical_observed(fit$categorical_data)[1:4],
               c(1, 1, 0, 3))
  expect_identical(.multilpa_categorical_observed(NULL), 0L)
  ordinary <- matrix(c(1L, NA_integer_, 2L, NA_integer_), 2L)
  expect_equal(.multilpa_categorical_observed(ordinary), c(2, 0))
  # Reporting conventions must leave the likelihood and posterior unchanged.
  legacy <- quietly(lca(data, names(data), n_classes = 2,
    local_dependence = "a:b", missing = "fiml", n_starts = 2,
    seed = 19, max_iter = 300, drop_empty = FALSE))
  expect_identical(fit$log_likelihood, legacy$log_likelihood)
  expect_identical(fit$subject_posteriors, legacy$subject_posteriors)
  expect_false(anyNA(get_results(legacy, "assignments")$profile))
})

test_that("ordinal and count observations keep LTA rows informative", {
  withr::local_seed(13)
  data <- data.frame(id = rep(seq_len(24L), each = 3L), time = rep(1:3, 24L))
  data$o <- ordered(sample(1:3, nrow(data), replace = TRUE))
  data$k <- stats::rpois(nrow(data), 3)
  data[1, "o"] <- NA
  data[2, "k"] <- NA
  data[3, c("o", "k")] <- NA
  fit <- quietly(lta(data, c("o", "k"), "id", 2L, time = "time",
    ordinal = "o", count = "k", missing = "fiml", n_starts = 1L,
    max_iter = 0L, seed = 51L))
  expect_identical(unname(fit$informative[1:4]), c(TRUE, TRUE, FALSE, TRUE))
  expect_equal(fit$n_informative, nrow(data) - 1L)
  assignments <- get_results(fit, "assignments")
  expect_equal(sum(!is.na(assignments$profile)), nrow(data) - 1L)
})
