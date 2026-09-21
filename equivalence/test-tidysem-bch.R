.step_data <- function(seed = 11L, n_groups = 50L, per = 10L, separation = 1.1,
                       difference = 10) {
  set.seed(seed)
  group <- rep(seq_len(n_groups), each = per)
  n <- length(group)
  truth <- 1L + as.integer(stats::runif(n) > 0.5)
  data.frame(g = group, truth = truth,
             a = stats::rnorm(n, ifelse(truth == 2L, separation, -separation)),
             b = stats::rnorm(n, ifelse(truth == 2L, separation, -separation)),
             y = stats::rnorm(n, ifelse(truth == 2L, difference, 0)))
}

.step_fit <- function(data, ...) {
  multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
           n_starts = 6, seed = 1, ...)
}

test_that("the correction agrees with an independent implementation", {
  skip_if_not_installed("tidySEM")
  data <- .step_data()
  fit <- .step_fit(data)
  posteriors <- fit$subject_posteriors

  # tidySEM builds the same matrix from the posteriors alone, so the two
  # implementations can be compared exactly rather than approximately.
  mine <- .multilpa_error_matrix(.multilpa_level_assignments(fit, "individuals"))
  theirs <- as.matrix(tidySEM:::classification_probs_mostlikely(posteriors))
  expect_equal(unname(mine), unname(theirs))

  weights_theirs <- solve(theirs)[apply(posteriors, 1L, which.max), ]
  estimate_theirs <- vapply(seq_len(2), function(class) {
    sum(weights_theirs[, class] * data$y) / sum(weights_theirs[, class])
  }, numeric(1))
  expect_equal(three_step(fit, data, "y", method = "bch")$estimate,
               estimate_theirs)
})

