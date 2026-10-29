# Comparisons with independent references moved from tests/testthat/test-transitions.R,
# so that they run in CI's equivalence job rather than under R CMD check.

.transition_fixture <- function(n_groups = 9L, occasions = 4L, seed = 99L) {
  set.seed(seed)
  frame <- data.frame(g = rep(seq_len(n_groups), each = occasions),
                      t = rep(seq_len(occasions), times = n_groups))
  frame$y1 <- stats::rnorm(nrow(frame))
  frame$y2 <- stats::rnorm(nrow(frame))
  frame$c1 <- sample(c("lo", "hi"), nrow(frame), replace = TRUE)
  frame
}

# Sums the complete-data likelihood over every state path of every group and
# every group class. Shares only the model definition with the recursion.
.transition_enumerate <- function(continuous, codes, group_index, time, fit) {
  n_profiles <- fit$n_profiles
  n_rows <- length(group_index)
  density <- vapply(seq_len(n_profiles), function(profile) {
    gaussian <- if (ncol(continuous) == 0L) numeric(n_rows) else
      rowSums(vapply(seq_len(ncol(continuous)), function(indicator) {
        value <- stats::dnorm(continuous[, indicator], fit$means[profile, indicator],
                              sqrt(fit$variances[profile, indicator]), log = TRUE)
        ifelse(is.na(value), 0, value)
      }, numeric(n_rows)))
    categorical <- if (is.null(codes)) numeric(n_rows) else
      rowSums(vapply(seq_along(fit$response_probabilities), function(indicator) {
        value <- log(fit$response_probabilities[[indicator]][profile, codes[, indicator]])
        ifelse(is.na(value), 0, value)
      }, numeric(n_rows)))
    gaussian + categorical
  }, numeric(n_rows))
  sum(vapply(split(seq_along(group_index), group_index), function(rows) {
    rows <- rows[order(time[rows])]
    paths <- as.matrix(expand.grid(rep(list(seq_len(n_profiles)), length(rows))))
    log(sum(vapply(seq_len(fit$n_group_classes), function(type) {
      initial <- fit$initial_probabilities[type, ]
      transition <- matrix(fit$transition_probabilities[, , type],
                           n_profiles, n_profiles)
      fit$group_probabilities[type] * sum(apply(paths, 1L, function(path) {
        moves <- if (length(path) < 2L) 1 else
          prod(transition[cbind(path[-length(path)], path[-1L])])
        initial[path[1L]] * moves * exp(sum(density[cbind(rows, path)]))
      }))
    }, numeric(1))))
  }, numeric(1)))
}

test_that("the recursion reproduces the likelihood enumeration gives", {
  data <- .transition_fixture()
  continuous <- as.matrix(data[, c("y1", "y2")])
  cases <- list(
    list(profiles = 2L, types = 1L), list(profiles = 3L, types = 1L),
    list(profiles = 2L, types = 2L), list(profiles = 3L, types = 2L))
  invisible(lapply(cases, function(case) {
    fit <- quietly(lta(
      data, c("y1", "y2"), "g", n_profiles = case$profiles, time = "t",
      n_group_classes = case$types, n_starts = 1, seed = 4, max_iter = 25))
    expect_equal(fit$log_likelihood,
                 .transition_enumerate(continuous, NULL, fit$group_index,
                                       data$t, fit),
                 tolerance = 1e-10,
                 label = sprintf("%d profiles, %d group classes",
                                 case$profiles, case$types))
  }))
})

test_that("enumeration agrees for categorical, mixed and incomplete indicators", {
  data <- .transition_fixture()
  mixed <- quietly(lta(
    data, c("y1", "y2", "c1"), "g", n_profiles = 2L, time = "t",
    categorical = "c1", n_starts = 1, seed = 8, max_iter = 25))
  expect_equal(mixed$log_likelihood,
               .transition_enumerate(as.matrix(data[, c("y1", "y2")]),
                                     mixed$categorical_data, mixed$group_index,
                                     data$t, mixed), tolerance = 1e-10)

  pure <- quietly(lta(
    data, "c1", "g", n_profiles = 2L, time = "t", categorical = "c1",
    n_starts = 1, seed = 8, max_iter = 25))
  expect_equal(pure$log_likelihood,
               .transition_enumerate(matrix(numeric(0), nrow(data), 0L),
                                     pure$categorical_data, pure$group_index,
                                     data$t, pure), tolerance = 1e-10)
  expect_identical(pure$measurement_model, "categorical")

  incomplete <- data
  incomplete$y1[c(2L, 9L, 20L)] <- NA
  incomplete$y2[c(5L, 9L)] <- NA
  fiml <- quietly(lta(
    incomplete, c("y1", "y2", "c1"), "g", n_profiles = 2L, time = "t",
    categorical = "c1", missing = "fiml", n_starts = 1, seed = 8, max_iter = 25))
  expect_equal(fiml$log_likelihood,
               .transition_enumerate(as.matrix(incomplete[, c("y1", "y2")]),
                                     fiml$categorical_data, fiml$group_index,
                                     incomplete$t, fiml), tolerance = 1e-10)
})

test_that("groups of one occasion contribute an initial state and no transition", {
  data <- .transition_fixture()
  short <- rbind(data, data.frame(g = 10L, t = 1L, y1 = 0.4, y2 = -0.2,
                                  c1 = "lo"))
  fit <- quietly(lta(
    short, c("y1", "y2"), "g", n_profiles = 2L, time = "t", n_starts = 1,
    seed = 7, max_iter = 25))
  expect_identical(unname(fit$sequence_lengths[10L]), 1L)
  expect_equal(sum(fit$transition_counts), sum(fit$sequence_lengths - 1L))
  expect_equal(fit$log_likelihood,
               .transition_enumerate(as.matrix(short[, c("y1", "y2")]), NULL,
                                     fit$group_index, short$t, fit),
               tolerance = 1e-10)
})
