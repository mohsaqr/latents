# Comparisons with independent references moved from tests/testthat/test-transitions-numerics.R,
# so that they run in CI's equivalence job rather than under R CMD check.

# Expected transition counts from first principles: every state path of every
# group, weighted by its own complete-data likelihood. Shares only the model
# definition with the recursion, so the recursion cannot verify itself.
.enumerate_transition_counts <- function(continuous, group_index, time, fit) {
  n_profiles <- fit$n_profiles
  n_rows <- length(group_index)
  density <- vapply(seq_len(n_profiles), function(profile) {
    rowSums(vapply(seq_len(ncol(continuous)), function(indicator) {
      stats::dnorm(continuous[, indicator], fit$means[profile, indicator],
                   sqrt(fit$variances[profile, indicator]), log = TRUE)
    }, numeric(n_rows)))
  }, numeric(n_rows))
  per_group <- lapply(split(seq_along(group_index), group_index), function(rows) {
    rows <- rows[order(time[rows])]
    paths <- as.matrix(expand.grid(rep(list(seq_len(n_profiles)), length(rows))))
    moves <- lapply(seq_len(nrow(paths)), function(path_index) {
      path <- paths[path_index, ]
      pairs <- cbind(path[-length(path)], path[-1L])
      tabulate((pairs[, 2L] - 1L) * n_profiles + pairs[, 1L],
               nbins = n_profiles * n_profiles)
    })
    weights <- lapply(seq_len(fit$n_group_classes), function(type) {
      initial <- fit$initial_probabilities[type, ]
      transition <- matrix(fit$transition_probabilities[, , type],
                           n_profiles, n_profiles)
      vapply(seq_len(nrow(paths)), function(path_index) {
        path <- paths[path_index, ]
        initial[path[1L]] *
          prod(transition[cbind(path[-length(path)], path[-1L])]) *
          exp(sum(density[cbind(rows, path)])) * fit$group_probabilities[type]
      }, numeric(1))
    })
    total <- sum(unlist(weights, use.names = FALSE))
    lapply(weights, function(weight) {
      matrix(Reduce(`+`, Map(function(share, move) (share / total) * move,
                             weight, moves)),
             n_profiles, n_profiles)
    })
  })
  lapply(seq_len(fit$n_group_classes), function(type) {
    Reduce(`+`, lapply(per_group, `[[`, type))
  })
}

.numerics_fixture <- function(n_groups = 9L, occasions = 4L, seed = 99L) {
  set.seed(seed)
  frame <- data.frame(g = rep(seq_len(n_groups), each = occasions),
                      t = rep(seq_len(occasions), times = n_groups))
  frame$y1 <- stats::rnorm(nrow(frame))
  frame$y2 <- stats::rnorm(nrow(frame))
  frame
}

test_that("expected counts match explicit path enumeration", {
  data <- .numerics_fixture()
  continuous <- as.matrix(data[, c("y1", "y2")])
  cases <- list(list(profiles = 2L, types = 1L),
                list(profiles = 3L, types = 1L),
                list(profiles = 2L, types = 2L))
  deviations <- vapply(cases, function(case) {
    fit <- quietly(lta(
      data, c("y1", "y2"), "g", n_profiles = case$profiles, time = "t",
      n_group_classes = case$types, n_starts = 1, seed = 4, max_iter = 25))
    reference <- .enumerate_transition_counts(continuous, fit$group_index,
                                              data$t, fit)
    recursion <- lapply(seq_len(case$types), function(type) {
      matrix(fit$transition_counts[, , type], case$profiles, case$profiles)
    })
    max(abs(unlist(recursion, use.names = FALSE) -
              unlist(reference, use.names = FALSE)))
  }, numeric(1))
  # Measured at 1.8e-15 on counts whose largest entry is about twenty, so the
  # bound below is the observed deviation with an order of magnitude of room.
  expect_lt(max(deviations), 2e-14)
})
