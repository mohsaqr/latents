# Numerics of the latent transition recursions.
#
# The expected transition counts are a second reported quantity of the E-step,
# and the log likelihood does not police them: the likelihood is read off the
# scaled forward pass, while the counts also read the unscaled backward pass.
# A long sequence drives the two apart, so these tests check the counts on
# their own -- against explicit path enumeration, against their own invariants,
# and at a length that used to overflow.

# One group, one indicator, perfectly alternating far-apart states measured
# with tiny variance. Every log density is either zero or -5000, which is what
# drives the two recursions apart; the deterministic path makes the answer
# known in advance.
.alternating_expectation <- function(n_occasions) {
  x <- matrix(rep(c(-5, 5), length.out = n_occasions), n_occasions, 1L)
  layout <- .multilpa_sequence_layout(rep(1L, n_occasions),
                                      seq_len(n_occasions), 1L)
  parameters <- list(
    means = matrix(c(-5, 5), 2L, 1L), variances = matrix(0.01, 2L, 1L),
    group_probabilities = 1, initial_probabilities = matrix(c(0.5, 0.5), 1L, 2L),
    transition_probabilities = array(c(0.85, 0.15, 0.15, 0.85), c(2L, 2L, 1L)))
  .multilpa_transition_expectation(x, layout, parameters)
}

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

test_that("a long sequence returns finite, nonnegative expected counts", {
  # 600 occasions of this fixture overflowed the factored form, which
  # exponentiated the forward and backward terms apart: the forward term
  # reached 1e248 where the backward term had underflowed to zero, so every
  # count came back NaN from a likelihood that was itself correct. The
  # factored form first loses a count near 380 occasions on this fixture.
  expectation <- .alternating_expectation(600L)
  counts <- expectation$sequence[[1L]]$transition
  expect_true(all(is.finite(counts)))
  expect_true(all(counts >= 0))
  # The log likelihood was never wrong and must not move.
  expect_equal(expectation$log_likelihood, -306.8801, tolerance = 1e-6)
  # The states alternate deterministically, so 300 of the 599 moves leave the
  # first state and 299 leave the second, and none is a stay.
  expect_equal(counts, matrix(c(0, 299, 300, 0), 2L, 2L), tolerance = 1e-8)
})

test_that("expected counts and posteriors satisfy their own invariants", {
  lengths <- c(2L, 7L, 50L, 400L)
  invisible(lapply(lengths, function(n_occasions) {
    expectation <- .alternating_expectation(n_occasions)
    counts <- expectation$sequence[[1L]]$transition
    # One group of `n_occasions` occasions makes exactly `n_occasions - 1`
    # within-sequence transitions, and the pair posteriors sum to one at each.
    expect_equal(sum(counts), n_occasions - 1,
                 tolerance = 1e-8, label = sprintf("total at %d", n_occasions))
    expect_equal(sum(expectation$sequence[[1L]]$initial), 1, tolerance = 1e-10)
    expect_equal(unname(rowSums(expectation$subject_posteriors)),
                 rep(1, n_occasions), tolerance = 1e-10,
                 label = sprintf("posterior rows at %d", n_occasions))
  }))
})

test_that("a fitted model's counts total its within-sequence transitions", {
  data <- .numerics_fixture(n_groups = 12L, occasions = 5L, seed = 11L)
  fit <- quietly(fit_transitions(
    data, c("y1", "y2"), "g", n_profiles = 2L, time = "t",
    n_group_classes = 2L, n_starts = 2, seed = 3, max_iter = 40))
  expect_true(all(is.finite(fit$transition_counts)))
  expect_true(all(fit$transition_counts >= 0))
  expect_equal(sum(fit$transition_counts), 12 * (5 - 1), tolerance = 1e-8)
  expect_equal(unname(rowSums(fit$subject_posteriors)),
               rep(1, nrow(data)), tolerance = 1e-10)
  expect_equal(unname(rowSums(fit$group_posteriors)), rep(1, 12),
               tolerance = 1e-10)
})

test_that("expected counts match explicit path enumeration", {
  data <- .numerics_fixture()
  continuous <- as.matrix(data[, c("y1", "y2")])
  cases <- list(list(profiles = 2L, types = 1L),
                list(profiles = 3L, types = 1L),
                list(profiles = 2L, types = 2L))
  deviations <- vapply(cases, function(case) {
    fit <- quietly(fit_transitions(
      data, c("y1", "y2"), "g", n_profiles = case$profiles, time = "t",
      n_group_classes = case$types, n_starts = 2, seed = 4, max_iter = 25))
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

test_that("zero iterations evaluate the start instead of failing to assemble", {
  set.seed(1)
  data <- data.frame(id = rep(seq_len(3), each = 8),
                     time = rep(seq_len(8), times = 3),
                     y = stats::rnorm(24))
  fit <- quietly(fit_transitions(
    data, "y", "id", n_profiles = 2L, time = "time", n_starts = 1,
    max_iter = 0, seed = 1))
  expect_s3_class(fit, "multilpa_transitions")
  expect_identical(fit$iterations, 0L)
  expect_false(fit$converged)
  # The prevalence is a summary of the reported expectation, so it exists and
  # is a distribution even when no M-step ever ran.
  prevalence <- as.data.frame(fit, what = "initial")
  expect_true(all(is.finite(prevalence$prevalence)))
  expect_equal(sum(prevalence$prevalence), 1, tolerance = 1e-10)
  expect_true(all(is.finite(transitions(fit)$expected_count)))
  expect_s3_class(summary(fit), "summary_multilpa_transitions")
})

test_that("the prevalence reported is the one the final expectation implies", {
  data <- .numerics_fixture(n_groups = 10L, occasions = 4L, seed = 21L)
  # With one group class every group belongs to it with probability one, so
  # the class's joint posteriors are the individual posteriors themselves and
  # the implied prevalence is their normalised column sums.
  single <- quietly(fit_transitions(
    data, c("y1", "y2"), "g", n_profiles = 2L, time = "t",
    n_starts = 2, seed = 5, max_iter = 60))
  shares <- colSums(single$subject_posteriors)
  expect_equal(unname(drop(single$profile_prevalence)),
               unname(shares / sum(shares)), tolerance = 1e-10)

  nested <- quietly(fit_transitions(
    data, c("y1", "y2"), "g", n_profiles = 2L, time = "t",
    n_group_classes = 2L, n_starts = 2, seed = 5, max_iter = 60))
  expect_true(all(is.finite(nested$profile_prevalence)))
  expect_equal(unname(rowSums(nested$profile_prevalence)), rep(1, 2),
               tolerance = 1e-10)
  # The prevalence is a class-conditional distribution, so on these equal-length
  # groups the class-weighted average of it is the overall profile share. The
  # tolerance is loose because only the prevalence is read from the reported
  # expectation: `group_probabilities` still comes from the last M-step, which
  # saw the previous one, so the two agree only to the EM tolerance.
  overall <- colSums(nested$subject_posteriors) / nrow(data)
  expect_equal(unname(drop(nested$group_probabilities %*%
                             nested$profile_prevalence)),
               unname(overall), tolerance = 1e-4)
})

test_that("two very long sequences fit rather than overflowing", {
  skip_on_cran()
  set.seed(1)
  n_occasions <- 1800L
  data <- data.frame(
    id = rep(seq_len(2), each = n_occasions),
    time = rep(seq_len(n_occasions), times = 2),
    y = rep(rep(c(-5, 5), length.out = n_occasions), times = 2) +
      stats::rnorm(2 * n_occasions, sd = 0.1))
  fit <- quietly(fit_transitions(
    data, "y", "id", n_profiles = 2L, time = "time", n_starts = 1,
    max_iter = 2, seed = 1))
  expect_s3_class(fit, "multilpa_transitions")
  expect_true(is.finite(fit$log_likelihood))
  moves <- transitions(fit)
  expect_true(all(is.finite(moves$expected_count)))
  expect_true(all(moves$expected_count >= 0))
  expect_equal(sum(moves$expected_count), 2 * (n_occasions - 1),
               tolerance = 1e-6)
})

test_that("broken contracts of the moment step raise by class", {
  data <- .numerics_fixture(n_groups = 6L, occasions = 3L, seed = 31L)
  fit <- quietly(fit_transitions(
    data, c("y1", "y2"), "g", n_profiles = 2L, time = "t", n_starts = 2,
    seed = 6, max_iter = 20))
  # A single profile has nothing to move between, and a fit that carries no
  # ordering has no sequence to read: both are classed, not message-matched.
  expect_error(fit_transitions(data, c("y1", "y2"), "g", n_profiles = 1L,
                               time = "t", n_starts = 1, seed = 6),
               class = "multilpa_bad_transition")
  expect_error(parameter_inference(fit), class = "multilpa_no_inference")
  expect_error(plot(fit), class = "multilpa_no_plot")
  # A group class with no effective membership leaves its profile prevalence
  # undefined, and is refused rather than divided by zero.
  expect_error(
    .multilpa_transition_prevalence(list(joint = list(matrix(0, 4L, 2L)))),
    class = "multilpa_empty_profile")
  # The moment step states its own contract rather than trusting the caller.
  expect_error(.multilpa_sequence_moments(list(alpha = list(), beta = list(),
                                               log_scaled = 0),
                                          list(), list(slot = matrix(1, 2L, 1L)),
                                          matrix(c(0.5, 0.5, 0.5, 0.5), 2L, 2L),
                                          weights = 1, n_observations = 2L),
               "one group-class posterior per group")
})
