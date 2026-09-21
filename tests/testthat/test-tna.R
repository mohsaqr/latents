# Handing a transition fit to tna. What is tested is the quantity handed over --
# that the aggregate is the marginal transition matrix and not an average of the
# class matrices -- and that the objects are the ones tna's own verbs accept.

transition_fit <- function() {
  activity <- c("browse", "lectures", "forum_read")
  quietly(fit_transitions(course_engagement, vars = activity, id = "student",
                          time = "sequence", n_profiles = 3,
                          n_group_classes = 2, n_starts = 3, max_iter = 1000,
                          seed = 1))
}

test_that("the aggregate is the marginal, not the mean of the class matrices", {
  fit <- transition_fit()
  aggregate <- multilpa:::.multilpa_aggregate_transitions(fit)

  # Every row of a transition matrix is a distribution.
  expect_equal(unname(rowSums(aggregate$probabilities)),
               rep(1, fit$n_profiles))
  expect_equal(sum(aggregate$initial), 1)

  # The marginal, computed independently here: pool the expected counts over
  # classes, then normalise. A class holding few students but many moves counts
  # for its moves.
  counts <- apply(fit$transition_counts, seq_len(2L), sum)
  expect_equal(unname(aggregate$probabilities),
               unname(counts / rowSums(counts)))

  # And it is not the class-probability-weighted average of the matrices, which
  # is the quantity it would be easy to hand over by mistake.
  averaged <- Reduce(`+`, lapply(seq_len(fit$n_group_classes), function(class) {
    fit$group_probabilities[class] * fit$transition_probabilities[, , class]
  }))
  expect_false(isTRUE(all.equal(unname(aggregate$probabilities), unname(averaged))))
})

test_that("get_tna returns a tna model of the whole sample", {
  skip_if_not_installed("tna")
  fit <- transition_fit()
  model <- get_tna(fit)
  expect_s3_class(model, "tna")
  expect_identical(dim(model$weights), c(3L, 3L))
  expect_equal(unname(rowSums(model$weights)), rep(1, 3))
  expect_equal(sum(model$inits), 1)
  # The states keep the fit's own names, so a network and a profile table can be
  # read against each other.
  expect_identical(model$labels, multilpa:::.multilpa_state_labels(fit))
  # tna's verbs accept it.
  expect_identical(nrow(tna::centralities(model)), 3L)
})

test_that("get_group_tna returns one model per latent class", {
  skip_if_not_installed("tna")
  fit <- transition_fit()
  grouped <- get_group_tna(fit)
  expect_s3_class(grouped, "group_tna")
  expect_length(grouped, fit$n_group_classes)
  expect_identical(names(grouped), sprintf("Group class %d", 1:2))
  invisible(lapply(grouped, function(model) expect_s3_class(model, "tna")))

  # Each class carries its own matrix, unaveraged.
  expect_equal(unname(grouped[[1L]]$weights),
               unname(fit$transition_probabilities[, , 1L]))
  expect_equal(unname(grouped[[2L]]$weights),
               unname(fit$transition_probabilities[, , 2L]))

  # tna's grouped verbs read it: one tidy table with a group column.
  centrality <- tna::centralities(grouped)
  expect_true("group" %in% names(centrality))
  expect_identical(nrow(centrality), 6L)
  expect_identical(length(unique(centrality$group)), 2L)
})

test_that("the class label is the caller's", {
  skip_if_not_installed("tna")
  grouped <- get_group_tna(transition_fit(), label = "Trajectory")
  expect_identical(names(grouped), sprintf("Trajectory %d", 1:2))
  expect_identical(attr(grouped, "label"), "Trajectory")
  expect_error(get_group_tna(transition_fit(), label = c("a", "b")),
               "single string")
})

test_that("a state nothing ever leaves is a self-transition, not NaN", {
  # Normalising an all-zero row divides by zero. "Never left" is a probability
  # one of staying, which every centrality can read; NaN is not.
  fit <- transition_fit()
  fit$transition_counts[2L, , ] <- 0
  aggregate <- multilpa:::.multilpa_aggregate_transitions(fit)
  expect_false(anyNA(aggregate$probabilities))
  expect_equal(unname(aggregate$probabilities[2L, ]), c(0, 1, 0))
  expect_equal(unname(rowSums(aggregate$probabilities)), rep(1, 3))
})
