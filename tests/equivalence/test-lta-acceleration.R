# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Comparisons with independent references moved from tests/testthat/test-lta-acceleration.R,
# so that they run in CI's equivalence job rather than under R CMD check.
# Reference: the original pre-scaling log-space recursion and its moments,
# kept in helper-lta-log-reference.R (base R). Package side: the internal
# scaled forward-backward kernels and expected counts behind lta(); they have
# no exported accessor.

# One row per number of a matrix (or vector): its row, its column, its value,
# with the labels passed in `...` as further key columns.
acceleration_cells <- function(values, ...) {
  values <- as.matrix(values)
  data.frame(..., row = as.vector(row(values)), column = as.vector(col(values)),
             value = as.vector(values))
}

# The forward (alpha) and backward (beta) log messages of a pass, one row per
# case, message, occasion, group and profile.
acceleration_messages <- function(pass, case) {
  do.call(rbind, c(
    lapply(seq_along(pass$alpha), \(t) {
      acceleration_cells(pass$alpha[[t]], case = case, message = "alpha", occasion = t)
    }),
    lapply(seq_along(pass$beta), \(t) {
      acceleration_cells(pass$beta[[t]], case = case, message = "beta", occasion = t)
    })))
}

# The expected occupancies and counts, one row per case, moment and cell.
acceleration_moments <- function(moments, case) {
  do.call(rbind, lapply(c("posterior", "initial", "transition"), \(moment) {
    acceleration_cells(moments[[moment]], case = case, moment = moment)
  }))
}

acceleration_note <- paste(
  "Package side: the internal scaled forward-backward kernels",
  "(.multilpa_forward_backward, .multilpa_sequence_moments) behind lta();",
  "reference: the pre-scaling log-space recursion of helper-lta-log-reference.R.")

# Scaled recursions are compared with the original log-space implementation,
# including its messages: inference still consumes alpha and beta in log space.
test_that("scaled homogeneous messages and counts match the log recursion", {
  set.seed(912)
  index <- rep(seq_len(4L), c(1L, 3L, 5L, 4L))
  time <- c(1, 1, 3, 5, 1:5, 1, 2, 4, 5)
  layout <- latents:::.multilpa_sequence_layout(index, time, 4L, "grid")
  emission <- latents:::.multilpa_sequence_emission(
    matrix(rnorm(length(index) * 3), ncol = 3), layout)$log_density
  weights <- cbind(c(0.2, 0.6, 0.7, 0.3), c(0.8, 0.4, 0.3, 0.7))
  # Two random parameter sets, drawn in the original order.
  cases <- lapply(seq_len(2L), function(h) {
    initial <- prop.table(runif(3))
    transition <- matrix(runif(9), 3)
    transition <- transition / rowSums(transition)
    log_initial <- matrix(log(initial), 4, 3, byrow = TRUE)
    arrays <- rep(list(array(rep(log(transition), each = 4), c(4, 3, 3))), 4)
    old <- .lta_log_reference(emission, layout, log_initial, arrays)
    new <- latents:::.multilpa_forward_backward(emission, layout, initial, transition)
    # Also exercise the shared kernel's homogeneous matrix interface.
    shared <- latents:::.latents_forward_backward(emission, layout, log_initial, log(transition))
    list(old = old, new = new, shared = shared,
         expected = .lta_moments_reference(old, emission, layout, transition,
                                           weights[, h], length(index)),
         actual = latents:::.multilpa_sequence_moments(new, emission, layout, transition,
                                                       weights[, h], length(index)))
  })

  # Unit checks: the two package interfaces agree, and the moments add up.
  invisible(lapply(seq_len(2L), function(h) {
    expect_equal(cases[[h]]$shared, cases[[h]]$new, tolerance = 1e-12)
    actual <- cases[[h]]$actual
    expect_equal(rowSums(actual$posterior), weights[index, h], tolerance = 1e-12)
    expect_equal(sum(actual$initial), sum(weights[, h]), tolerance = 1e-12)
    expect_equal(sum(actual$transition), sum(weights[, h] * (layout$span - 1)),
                 tolerance = 1e-12)
  }))
  expect_equal(rowSums(Reduce(`+`, lapply(cases, \(case) case$actual$posterior))),
               rep(1, length(index)), tolerance = 1e-12)

  # The log likelihood of each group, both parameter sets.
  expect_reference(
    do.call(rbind, lapply(seq_len(2L), \(h) acceleration_cells(cases[[h]]$new$log_scaled, case = h))),
    do.call(rbind, lapply(seq_len(2L), \(h) acceleration_cells(cases[[h]]$old$log_scaled, case = h))),
    reference = "independent: log-space forward-backward recursion (helper-lta-log-reference.R)",
    covers = "lta()", quantity = "group log likelihood, homogeneous kernel",
    tolerance = 1e-12, measure = "relative", key = c("case", "row", "column"),
    note = acceleration_note)
  # Log messages are exactly 0 where a group's sequence has ended (beta at and
  # after its last occasion) and pass near zero elsewhere: a relative measure
  # would judge rounding noise of those entries, so the difference is scaled
  # by max(1, |reference|).
  expect_reference(
    do.call(rbind, lapply(seq_len(2L), \(h) acceleration_messages(cases[[h]]$new, h))),
    do.call(rbind, lapply(seq_len(2L), \(h) acceleration_messages(cases[[h]]$old, h))),
    reference = "independent: log-space forward-backward recursion (helper-lta-log-reference.R)",
    covers = "lta()", quantity = "log messages, homogeneous kernel",
    tolerance = 1e-12, measure = "scaled",
    key = c("case", "message", "occasion", "row", "column"), note = acceleration_note)
  # Posterior occupancies, initial and transition expected counts.
  expect_reference(
    do.call(rbind, lapply(seq_len(2L), \(h) acceleration_moments(cases[[h]]$actual, h))),
    do.call(rbind, lapply(seq_len(2L), \(h) acceleration_moments(cases[[h]]$expected, h))),
    reference = "independent: log-space forward-backward recursion (helper-lta-log-reference.R)",
    covers = "lta()", quantity = "posteriors and expected counts, homogeneous kernel",
    tolerance = 1e-12, measure = "relative", key = c("case", "moment", "row", "column"),
    note = acceleration_note)
})

test_that("the general per-group per-occasion kernel preserves log messages", {
  set.seed(715)
  layout <- latents:::.multilpa_sequence_layout(rep(1:3, c(2, 4, 3)),
                                                c(1:2, 1:4, 1:3), 3L)
  emission <- lapply(1:4, function(t) matrix(-rexp(9), 3, 3))
  initial <- matrix(runif(9), 3)
  initial <- log(initial / rowSums(initial))
  transition <- lapply(1:3, function(t) {
    p <- array(runif(27), c(3, 3, 3))
    log(sweep(p, c(1, 2), apply(p, c(1, 2), sum), "/"))
  })
  old <- .lta_log_reference(emission, layout, initial, transition)
  new <- latents:::.latents_forward_backward(emission, layout, initial, transition)
  general_note <- paste(
    "Package side: the internal general forward-backward kernel (.latents_forward_backward)",
    "behind lta() with occasion-varying or covariate transitions; reference: the pre-scaling",
    "log-space recursion of helper-lta-log-reference.R.")
  expect_reference(acceleration_cells(new$log_scaled), acceleration_cells(old$log_scaled),
                   reference = "independent: log-space forward-backward recursion (helper-lta-log-reference.R)",
                   covers = "lta()", quantity = "group log likelihood, general kernel",
                   tolerance = 1e-12, measure = "relative", key = c("row", "column"),
                   note = general_note)
  # Scaled for the reason given above: log messages are 0 where a sequence has
  # ended and can sit near zero.
  expect_reference(acceleration_messages(new, 1L), acceleration_messages(old, 1L),
                   reference = "independent: log-space forward-backward recursion (helper-lta-log-reference.R)",
                   covers = "lta()", quantity = "log messages, general kernel",
                   tolerance = 1e-12, measure = "scaled",
                   key = c("case", "message", "occasion", "row", "column"), note = general_note)
  # Unit check: the occupancies of every occasion are distributions.
  invisible(lapply(seq_len(4), function(t) {
    expect_equal(rowSums(exp(new$alpha[[t]] + new$beta[[t]] - new$log_scaled)),
                 rep(1, 3), tolerance = 1e-12)
  }))
})
