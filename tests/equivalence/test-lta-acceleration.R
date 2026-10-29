# Comparisons with independent references moved from tests/testthat/test-lta-acceleration.R,
# so that they run in CI's equivalence job rather than under R CMD check.

# Scaled recursions are compared with the original log-space implementation,
# including its messages: inference still consumes alpha and beta in log space.
test_that("scaled homogeneous messages and counts match the log recursion", {
  set.seed(912)
  index <- rep(seq_len(4L), c(1L, 3L, 5L, 4L))
  time <- c(1, 1, 3, 5, 1:5, 1, 2, 4, 5)
  layout <- .multilpa_sequence_layout(index, time, 4L, "grid")
  emission <- .multilpa_sequence_emission(matrix(rnorm(length(index) * 3), ncol = 3),
                                          layout)$log_density
  weights <- cbind(c(0.2, 0.6, 0.7, 0.3), c(0.8, 0.4, 0.3, 0.7))
  moments <- lapply(seq_len(2L), function(h) {
    initial <- prop.table(runif(3))
    transition <- matrix(runif(9), 3)
    transition <- transition / rowSums(transition)
    log_initial <- matrix(log(initial), 4, 3, byrow = TRUE)
    arrays <- rep(list(array(rep(log(transition), each = 4), c(4, 3, 3))), 4)
    old <- .lta_log_reference(emission, layout, log_initial, arrays)
    new <- .multilpa_forward_backward(emission, layout, initial, transition)
    # Also exercise the shared kernel's homogeneous matrix interface.
    shared <- .latents_forward_backward(emission, layout, log_initial, log(transition))
    expect_equal(shared, new, tolerance = 1e-12)
    expect_equal(new$log_scaled, old$log_scaled, tolerance = 1e-12)
    expect_equal(new$alpha, old$alpha, tolerance = 1e-12)
    expect_equal(new$beta, old$beta, tolerance = 1e-12)
    expected <- .lta_moments_reference(old, emission, layout, transition,
                                       weights[, h], length(index))
    actual <- .multilpa_sequence_moments(new, emission, layout, transition,
                                         weights[, h], length(index))
    expect_equal(actual, expected, tolerance = 1e-12)
    expect_equal(rowSums(actual$posterior), weights[index, h], tolerance = 1e-12)
    expect_equal(sum(actual$initial), sum(weights[, h]), tolerance = 1e-12)
    expect_equal(sum(actual$transition), sum(weights[, h] * (layout$span - 1)),
                 tolerance = 1e-12)
    actual
  })
  expect_equal(rowSums(Reduce(`+`, lapply(moments, `[[`, "posterior"))),
               rep(1, length(index)), tolerance = 1e-12)
})

test_that("the general per-group per-occasion kernel preserves log messages", {
  set.seed(715)
  layout <- .multilpa_sequence_layout(rep(1:3, c(2, 4, 3)),
                                      c(1:2, 1:4, 1:3), 3L)
  emission <- lapply(1:4, function(t) matrix(-rexp(9), 3, 3))
  initial <- matrix(runif(9), 3)
  initial <- log(initial / rowSums(initial))
  transition <- lapply(1:3, function(t) {
    p <- array(runif(27), c(3, 3, 3))
    log(sweep(p, c(1, 2), apply(p, c(1, 2), sum), "/"))
  })
  old <- .lta_log_reference(emission, layout, initial, transition)
  new <- .latents_forward_backward(emission, layout, initial, transition)
  expect_equal(new, old, tolerance = 1e-12)
  invisible(lapply(seq_len(4), function(t) {
    expect_equal(rowSums(exp(new$alpha[[t]] + new$beta[[t]] - new$log_scaled)),
                 rep(1, 3), tolerance = 1e-12)
  }))
})
