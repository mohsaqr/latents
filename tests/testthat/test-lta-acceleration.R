# The comparisons with independent references are in tests/equivalence/test-lta-acceleration.R.


test_that("singleton sequences and one-profile kernels retain their dimensions", {
  invisible(lapply(c(1L, 3L), function(n_groups) {
    invisible(lapply(c(1L, 2L), function(k) {
      layout <- .multilpa_sequence_layout(seq_len(n_groups), rep(1, n_groups), n_groups)
      emission <- list(matrix(0, n_groups, k))
      pass <- .multilpa_forward_backward(emission, layout, rep(1 / k, k),
                                          matrix(1 / k, k, k))
      expect_equal(pass$log_scaled, rep(0, n_groups), tolerance = 1e-12)
      expect_length(pass$alpha, 1L)
      expect_length(pass$beta, 1L)
      counts <- .multilpa_sequence_moments(pass, emission, layout,
                                           matrix(1 / k, k, k), rep(1, n_groups), n_groups)
      expect_equal(rowSums(counts$posterior), rep(1, n_groups), tolerance = 1e-12)
      expect_equal(sum(counts$transition), 0)
    }))
  }))
})

test_that("transition SQUAREM maps round-trip and reject invalid probabilities", {
  p <- list(means = matrix(c(-1, 1), 2), variances = matrix(1, 2),
            group_probabilities = c(0.4, 0.6),
            initial_probabilities = matrix(c(0.3, 0.7, 0.7, 0.3), 2),
            transition_probabilities = array(c(0.8, 0.3, 0.2, 0.7), c(2, 2, 2)))
  x <- matrix(c(-3, -1, 1, 3), ncol = 1L)
  accelerator <- .lta_accelerator("squarem", x, p)
  expect_equal(accelerator$restore(p, accelerator$flatten(p)), p, tolerance = 1e-15)
  expect_true(accelerator$valid(p))
  # The coordinates are in standard deviations of the data, so the same
  # parameters in other units give the same point.
  scaled <- p
  scaled$means <- p$means * 1000
  scaled$variances <- p$variances * 1000^2
  expect_equal(.lta_accelerator("squarem", x * 1000, scaled)$flatten(scaled),
               accelerator$flatten(p), tolerance = 1e-14)
  bad <- p
  bad$initial_probabilities[1, 1] <- 0
  expect_false(accelerator$valid(bad))
  bad <- p
  bad$initial_probabilities[1, 1] <- 0.5
  expect_false(accelerator$valid(bad))
  bad <- p
  bad$transition_probabilities[1, 1, 1] <- -0.1
  expect_false(accelerator$valid(bad))
  bad <- p
  bad$transition_probabilities[1, 1, 1] <- 0.5
  expect_false(accelerator$valid(bad))
})

test_that("SQUAREM and ordinary lta reach the same seeded likelihood", {
  set.seed(519)
  panel <- data.frame(id = rep(1:35, each = 6), time = rep(1:6, 35))
  state <- sample(c(-2, 2), nrow(panel), replace = TRUE)
  panel$y <- rnorm(nrow(panel), state, 0.5)
  panel$z <- rnorm(nrow(panel), state, 0.7)
  panel$y[c(5, 31, 88)] <- NA
  fits <- lapply(c("none", "squarem"), function(acceleration) {
    quietly(lta(panel, c("y", "z"), "id", 2, "time", n_starts = 1, seed = 1,
                missing = "fiml", acceleration = acceleration))
  })
  expect_true(all(vapply(fits, `[[`, logical(1), "converged")))
  expect_lt(abs(fits[[1]]$log_likelihood - fits[[2]]$log_likelihood), 1e-6)
  expect_true(all(diff(fits[[2]]$log_likelihood_history) >= -1e-8))
  # Enumeration accepts and forwards the new argument through its dots.
  grid <- quietly(enumerate_classes(panel, c("y", "z"), id = "id", time = "time",
                                    n_profiles = 2, n_starts = 1, seed = 1,
                                    missing = "fiml", acceleration = "none"))
  refit <- quietly(.lta_refit(fits[[1]], panel, n_starts = 1, max_iter = 0, seed = 1))
  expect_identical(refit$acceleration, "none")
  expect_identical(grid$fits[[1]]$acceleration, "none")
  expect_equal(grid$fits[[1]]$log_likelihood, fits[[1]]$log_likelihood)
})

test_that("general SQUAREM maps preserve every updated block and the stayer mask", {
  parameters <- list(
    measurement = rep(list(list(means = matrix(c(-1, 1), 2),
                                 variances = matrix(1, 2),
                                 response_probabilities = list(matrix(0.5, 2, 2)))), 2),
    initial = array(0, c(2, 1, 2)),
    transition = array(0, c(3, 1, 2, 2)),
    transition2 = array(0, c(3, 1, 4, 2)),
    group_probabilities = c(0.6, 0.4), stayer = c(FALSE, TRUE))
  attr(parameters$transition, "n_occasion_intercepts") <- 2L
  x <- matrix(c(-3, -1, 1, 3), ncol = 1L)
  accelerator <- .lta_accelerator("squarem", x, parameters, general = TRUE)
  values <- accelerator$flatten(parameters)
  expect_length(values, 2 * (2 + 2 + 4) + 4 + 12 + 24 + 2)
  expect_equal(accelerator$restore(parameters, values), parameters, tolerance = 1e-15)
  scaled <- parameters
  scaled$measurement <- lapply(parameters$measurement, function(block) {
    block$means <- block$means / 1000
    block$variances <- block$variances / 1000^2
    block
  })
  expect_equal(.lta_accelerator("squarem", x / 1000, scaled, general = TRUE)$flatten(scaled),
               values, tolerance = 1e-14)
  expect_true(accelerator$valid(parameters))
  bad <- parameters
  bad$transition2[1] <- Inf
  expect_false(accelerator$valid(bad))
  bad <- parameters
  bad$measurement[[1]]$response_probabilities[[1]][1, 1] <- 0.3
  expect_false(accelerator$valid(bad))
})
