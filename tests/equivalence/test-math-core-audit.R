# Comparisons with independent references moved from tests/testthat/test-math-core-audit.R,
# so that they run in CI's equivalence job rather than under R CMD check.

test_that("mixed full-covariance FIML agrees with exhaustive latent assignments", {
  x <- rbind(c(.2, -1, .6), c(NA, .3, -.8), c(NA, NA, NA),
             c(.4, NA, -.5), c(.9, .1, NA))
  codes <- rbind(c(1L, 2L), c(2L, NA), c(NA, NA), c(3L, 1L), c(1L, 2L))
  group <- c(1L, 2L, 1L, 3L, 2L)
  covariance <- array(c(1, .3, -.2, .3, 2, .4, -.2, .4, 1.5,
                         2, -.4, .1, -.4, 1.5, .2, .1, .2, .8), c(3L, 3L, 2L))
  parameters <- list(means = rbind(c(-.4, .6, -.2), c(.7, -.3, .8)),
    variances = t(vapply(1:2, function(k) diag(covariance[, , k]), numeric(3))),
    covariances = covariance,
    profile_probabilities = rbind(c(.8, .2), c(.35, .65)),
    group_probabilities = c(.4, .6),
    response_probabilities = list(rbind(c(.2, .3, .5), c(.6, .3, .1)),
                                  rbind(c(.7, .3), c(.1, .9))))
  # Direct densities use determinants/inverses, and the joint latent states are
  # summed exhaustively. This shares neither log-domain arithmetic nor the E-step.
  density <- vapply(1:2, function(k) {
    vapply(seq_len(nrow(x)), function(i) {
      observed <- which(!is.na(x[i, ]))
      gaussian <- if (!length(observed)) 1 else {
        residual <- x[i, observed] - parameters$means[k, observed]
        sigma <- covariance[observed, observed, k, drop = FALSE]
        dim(sigma) <- c(length(observed), length(observed))
        exp(-sum(residual * solve(sigma, residual)) / 2) /
          sqrt((2 * pi)^length(observed) * det(sigma))
      }
      categorical <- prod(vapply(seq_len(ncol(codes)), function(j) {
        if (is.na(codes[i, j])) 1 else parameters$response_probabilities[[j]][k, codes[i, j]]
      }, numeric(1)))
      gaussian * categorical
    }, numeric(1))
  }, numeric(nrow(x)))
  references <- lapply(unique(group), function(g) {
    rows <- which(group == g)
    assignments <- as.matrix(expand.grid(rep(list(1:2), length(rows))))
    masses <- vapply(1:2, function(h) {
      apply(assignments, 1L, function(z) {
        parameters$group_probabilities[h] *
          prod(parameters$profile_probabilities[h, z]) *
          prod(density[cbind(rows, z)])
      })
    }, numeric(nrow(assignments)))
    subject <- vapply(seq_along(rows), function(i) {
      vapply(1:2, function(k) sum(masses[assignments[, i] == k, , drop = FALSE]), numeric(1))
    }, numeric(2))
    list(rows = rows, likelihood = sum(masses),
         group = colSums(masses) / sum(masses), subject = t(subject) / sum(masses))
  })
  expectation <- .multilpa_expectation(x, group, parameters, codes)
  expect_equal(expectation$log_likelihood,
               sum(vapply(references, function(r) log(r$likelihood), numeric(1))),
               tolerance = 1e-12)
  expect_equal(expectation$group_posteriors,
               do.call(rbind, lapply(references, `[[`, "group")), tolerance = 1e-12)
  invisible(lapply(references, function(r) {
    expect_equal(expectation$subject_posteriors[r$rows, , drop = FALSE], r$subject,
                 tolerance = 1e-12)
  }))
})

test_that("bounded multinomial M-steps satisfy the constrained optimum", {
  # Enumerating all feasible three-category grid points supplies an independent
  # lower bound on the optimum, including cases with multiple active bounds.
  grid <- expand.grid(first = seq(.15, .7, by = .005),
                       second = seq(.15, .7, by = .005))
  grid$third <- 1 - grid$first - grid$second
  feasible <- as.matrix(grid[grid$third >= .15, ])
  invisible(lapply(list(c(0, 1, 19), c(3, 7, 10), c(1, 1, 1), c(0, 0, 0)), function(counts) {
    probabilities <- .multilpa_bound_probabilities(counts, .15)
    expect_equal(sum(probabilities), 1, tolerance = 1e-14)
    expect_true(all(probabilities >= .15))
    optimum <- sum(counts * log(probabilities))
    expect_gte(optimum + 1e-12, max(as.vector(log(feasible) %*% counts)))
    free <- probabilities > .15 + 1e-12
    if (any(free) && sum(counts) > 0) {
      multiplier <- counts[free] / probabilities[free]
      expect_lt(diff(range(multiplier)), 1e-12)
      expect_true(all(counts[!free] / probabilities[!free] <= multiplier[1] + 1e-12))
    }
  }))
})
