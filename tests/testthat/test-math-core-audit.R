test_that("extreme common Gaussian densities retain profile and group priors", {
  parameters <- list(means = matrix(0, 2L, 1L), variances = matrix(1, 2L, 1L),
    profile_probabilities = rbind(c(.2, .8), c(.7, .3)),
    group_probabilities = c(.3, .7))
  x <- matrix(c(1e8, -1e8, 2e8), 3L, 1L)
  groups <- c(1L, 2L, 1L)
  expectation <- .multilpa_expectation(x, groups, parameters)
  # Identical measurement distributions supply no class information, regardless
  # of how surprising the observations are under those distributions.
  expect_equal(expectation$group_posteriors,
               matrix(rep(c(.3, .7), each = 2L), 2L), tolerance = 1e-14)
  expect_equal(expectation$subject_posteriors,
               matrix(rep(c(.55, .45), each = 3L), 3L), tolerance = 1e-14)
  expect_equal(rowSums(expectation$subject_posteriors), rep(1, 3L))
  expect_equal(rowSums(expectation$group_posteriors), rep(1, 2L))
  expect_equal(expectation$log_likelihood, sum(dnorm(x, log = TRUE)))

  # A categorical observation must still provide information when the Gaussian
  # block is extremely unlikely. Its reference is the categorical-only model.
  codes <- matrix(c(1L, 2L, 1L), 3L, 1L)
  parameters$response_probabilities <- list(rbind(c(.8, .2), c(.1, .9)))
  mixed <- .multilpa_expectation(x, groups, parameters, codes)
  categorical <- parameters
  categorical$means <- categorical$variances <- matrix(numeric(0), 2L, 0L)
  reference <- .multilpa_expectation(matrix(numeric(0), 3L, 0L), groups,
                                    categorical, codes)
  expect_equal(mixed$subject_posteriors, reference$subject_posteriors,
               tolerance = 1e-14)
  expect_equal(mixed$group_posteriors, reference$group_posteriors,
               tolerance = 1e-14)
})

test_that("all-categorical starts have no Gaussian covariance requirement", {
  data <- data.frame(group = rep(1:4, each = 3), item = rep(c(1, 2, 2), 4))
  start <- list(response_probabilities = list(matrix(c(1 / 3, 2 / 3), 1L)),
                profile_probabilities = matrix(1), group_probabilities = 1)
  fits <- lapply(c("diagonal", "full"), function(covariance) {
    multilpa(data, "item", "group", 1, 1, categorical = "item",
             covariance_model = covariance, start = start, n_starts = 1)
  })
  expect_equal(fits[[1]]$log_likelihood, 4 * log(1 / 3) + 8 * log(2 / 3))
  expect_equal(fits[[2]]$log_likelihood, fits[[1]]$log_likelihood)
  expect_equal(fits[[2]]$response_probabilities, fits[[1]]$response_probabilities)
  expect_null(fits[[2]]$covariances)
  expect_equal(fits[[2]]$n_parameters, 1)
})

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

test_that("time uniqueness respects distinct native numeric group identifiers", {
  data <- data.frame(group = rep(c(1, 1 + 1e-15), each = 3),
                     time = rep(1:3, 2), score = c(-2, -1, -3, 2, 1, 3))
  expect_equal(length(unique(data$group)), 2L)
  fit <- multilpa(data, "score", "group", 1, 1, n_starts = 1, time = "time")
  expect_equal(fit$n_groups, 2L)
  expect_equal(fit$group_index, rep(1:2, each = 3))
  expect_equal(fit$time_values, data$time)
  data$time[2] <- data$time[1]
  expect_error(multilpa(data, "score", "group", 1, 1, n_starts = 1, time = "time"),
               class = "multilpa_bad_time")
})
