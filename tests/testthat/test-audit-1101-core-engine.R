# Audit 2026-11-01, slice core-engine: regression tests for the defects fixed
# in R/kernel-structures.R and R/acceleration.R, and the checks of the
# covariance-structure M-steps against their definitions and mclust.
# Log: validation/audit-2026-11-01/core-engine.md.

skip_on_cran()

# A brute-force two-level likelihood and group posteriors for the nested
# structure: L_j = sum_h omega_h prod_{i in j} sum_k pi_{k|h} f_k(y_i).
audit_nested_reference <- function(f, group_index, profile_probabilities,
                                   group_probabilities) {
  groups <- seq_len(max(group_index))
  terms <- t(vapply(groups, function(j) {
    rows <- which(group_index == j)
    vapply(seq_along(group_probabilities), function(h) {
      group_probabilities[h] *
        prod(f[rows, , drop = FALSE] %*% profile_probabilities[h, ])
    }, numeric(1))
  }, numeric(length(group_probabilities))))
  list(log_likelihood = sum(log(rowSums(terms))),
       group_posteriors = terms / rowSums(terms))
}

test_that("the nested structure puts group j in row j for any group numbering", {
  set.seed(11)
  n <- 30L
  f <- matrix(stats::runif(n * 2L, 0.05, 1), n, 2L)
  # Groups numbered 1..5, but not in order of first appearance.
  group_index <- c(3L, 1L, 5L, 2L, 4L)[sample(5L, n, TRUE)]
  parameters <- list(profile_probabilities = rbind(c(.3, .7), c(.8, .2)),
                     group_probabilities = c(.4, .6))
  reference <- audit_nested_reference(f, group_index,
                                      parameters$profile_probabilities,
                                      parameters$group_probabilities)
  # Split each density into a removed offset and the rest, as the E-step does.
  offset <- stats::rnorm(n)
  found <- .multilpa_nested_structure(log(f) - offset, offset, group_index,
                                      parameters)
  expect_equal(found$log_likelihood, reference$log_likelihood, tolerance = 1e-12)
  # Before the fix the rows came in first-appearance order: off by 0.73.
  expect_equal(found$group_posteriors, reference$group_posteriors,
               tolerance = 1e-12)
  # Weights are per group index, so they must meet the same rows.
  weights <- c(2, 1, 3, 1, 2)
  weighted <- .multilpa_weigh(found, weights, group_index)
  expect_equal(weighted$log_likelihood,
               sum(weights * found$group_log_likelihood), tolerance = 1e-12)
  by_group <- vapply(seq_len(5L), function(j) {
    log(sum(vapply(1:2, function(h) {
      parameters$group_probabilities[h] *
        prod(f[group_index == j, , drop = FALSE] %*% parameters$profile_probabilities[h, ])
    }, numeric(1))))
  }, numeric(1))
  expect_equal(weighted$log_likelihood, sum(weights * by_group), tolerance = 1e-12)
})

test_that("renumbering the groups only permutes the nested structure's rows", {
  set.seed(5)
  n <- 24L
  log_density <- matrix(stats::rnorm(n * 3L, -2), n, 3L)
  first <- rep(seq_len(6L), each = 4L)
  parameters <- list(profile_probabilities = rbind(c(.2, .3, .5), c(.6, .3, .1)),
                     group_probabilities = c(.5, .5))
  permutation <- c(4L, 6L, 1L, 3L, 2L, 5L)
  ordered <- .multilpa_nested_structure(log_density, numeric(n), first, parameters)
  renumbered <- .multilpa_nested_structure(log_density, numeric(n),
                                           permutation[first], parameters)
  expect_equal(renumbered$group_posteriors[permutation, ], ordered$group_posteriors,
               tolerance = 1e-14)
  expect_equal(renumbered$subject_posteriors, ordered$subject_posteriors,
               tolerance = 1e-14)
  expect_equal(renumbered$log_likelihood, ordered$log_likelihood, tolerance = 1e-14)
})

test_that("the logit structure's likelihood does not depend on the group numbering", {
  set.seed(8)
  n <- 25L
  group <- sample(c("q", "a", "m", "z"), n, TRUE)
  first <- match(group, unique(group))
  log_density <- matrix(stats::rnorm(n * 3L, -2), n, 3L)
  offset <- stats::rnorm(n)
  covariate <- stats::rnorm(n)
  profile_design <- lapply(1:2, function(h) {
    cbind(diag(2L)[rep(h, n), , drop = FALSE], covariate)
  })
  group_design <- cbind(1, stats::rnorm(4L))
  beta <- matrix(stats::rnorm(6L), 3L)
  gamma <- matrix(stats::rnorm(2L), 2L)
  softmax <- function(scores) {
    scores <- c(scores, 0)
    exp(scores - max(scores)) / sum(exp(scores - max(scores)))
  }
  by_hand <- sum(vapply(seq_len(4L), function(j) {
    rows <- which(first == j)
    eta <- softmax(as.vector(group_design[j, ] %*% gamma))
    log(sum(vapply(1:2, function(h) {
      eta[h] * prod(vapply(rows, function(i) {
        sum(softmax(as.vector(profile_design[[h]][i, ] %*% beta)) *
              exp(log_density[i, ] + offset[i]))
      }, numeric(1)))
    }, numeric(1))))
  }, numeric(1)))
  ordered <- .multilpa_logit_structure(log_density, offset, first, profile_design,
                                       group_design, beta, gamma)
  expect_equal(ordered$log_likelihood, by_hand, tolerance = 1e-12)
  # Group j's design row sits in row j under any numbering. Before the fix the
  # likelihood came out 0.033 lower, with group posteriors off by 0.11.
  permutation <- c(3L, 1L, 4L, 2L)
  renumbered_design <- group_design
  renumbered_design[permutation, ] <- group_design
  renumbered <- .multilpa_logit_structure(log_density, offset, permutation[first],
                                          profile_design, renumbered_design,
                                          beta, gamma)
  expect_equal(renumbered$log_likelihood, by_hand, tolerance = 1e-12)
  expect_equal(renumbered$group_posteriors[permutation, ], ordered$group_posteriors,
               tolerance = 1e-12)
  expect_equal(renumbered$subject_posteriors, ordered$subject_posteriors,
               tolerance = 1e-12)
})

test_that("the homogeneous transition accelerator measures rows of positive weight", {
  set.seed(2)
  x <- cbind(stats::rnorm(50, 0, 2), stats::rnorm(50, 0, 10))
  parameters <- list(means = matrix(0, 2L, 2L), variances = matrix(1, 2L, 2L),
                     profile_probabilities = matrix(c(.5, .5), 1L),
                     group_probabilities = 1)
  row_weights <- c(rep(0, 10L), rep(1, 40L))
  # The excluded rows take absurd values: they must not set the chart.
  outlying <- x
  outlying[1:10, ] <- 1e6
  accelerator <- .lta_accelerator("squarem", outlying, parameters,
                                  general = FALSE, row_weights = row_weights)
  # Means are charted in reference standard deviations: with every mean 1
  # the flattened means are 1 / sd over the 40 rows of positive weight.
  unit_means <- utils::modifyList(parameters, list(means = matrix(1, 2L, 2L)))
  spread <- apply(x[11:50, ], 2L, stats::sd)
  expect_equal(accelerator$flatten(unit_means)[1:4], rep(1 / spread, each = 2L),
               tolerance = 1e-12)
  expect_identical(accelerator$flatten(parameters),
                   .multilpa_accelerator(outlying, parameters,
                                         row_weights)$flatten(parameters))
})

# The covariance-structure M-steps ------------------------------------------

audit_mstep_data <- function() {
  set.seed(3)
  n <- 300L
  x <- matrix(stats::rnorm(n * 3L), n, 3L)
  x[1:100, ] <- x[1:100, ] %*% matrix(c(2, .8, 0, 0, 1, .3, 0, 0, .5), 3L) + 2
  x[101:200, ] <- x[101:200, ] %*% diag(c(.5, 2, 1)) - 1
  z <- matrix(stats::runif(n * 3L), n, 3L)
  list(x = x, z = z / rowSums(z))
}

# The M-step through the path EM uses: posteriors and conditional moments in,
# means and covariances out, for one structure.
audit_mstep <- function(data, code, inner_steps = TRUE) {
  d <- ncol(data$x)
  k <- ncol(data$z)
  moments <- .multilpa_gaussian_moments(data$x, list(
    means = matrix(0, k, d), variances = matrix(1, k, d),
    covariances = array(diag(d), c(d, d, k))))
  full <- .multilpa_is_ellipsoidal(code)
  equal <- startsWith(code, "EE")
  out <- .multilpa_maximize_moments(
    data$x, list(subject_posteriors = data$z, gaussian_moments = moments$moments),
    if (equal) "equal" else "varying", 1e-12,
    if (full) "full" else "diagonal", structure = code)
  sigma <- if (full) out$covariances else array(unlist(lapply(seq_len(k), function(p) {
    diag(out$variances[p, ], d)
  })), c(d, d, k))
  list(means = out$means, sigma = sigma)
}

# What the M-step minimises: sum_k n_k log|Sigma_k| + tr(Sigma_k^-1 W_k).
audit_mstep_objective <- function(data, sigma, means) {
  sum(vapply(seq_len(ncol(data$z)), function(p) {
    block <- matrix(sigma[, , p], ncol(data$x))
    residuals <- sweep(data$x, 2L, means[p, ], "-")
    scatter <- crossprod(residuals, residuals * data$z[, p])
    sum(data$z[, p]) * as.numeric(determinant(block, logarithm = TRUE)$modulus) +
      sum(chol2inv(chol(block)) * scatter)
  }, numeric(1)))
}

test_that("every structure's M-step is mclust's, or no worse where it iterates", {
  skip_if_not_installed("mclust")
  data <- audit_mstep_data()
  invisible(lapply(.multilpa_structures(), function(code) {
    reference <- do.call(mclust::mstep, list(modelName = code, data = data$x, z = data$z),
                         envir = asNamespace("mclust"))
    theirs <- reference$parameters$variance$sigma
    ours <- audit_mstep(data, code)
    expect_equal(ours$means, t(reference$parameters$mean), tolerance = 1e-12,
                 ignore_attr = TRUE, info = code)
    if (code %in% c("EVE", "VVE")) {
      # Inside EM these take five warm-started orientation steps (a
      # generalized M-step); solved to convergence they reach mclust's
      # minimum or a lower one.
      w <- colSums(data$z)
      scatter <- lapply(seq_len(3L), function(p) {
        residuals <- sweep(data$x, 2L, ours$means[p, ], "-")
        crossprod(residuals, residuals * data$z[, p])
      })
      solved <- .multilpa_structure_covariances(scatter, w, code, 1e-12)
      expect_lte(audit_mstep_objective(data, solved, ours$means),
                 audit_mstep_objective(data, theirs, ours$means) + 1e-8,
                 label = code)
    } else {
      expect_equal(ours$sigma, theirs, tolerance = 1e-8, ignore_attr = TRUE,
                   info = code)
    }
  }))
})

test_that("a structure's M-step objective is no worse than any structure it nests", {
  data <- audit_mstep_data()
  objective <- vapply(.multilpa_structures(), function(code) {
    fitted <- audit_mstep(data, code)
    audit_mstep_objective(data, fitted$sigma, fitted$means)
  }, numeric(1))
  nests <- list(VII = "EII", EEI = "EII", VEI = c("VII", "EEI"), EVI = "EEI",
                VVI = c("VEI", "EVI"), EEE = "EEI", VEE = c("EEE", "VEI"),
                EVE = c("EEE", "EVI"), VVE = c("VEE", "EVE", "VVI"),
                EEV = "EEE", VEV = c("EEV", "VEE"), EVV = c("EEV", "EVE"),
                VVV = c("VEV", "EVV", "VVE"))
  invisible(lapply(names(nests), function(looser) {
    lapply(nests[[looser]], function(tighter) {
      expect_lte(objective[[looser]], objective[[tighter]] + 1e-8,
                 label = paste(looser, "over", tighter))
    })
  }))
})

test_that("a looser structure started at a nested fit ends no lower", {
  data <- audit_mstep_data()
  frame <- stats::setNames(as.data.frame(data$x), c("a", "b", "c"))
  fit <- function(code, start = NULL) {
    quietly(multilpa(frame, c("a", "b", "c"), id = NULL, n_profiles = 3L,
                     n_starts = 1L, seed = 1, tol = 1e-10, max_iter = 3000L,
                     model = code, start = start))
  }
  start_from <- function(nested, full) {
    start <- list(means = unname(nested$means), variances = unname(nested$variances),
                  profile_probabilities = unname(nested$profile_probabilities),
                  group_probabilities = 1)
    if (full) {
      start$covariances <- if (!is.null(nested$covariances)) unname(nested$covariances) else
        array(unlist(lapply(seq_len(3L), function(p) diag(start$variances[p, ]))),
              c(3L, 3L, 3L))
    }
    start
  }
  pairs <- list(c("EVI", "EEI"), c("VVI", "EVI"), c("EVE", "EVI"),
                c("VVE", "EVE"), c("VEV", "VEE"), c("EVV", "EVE"),
                c("VVV", "VVE"))
  invisible(lapply(pairs, function(pair) {
    nested <- fit(pair[[2L]])
    looser <- fit(pair[[1L]], start_from(nested, .multilpa_is_ellipsoidal(pair[[1L]])))
    expect_gte(looser$log_likelihood, nested$log_likelihood - 1e-6,
               label = paste(pair, collapse = " from "))
  }))
})

test_that("the free-parameter count of every structure is mclust's", {
  skip_if_not_installed("mclust")
  grid <- expand.grid(code = .multilpa_structures(), d = 1:4, k = 1:4,
                      stringsAsFactors = FALSE)
  ours <- mapply(.multilpa_structure_parameters, grid$code, grid$k, grid$d)
  theirs <- mapply(mclust::nVarParams, grid$code, grid$d, grid$k)
  expect_equal(unname(ours), unname(theirs))
})

test_that("the full-covariance likelihood with missing values is the marginal one", {
  set.seed(11)
  n <- 40L
  x <- matrix(stats::rnorm(n * 3L), n, 3L)
  x[sample(n * 3L, 15L)] <- NA
  x[5L, ] <- NA
  group_index <- match(sample(c(7, 3, 9, 1, 5), n, TRUE), c(7, 3, 9, 1, 5))
  factor <- matrix(c(1, .4, .2, 0, 1, .3, 0, 0, 1), 3L)
  covariances <- array(c(crossprod(factor), diag(c(2, .5, 1))), c(3L, 3L, 2L))
  parameters <- list(means = rbind(c(0, 1, -1), c(1, -1, 0)),
                     variances = rbind(diag(covariances[, , 1]), diag(covariances[, , 2])),
                     covariances = covariances,
                     profile_probabilities = rbind(c(.3, .7), c(.8, .2)),
                     group_probabilities = c(.4, .6))
  density <- function(row, k) {
    observed <- which(!is.na(row))
    if (length(observed) == 0L) return(1)
    sigma <- matrix(covariances[observed, observed, k], length(observed))
    residual <- row[observed] - parameters$means[k, observed]
    exp(-0.5 * (length(observed) * log(2 * pi) +
                  as.numeric(determinant(sigma)$modulus) +
                  sum(residual * solve(sigma, residual))))
  }
  f <- outer(seq_len(n), 1:2, Vectorize(function(i, k) density(x[i, ], k)))
  reference <- audit_nested_reference(f, group_index, parameters$profile_probabilities,
                                      parameters$group_probabilities)
  found <- .multilpa_expectation(x, group_index, parameters)
  expect_equal(found$log_likelihood, reference$log_likelihood, tolerance = 1e-10)
  expect_equal(found$group_posteriors, reference$group_posteriors, tolerance = 1e-12)
})
