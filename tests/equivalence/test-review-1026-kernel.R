# Comparisons with independent references moved from tests/testthat/test-review-1026-kernel.R,
# so that they run in CI's equivalence job rather than under R CMD check.

kernel_review_softmax <- function(scores) exp(scores) / rowSums(exp(scores))

test_that("structure posteriors match brute force, and survive a -5e15 offset per row", {
  set.seed(5)
  n <- 30L
  k <- 3L
  n_group_classes <- 2L
  n_groups <- 6L
  group_index <- rep(seq_len(n_groups), each = 5L)
  w <- cbind(1, stats::rnorm(n))
  v <- cbind(1, stats::rnorm(n_groups))
  w_group <- cbind(1, stats::rnorm(n_groups))
  gamma <- cbind(0, matrix(stats::rnorm(2L * (k - 1L)), 2L))
  delta <- cbind(0, matrix(stats::rnorm(2L * (n_group_classes - 1L)), 2L))
  class_logits <- cbind(0, matrix(stats::rnorm((n_group_classes + 2L) * (k - 1L)),
                                  n_group_classes + 2L))
  observation <- list(nesting = "observation", w = w)
  group <- list(nesting = "group", w = w_group, group_index = group_index)
  two_level <- list(nesting = "two-level", w = w, v = v, n = n,
                    n_group_classes = n_group_classes, group_index = group_index,
                    n_groups = n_groups)
  brute <- function(relative) {
    density <- exp(relative)
    prior <- kernel_review_softmax(w %*% gamma)
    group_prior <- kernel_review_softmax(w_group %*% gamma)
    group_density <- exp(rowsum(relative, group_index))
    eta <- kernel_review_softmax(v %*% delta)
    conditional <- lapply(seq_len(n_group_classes), function(h) {
      kernel_review_softmax(.mixture_two_level_design(two_level, h) %*% class_logits)
    })
    by_class <- vapply(seq_len(n_group_classes), function(h) {
      as.vector(exp(rowsum(log(rowSums(density * conditional[[h]])), group_index)))
    }, numeric(n_groups))
    rho <- eta * by_class / rowSums(eta * by_class)
    list(tau = density * prior / rowSums(density * prior),
         log_likelihood = sum(log(rowSums(density * prior))),
         group_tau = group_density * group_prior / rowSums(group_density * group_prior),
         rho = rho,
         two_level_tau = Reduce(`+`, lapply(seq_len(n_group_classes), function(h) {
           rho[group_index, h] * density * conditional[[h]] /
             rowSums(density * conditional[[h]])
         })))
  }
  log_density <- matrix(stats::rnorm(n * k, -3), n, k)
  reference <- brute(log_density)
  first <- .mixture_structure_expectation(observation, list(gamma = gamma), log_density)
  expect_equal(first$tau, reference$tau, tolerance = 1e-13)
  expect_equal(first$log_likelihood, reference$log_likelihood, tolerance = 1e-13)
  second <- .mixture_structure_expectation(group, list(gamma = gamma), log_density)
  expect_equal(second$group_tau, reference$group_tau, tolerance = 1e-13)
  third <- .mixture_structure_expectation(
    two_level, list(delta = delta, class_logits = class_logits), log_density)
  expect_equal(third$rho, reference$rho, tolerance = 1e-13)
  expect_equal(third$tau, reference$two_level_tau, tolerance = 1e-13)
  # A large per-row constant changes no posterior. Before 8ac3576 the
  # observation posteriors were off by about 0.44 here.
  extreme <- log_density + (-5e15 * stats::runif(n, 0.5, 1.5))
  relative <- extreme - .multilpa_row_max(extreme)
  truth <- brute(relative)
  shifted <- .mixture_structure_expectation(observation, list(gamma = gamma), extreme)
  expect_equal(shifted$tau, truth$tau, tolerance = 1e-13)
  shifted_group <- .mixture_structure_expectation(group, list(gamma = gamma), extreme)
  expect_equal(shifted_group$group_tau, truth$group_tau, tolerance = 1e-13)
  shifted_two <- .mixture_structure_expectation(
    two_level, list(delta = delta, class_logits = class_logits), extreme)
  expect_equal(shifted_two$rho, truth$rho, tolerance = 1e-13)
  expect_equal(shifted_two$tau, truth$two_level_tau, tolerance = 1e-13)
  expect_equal(rowSums(shifted_two$tau), rep(1, n), tolerance = 1e-14)
})

test_that("the forward-backward recursions equal path enumeration", {
  set.seed(9)
  n_groups <- 4L
  k <- 3L
  n_occasions <- 4L
  lengths <- c(4L, 2L, 3L, 1L)
  layout <- list(slot = matrix(1, n_groups, n_occasions), n_occasions = n_occasions,
                 within = outer(lengths, seq_len(n_occasions), ">="))
  emission <- lapply(seq_len(n_occasions), function(t) {
    matrix(stats::rnorm(n_groups * k, -1), n_groups, k)
  })
  initial <- matrix(c(.5, .3, .2), n_groups, k, byrow = TRUE)
  initial[2L, ] <- c(.1, .1, .8)
  transition <- matrix(c(.7, .2, .1, .1, .8, .1, .3, .3, .4), k, k, byrow = TRUE)
  paths <- as.matrix(expand.grid(rep(list(seq_len(k)), n_occasions)))
  brute <- lapply(seq_len(n_groups), function(g) {
    prefixes <- unique(paths[, seq_len(lengths[g]), drop = FALSE])
    log_path <- apply(prefixes, 1L, function(state) {
      log(initial[g, state[1L]]) + emission[[1L]][g, state[1L]] +
        sum(vapply(seq_along(state)[-1L], function(t) {
          log(transition[state[t - 1L], state[t]]) + emission[[t]][g, state[t]]
        }, numeric(1)))
    })
    total <- log(sum(exp(log_path)))
    list(log_likelihood = total, posterior = vapply(seq_len(lengths[g]), function(t) {
      vapply(seq_len(k), function(state) {
        sum(exp(log_path[prefixes[, t] == state] - total))
      }, numeric(1))
    }, numeric(k)))
  })
  scaled <- .latents_forward_backward(emission, layout, log(initial), log(transition))
  logged <- .latents_forward_backward(emission, layout, log(initial),
    lapply(seq_len(n_occasions - 1L), function(t) {
      array(rep(log(transition), each = n_groups), c(n_groups, k, k))
    }))
  expected <- vapply(brute, `[[`, numeric(1), "log_likelihood")
  expect_equal(scaled$log_scaled, expected, tolerance = 1e-13)
  expect_equal(logged$log_scaled, expected, tolerance = 1e-13)
  invisible(lapply(seq_len(n_groups), function(g) {
    invisible(lapply(seq_len(lengths[g]), function(t) {
      expect_equal(exp(scaled$alpha[[t]][g, ] + scaled$beta[[t]][g, ] - expected[g]),
                   brute[[g]]$posterior[, t], tolerance = 1e-12)
      expect_equal(exp(logged$alpha[[t]][g, ] + logged$beta[[t]][g, ] - expected[g]),
                   brute[[g]]$posterior[, t], tolerance = 1e-12)
    }))
  }))
})

test_that("the analytic FIML score is the numerical gradient of the likelihood", {
  skip_if_not_installed("numDeriv")
  skip_if_not_installed("MASS")
  set.seed(4)
  y <- rbind(MASS::mvrnorm(60L, c(-1.5, -1), matrix(c(1, .5, .5, 1), 2L)),
             MASS::mvrnorm(60L, c(1.5, 1), matrix(c(1, -.3, -.3, 1.5), 2L)))
  y <- y[sample(nrow(y)), ]
  y[sample(length(y), 30L)] <- NA
  frame <- data.frame(g = rep(seq_len(20L), each = 6L), a = y[, 1L], b = y[, 2L])
  centre <- colMeans(y, na.rm = TRUE)
  x <- sweep(y, 2L, centre, "-")
  invisible(lapply(c("varying", "equal"), function(variance_model) {
    fit <- suppressMessages(multilpa(frame, c("a", "b"), "g", 2L, 2L,
      n_starts = 2L, seed = 3L, missing = "fiml", covariance_model = "full",
      variance_model = variance_model))
    centred <- fit
    centred$means <- sweep(fit$means, 2L, centre, "-")
    theta <- .multilpa_coefficients(centred, "unconstrained") + 0.03
    log_likelihood <- function(point) {
      .multilpa_expectation(x, fit$group_index, .multilpa_decode(point, fit))$log_likelihood
    }
    numerical <- numDeriv::grad(log_likelihood, theta,
                                method.args = list(eps = 1e-4, d = 1e-3, r = 6))
    expect_equal(unname(-.multilpa_score(theta, x, fit)), numerical, tolerance = 1e-6)
    expect_equal(unname(colSums(.multilpa_group_scores(theta, x, fit))), numerical,
                 tolerance = 1e-6)
  }))
})

test_that("each covariance M-step reproduces mclust's, or attains a lower objective", {
  skip_if_not_installed("mclust")
  skip_if_not_installed("MASS")
  set.seed(11)
  d <- 3L
  k <- 3L
  x <- rbind(MASS::mvrnorm(100L, c(0, 0, 0), matrix(c(1, .5, .2, .5, 2, .3, .2, .3, 1.5), 3L)),
             MASS::mvrnorm(100L, c(3, 1, -1), matrix(c(2, -.4, 0, -.4, 1, .1, 0, .1, .7), 3L)),
             MASS::mvrnorm(100L, c(-2, 3, 2), diag(c(.5, 1.5, 3))))
  z <- matrix(stats::runif(300L * k), 300L, k) +
    3 * outer(rep(seq_len(k), each = 100L), seq_len(k), "==")
  z <- z / rowSums(z)
  weights <- colSums(z)
  means <- t(crossprod(x, z)) / weights
  scatter <- lapply(seq_len(k), function(profile) {
    residuals <- sweep(x, 2L, means[profile, ], "-")
    crossprod(residuals, residuals * z[, profile])
  })
  objective <- function(sigma) sum(vapply(seq_len(k), function(profile) {
    weights[profile] * as.numeric(determinant(sigma[, , profile])$modulus) +
      sum(solve(sigma[, , profile]) * scatter[[profile]])
  }, numeric(1)))
  invisible(lapply(.multilpa_structures(), function(structure) {
    ours <- .multilpa_structure_covariances(scatter, weights, structure, 1e-8)
    reference <- do.call(mclust::mstep, list(data = x, modelName = structure, z = z),
                         envir = asNamespace("mclust"))$parameters$variance$sigma
    expect_lte(objective(ours), objective(reference) + 1e-9 * abs(objective(reference)))
    if (!structure %in% c("VEE", "EVE", "VVE")) {
      # Closed forms, and fixed points solved to tolerance.
      expect_equal(unname(ours), unname(reference), tolerance = 1e-8, label = structure)
    }
    expect_identical(.multilpa_structure_parameters(structure, k, d),
                     as.integer(mclust::nVarParams(structure, d, k)), label = structure)
  }))
})

test_that("a floored axis-parallel M-step stays the constrained maximizer", {
  objective <- function(variances, scatter, weights) {
    sum(weights * rowSums(log(variances))) + sum(scatter / variances)
  }
  # Profile 2 has an effective membership below one: its scatter diagonal is
  # below the floor while its variance is not. The floor does not bind at the
  # maximum, so the closed form must equal the exact unconstrained one.
  weights <- c(50, 0.5)
  scatter <- rbind(c(40, 60), c(8e-7, 1))
  floor_value <- 1e-6
  evi <- .multilpa_structure_variances(scatter, weights, "EVI", floor_value)
  log_det <- rowSums(log(scatter)) / 2
  exact <- (scatter / exp(log_det)) * sum(exp(log_det)) / sum(weights)
  expect_gt(min(exact), floor_value)
  expect_equal(evi, exact, tolerance = 1e-12)
  expect_equal(objective(evi, scatter, weights),
               objective(.multilpa_bounded_variances(scatter, weights, "EVI", floor_value),
                         scatter, weights), tolerance = 1e-13)
  vei <- .multilpa_structure_variances(scatter, weights, "VEI", floor_value)
  expect_equal(vei, .multilpa_bounded_variances(scatter, weights, "VEI", floor_value),
               tolerance = 1e-7)
  # A binding floor: the KKT solve is no worse than constrained Nelder-Mead
  # in the volume/shape chart, and the result stays inside the family.
  weights <- c(40, 30, 0.6)
  scatter <- rbind(c(40, 20, 60), c(30, 60, 3), c(1e-7, 0.3, 0.9))
  floor_value <- 1e-3
  ours <- .multilpa_structure_variances(scatter, weights, "EVI", floor_value)
  design <- .multilpa_eigenvalue_design(3L, 3L, "EVI")
  chart <- function(theta) matrix(exp(design %*% theta), 3L, 3L, byrow = TRUE)
  found <- stats::constrOptim(c(log(2), rep(0, ncol(design) - 1L)),
    function(theta) objective(chart(theta), scatter, weights), NULL,
    ui = design, ci = rep(log(floor_value), nrow(design)) - 1e-12,
    method = "Nelder-Mead", control = list(maxit = 20000L, reltol = 1e-14))
  expect_gte(min(ours), floor_value * (1 - 1e-12))
  expect_lte(objective(ours, scatter, weights), found$value + 1e-9)
  expect_lt(diff(range(rowMeans(log(ours)))), 1e-12)
})

test_that("the noise hypervolume is mclust's, in logs", {
  skip_if_not_installed("mclust")
  set.seed(7)
  x <- matrix(stats::rnorm(150L), 50L, 3L) %*% matrix(c(1, .3, 0, 0, 1, .5, 0, 0, 2), 3L)
  expect_equal(.multilpa_log_hypervolume(x), log(mclust::hypvol(x)), tolerance = 1e-12)
  expect_equal(.multilpa_log_hypervolume(x[, 1L, drop = FALSE]),
               log(diff(range(x[, 1L]))), tolerance = 1e-14)
})
