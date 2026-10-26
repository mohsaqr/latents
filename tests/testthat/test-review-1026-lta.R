# Calibration and invariance checks of the latent transition engines against
# a brute-force sum over every state path (Rabiner, 1989): for a person with
# T occasions and K profiles, L = sum over the K^T paths s of
# pi(s_1) prod_t A_t(s_{t-1}, s_t) prod_t f(y_t | s_t), with an occasion that
# has no observed indicator contributing f = 1.

review_lta_data <- function(n = 30L, n_occasions = 3L, seed = 1L, stay = 0.75) {
  withr::with_seed(seed, {
    do.call(rbind, lapply(seq_len(n), function(i) {
      states <- Reduce(function(previous, t) {
        if (stats::runif(1L) < stay) previous else 3L - previous
      }, seq_len(n_occasions - 1L), init = sample.int(2L, 1L), accumulate = TRUE)
      data.frame(id = i, time = seq_len(n_occasions), z = stats::rnorm(1L),
                 w = stats::rnorm(n_occasions),
                 y1 = stats::rnorm(n_occasions, c(0, 2)[states], 0.7),
                 y2 = stats::rnorm(n_occasions, c(0, 1.5)[states], 0.7))
    }))
  })
}

# Diagonal Gaussian log densities, rows x K; an all-missing row gives 0.
review_gauss <- function(data, means, variances) {
  y <- as.matrix(data[c("y1", "y2")])
  vapply(seq_len(nrow(means)), function(k) {
    dens <- stats::dnorm(y, rep(means[k, ], each = nrow(y)),
                         rep(sqrt(variances[k, ]), each = nrow(y)), log = TRUE)
    rowSums(matrix(dens, nrow(y)), na.rm = TRUE)
  }, numeric(nrow(y)))
}

review_rows <- function(data) {
  lapply(split(seq_len(nrow(data)), data$id), function(r) r[order(data$time[r])])
}

# Enumerate every path of each person: per-person log likelihood, occupancy
# posteriors (occasions x K) and expected move counts (K x K).
review_brute <- function(log_density, initial, transition) {
  n_profiles <- ncol(log_density[[1L]])
  lapply(seq_along(log_density), function(g) {
    dens <- log_density[[g]]
    n_t <- nrow(dens)
    paths <- as.matrix(expand.grid(rep(list(seq_len(n_profiles)), n_t)))
    log_path <- vapply(seq_len(nrow(paths)), function(r) {
      s <- paths[r, ]
      Reduce(function(total, t) {
        total + log(transition(g, t)[s[t - 1L], s[t]]) + dens[t, s[t]]
      }, seq_len(n_t)[-1L], init = log(initial(g)[s[1L]]) + dens[1L, s[1L]])
    }, numeric(1))
    top <- max(log_path)
    weight <- exp(log_path - top) / sum(exp(log_path - top))
    occupancy <- t(vapply(seq_len(n_t), function(t) {
      vapply(seq_len(n_profiles), function(k) sum(weight[paths[, t] == k]), numeric(1))
    }, numeric(n_profiles)))
    moves <- outer(seq_len(n_profiles), seq_len(n_profiles), Vectorize(function(k, l) {
      sum(vapply(seq_len(n_t)[-1L], function(t) {
        sum(weight[paths[, t - 1L] == k & paths[, t] == l])
      }, numeric(1)))
    }))
    list(log_likelihood = top + log(sum(exp(log_path - top))),
         occupancy = matrix(occupancy, n_t), moves = moves)
  })
}

review_softmax <- function(eta) exp(eta - max(eta)) / sum(exp(eta - max(eta)))

# Initial logits against the last profile; transition logits against staying.
review_initial <- function(coefficients, design_row) {
  review_softmax(c(drop(design_row %*% coefficients), 0))
}
review_transition <- function(coefficients, design_row) {
  n_profiles <- dim(coefficients)[3L]
  t(vapply(seq_len(n_profiles), function(k) {
    eta <- numeric(n_profiles)
    eta[-k] <- drop(design_row %*% matrix(coefficients[, , k], length(design_row)))
    review_softmax(eta)
  }, numeric(n_profiles)))
}

review_place <- function(brute, rows, n_rows) {
  out <- matrix(0, n_rows, ncol(brute[[1L]]$occupancy))
  out[unlist(rows), ] <- do.call(rbind, lapply(brute, `[[`, "occupancy"))
  out
}

test_that("homogeneous lta() matches the brute-force path sum on unequal sequences", {
  data <- review_lta_data(n = 30L, n_occasions = 4L, seed = 3L)
  data <- subset(data, !(id %in% 1:5 & time == 4) & !(id %in% 6:8 & time >= 3))
  data <- subset(data, !(id %in% 9:10 & time >= 2))  # length-one sequences
  fit <- lta(data, vars = c("y1", "y2"), id = "id", n_profiles = 2, time = "time",
             n_starts = 2, seed = 1)
  rows <- review_rows(data)
  density <- review_gauss(data, fit$means, fit$variances)
  brute <- review_brute(lapply(rows, function(r) density[r, , drop = FALSE]),
                        function(g) fit$initial_probabilities[1L, ],
                        function(g, t) fit$transition_probabilities[, , 1L])
  per_person <- vapply(brute, `[[`, numeric(1), "log_likelihood")
  expect_equal(fit$log_likelihood, sum(per_person), tolerance = 1e-10)
  expect_equal(unname(fit$group_log_likelihood), per_person, tolerance = 1e-10)
  expect_equal(unname(fit$subject_posteriors),
               review_place(brute, rows, nrow(data)), tolerance = 1e-10)
  expect_equal(unname(fit$transition_counts[, , 1L]),
               Reduce(`+`, lapply(brute, `[[`, "moves")), tolerance = 1e-10)
  expect_equal(unname(rowSums(fit$subject_posteriors)), rep(1, nrow(data)),
               tolerance = 1e-12)
  expect_equal(unname(rowSums(fit$transition_probabilities[, , 1L])), c(1, 1),
               tolerance = 1e-12)
})

test_that("covariate transitions and initial states match the brute-force path sum", {
  data <- review_lta_data(n = 30L, n_occasions = 3L, seed = 5L)
  data <- subset(data, !(id %in% 1:6 & time == 3))
  fit <- suppressWarnings(lta(
    data, vars = c("y1", "y2"), id = "id", n_profiles = 2, time = "time",
    n_starts = 2, seed = 1, transitions = "occasion", transition_covariates = "w",
    initial_covariates = "z"))
  rows <- review_rows(data)
  block <- fit$measurement[[1L]]
  density <- review_gauss(data, block$means, block$variances)
  transition <- array(fit$transition_coefficients[, , , 1L],
                      dim(fit$transition_coefficients)[1:3])
  initial <- matrix(fit$initial_coefficients[, , 1L], 2L)
  brute <- review_brute(
    lapply(rows, function(r) density[r, , drop = FALSE]),
    function(g) review_initial(initial, c(1, data$z[rows[[g]][1L]])),
    function(g, t) review_transition(transition, c(as.numeric(2:3 == t),
                                                   data$w[rows[[g]][t]])))
  expect_equal(fit$log_likelihood, sum(vapply(brute, `[[`, numeric(1), "log_likelihood")),
               tolerance = 1e-10)
  expect_equal(unname(fit$subject_posteriors), review_place(brute, rows, nrow(data)),
               tolerance = 1e-10)
  pairs <- Reduce(`+`, lapply(fit$expectation$moments[[1L]]$pairs,
                              function(p) apply(p, c(2L, 3L), sum)))
  expect_equal(pairs, Reduce(`+`, lapply(brute, `[[`, "moves")), tolerance = 1e-10)
  # Every model-implied transition row is a distribution.
  probabilities <- get_results(fit, "transitions")
  totals <- aggregate(probability ~ occasion + from, data = probabilities, FUN = sum)
  expect_equal(totals$probability, rep(1, nrow(totals)), tolerance = 1e-12)
})

test_that("second-order and mover-stayer likelihoods match the brute-force path sum", {
  data <- review_lta_data(n = 30L, n_occasions = 4L, seed = 7L, stay = 0.65)
  data <- subset(data, !(id %in% 1:5 & time == 4))
  rows <- review_rows(data)
  second <- suppressWarnings(lta(data, vars = c("y1", "y2"), id = "id",
                                 n_profiles = 2, time = "time", n_starts = 2,
                                 seed = 1, order = 2))
  block <- second$measurement[[1L]]
  density <- lapply(rows, function(r) {
    review_gauss(data, block$means, block$variances)[r, , drop = FALSE]
  })
  first_order <- review_transition(array(second$transition_coefficients[, , , 1L],
                                         dim(second$transition_coefficients)[1:3]), 1)
  later <- second$second_order_coefficients
  # A path's log probability under the second-order chain: pair (i, j) is
  # index (i - 1) K + j, its logits against staying in j.
  path_log <- function(dens, s) {
    start <- log(review_initial(matrix(second$initial_coefficients[, , 1L], 1L), 1)[s[1L]]) +
      dens[1L, s[1L]]
    if (length(s) < 2L) return(start)
    start <- start + log(first_order[s[1L], s[2L]]) + dens[2L, s[2L]]
    Reduce(function(total, t) {
      eta <- numeric(2L)
      eta[-s[t - 1L]] <- later[1L, , (s[t - 2L] - 1L) * 2L + s[t - 1L], 1L]
      total + log(review_softmax(eta)[s[t]]) + dens[t, s[t]]
    }, seq_len(length(s))[-(1:2)], init = start)
  }
  brute_second <- vapply(density, function(dens) {
    paths <- as.matrix(expand.grid(rep(list(1:2), nrow(dens))))
    values <- vapply(seq_len(nrow(paths)), function(r) path_log(dens, paths[r, ]), numeric(1))
    max(values) + log(sum(exp(values - max(values))))
  }, numeric(1))
  expect_equal(second$log_likelihood, sum(brute_second), tolerance = 1e-10)

  stayers <- suppressWarnings(lta(data, vars = c("y1", "y2"), id = "id",
                                  n_profiles = 2, time = "time", n_starts = 2,
                                  seed = 1, mover_stayer = TRUE))
  block <- stayers$measurement[[1L]]
  density <- lapply(rows, function(r) {
    review_gauss(data, block$means, block$variances)[r, , drop = FALSE]
  })
  by_class <- vapply(seq_len(2L), function(h) {
    moves <- if (stayers$stayer[[h]]) diag(2L) else
      review_transition(array(stayers$transition_coefficients[, , , h],
                              dim(stayers$transition_coefficients)[1:3]), 1)
    vapply(review_brute(density,
                        function(g) review_initial(matrix(stayers$initial_coefficients[, , h], 1L), 1),
                        function(g, t) moves), `[[`, numeric(1), "log_likelihood")
  }, numeric(length(rows)))
  mixed <- log(exp(by_class) %*% stayers$group_probabilities)
  expect_equal(stayers$log_likelihood, sum(mixed), tolerance = 1e-10)
})

test_that("the brute-force likelihood is invariant to relabelling the profiles", {
  data <- review_lta_data(n = 25L, n_occasions = 3L, seed = 11L)
  fit <- lta(data, vars = c("y1", "y2"), id = "id", n_profiles = 2, time = "time",
             n_starts = 1, seed = 1)
  permuted <- .lta_permute_transitions(fit, c(2L, 1L), 1L)
  rows <- review_rows(data)
  total <- function(x) {
    density <- review_gauss(data, x$means, x$variances)
    sum(vapply(review_brute(lapply(rows, function(r) density[r, , drop = FALSE]),
                            function(g) x$initial_probabilities[1L, ],
                            function(g, t) x$transition_probabilities[, , 1L]),
               `[[`, numeric(1), "log_likelihood"))
  }
  expect_equal(total(permuted), total(fit), tolerance = 1e-12)
  expect_equal(total(fit), fit$log_likelihood, tolerance = 1e-10)
})

test_that("relabelled general-fit logits imply the permuted probabilities", {
  data <- review_lta_data(n = 30L, n_occasions = 3L, seed = 13L)
  fit <- suppressWarnings(lta(data, vars = c("y1", "y2"), id = "id", n_profiles = 2,
                              time = "time", n_starts = 1, seed = 1,
                              transition_covariates = "w", initial_covariates = "z"))
  estimates <- .lta_bootstrap_estimates(fit)
  implied <- function(values, w, z) {
    pick <- function(name) values[[name]]
    moves <- t(vapply(1:2, function(k) {
      l <- 3L - k
      prefix <- sprintf("transition.group_class_1.profile_%d->profile_%d", k, l)
      eta <- numeric(2L)
      eta[l] <- pick(paste0(prefix, ".(Intercept)")) + w * pick(paste0(prefix, ".w"))
      review_softmax(eta)
    }, numeric(2)))
    start <- review_softmax(c(pick("initial.group_class_1.profile_1.(Intercept)") +
                                z * pick("initial.group_class_1.profile_1.z"), 0))
    list(moves = moves, start = start)
  }
  original <- implied(estimates, 0.4, -1.1)
  relabelled <- implied(.lta_relabel_profiles(estimates, c(2L, 1L)), 0.4, -1.1)
  expect_equal(relabelled$moves, original$moves[2:1, 2:1], tolerance = 1e-12)
  expect_equal(relabelled$start, original$start[2:1], tolerance = 1e-12)
})

test_that("an occasion with no observed indicator leaves the likelihood unchanged", {
  data <- review_lta_data(n = 30L, n_occasions = 3L, seed = 17L)
  trailing <- rbind(data, data.frame(id = 1:8, time = 4L, z = 0, w = 0,
                                     y1 = NA_real_, y2 = NA_real_))
  fit <- lta(data, vars = c("y1", "y2"), id = "id", n_profiles = 2, time = "time",
             n_starts = 1, seed = 1, missing = "fiml")
  padded <- lta(trailing, vars = c("y1", "y2"), id = "id", n_profiles = 2,
                time = "time", n_starts = 1, seed = 1, missing = "fiml")
  expect_equal(padded$log_likelihood, fit$log_likelihood, tolerance = 1e-8)
  expect_equal(padded$transition_probabilities, fit$transition_probabilities,
               tolerance = 1e-5)
  expect_identical(nobs(padded), nobs(fit))
  # Under the grid, a skipped position and an explicit all-missing row are
  # the same missing emission.
  skipped <- subset(data, !(id %in% 1:6 & time == 2))
  explicit <- data
  explicit[explicit$id %in% 1:6 & explicit$time == 2, c("y1", "y2")] <- NA
  grid_skip <- lta(skipped, vars = c("y1", "y2"), id = "id", n_profiles = 2,
                   time = "time", n_starts = 1, seed = 1, missing = "fiml",
                   occasions = "grid")
  grid_na <- lta(explicit, vars = c("y1", "y2"), id = "id", n_profiles = 2,
                 time = "time", n_starts = 1, seed = 1, missing = "fiml",
                 occasions = "grid")
  expect_equal(grid_na$log_likelihood, grid_skip$log_likelihood, tolerance = 1e-8)
  expect_equal(BIC(grid_na), grid_na$bic, tolerance = 1e-12)
  expect_equal(BIC(grid_na), BIC(grid_skip), tolerance = 1e-6)
})

test_that("the step-three chain matches its brute-force likelihood and numerical scores", {
  data <- review_lta_data(n = 40L, n_occasions = 3L, seed = 19L)
  # An all-missing middle occasion and a wholly missing person.
  data[data$id %in% 1:4 & data$time == 2, c("y1", "y2")] <- NA
  data[data$id == 40, c("y1", "y2")] <- NA
  fit <- lta(data, vars = c("y1", "y2"), id = "id", n_profiles = 2, time = "time",
             n_starts = 1, seed = 1, missing = "fiml")
  inputs <- .lta_three_step_inputs(fit, data, "z", "w")
  observed <- rowSums(!is.na(data[c("y1", "y2")])) > 0
  expect_true(all(is.na(inputs$modal[!observed])))
  expect_false(anyNA(inputs$modal[observed]))
  expect_identical(sum(inputs$informative_groups), 39L)
  # Step two: pooled occasion posteriors and the classification-error matrix
  # D[k, j] = sum_i p_ik 1(m_i = j) / sum_i p_ik over observed occasions.
  shares <- colMeans(fit$subject_posteriors[observed, ])
  joint <- sweep(exp(review_gauss(data, fit$means, fit$variances)), 2L, shares, "*")
  posterior <- joint / rowSums(joint)
  expect_equal(.lta_occasion_posteriors(fit), posterior, tolerance = 1e-12)
  modal <- max.col(posterior[observed, ], ties.method = "first")
  errors <- outer(1:2, 1:2, Vectorize(function(k, j) {
    sum(posterior[observed, k] * (modal == j)) / sum(posterior[observed, k])
  }))
  expect_equal(unname(inputs$errors), errors, tolerance = 1e-12)
  expect_equal(unname(rowSums(inputs$errors)), c(1, 1), tolerance = 1e-12)
  # Step three at an arbitrary point.
  theta <- c(0.3, -0.4, -1.2, 0.5, -0.8, -0.6)
  evaluated <- .lta_three_step_evaluate(theta, inputs)
  rows <- review_rows(data)
  emission <- lapply(rows, function(r) {
    t(vapply(r, function(i) if (is.na(inputs$modal[i])) c(0, 0) else
      log(errors[, inputs$modal[i]]), numeric(2)))
  })
  brute <- vapply(review_brute(
    emission, function(g) review_initial(matrix(theta[1:2], 2L), c(1, data$z[rows[[g]][1L]])),
    function(g, t) review_transition(array(theta[3:6], c(2L, 1L, 2L)),
                                     c(1, data$w[rows[[g]][t]]))),
    `[[`, numeric(1), "log_likelihood")
  expect_equal(evaluated$log_likelihood, brute, tolerance = 1e-10)
  # The wholly missing person contributes neither likelihood nor score.
  expect_equal(evaluated$log_likelihood[40L], 0, tolerance = 1e-12)
  expect_equal(unname(evaluated$scores[40L, ]), rep(0, 6L), tolerance = 1e-12)
  skip_if_not_installed("numDeriv")
  numerical <- numDeriv::jacobian(function(v) {
    .lta_three_step_evaluate(v, inputs)$log_likelihood
  }, theta)
  expect_equal(unname(evaluated$scores), numerical, tolerance = 1e-6)
})

test_that("nested transition comparisons refuse different fitting bounds", {
  data <- review_lta_data(n = 25L, n_occasions = 3L, seed = 29L)
  null <- lta(data, vars = c("y1", "y2"), id = "id", n_profiles = 2, time = "time",
              n_starts = 1, seed = 1, variance_model = "equal")
  alternative <- lta(data, vars = c("y1", "y2"), id = "id", n_profiles = 2,
                     time = "time", n_starts = 1, seed = 1, min_variance = 1e-4)
  same_bounds <- lta(data, vars = c("y1", "y2"), id = "id", n_profiles = 2,
                     time = "time", n_starts = 1, seed = 1)
  expect_true(.lta_nested(null, same_bounds))
  expect_false(.lta_nested(null, alternative))
  expect_error(bootstrap_lrt(null, alternative, data = data, iter = 2L, seed = 1),
               class = "latents_bad_argument")
})
