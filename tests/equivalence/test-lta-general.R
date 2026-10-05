# General latent transition model: occasion-varying and covariate-dependent
# transitions, covariate initial distribution, occasion-specific measurement,
# second order, covariance structures. References: an exact sum over every
# profile path (independent of the forward-backward code), numerical
# derivatives, external fits (depmixS4, LMest, Mplus) bundled as fixtures.

lta_references <- readRDS(test_path("..", "testthat", "fixtures", "lta-references.rds"))

# Small two-profile panel with a time-varying covariate and a group-level one.
lta_small <- local({
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(8)
  n <- 60L
  occasions <- 4L
  w <- stats::rnorm(n)
  rows <- do.call(rbind, lapply(seq_len(n), function(j) {
    z <- stats::rnorm(occasions)
    state <- integer(occasions)
    state[1L] <- 1L + (stats::runif(1L) < stats::plogis(0.3 + 0.8 * w[j]))
    for (t in seq_len(occasions)[-1L]) {  # simulation only
      move <- stats::plogis(-1 + 0.9 * z[t])
      state[t] <- if (stats::runif(1L) < move) 3L - state[t - 1L] else state[t - 1L]
    }
    data.frame(id = j, time = seq_len(occasions), z = z, w = w[j],
               y = stats::rnorm(occasions, c(-1, 1)[state], 0.7))
  }))
  rows
})

# Exact likelihood by summing over every profile path of every group.
lta_brute_force <- function(fit, data) {
  n_states <- fit$n_profiles
  occasions <- max(data$time)
  paths <- as.matrix(expand.grid(rep(list(seq_len(n_states)), occasions)))
  parameters <- .lta_parameters(fit)
  sum(vapply(split(data, data$id), function(rows) {
    rows <- rows[order(rows$time), ]
    j <- match(rows$id[1L], fit$group_ids)
    initial <- exp(.lta_log_initial(fit$designs$initial[j, , drop = FALSE],
                                    matrix(parameters$initial[, , 1L],
                                           ncol(fit$designs$initial))))
    move <- lapply(seq_len(occasions)[-1L], function(t) {
      exp(.lta_log_transition(fit$designs$transition[[t - 1L]][j, , drop = FALSE],
                              array(parameters$transition[, , , 1L],
                                    dim(parameters$transition)[1:3])))[1L, , ]
    })
    block <- function(t) fit$measurement[[if (length(fit$measurement) == 1L) 1L else t]]
    log(sum(apply(paths, 1L, function(s) {
      probability <- initial[s[1L]] *
        prod(vapply(seq_len(occasions)[-1L], function(t) move[[t - 1L]][s[t - 1L], s[t]],
                    numeric(1)))
      density <- prod(vapply(seq_len(occasions), function(t) {
        stats::dnorm(rows$y[t], block(t)$means[s[t], 1L], sqrt(block(t)$variances[s[t], 1L]))
      }, numeric(1)))
      probability * density
    })))
  }, numeric(1)))
}

test_that("analytic scores equal numerical derivatives in every extension", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  fits <- list(
    # Two group classes on 60 groups: a class can have a move that never
    # occurs. The scores are checked at a perturbed point regardless.
    covariates = quietly(lta(lta_small, "y", "id", n_profiles = 2, time = "time",
                             n_group_classes = 2, transition_covariates = "z",
                             initial_covariates = "w", n_starts = 2, seed = 1,
                             max_iter = 20), "latents_boundary"),
    occasion = quietly(lta(lta_small, "y", "id", n_profiles = 2, time = "time",
                           transitions = "occasion", measurement = "occasion",
                           n_starts = 2, seed = 1, max_iter = 20), "latents_boundary"),
    second = lta(lta_small, "y", "id", n_profiles = 2, time = "time", order = 2,
                 n_starts = 2, seed = 1, max_iter = 20))
  lapply(fits, \(fit) {
    set.seed(2)
    theta <- coef(fit) + stats::rnorm(length(coef(fit)), 0, 0.05)
    analytic <- colSums(.lta_group_scores(theta, fit))
    numeric_gradient <- numDeriv::grad(\(v) .lta_expectation(
      fit$x, fit$codes, fit$layout, fit$designs, .lta_unpack(v, fit),
      fit$occasion_of_row)$log_likelihood, theta)
    expect_lt(max(abs(analytic - numeric_gradient)), 1e-5 * (1 + max(abs(numeric_gradient))))
  })
})

test_that("occasion-specific measurement agrees with Mplus (bundled fixture)", {
  skip_on_cran()
  ref <- lta_references$mplus_measurement
  long <- do.call(rbind, lapply(1:2, \(t) {
    frame <- as.data.frame(ref$data[, (t - 1L) * 3L + 1:3])
    names(frame) <- paste0("u", 1:3)
    cbind(person = seq_len(nrow(frame)), occasion = t, frame)
  }))
  long[paste0("u", 1:3)] <- lapply(long[paste0("u", 1:3)], factor)
  fit <- lta(long, paste0("u", 1:3), "person", n_profiles = 2, time = "occasion",
             categorical = paste0("u", 1:3), measurement = "occasion", n_starts = 5,
             seed = 1, tol = 1e-10)
  # Mplus prints the likelihood to three decimals.
  expect_lt(abs(fit$log_likelihood - ref$log_likelihood), 5e-3)
  expect_identical(fit$n_parameters, ref$n_parameters)
})

test_that("missing indicators (FIML) in the extended model match the exact path sum", {
  data <- lta_small
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(2)
  data$y[sample(nrow(data), 30L)] <- NA
  fit <- lta(data, "y", "id", n_profiles = 2, time = "time",
             transition_covariates = "z", initial_covariates = "w", missing = "fiml",
             n_starts = 3, seed = 1, tol = 1e-10)
  p <- .lta_parameters(fit)
  m <- fit$measurement[[1L]]
  paths <- as.matrix(expand.grid(rep(list(1:2), 4L)))
  exact <- sum(vapply(split(data, data$id), \(rows) {
    rows <- rows[order(rows$time), ]
    j <- match(rows$id[1L], fit$group_ids)
    initial <- exp(.lta_log_initial(fit$designs$initial[j, , drop = FALSE],
                                    matrix(p$initial[, , 1L], ncol(fit$designs$initial))))
    moves <- lapply(2:4, \(t) exp(.lta_log_transition(
      fit$designs$transition[[t - 1L]][j, , drop = FALSE],
      array(p$transition[, , , 1L], dim(p$transition)[1:3])))[1L, , ])
    observed <- !is.na(rows$y)
    log(sum(apply(paths, 1L, \(s) {
      initial[s[1L]] * prod(vapply(2:4, \(t) moves[[t - 1L]][s[t - 1L], s[t]], numeric(1))) *
        prod(stats::dnorm(rows$y[observed], m$means[s[observed], 1L],
                          sqrt(m$variances[s[observed], 1L])))
    })))
  }, numeric(1)))
  expect_equal(fit$log_likelihood, exact, tolerance = 1e-9)
  skip_if_not_installed("numDeriv")
  theta <- coef(fit) + 0.05
  analytic <- colSums(.lta_group_scores(theta, fit))
  numeric_gradient <- numDeriv::grad(\(v) .lta_expectation(
    fit$x, fit$codes, fit$layout, fit$designs, .lta_unpack(v, fit),
    fit$occasion_of_row)$log_likelihood, theta)
  expect_lt(max(abs(analytic - numeric_gradient)), 1e-5 * (1 + max(abs(numeric_gradient))))
})
