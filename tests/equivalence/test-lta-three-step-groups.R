# r3step(initial_covariates =, transition_covariates =) on an lta() fit with
# group classes: the step-three model is a mixture of hidden Markov chains for
# the modal assignments, one chain per group class, classification errors
# fixed. Two independent references:
#   1. the likelihood written out as a sum over every group class and every
#      latent path, in base R, at the reported estimate, and its numerical
#      gradient (the estimate is a maximum);
#   2. depmixS4 (Visser and Speekenbrink 2010) on the expanded state space
#      (group class x profile): a block-diagonal chain whose cross-class
#      moves are fixed at probability ~ 0, with the responses fixed at the
#      step-two error rates. Its log likelihood at latents' estimates and its
#      own maximum are compared, and its maximum is mapped back to latents'
#      coefficients.
# The classification inputs (each occasion's modal profile and the error
# matrix) are rebuilt here from the fitted normal densities, not read from
# latents.

group_markov_fixture <- function(n = 600L, occasions = 4L, seed = 41L) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(seed)
  g <- 1L + (stats::runif(n) > 0.5)
  w <- stats::rnorm(n)
  z <- matrix(stats::rnorm(n * occasions), n)
  initial <- rbind(c(1, 0.8), c(-1, 0.8))
  # Per group class and origin: intercept and slope on z of the log odds of
  # moving. Class 1 is stable, class 2 moves often.
  moves <- array(c(-2.0, 0.8, -1.8, -0.6, -0.4, 0.8, -0.2, -0.6), c(2L, 2L, 2L))
  state <- matrix(NA_integer_, n, occasions)
  state[, 1L] <- ifelse(stats::runif(n) < stats::plogis(initial[g, 1L] + initial[g, 2L] * w),
                        1L, 2L)
  # Sequential by construction: each occasion's state depends on the last.
  for (t in seq_len(occasions)[-1L]) {  # simulation only
    from <- state[, t - 1L]
    logit <- moves[cbind(1L, from, g)] + moves[cbind(2L, from, g)] * z[, t]
    moved <- stats::runif(n) < stats::plogis(logit)
    state[, t] <- ifelse(moved, 3L - from, from)
  }
  data <- data.frame(id = rep(seq_len(n), each = occasions),
                     time = rep(seq_len(occasions), n), w = rep(w, each = occasions),
                     z = as.vector(t(z)))
  s <- as.vector(t(state))
  data$y1 <- stats::rnorm(length(s), c(-1.5, 1.5)[s], 1)
  data$y2 <- stats::rnorm(length(s), c(-1.3, 1.3)[s], 1)
  data
}

group_markov_fit <- function(data) {
  lta(data, c("y1", "y2"), "id", n_profiles = 2, time = "time", n_group_classes = 2,
      n_starts = 3, seed = 1)
}

# Modal assignments and error matrix from the fit's own normal densities.
step_two <- function(fit, data) {
  shares <- colMeans(fit$subject_posteriors)
  joint <- vapply(1:2, function(k) {
    shares[k] * stats::dnorm(data$y1, fit$means[k, "y1"], sqrt(fit$variances[k, "y1"])) *
      stats::dnorm(data$y2, fit$means[k, "y2"], sqrt(fit$variances[k, "y2"]))
  }, numeric(nrow(data)))
  posteriors <- joint / rowSums(joint)
  modal <- max.col(posteriors, ties.method = "first")
  errors <- t(vapply(1:2, function(k) {
    vapply(1:2, function(j) sum(posteriors[modal == j, k]), numeric(1)) / sum(posteriors[, k])
  }, numeric(2)))
  list(modal = modal, errors = errors)
}

coefficient <- function(result, theta, block, class, from, term) {
  theta[result$block == block & result$group_class == class &
          (is.na(result$from) | result$from == from) & result$term == term]
}

test_that("the group-class step-three likelihood equals a sum over classes and paths", {
  data <- group_markov_fixture()
  fit <- group_markov_fit(data)
  result <- r3step(fit, data, initial_covariates = "w", transition_covariates = "z")
  expect_named(result, c("block", "group_class", "from", "to", "term", "estimate",
                         "standard_error", "statistic", "p_value", "p_adjusted",
                         "conf_low", "conf_high"))
  expect_identical(nrow(result), 1L + 2L * (2L + 4L))
  two <- step_two(fit, data)
  expect_equal(unname(attr(result, "classification_errors")), unname(two$errors),
               tolerance = 1e-10)
  paths <- as.matrix(expand.grid(rep(list(1:2), 4L)))
  persons <- split(seq_len(nrow(data)), data$id)
  path_likelihood <- function(theta) {
    gamma <- theta[result$block == "group"]
    shares <- c(stats::plogis(gamma), 1 - stats::plogis(gamma))
    sum(vapply(persons, function(rows) {
      log(sum(vapply(1:2, function(h) {
        class <- paste0("group_class_", h)
        first <- stats::plogis(
          coefficient(result, theta, "initial", class, NA, "(Intercept)") +
            coefficient(result, theta, "initial", class, NA, "w") * data$w[rows[1L]])
        initial <- c(first, 1 - first)
        move <- function(from, t) stats::plogis(
          coefficient(result, theta, "transition", class, paste0("profile_", from),
                      "(Intercept)") +
            coefficient(result, theta, "transition", class, paste0("profile_", from), "z") *
            data$z[rows[t]])
        shares[h] * sum(apply(paths, 1L, function(path) {
          steps <- vapply(2:4, function(t) {
            p <- move(path[t - 1L], t)
            if (path[t] == path[t - 1L]) 1 - p else p
          }, numeric(1))
          initial[path[1L]] * prod(steps) * prod(two$errors[cbind(path, two$modal[rows])])
        }))
      }, numeric(1))))
    }, numeric(1)))
  }
  expect_equal(attr(result, "log_likelihood"), path_likelihood(result$estimate),
               tolerance = 1e-9)
  gradient <- vapply(seq_along(result$estimate), function(j) {
    step <- replace(numeric(nrow(result)), j, 1e-5)
    (path_likelihood(result$estimate + step) - path_likelihood(result$estimate - step)) / 2e-5
  }, numeric(1))
  expect_lt(max(abs(gradient)), 1e-3)
  # The step-three group posteriors are the classes' shares of each person's
  # likelihood: rows sum to one, one row per person.
  posteriors <- attr(result, "group_posteriors")
  expect_identical(dim(posteriors), c(600L, 2L))
  expect_equal(unname(rowSums(posteriors)), rep(1, 600L), tolerance = 1e-12)
})

test_that("the group-class step-three chain matches depmixS4 on the expanded states", {
  skip_if_not_installed("depmixS4")
  data <- group_markov_fixture(seed = 42L)
  fit <- group_markov_fit(data)
  ours <- r3step(fit, data, transition_covariates = "z", vcov_type = "observed")
  two <- step_two(fit, data)
  frame <- data.frame(a = factor(two$modal, levels = 1:2), z = data$z)
  # State s = (h - 1) * 2 + k. Cross-class moves are held at log odds -40
  # (probability about 4e-18) against staying.
  n_states <- 4L
  class_of <- (seq_len(n_states) - 1L) %/% 2L + 1L
  profile_of <- (seq_len(n_states) - 1L) %% 2L + 1L
  build <- function() {
    responses <- lapply(seq_len(n_states), function(s) list(depmixS4::GLMresponse(
      a ~ 1, data = frame, family = depmixS4::multinomial("identity"),
      pstart = two$errors[profile_of[s], ])))
    transitions <- lapply(seq_len(n_states), function(s) depmixS4::transInit(
      ~ z, nstates = n_states, data = frame,
      family = depmixS4::multinomial("mlogit", base = s)))
    prior <- depmixS4::transInit(~ 1, nstates = n_states,
                                 family = depmixS4::multinomial("identity"),
                                 data = data.frame(one = rep(1, 600L)),
                                 pstart = rep(0.25, n_states))
    depmixS4::makeDepmix(response = responses, transition = transitions, prior = prior,
                         ntimes = rep(4L, 600L), homogeneous = FALSE)
  }
  model <- build()
  # Parameter layout: 4 prior probabilities; per origin state the intercepts
  # of the 4 destinations, then their z slopes; per state 2 response
  # probabilities.
  transition_index <- function(s, destination, term) {
    n_states + (s - 1L) * 2L * n_states + (term - 1L) * n_states + destination
  }
  # Our coefficients on the expanded states.
  expand <- function(result) {
    theta <- result$estimate
    gamma <- theta[result$block == "group"]
    shares <- c(stats::plogis(gamma), 1 - stats::plogis(gamma))
    parameters <- depmixS4::getpars(model)
    prior <- unlist(lapply(1:2, function(h) {
      first <- stats::plogis(coefficient(result, theta, "initial", paste0("group_class_", h),
                                         NA, "(Intercept)"))
      shares[h] * c(first, 1 - first)
    }))
    parameters[seq_len(n_states)] <- prior
    invisible(lapply(seq_len(n_states), function(s) {
      other <- (class_of[s] - 1L) * 2L + 3L - profile_of[s]
      class <- paste0("group_class_", class_of[s])
      from <- paste0("profile_", profile_of[s])
      across <- setdiff(seq_len(n_states), c(s, other))
      parameters[transition_index(s, other, 1L)] <<-
        coefficient(result, theta, "transition", class, from, "(Intercept)")
      parameters[transition_index(s, other, 2L)] <<-
        coefficient(result, theta, "transition", class, from, "z")
      parameters[transition_index(s, across, 1L)] <<- -40
      parameters[transition_index(s, across, 2L)] <<- 0
      parameters[transition_index(s, s, 1:2)] <<- 0
    }))
    parameters
  }
  at_ours <- depmixS4::setpars(model, expand(ours))
  # depmixS4 reads a transition's covariates on the origin row; latents on
  # the destination row. The fixture's z is shifted accordingly below, so
  # first establish which row convention reproduces the likelihood.
  depmix_ll <- as.numeric(depmixS4::logLik(at_ours))
  if (abs(depmix_ll - attr(ours, "log_likelihood")) > 1e-6) {
    shifted <- stats::ave(data$z, data$id, FUN = function(v) c(v[-1L], 0))
    frame$z <- shifted
    model <- build()
    at_ours <- depmixS4::setpars(model, expand(ours))
    depmix_ll <- as.numeric(depmixS4::logLik(at_ours))
  }
  expect_lt(abs(depmix_ll - attr(ours, "log_likelihood")), 1e-6)
  # The maximum of depmixS4's likelihood: the 11 free coordinates (prior
  # logits against the last state, and each origin's within-class move
  # intercept and slope) by a generic optimizer, everything else held,
  # started away from ours. depmixS4's own fit() with `fixed` runs Rsolnp,
  # which stops short here ("solution not reliable").
  moving <- unlist(lapply(seq_len(n_states), function(s) {
    transition_index(s, (class_of[s] - 1L) * 2L + 3L - profile_of[s], 1:2)
  }))
  held <- expand(ours)
  depmix_parameters <- function(phi) {
    prior <- exp(c(phi[1:3], 0))
    replace(held, c(seq_len(n_states), moving), c(prior / sum(prior), phi[-(1:3)]))
  }
  negative_ll <- function(phi) {
    -as.numeric(depmixS4::logLik(depmixS4::setpars(model, depmix_parameters(phi))))
  }
  at_ours_phi <- c(log(held[1:3] / held[4]), held[moving])
  start <- at_ours_phi + rep(c(0.3, -0.2), length.out = length(at_ours_phi))
  searched <- stats::optim(start, negative_ll, method = "BFGS",
                           control = list(maxit = 2000L, reltol = 1e-14))
  searched <- stats::optim(searched$par, negative_ll, method = "Nelder-Mead",
                           control = list(maxit = 20000L, reltol = 1e-14))
  searched <- stats::optim(searched$par, negative_ll, method = "BFGS",
                           control = list(maxit = 2000L, reltol = 1e-14, ndeps = rep(1e-5, 11L)))
  theirs <- depmix_parameters(searched$par)
  expect_lt(abs(-searched$value - attr(ours, "log_likelihood")), 1e-5)
  expect_gte(attr(ours, "log_likelihood"), -searched$value - 1e-6)
  # Map depmixS4's maximum back: group shares and initial logits from the
  # prior, transition logits directly.
  prior <- theirs[seq_len(n_states)]
  gap <- c(
    log(sum(prior[1:2]) / sum(prior[3:4])) - ours$estimate[ours$block == "group"],
    vapply(1:2, function(h) {
      log(prior[2L * h - 1L] / prior[2L * h]) -
        coefficient(ours, ours$estimate, "initial", paste0("group_class_", h), NA,
                    "(Intercept)")
    }, numeric(1)),
    unlist(lapply(seq_len(n_states), function(s) {
      other <- (class_of[s] - 1L) * 2L + 3L - profile_of[s]
      class <- paste0("group_class_", class_of[s])
      from <- paste0("profile_", profile_of[s])
      c(theirs[transition_index(s, other, 1L)] -
          coefficient(ours, ours$estimate, "transition", class, from, "(Intercept)"),
        theirs[transition_index(s, other, 2L)] -
          coefficient(ours, ours$estimate, "transition", class, from, "z"))
    })))
  expect_lt(max(abs(gap)), 2e-3)
  # The generic search is a cross-check on the maximum; the decisive check is
  # that depmixS4's likelihood has no ascent direction at latents' estimate.
  gradient <- vapply(seq_along(at_ours_phi), function(j) {
    step <- replace(numeric(length(at_ours_phi)), j, 1e-5)
    (negative_ll(at_ours_phi + step) - negative_ll(at_ours_phi - step)) / 2e-5
  }, numeric(1))
  expect_lt(max(abs(gradient)), 1e-3)
})
