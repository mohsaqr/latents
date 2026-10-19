# Three-step latent Markov model with covariates (Bolck, Croon and Hagenaars
# 2004; Vermunt 2010; Di Mari, Oberski and Vermunt 2016). Step one is a fitted
# lta(); step two assigns every occasion to its modal profile and estimates
# the classification-error matrix D[k, j] = P(assigned j | profile k) from the
# step-one posteriors; step three is a hidden Markov model whose single
# indicator is the assigned profile with response probabilities held at D,
# and whose initial and transition logits depend on covariates. The logits
# use lta()'s coding: initial logits against the last profile, transition
# logits against staying, so a coefficient is the log odds of moving to that
# profile rather than staying (Latent GOLD's `~tra` coding).

#' The data of a three-step latent Markov model
#'
#' @param x A covariate-free `lta()` fit with one group class.
#' @param data The data frame it was fitted to.
#' @param initial_covariates,transition_covariates Covariate names.
#' @return A list with `layout`, `designs`, `emission`, `errors`, `n_profiles`
#'   and `n_observations`.
#' @noRd
.lta_three_step_inputs <- function(x, data, initial_covariates, transition_covariates) {
  unsupported <- function(message) {
    stop(errorCondition(message, class = "latents_unsupported_three_step", call = NULL))
  }
  if (x$n_group_classes > 1L) {
    unsupported(paste(
      "A three-step transition model is implemented for one group class: the",
      "step-three chain would also need the group classes."))
  }
  if (inherits(x, "multilpa_lta") &&
      (length(x$arguments$initial_covariates) > 0L ||
       length(x$arguments$transition_covariates) > 0L || isTRUE(x$mover_stayer))) {
    unsupported(paste(
      "The step-one transition model already has covariates (or a stayer",
      "class): its classification errors are conditional on them, which the",
      "fixed error matrix of the three-step correction does not allow. Fit",
      "the step-one model without them."))
  }
  .lta_check_data(x, data)
  a <- .lta_arguments(x)
  reused <- intersect(c(initial_covariates, transition_covariates),
                      unique(c(a$vars, .multilpa_measured_items(x))))
  if (length(reused) > 0L) {
    stop(errorCondition(sprintf(
      "Measurement indicator(s) %s already helped define the profiles and cannot be used as covariates.",
      paste(sprintf("`%s`", reused), collapse = ", ")),
      class = c("multilpa_indicator_reused", "latents_bad_covariate"), call = NULL))
  }
  groups <- .multilpa_prepare_groups(data[[a$id]])
  time_values <- .multilpa_time_values(data, a$time, a$id, a$vars)
  layout <- .multilpa_sequence_layout(groups$index, time_values, groups$n, a$occasions)
  .lta_check_covariate_gaps(layout, groups, time_values, transition_covariates)
  initial_rows <- .lta_covariate_matrix(data, initial_covariates)
  if (!is.null(initial_rows)) {
    by_group <- split(as.data.frame(initial_rows), groups$index)
    constant <- vapply(by_group, function(rows) {
      all(vapply(rows, function(v) length(unique(v)) == 1L, logical(1)))
    }, logical(1))
    if (!all(constant)) {
      stop(errorCondition(paste(
        "`initial_covariates` must be constant within each person: they",
        "predict the profile at the first occasion."),
        class = "latents_bad_covariate", call = NULL))
    }
  }
  designs <- .lta_designs(layout, .lta_covariate_matrix(data, transition_covariates),
                          initial_rows, "homogeneous")
  posteriors <- .lta_occasion_posteriors(x)
  pieces <- list(posteriors = posteriors,
                 modal = max.col(posteriors, ties.method = "first"),
                 n_classes = ncol(posteriors))
  errors <- .multilpa_error_matrix(pieces)
  # Each occasion's indicator is its assigned profile, observed with the
  # step-two error rates: the log density under profile k is log D[k, m].
  log_density <- log(pmax(t(errors)[pieces$modal, , drop = FALSE], 1e-300))
  list(layout = layout, designs = designs,
       emission = .multilpa_sequence_emission(log_density, layout),
       errors = errors, modal = pieces$modal, n_profiles = ncol(posteriors),
       n_observations = nrow(posteriors))
}

#' Each occasion's profile posteriors from its own indicators
#'
#' The step-two classification of the three-step latent Markov model reads
#' one occasion at a time: p(k | y_it) is proportional to w_k f_k(y_it), with
#' the transition fit's measurement densities f_k and its profile shares w_k
#' pooled over occasions, the latent class model of the pooled data (Di Mari,
#' Oberski and Vermunt 2016). The transition fit's own posteriors would not
#' do: they borrow from the person's other occasions through the chain, so
#' their errors are not independent across occasions given the profiles, and
#' the step-three model, which assumes they are, then removes the apparent
#' moves the smoothing already removed.
#'
#' @param x A transition fit.
#' @return An observations-by-profiles matrix whose rows sum to one.
#' @noRd
.lta_occasion_posteriors <- function(x) {
  general <- inherits(x, "multilpa_lta")
  posteriors <- if (general) x$expectation$subject_posteriors else x$subject_posteriors
  log_density <- if (general) {
    .lta_log_density(x$x, x$codes, .lta_parameters(x), x$occasion_of_row, x$extra_data)
  } else {
    .latents_measurement_log_density(x$indicator_data, x, x$categorical_data,
                                     extra = x$extra_data, offset = FALSE)$log_density
  }
  joint <- sweep(log_density, 2L, log(colMeans(posteriors)), "+")
  joint <- exp(joint - .multilpa_row_max(joint))
  joint / rowSums(joint)
}

#' Centre and scale a logit design's covariate columns
#' @param design A design matrix whose first column is the intercept.
#' @param rows Which rows carry information (inside a person's span).
#' @return The `p x p` matrix `T` with `design %*% T` standardized: natural
#'   coefficients are `T %*% standardized`.
#' @noRd
.lta_three_step_basis <- function(design, rows = seq_len(nrow(design))) {
  p <- ncol(design)
  basis <- diag(p)
  if (p == 1L) return(basis)
  values <- design[rows, -1L, drop = FALSE]
  centre <- colMeans(values)
  spread <- sqrt(colMeans(sweep(values, 2L, centre)^2))
  if (any(!is.finite(spread)) || any(spread <= sqrt(.Machine$double.eps) *
                                       pmax(1, abs(centre)))) {
    stop(errorCondition(
      "A covariate is constant over the occasions it applies to; remove it.",
      class = "latents_bad_covariate", call = NULL))
  }
  basis[1L, -1L] <- -centre / spread
  basis[-1L, -1L] <- diag(1 / spread, p - 1L)
  basis
}

#' Log likelihood and per-person scores of the step-three chain
#'
#' The score of a multinomial logit with known emissions is the design times
#' the posterior count minus its expectation: for the initial logits,
#' P(s_1 = l | data) - P(s_1 = l); for a move out of k,
#' P(s_{t-1} = k, s_t = l | data) - P(s_{t-1} = k | data) P(l | k).
#'
#' @param theta Coefficients on the standardized designs.
#' @param inputs From `.lta_three_step_inputs()`, with standardized designs.
#' @return A list with `log_likelihood` (per person) and `scores` (persons by
#'   coefficients).
#' @noRd
.lta_three_step_evaluate <- function(theta, inputs) {
  designs <- inputs$designs
  layout <- inputs$layout
  k_profiles <- inputs$n_profiles
  p0 <- ncol(designs$initial)
  p <- ncol(designs$transition[[1L]])
  initial <- matrix(theta[seq_len(p0 * (k_profiles - 1L))], p0, k_profiles - 1L)
  transition <- array(theta[-seq_len(p0 * (k_profiles - 1L))],
                      c(p, k_profiles - 1L, k_profiles))
  log_initial <- .lta_log_initial(designs$initial, initial)
  log_transition <- lapply(designs$transition, .lta_log_transition,
                           coefficients = transition)
  pass <- .latents_forward_backward(inputs$emission$log_density, layout, log_initial,
                                    log_transition)
  n_groups <- nrow(layout$slot)
  moments <- .lta_moments(pass, inputs$emission$log_density, layout, log_transition,
                          rep(1, n_groups), inputs$n_observations)
  residual <- moments$initial - exp(log_initial)
  initial_scores <- do.call(cbind, lapply(seq_len(k_profiles - 1L), function(l) {
    designs$initial * residual[, l]
  }))
  transition_scores <- do.call(cbind, lapply(seq_len(k_profiles), function(k) {
    others <- setdiff(seq_len(k_profiles), k)
    do.call(cbind, lapply(others, function(l) {
      Reduce(`+`, lapply(seq_along(moments$pairs), function(t) {
        counts <- matrix(moments$pairs[[t]][, k, ], n_groups)
        expected <- rowSums(counts) * exp(log_transition[[t]][, k, l])
        designs$transition[[t]] * (counts[, l] - expected)
      }))
    }))
  }))
  list(log_likelihood = pass$log_scaled + inputs$emission$offset,
       scores = cbind(initial_scores, transition_scores))
}

#' Fit the step-three chain and its covariance
#'
#' @param inputs From `.lta_three_step_inputs()`.
#' @param vcov_type `"robust"` (the sandwich over persons) or `"observed"`.
#' @return A list with `estimates`, `covariance`, `log_likelihood` and `labels`
#'   (a data frame: `block`, `from`, `to`, `term`).
#' @noRd
.lta_three_step_fit <- function(inputs, vcov_type) {
  designs <- inputs$designs
  layout <- inputs$layout
  k_profiles <- inputs$n_profiles
  inside <- lapply(seq_along(designs$transition), function(t) layout$within[, t + 1L])
  stacked <- do.call(rbind, designs$transition)
  initial_basis <- .lta_three_step_basis(designs$initial)
  transition_basis <- .lta_three_step_basis(stacked, unlist(inside))
  standardized <- inputs
  standardized$designs <- list(
    initial = designs$initial %*% initial_basis,
    transition = lapply(designs$transition, function(d) d %*% transition_basis))
  p0 <- ncol(designs$initial)
  p <- ncol(stacked)
  # Starts: intercepts from the assigned profiles' own counts (half a count
  # added to every cell), covariate effects at zero.
  first <- inputs$modal[layout$slot[, 1L]]
  first_counts <- tabulate(first, k_profiles) + 0.5
  initial_start <- matrix(0, p0, k_profiles - 1L)
  initial_start[1L, ] <- log(first_counts[-k_profiles] / first_counts[k_profiles])
  moves <- do.call(rbind, lapply(seq_len(layout$n_occasions)[-1L], function(t) {
    observed <- !is.na(layout$slot[, t - 1L]) & !is.na(layout$slot[, t])
    cbind(inputs$modal[layout$slot[observed, t - 1L]],
          inputs$modal[layout$slot[observed, t]])
  }))
  move_counts <- matrix(tabulate((moves[, 1L] - 1L) * k_profiles + moves[, 2L],
                                 k_profiles^2), k_profiles, byrow = TRUE) + 0.5
  transition_start <- array(0, c(p, k_profiles - 1L, k_profiles))
  invisible(lapply(seq_len(k_profiles), function(k) {
    others <- setdiff(seq_len(k_profiles), k)
    transition_start[1L, , k] <<- log(move_counts[k, others] / move_counts[k, k])
  }))
  # One basis block per logit column: initial columns first, then every
  # origin's destinations, in the order the coefficients are vectorized.
  initial_block <- kronecker(diag(k_profiles - 1L), initial_basis)
  transition_block <- kronecker(diag((k_profiles - 1L) * k_profiles), transition_basis)
  natural_basis <- rbind(
    cbind(initial_block, matrix(0, nrow(initial_block), ncol(transition_block))),
    cbind(matrix(0, nrow(transition_block), ncol(initial_block)), transition_block))
  start <- solve(natural_basis, c(as.vector(initial_start), as.vector(transition_start)))
  objective <- function(theta) -sum(.lta_three_step_evaluate(theta, standardized)$log_likelihood)
  gradient <- function(theta) -colSums(.lta_three_step_evaluate(theta, standardized)$scores)
  fitted <- stats::optim(start, objective, gradient, method = "BFGS",
                         control = list(maxit = 2000L, reltol = 1e-14))
  scores <- .lta_three_step_evaluate(fitted$par, standardized)$scores
  if (!identical(fitted$convergence, 0L) ||
      max(abs(colSums(scores))) > 1e-3 * max(1, sqrt(nrow(scores)))) {
    stop(errorCondition(sprintf(paste(
      "The step-three transition model did not converge (optim code %d, largest",
      "score %.3g)."), fitted$convergence, max(abs(colSums(scores)))),
      class = "latents_no_converge", call = NULL))
  }
  total <- function(theta) colSums(.lta_three_step_evaluate(theta, standardized)$scores)
  information <- .inference_score_information(total, fitted$par, 1e-5, scheme = "central",
                                              h_rule = "plus1")
  inverse <- .mixture_invert(information)
  covariance <- if (identical(vcov_type, "robust")) {
    .inference_sandwich(inverse, crossprod(scores), nrow(scores))
  } else inverse
  profiles <- paste0("profile_", seq_len(k_profiles))
  initial_labels <- expand.grid(term = colnames(designs$initial),
                                to = profiles[-k_profiles], stringsAsFactors = FALSE)
  transition_labels <- do.call(rbind, lapply(seq_len(k_profiles), function(k) {
    expand.grid(term = colnames(stacked), to = profiles[-k], from = profiles[k],
                stringsAsFactors = FALSE)
  }))
  labels <- rbind(
    data.frame(block = "initial", from = NA_character_, to = initial_labels$to,
               term = initial_labels$term, stringsAsFactors = FALSE),
    data.frame(block = "transition", from = transition_labels$from,
               to = transition_labels$to, term = transition_labels$term,
               stringsAsFactors = FALSE))
  list(estimates = as.vector(natural_basis %*% fitted$par),
       covariance = natural_basis %*% covariance %*% t(natural_basis),
       log_likelihood = -fitted$value, labels = labels)
}

#' The three-step transition table of r3step()
#' @return A base data frame, one row per initial and transition coefficient.
#' @noRd
.lta_three_step <- function(x, data, initial_covariates, transition_covariates,
                            ci_level, vcov_type, adjust) {
  inputs <- .lta_three_step_inputs(x, data, initial_covariates, transition_covariates)
  if (identical(vcov_type, "robust")) {
    n_coefficients <- (inputs$n_profiles - 1L) *
      (ncol(inputs$designs$initial) + inputs$n_profiles *
         ncol(inputs$designs$transition[[1L]]))
    .multilpa_require_clusters(
      nrow(inputs$layout$slot), n_coefficients,
      what = "A cluster-robust covariance for the step-three transition model",
      alternative = paste("`vcov_type = \"observed\"` uses the model-based",
                          "information instead."))
  }
  fitted <- .lta_three_step_fit(inputs, vcov_type)
  errors <- sqrt(pmax(diag(fitted$covariance), 0))
  statistic <- ifelse(errors > 0, fitted$estimates / errors, NA_real_)
  p_value <- 2 * stats::pnorm(-abs(statistic))
  tested <- fitted$labels$term != "(Intercept)"
  adjusted <- rep(NA_real_, length(p_value))
  adjusted[tested] <- stats::p.adjust(p_value[tested], method = adjust)
  quantile <- .inference_critical(ci_level, "one_minus")
  result <- data.frame(fitted$labels, estimate = fitted$estimates,
                       standard_error = errors, statistic = statistic,
                       p_value = p_value, p_adjusted = adjusted,
                       conf_low = fitted$estimates - quantile * errors,
                       conf_high = fitted$estimates + quantile * errors,
                       row.names = NULL, stringsAsFactors = FALSE)
  structure(result, vcov_type = vcov_type, adjust = adjust,
            log_likelihood = fitted$log_likelihood,
            classification_errors = inputs$errors)
}
