# Maximum-likelihood step three for a distal outcome (Vermunt 2010; Bakk,
# Tekle and Vermunt 2013). The assigned class is an indicator of the true
# class with the classification errors of step one held fixed, and the
# outcome has a normal distribution in each true class:
#
#   L = sum_i log sum_t pi_t D[t, W_i] phi(z_i; mu_t, sigma_t),
#
# with pi_t the average posterior probability of class t and D[t, s] the
# probability that a unit of class t is assigned to class s. As in Latent
# GOLD's step-three model for a dependent variable, pi and D are fixed and
# only the class means and variances are estimated.

#' Fit the step-three model of a continuous distal outcome
#'
#' @param pieces The list from `.multilpa_level_assignments()` (or the
#'   trajectory equivalent): `posteriors`, `modal`, `n_classes`.
#' @param values The outcome, one value per unit.
#' @param variance `"varying"` (a variance per class) or `"equal"`.
#' @param tol Relative change in the log likelihood at which EM stops.
#' @param max_iter Maximum number of EM iterations.
#' @return A list: `means`, `variances` (one per class), `responsibilities`
#'   (units by classes), `log_likelihood`, `converged`, `iterations`, and the
#'   fixed `weights` (units by classes: the prior times the probability of the
#'   unit's assigned class).
#' @noRd
.multilpa_step_ml_fit <- function(pieces, values, variance, tol = 1e-12,
                                  max_iter = 5000L) {
  stopifnot("`variance` must be \"varying\" or \"equal\"" =
              variance %in% c("varying", "equal"))
  if (stats::var(values) <= 0) {
    stop(errorCondition(
      "The outcome takes one value, so it has no distribution to estimate by class.",
      class = "latents_bad_outcome", call = NULL))
  }
  errors <- .multilpa_error_matrix(pieces)
  prior <- colMeans(pieces$posteriors)
  # Each unit's fixed weight on each class: the class's share times the
  # probability that a unit of that class is assigned where this one was.
  weights <- t(errors[, pieces$modal, drop = FALSE]) *
    matrix(prior, length(values), pieces$n_classes, byrow = TRUE)
  log_weights <- log(weights)
  # The proportional (posterior-weighted) moments start EM, a non-random start
  # as in Latent GOLD's step three.
  totals <- colSums(pieces$posteriors)
  means <- as.vector(crossprod(pieces$posteriors, values)) / totals
  smallest <- 1e-10 * stats::var(values)
  moments <- function(responsibilities, means) {
    squares <- (values - matrix(means, length(values), length(means), byrow = TRUE))^2
    by_class <- colSums(responsibilities * squares)
    if (identical(variance, "equal")) {
      rep(sum(by_class) / length(values), length(means))
    } else by_class / colSums(responsibilities)
  }
  variances <- pmax(moments(pieces$posteriors, means), smallest)
  e_step <- function(means, variances) {
    log_density <- vapply(seq_along(means), function(class) {
      stats::dnorm(values, means[class], sqrt(variances[class]), log = TRUE)
    }, numeric(length(values)))
    joint <- log_weights + log_density
    largest <- apply(joint, 1L, max)
    log_marginal <- largest + log(rowSums(exp(joint - largest)))
    list(responsibilities = exp(joint - log_marginal),
         log_likelihood = sum(log_marginal))
  }
  state <- e_step(means, variances)
  converged <- FALSE
  iterations <- 0L
  # EM is iterative by nature: each step needs the previous one's estimates,
  # so this loop cannot be vectorised.
  while (iterations < max_iter) {
    iterations <- iterations + 1L
    responsibilities <- state$responsibilities
    means <- as.vector(crossprod(responsibilities, values)) / colSums(responsibilities)
    variances <- pmax(moments(responsibilities, means), smallest)
    updated <- e_step(means, variances)
    change <- abs(updated$log_likelihood - state$log_likelihood)
    state <- updated
    if (change <= tol * (1 + abs(state$log_likelihood))) {
      converged <- TRUE
      break
    }
  }
  fit <- list(means = unname(means), variances = unname(variances),
              responsibilities = state$responsibilities,
              log_likelihood = state$log_likelihood, converged = converged,
              iterations = iterations, weights = weights,
              at_floor = any(variances <= smallest * (1 + 1e-8)))
  if (fit$at_floor) return(fit)
  .multilpa_step_ml_polish(fit, values, variance)
}

#' Newton steps from the EM estimate to the stationary point
#'
#' EM approaches the maximum slowly near it, leaving a score that is small in
#' standard-error units but visible against another program's printed
#' estimates. A few Newton steps on the means and log standard deviations
#' finish the job; each step is halved until the likelihood does not fall, so
#' the polish never moves away from the EM estimate's likelihood.
#' @param fit From `.multilpa_step_ml_fit()`'s EM.
#' @param values The outcome.
#' @param variance `"varying"` or `"equal"`.
#' @return The fit at the stationary point, with `converged` set when the
#'   score vanished.
#' @noRd
.multilpa_step_ml_polish <- function(fit, values, variance) {
  objective <- function(parameters) {
    -.multilpa_step_ml_log_likelihood(fit, values, variance, parameters)
  }
  gradient <- function(parameters) {
    -colSums(.multilpa_step_ml_scores(fit, values, variance, parameters))
  }
  newton <- function(state, step) {
    if (isTRUE(state$done)) return(state)
    score <- gradient(state$parameters)
    if (max(abs(score)) < 1e-9 * (1 + abs(state$value))) {
      return(c(state[c("parameters", "value")], list(done = TRUE)))
    }
    hessian <- stats::optimHess(state$parameters, objective, gradient)
    direction <- tryCatch(solve(hessian, score), error = function(condition) NULL)
    if (is.null(direction)) return(c(state[c("parameters", "value")], list(done = FALSE)))
    # Halve the step until the likelihood does not fall.
    candidates <- lapply(2^-(0:20), \(size) state$parameters - size * direction)
    values_at <- vapply(candidates, objective, numeric(1))
    accepted <- which(values_at <= state$value)[1L]
    if (is.na(accepted)) return(c(state[c("parameters", "value")], list(done = FALSE)))
    list(parameters = candidates[[accepted]], value = values_at[accepted], done = NA)
  }
  start <- .multilpa_step_ml_coordinates(fit, variance)
  state <- Reduce(newton, seq_len(25L),
                  init = list(parameters = start, value = objective(start), done = NA))
  classes <- length(fit$means)
  sds <- exp(state$parameters[-seq_len(classes)])
  fit$means <- unname(state$parameters[seq_len(classes)])
  fit$variances <- unname(if (identical(variance, "equal")) rep(sds^2, classes) else sds^2)
  fit$log_likelihood <- -state$value
  # Converged when EM met its tolerance or Newton found the score zero; a
  # Newton step that cannot improve on a converged EM estimate (the score
  # left is rounding) does not undo EM's convergence.
  fit$converged <- isTRUE(fit$converged) || isTRUE(state$done)
  fit$responsibilities <- .multilpa_step_ml_responsibilities(fit, values)
  fit
}

#' Each unit's posterior class probabilities under the step-three model
#' @noRd
.multilpa_step_ml_responsibilities <- function(fit, values) {
  joint <- log(fit$weights) + vapply(seq_along(fit$means), function(class) {
    stats::dnorm(values, fit$means[class], sqrt(fit$variances[class]), log = TRUE)
  }, numeric(length(values)))
  largest <- apply(joint, 1L, max)
  exp(joint - largest - log(rowSums(exp(joint - largest))))
}

#' Per-unit scores of the step-three model
#'
#' Coordinates: the class means, then the log standard deviations (one per
#' class, or one shared).
#' @param fit From `.multilpa_step_ml_fit()`.
#' @param values The outcome.
#' @param variance `"varying"` or `"equal"`.
#' @param parameters Optional coordinates at which to evaluate; the fit's own
#'   by default.
#' @return A units-by-coordinates matrix.
#' @noRd
.multilpa_step_ml_scores <- function(fit, values, variance, parameters = NULL) {
  classes <- length(fit$means)
  parameters <- parameters %||% .multilpa_step_ml_coordinates(fit, variance)
  means <- parameters[seq_len(classes)]
  sds <- exp(parameters[-seq_len(classes)])
  sds <- if (identical(variance, "equal")) rep(sds, classes) else sds
  log_density <- vapply(seq_len(classes), function(class) {
    stats::dnorm(values, means[class], sds[class], log = TRUE)
  }, numeric(length(values)))
  joint <- log(fit$weights) + log_density
  largest <- apply(joint, 1L, max)
  responsibilities <- exp(joint - largest - log(rowSums(exp(joint - largest))))
  residuals <- values - matrix(means, length(values), classes, byrow = TRUE)
  scaled <- residuals / matrix(sds, length(values), classes, byrow = TRUE)
  mean_scores <- responsibilities * scaled / matrix(sds, length(values), classes,
                                                    byrow = TRUE)
  sd_scores <- responsibilities * (scaled^2 - 1)
  if (identical(variance, "equal")) sd_scores <- matrix(rowSums(sd_scores), ncol = 1L)
  unname(cbind(mean_scores, sd_scores))
}

#' The step-three model's coordinates at its estimate
#' @noRd
.multilpa_step_ml_coordinates <- function(fit, variance) {
  sds <- sqrt(fit$variances)
  c(fit$means, log(if (identical(variance, "equal")) sds[1L] else sds))
}

#' The step-three log likelihood at given coordinates
#' @noRd
.multilpa_step_ml_log_likelihood <- function(fit, values, variance, parameters) {
  classes <- length(fit$means)
  means <- parameters[seq_len(classes)]
  sds <- exp(parameters[-seq_len(classes)])
  sds <- if (identical(variance, "equal")) rep(sds, classes) else sds
  joint <- log(fit$weights) + vapply(seq_len(classes), function(class) {
    stats::dnorm(values, means[class], sds[class], log = TRUE)
  }, numeric(length(values)))
  largest <- apply(joint, 1L, max)
  sum(largest + log(rowSums(exp(joint - largest))))
}

#' Influence contributions of the class means, one row per unit
#'
#' The sandwich A^-1 B A^-1 is the sum of the outer products of these rows,
#' where A is the observed information of the step-three likelihood and B
#' the cross-product of the unit scores. Rows for the means only.
#' @param fit From `.multilpa_step_ml_fit()`.
#' @param values The outcome.
#' @param variance `"varying"` or `"equal"`.
#' @return A units-by-classes matrix.
#' @noRd
.multilpa_step_ml_influence <- function(fit, values, variance) {
  coordinates <- .multilpa_step_ml_coordinates(fit, variance)
  information <- stats::optimHess(
    coordinates,
    function(parameters) -.multilpa_step_ml_log_likelihood(fit, values, variance,
                                                          parameters),
    function(parameters) -colSums(.multilpa_step_ml_scores(fit, values, variance,
                                                          parameters)))
  inverse <- tryCatch(solve(information), error = function(condition) {
    stop(errorCondition(
      "The information of the step-three model is singular; the class means are not identified.",
      class = "latents_singular_information", call = NULL))
  })
  scores <- .multilpa_step_ml_scores(fit, values, variance)
  (scores %*% inverse)[, seq_along(fit$means), drop = FALSE]
}

#' The ML three-step table of a distal outcome
#'
#' @param pieces,values The units' assignments and outcome values.
#' @param level,variance,ci_level,contrast,adjust,vcov_type As in
#'   [three_step()].
#' @return The table [three_step()] returns.
#' @noRd
.multilpa_three_step_ml <- function(pieces, values, level, variance, ci_level,
                                    contrast, adjust, vcov_type) {
  fit <- .multilpa_step_ml_fit(pieces, values, variance)
  if (!fit$converged) {
    warning(warningCondition(sprintf(paste(
      "The step-three model did not converge in %d iterations, so its class",
      "means are not a maximum."), fit$iterations),
      class = "latents_unconverged", call = NULL))
  }
  if (fit$at_floor) {
    warning(warningCondition(paste(
      "A class variance of the outcome reached its lower bound: the class has",
      "almost no spread in the outcome, and its standard error is unreliable."),
      class = "latents_boundary", call = NULL))
  }
  quantile <- stats::qnorm(1 - (1 - ci_level) / 2)
  clustered <- .multilpa_step_clusters(.multilpa_step_ml_influence(fit, values, variance),
                                       pieces, vcov_type)
  if (identical(contrast, "pairs")) {
    return(.multilpa_step_pairs(fit$means, clustered, level, "ml", quantile,
                                adjust, vcov_type))
  }
  errors <- sqrt(colSums(clustered^2))
  result <- data.frame(level = level, method = "ml",
                       class = seq_len(pieces$n_classes), estimate = fit$means,
                       standard_error = errors,
                       conf_low = fit$means - quantile * errors,
                       conf_high = fit$means + quantile * errors,
                       variance = fit$variances,
                       effective_n = colSums(fit$responsibilities),
                       row.names = NULL, stringsAsFactors = FALSE)
  attr(result, "vcov_type") <- vcov_type
  attr(result, "variance") <- variance
  column <- .latents_column
  .latents_table(result, ci_level,
                 sprintf("Distal outcome by class (ml, %s variances, %s%% CI)",
                         variance, format(100 * ci_level)),
                 list(column("Class", "class", "integer"), column("Mean", "estimate"),
                      column(sprintf("%s%% CI", format(100 * ci_level)),
                             c("conf_low", "conf_high"), "ci", "estimate"),
                      column("SE", "standard_error"), column("Variance", "variance"),
                      column("Effective n", "effective_n")))
}
