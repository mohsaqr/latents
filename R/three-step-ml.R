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
  # Fit and differentiate in a chart whose origin and units are determined
  # by the outcome spread. A finite-difference step around a raw mean near
  # 1e12 can round to zero even though the outcome varies by whole units.
  outcome_center <- mean(values)
  outcome_scale <- stats::sd(values)
  values <- (values - outcome_center) / outcome_scale
  errors <- .multilpa_error_matrix(pieces)
  prior <- colMeans(pieces$posteriors)
  # Each unit's fixed weight on each class: the class's share times the
  # probability that a unit of that class is assigned where this one was.
  weights <- t(errors[, pieces$modal, drop = FALSE]) *
    matrix(prior, length(values), pieces$n_classes, byrow = TRUE)
  log_weights <- log(weights)
  smallest <- 1e-10 * stats::var(values)
  moments <- function(responsibilities, means) {
    squares <- (values - matrix(means, length(values), length(means), byrow = TRUE))^2
    by_class <- colSums(responsibilities * squares)
    if (identical(variance, "equal")) {
      rep(sum(by_class) / length(values), length(means))
    } else by_class / colSums(responsibilities)
  }
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
  run <- function(means, variances) {
    state <- e_step(means, variances)
    converged <- FALSE
    failed <- FALSE
    iterations <- 0L
    # EM is iterative by nature: each step needs the previous one's estimates,
    # so this loop cannot be vectorised.
    while (iterations < max_iter) {
      iterations <- iterations + 1L
      responsibilities <- state$responsibilities
      # A start that leaves a class with no weight has no estimate for it;
      # the start is dropped rather than carried on with an undefined mean.
      if (!all(colSums(responsibilities) > 0)) {
        failed <- TRUE
        break
      }
      means <- as.vector(crossprod(responsibilities, values)) / colSums(responsibilities)
      variances <- pmax(moments(responsibilities, means), smallest)
      updated <- e_step(means, variances)
      if (!is.finite(updated$log_likelihood)) {
        failed <- TRUE
        break
      }
      change <- abs(updated$log_likelihood - state$log_likelihood)
      state <- updated
      if (change <= tol * (1 + abs(state$log_likelihood))) {
        converged <- TRUE
        break
      }
    }
    list(means = unname(means), variances = unname(variances),
         responsibilities = state$responsibilities,
         log_likelihood = state$log_likelihood, converged = converged,
         iterations = iterations, weights = weights, failed = failed,
         at_floor = !failed && any(variances <= smallest * (1 + 1e-8)))
  }
  # The step-three likelihood can have several maxima (a heavy-tailed or
  # skewed outcome lets one class absorb the tail), so EM runs from a fixed,
  # reproducible set of starts: the proportional (posterior-weighted) moments,
  # Latent GOLD's start, first; then the modal-class means; then the class
  # means at the outcome's quantiles in every order (every cyclic order, both
  # ways, above five classes), from the outcome's own variance and from the
  # proportional spread; then the seeded draws below. The equal-variance
  # model also tries quasi-Newton optimization from the quantile placements.
  proportional <- as.vector(crossprod(pieces$posteriors, values)) / colSums(pieces$posteriors)
  spread <- pmax(moments(pieces$posteriors, proportional), smallest)
  modal <- vapply(seq_len(pieces$n_classes), function(class) {
    members <- pieces$modal == class
    if (any(members)) mean(values[members]) else proportional[class]
  }, numeric(1))
  k <- pieces$n_classes
  quantiles <- stats::quantile(values, (seq_len(k) - 0.5) / k, names = FALSE)
  orders <- if (k <= 5L) .multilpa_permutations(k) else {
    shifts <- t(vapply(seq_len(k) - 1L, function(shift) (seq_len(k) + shift - 1L) %% k + 1L,
                       integer(k)))
    rbind(shifts, shifts[, rev(seq_len(k)), drop = FALSE])
  }
  placed <- lapply(seq_len(nrow(orders)), function(row) quantiles[orders[row, ]])
  # A heavy tail can be the likelihood's best class: each class in turn starts
  # at the mean of the outcome's top and bottom 2%, the others proportional.
  tails <- c(mean(values[values >= stats::quantile(values, 0.98)]),
             mean(values[values <= stats::quantile(values, 0.02)]))
  tailed <- do.call(c, lapply(tails, function(tail) lapply(seq_len(k), function(class) {
    replace(proportional, class, tail)
  })))
  # Twenty more starts at random quantiles with random standard deviations
  # (0.3 to 1 of the outcome's) reach basins the ordered placements miss.
  # They are drawn from a fixed seed in a local stream, so the result is the
  # same whatever the caller's seed and the caller's .Random.seed is restored
  # (with the Box-Muller cache limitation documented by the seed service).
  drawn <- (function() {
    .latents_local_seed(20261021L)
    lapply(seq_len(20L), function(draw) {
      sds <- if (identical(variance, "equal")) rep(stats::runif(1L, 0.3, 1), k) else
        stats::runif(k, 0.3, 1)
      list(means = stats::quantile(values, stats::runif(k), names = FALSE),
           variances = sds^2)
    })
  })()
  # A start keeps only its estimates, not its observations-by-classes
  # responsibilities: with five classes there are 272 starts. The winner's
  # are recomputed below from its own estimates, by the same E-step.
  slim <- function(means, variances) {
    attempt <- run(means, variances)
    attempt$responsibilities <- NULL
    attempt
  }
  runs <- c(lapply(c(list(proportional, modal), tailed), slim, variances = spread),
            lapply(placed, slim, variances = rep(1, k)),
            lapply(placed, slim, variances = spread),
            lapply(drawn, function(start) slim(start$means, start$variances)))
  if (identical(variance, "equal")) {
    # EM and quasi-Newton iterations can enter different basins from the
    # same mean placement. For the shared-variance model, also optimize the
    # existing quantile starts directly. The chart variance = floor + exp(2a)
    # enforces the same lower bound as EM without changing the likelihood.
    template <- list(means = numeric(k), variances = rep(1, k), weights = weights)
    log_floor <- log(smallest)
    decode <- function(point) {
      log_excess <- 2 * point[k + 1L]
      log_variance <- max(log_floor, log_excess) +
        log1p(exp(-abs(log_floor - log_excess)))
      c(point[seq_len(k)], log_variance / 2)
    }
    objective <- function(point) {
      -.multilpa_step_ml_log_likelihood(template, values, variance, decode(point))
    }
    gradient <- function(point) {
      decoded <- decode(point)
      score <- -colSums(.multilpa_step_ml_scores(template, values, variance, decoded))
      score[k + 1L] <- score[k + 1L] * exp(2 * point[k + 1L] - 2 * decoded[k + 1L])
      score
    }
    optimized <- lapply(placed, function(means) tryCatch({
      start <- c(means, log(1 - smallest) / 2)
      found <- stats::optim(start, objective, gradient, method = "BFGS",
                            control = list(maxit = max_iter, reltol = 1e-15))
      if (!identical(found$convergence, 0L) || !is.finite(found$value) ||
          any(!is.finite(found$par))) return(NULL)
      point <- decode(found$par)
      candidate <- template
      candidate$means <- point[seq_len(k)]
      candidate$variances <- rep(exp(2 * point[k + 1L]), k)
      if (any(candidate$variances <= smallest * (1 + 1e-8))) return(NULL)
      scale <- .multilpa_step_ml_scale(candidate, variance)
      score <- colSums(.multilpa_step_ml_scores(candidate, values, variance, point))
      if (any(!is.finite(score)) || max(abs(score * scale)) >
          sqrt(.Machine$double.eps) * (1 + length(values))) return(NULL)
      information <- .multilpa_step_ml_hessian(function(theta) {
        -colSums(.multilpa_step_ml_scores(candidate, values, variance, theta))
      }, point, scale)
      chol(information * tcrossprod(scale))
      slim(candidate$means, candidate$variances)
    }, error = function(condition) NULL))
    # Appending preserves the original EM winner whenever the likelihoods
    # are tied under the existing selection rule below.
    runs <- c(runs, Filter(Negate(is.null), optimized))
  }
  failed <- vapply(runs, `[[`, logical(1), "failed")
  if (all(failed)) {
    stop(errorCondition(paste(
      "Every start of the step-three model left a class without weight; the",
      "class means are not identified."),
      class = "latents_singular_information", call = NULL))
  }
  scores <- vapply(runs, function(attempt) {
    if (attempt$failed || attempt$at_floor) -Inf else attempt$log_likelihood
  }, numeric(1))
  # A later start replaces the first only when it is clearly higher, so a
  # likelihood with one maximum keeps Latent GOLD's solution exactly. When
  # every start ends at the variance floor, the first that ran is kept.
  best <- if (all(is.infinite(scores))) which(!failed)[1L] else {
    top <- max(scores)
    which(scores >= top - 1e-10 * (1 + abs(top)))[1L]
  }
  fit <- runs[[best]]
  fit$responsibilities <- e_step(fit$means, fit$variances)$responsibilities
  if (!fit$at_floor) fit <- .multilpa_step_ml_polish(fit, values, variance)
  # Keep the conditioned means for inference: subtracting the origin from
  # already rounded reported means would throw away the precision just saved.
  fit$outcome_coordinates <- list(center = outcome_center, scale = outcome_scale,
                                  means = fit$means)
  fit$means <- outcome_center + outcome_scale * fit$means
  fit$variances <- outcome_scale^2 * fit$variances
  fit$log_likelihood <- fit$log_likelihood - length(values) * log(outcome_scale)
  fit
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
  scale <- .multilpa_step_ml_scale(fit, variance)
  newton <- function(state, step) {
    if (isTRUE(state$done)) return(state)
    score <- gradient(state$parameters)
    if (max(abs(score * scale)) < 1e-9 * (1 + abs(state$value))) {
      return(c(state[c("parameters", "value")], list(done = TRUE)))
    }
    hessian <- .multilpa_step_ml_hessian(gradient, state$parameters, scale)
    # Solved on the coordinates' own scale, like the sandwich's inverse.
    direction <- tryCatch(scale * solve(hessian * tcrossprod(scale), score * scale),
                          error = function(condition) NULL)
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

#' The scale of each step-three coordinate
#'
#' The means are in the outcome's units and the log standard deviations are
#' unit free. The finite differences of the Hessian and the stopping rule on
#' the score are taken on this scale, so neither depends on the outcome's
#' units: an absolute step of 1e-3 is far wider than the curvature of an
#' outcome measured in units of 1e-4.
#' @param fit A step-three fit with `means` and `variances`.
#' @param variance `"varying"` or `"equal"`.
#' @return One positive scale per coordinate of `.multilpa_step_ml_coordinates()`.
#' @noRd
.multilpa_step_ml_scale <- function(fit, variance) {
  spread <- sqrt(mean(fit$variances))
  c(rep(spread, length(fit$means)),
    rep(1, if (identical(variance, "equal")) 1L else length(fit$means)))
}

#' The Hessian of the step-three objective from its analytic gradient
#'
#' Central differences of the gradient, each coordinate stepped by a fixed
#' fraction of its own scale, then symmetrized. `stats::optimHess()` steps by
#' 1e-3 in the coordinates' units, which on an outcome measured in units of
#' 1e-4 is wider than the curvature it measures (an information of 0.15 for
#' 365, and standard errors 3000 times too large); with `parscale` it was
#' still wrong.
#' @param gradient A function of the coordinates.
#' @param parameters The point.
#' @param scale From `.multilpa_step_ml_scale()`.
#' @return A symmetric matrix.
#' @noRd
.multilpa_step_ml_hessian <- function(gradient, parameters, scale) {
  steps <- 1e-5 * scale
  columns <- vapply(seq_along(parameters), function(j) {
    up <- parameters
    up[j] <- up[j] + steps[j]
    down <- parameters
    down[j] <- down[j] - steps[j]
    (gradient(up) - gradient(down)) / (2 * steps[j])
  }, numeric(length(parameters)))
  hessian <- matrix(columns, length(parameters))
  (hessian + t(hessian)) / 2
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
  chart <- fit$outcome_coordinates
  outcome_scale <- if (is.null(chart)) 1 else chart$scale
  if (!is.null(chart)) {
    values <- (values - chart$center) / outcome_scale
    fit$means <- chart$means
    fit$variances <- fit$variances / outcome_scale^2
  }
  coordinates <- .multilpa_step_ml_coordinates(fit, variance)
  scale <- .multilpa_step_ml_scale(fit, variance)
  information <- .multilpa_step_ml_hessian(
    function(parameters) -colSums(.multilpa_step_ml_scores(fit, values, variance,
                                                          parameters)),
    coordinates, scale)
  # Inverted on the coordinates' own scale, D (D H D)^-1 D: in raw units the
  # means' and log standard deviations' information differ by the outcome's
  # units squared, and at 1e8 that alone made an identified model singular.
  inverse <- tryCatch({
    conditioned <- information * tcrossprod(scale)
    # A genuine saddle or unidentified class must not produce uncertainty
    # from an arbitrary invertible but nonpositive information matrix.
    chol(conditioned)
    solve(conditioned) * tcrossprod(scale)
  }, error = function(condition) {
    stop(errorCondition(
      "The information of the step-three model is singular; the class means are not identified.",
      class = "latents_singular_information", call = NULL))
  })
  scores <- .multilpa_step_ml_scores(fit, values, variance)
  outcome_scale * (scores %*% inverse)[, seq_along(fit$means), drop = FALSE]
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
  # A sandwich for a class carried by a single unit is built from that unit's
  # score, which is zero at its own mean: it would report a standard error of
  # 0 rather than none. The units carrying a class are counted by Kish's
  # effective number, (sum r)^2 / sum r^2: 1 for one unit, 2 for two equal.
  # The units are those the sandwich sums over: groups under "cluster".
  responsibilities <- if (identical(vcov_type, "cluster")) {
    rowsum(fit$responsibilities, pieces$group_index, reorder = FALSE)
  } else fit$responsibilities
  carriers <- colSums(responsibilities)^2 / colSums(responsibilities^2)
  thin <- carriers < 1.5
  if (any(thin)) {
    warning(warningCondition(sprintf(paste(
      "Class %s of the step-three model is carried by a single unit, so its",
      "standard error cannot be estimated and is reported as NA."),
      paste(which(thin), collapse = ", ")), class = "latents_small_classes", call = NULL))
    clustered[, thin] <- NA_real_
  }
  if (identical(contrast, "pairs")) {
    # The common origin cancels from contrasts; remove it before subtraction
    # so reported rounding at a large origin cannot erase a small difference.
    means <- fit$outcome_coordinates$means * fit$outcome_coordinates$scale
    return(.multilpa_step_pairs(means, clustered, level, "ml", quantile,
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
