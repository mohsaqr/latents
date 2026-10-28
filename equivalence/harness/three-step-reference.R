# Step-three references for the three_step and r3step verbs, written from
# their definitions in base R and checked against two external programs.
#
# A three-step analysis takes the step-one posteriors as given; the verbs under
# test are the step-two classification (modal assignment, the error matrix D)
# and the step-three estimators that use it. Everything here starts from the
# posterior matrix and recomputes the rest without calling the package:
#
#   - D[t, s] = sum_i p_it 1(W_i = s) / sum_i p_it, rows the true class
#     (Vermunt 2010, eq. 6), with W_i the modal class, ties to the first;
#   - weighted class means and their cluster-robust sandwich, written as an
#     estimating equation (score w_ic (z_i - m_c), bread sum_i w_ic, CR1 factor
#     G / (G - 1));
#   - the ML step-three likelihood sum_i log sum_t pi_t D[t, W_i] phi(z_i;
#     mu_t, sigma_t) (Bakk, Tekle & Vermunt 2013);
#   - the R3STEP likelihood sum_i log sum_t pi_t(x_i) D[t, W_i], pi a
#     multinomial logit against the last class (Vermunt 2010);
#   - the same with group-class intercepts, sum_g log sum_h eta_h prod_{i in g}
#     sum_k pi_hk(x_i) D[k, W_i];
#   - the three-step latent Markov chain, a forward recursion over the
#     assigned profiles with emissions D, initial logits against the last
#     profile and transition logits against staying (Di Mari, Oberski &
#     Vermunt 2016).
#
# Maxima are found by depmixS4 (mixture and hidden Markov models with D held
# fixed) where it can express the model, and then carried to the fixed point
# by Newton steps on the likelihood written here, after the two likelihoods
# have been shown to agree at depmixS4's estimates. The ML, R3STEP and chain
# likelihoods carry their score in closed form (the chain's by a scaled
# forward-backward pass), checked against numDeriv's gradient; a Newton step
# is accepted when the score's norm falls, so the iteration reaches a score
# of about 1e-14 instead of stalling on the rounding of the summed
# likelihood. Their Hessians are closed form (ML, R3STEP) or numDeriv's
# Jacobian of the closed-form score (chain); the group-class R3STEP, which no
# program fits, uses numDeriv throughout. Everything is done in a chart whose
# covariates and outcome are centred, so that a covariate held far from zero
# does not make the intercept's step meaningless, and mapped back.

#' Modal class of each unit, ties to the first class
#' @param posteriors Units-by-classes posterior matrix.
#' @return An integer vector.
three_step_modal <- function(posteriors) {
  max.col(posteriors, ties.method = "first")
}

#' Classification error matrix, rows the true class
#' @param posteriors Units-by-classes posterior matrix.
#' @param modal Modal class of each unit.
#' @return A square matrix whose rows sum to one.
three_step_error_matrix <- function(posteriors, modal) {
  classes <- seq_len(ncol(posteriors))
  assigned <- vapply(classes, function(s) as.numeric(modal == s), numeric(length(modal)))
  crossprod(posteriors, assigned) / colSums(posteriors)
}

#' Class means of a distal outcome under fixed unit weights, with their
#' cluster-robust covariance written as an estimating equation
#'
#' The mean of class c solves sum_i w_ic (z_i - m_c) = 0. The sandwich is
#' A^-1 B A^-1 with A = diag(sum_i w_ic) and B the cross-product of the
#' cluster sums of the unit scores, times G / (G - 1) for G clusters.
#' @param z The outcome, one value per unit.
#' @param weights Units-by-classes weights (negative allowed).
#' @param cluster Cluster of each unit.
#' @return A list with `estimates` and `covariance`.
three_step_estimating_sandwich <- function(z, weights, cluster) {
  totals <- colSums(weights)
  estimates <- as.vector(crossprod(weights, z)) / totals
  scores <- weights * outer(z, estimates, "-")
  clustered <- rowsum(scores, cluster, reorder = FALSE)
  n_clusters <- nrow(clustered)
  bread <- diag(1 / totals, length(totals))
  list(estimates = estimates,
       covariance = bread %*% crossprod(clustered) %*% bread *
         n_clusters / (n_clusters - 1))
}

#' The same class means by weighted least squares, with sandwich::vcovCL
#'
#' The outcome is stacked once per class with the unit's weight on that class
#' and regressed on the class indicators without an intercept; the
#' coefficients are the weighted class means. `vcovCL(type = "HC0", cadjust =
#' TRUE)` is the CR1 sandwich clustered on the unit's cluster. lm() refuses
#' negative weights, so this reference serves the proportional and modal
#' weights only.
#'
#' Rows of weight zero are dropped before the fit. They add nothing to the
#' estimating equations, but lm() keeps them in its residuals (so in
#' sandwich's scores and its n) and not in its residual degrees of freedom
#' (so not in sandwich's bread): with the modal weights, two thirds of the
#' stacked rows of a three-class model, vcovCL() then mis-scales the variance
#' by (n_positive / n)^2.
#' @inheritParams three_step_estimating_sandwich
#' @return A list with `estimates` and `covariance`.
three_step_lm_sandwich <- function(z, weights, cluster) {
  classes <- ncol(weights)
  stacked <- data.frame(z = rep(z, classes),
                        class = factor(rep(seq_len(classes), each = length(z))),
                        weight = as.vector(weights),
                        cluster = rep(cluster, classes))
  stacked <- stacked[stacked$weight != 0, , drop = FALSE]
  model <- stats::lm(z ~ 0 + class, data = stacked, weights = weight)
  list(estimates = unname(stats::coef(model)),
       covariance = unname(sandwich::vcovCL(model, cluster = ~cluster, type = "HC0",
                                            cadjust = TRUE)))
}

#' The table three_step() reports for class means
#' @param level,method Carried into the table.
#' @param estimates,covariance Class means and their covariance.
#' @param effective_n One value per class.
#' @param variance `NULL`, or the ML class variances.
#' @param ci_level Confidence level.
#' @return A `data.frame` shaped like three_step(contrast = "none").
three_step_mean_table <- function(level, method, estimates, covariance, effective_n,
                                  variance = NULL, ci_level = 0.95) {
  errors <- sqrt(diag(covariance))
  quantile <- stats::qnorm(1 - (1 - ci_level) / 2)
  table <- data.frame(level = level, method = method, class = seq_along(estimates),
                      estimate = estimates, standard_error = errors,
                      conf_low = estimates - quantile * errors,
                      conf_high = estimates + quantile * errors,
                      stringsAsFactors = FALSE)
  if (!is.null(variance)) table$variance <- variance
  table$effective_n <- effective_n
  table
}

#' The table three_step(contrast = "pairs") reports
#'
#' Every unordered pair (a, b), a < b: the difference m_b - m_a, its standard
#' error sqrt(V_aa + V_bb - 2 V_ab), the Wald test against zero and the
#' multiplicity correction across the pairs (stats::p.adjust).
#' @inheritParams three_step_mean_table
#' @param adjust A p.adjust() method.
#' @return A `data.frame` shaped like three_step(contrast = "pairs").
three_step_pairs_table <- function(level, method, estimates, covariance, adjust = "BH",
                                   ci_level = 0.95) {
  pairs <- utils::combn(length(estimates), 2L)
  first <- pairs[1L, ]
  second <- pairs[2L, ]
  difference <- estimates[second] - estimates[first]
  errors <- sqrt(covariance[cbind(second, second)] + covariance[cbind(first, first)] -
                   2 * covariance[cbind(first, second)])
  statistic <- difference / errors
  p_value <- 2 * stats::pnorm(-abs(statistic))
  quantile <- stats::qnorm(1 - (1 - ci_level) / 2)
  data.frame(level = level, method = method, class = second, reference_class = first,
             estimate = difference, standard_error = errors, statistic = statistic,
             p_value = p_value, p_adjusted = stats::p.adjust(p_value, method = adjust),
             conf_low = difference - quantile * errors,
             conf_high = difference + quantile * errors, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Newton steps on a likelihood written here
# ---------------------------------------------------------------------------

#' A likelihood with numDeriv derivatives
#'
#' The shape every reference below hands to [three_step_newton()] and
#' [three_step_covariances()]: the total, its score and Hessian, and the
#' per-unit scores. Here all derivatives are numDeriv's (Richardson
#' extrapolation); the ML step-three and R3STEP models supply them in closed
#' form instead ([three_step_ml_model()], [three_step_r3_model()]).
#' @param unit_log_likelihood A function of the parameters returning one log
#'   likelihood per unit.
#' @return A list of functions `units`, `total`, `score`, `unit_scores` and
#'   `hessian`, and `kind = "numerical"`.
three_step_numerical_model <- function(unit_log_likelihood) {
  total <- function(point) sum(unit_log_likelihood(point))
  list(units = unit_log_likelihood, total = total,
       score = function(point) numDeriv::grad(total, point),
       unit_scores = function(point) numDeriv::jacobian(unit_log_likelihood, point),
       hessian = function(point) numDeriv::hessian(total, point),
       kind = "numerical")
}

#' How far a closed-form score and Hessian are from numDeriv's
#'
#' A closed-form derivative is trusted only because it reproduces numDeriv's:
#' the score against the Richardson gradient of the likelihood, the Hessian
#' against the Richardson Jacobian of that score (a first derivative of a
#' smooth function, good to about 1e-10, where numDeriv's second derivative
#' of the likelihood is good to about 1e-6). Each is compared at one point,
#' as a fraction of its largest entry or of one, whichever is larger.
#' @param model A model with a closed-form score.
#' @param theta The point.
#' @return A named vector: `score` and `hessian` relative differences.
three_step_derivative_check <- function(model, theta) {
  relative <- function(a, b) max(abs(a - b)) / max(abs(b), 1)
  c(score = relative(model$score(theta), numDeriv::grad(model$total, theta)),
    hessian = relative(model$hessian(theta), numDeriv::jacobian(model$score, theta)))
}

#' Carry an estimate to the stationary point of a log likelihood
#'
#' Newton steps, each halved until the norm of the score falls. Near the
#' maximum the likelihood gains less per step (1e-13 on these models) than
#' the rounding of a sum over a thousand units, so a rule that accepts steps
#' on the likelihood stalls short of the fixed point with a score of 1e-6
#' (the review of 2026-10-27); the score keeps falling quadratically. Started
#' from another program's maximum, a handful of steps reach the fixed point.
#' The iteration stops when no halving lowers the score or a step moves the
#' estimate by less than 1e-13 relative. The result reports the score and
#' whether the Hessian is negative definite there, so the caller can refuse a
#' saddle.
#' @param model A model, as from [three_step_numerical_model()].
#' @param theta The starting point.
#' @param steps Maximum number of Newton steps.
#' @return A list with `theta`, `log_likelihood`, `gradient` (the score at
#'   `theta`), `concave`, `kind` (the model's: `"closed form"` or
#'   `"numerical"` derivatives), `steps` (taken) and `moved` (the largest
#'   absolute change from the start).
three_step_newton <- function(model, theta, steps = 30L) {
  start <- theta
  norm <- function(point) sqrt(sum(model$score(point)^2))
  step <- function(state, iteration) {
    if (isTRUE(state$done)) return(state)
    direction <- solve(-model$hessian(state$theta), model$score(state$theta))
    # Find() tries the sizes in order and stops at the first that lowers the
    # score, so a full Newton step costs one score.
    size <- Find(function(size) norm(state$theta + size * direction) < state$norm,
                 2^-(0:30))
    if (is.null(size)) return(modifyList(state, list(done = TRUE)))
    moved <- size * max(abs(direction))
    candidate <- state$theta + size * direction
    list(theta = candidate, norm = norm(candidate), taken = state$taken + 1L,
         done = moved <= 1e-13 * (1 + max(abs(state$theta))))
  }
  # Newton's iteration needs the previous point; Reduce() carries it.
  state <- Reduce(step, seq_len(steps),
                  init = list(theta = theta, norm = norm(theta), taken = 0L, done = FALSE))
  list(theta = state$theta, log_likelihood = model$total(state$theta),
       gradient = model$score(state$theta),
       concave = all(eigen(model$hessian(state$theta), symmetric = TRUE,
                           only.values = TRUE)$values < 0),
       kind = model$kind, steps = state$taken, moved = max(abs(state$theta - start)))
}

#' Observed and cluster-robust covariances from a likelihood written here
#'
#' The observed information is minus the Hessian of the total; the robust
#' covariance is the sandwich with the cross-product of the cluster-summed
#' unit scores, times G / (G - 1). The derivatives are the model's: closed
#' form for ML and R3STEP, numDeriv otherwise.
#' @param model A model, as from [three_step_numerical_model()].
#' @param theta The estimate.
#' @param cluster Cluster of each unit.
#' @return A list with `observed` and `robust` covariance matrices.
three_step_covariances <- function(model, theta, cluster) {
  bread <- solve(-model$hessian(theta))
  scores <- rowsum(model$unit_scores(theta), cluster, reorder = FALSE)
  n_clusters <- nrow(scores)
  list(observed = bread,
       robust = bread %*% crossprod(scores) %*% bread * n_clusters / (n_clusters - 1))
}

# ---------------------------------------------------------------------------
# ML step three for a distal outcome
# ---------------------------------------------------------------------------

#' Per-unit log likelihood of the ML step-three model
#' @param theta Class means, then log standard deviations (one per class, or
#'   one shared under `variance = "equal"`).
#' @param z The outcome.
#' @param log_weights Units-by-classes log(pi_t D[t, W_i]).
#' @return One log likelihood per unit.
three_step_ml_units <- function(theta, z, log_weights) {
  classes <- ncol(log_weights)
  means <- theta[seq_len(classes)]
  sds <- rep_len(exp(theta[-seq_len(classes)]), classes)
  joint <- log_weights + vapply(seq_len(classes), function(class) {
    stats::dnorm(z, means[[class]], sds[[class]], log = TRUE)
  }, numeric(length(z)))
  peak <- apply(joint, 1L, max)
  peak + log(rowSums(exp(joint - peak)))
}

#' The ML step-three likelihood with its score and Hessian in closed form
#'
#' Unit i's log likelihood is l_i = log sum_t exp(c_it), c_it = log(pi_t
#' D[t, W_i]) + log phi(z_i; mu_t, s_t), s_t = exp(eta_t). With r_it the
#' responsibilities, u = z_i - mu_t and v = u^2 / s_t^2, the component
#' derivatives are dc/dmu_t = u / s_t^2 and dc/deta_t = v - 1, and the second
#' derivatives -1 / s_t^2, -2 u / s_t^2 and -2 v. The unit score is sum_t r_it
#' dc_it, and the Hessian of a log-sum-exp is sum_i [sum_t r_it (d2c_it +
#' dc_it dc_it') - g_i g_i']: one 2 x 2 block per class less the cross-product
#' of the unit scores. Under `variance = "equal"` one log sd is shared, a
#' linear map of the per-class one.
#' @param z The outcome.
#' @param log_weights Units-by-classes log(pi_t D[t, W_i]).
#' @param variance `"varying"` or `"equal"`.
#' @return A model, shaped like [three_step_numerical_model()]'s, with `kind =
#'   "closed form"`.
three_step_ml_model <- function(z, log_weights, variance) {
  classes <- ncol(log_weights)
  # Parameters (means, log sds) -> the per-class chart (means, K log sds).
  chart <- if (identical(variance, "equal")) {
    rbind(cbind(diag(classes), 0), cbind(matrix(0, classes, classes), 1))
  } else diag(2L * classes)
  pieces <- function(theta) {
    means <- theta[seq_len(classes)]
    sds <- rep_len(exp(theta[-seq_len(classes)]), classes)
    units <- three_step_ml_units(theta, z, log_weights)
    joint <- log_weights + vapply(seq_len(classes), function(class) {
      stats::dnorm(z, means[[class]], sds[[class]], log = TRUE)
    }, numeric(length(z)))
    spread <- matrix(sds^2, length(z), classes, byrow = TRUE)
    deviation <- outer(z, means, "-")
    list(units = units, responsibilities = exp(joint - units), spread = spread,
         deviation = deviation, scaled = deviation^2 / spread)
  }
  unit_scores <- function(theta) {
    with(pieces(theta), cbind(responsibilities * deviation / spread,
                              responsibilities * (scaled - 1))) %*% chart
  }
  hessian <- function(theta) {
    parts <- pieces(theta)
    scores <- with(parts, cbind(responsibilities * deviation / spread,
                                responsibilities * (scaled - 1)))
    within_class <- with(parts, list(
      mean_mean = colSums(responsibilities * (deviation^2 / spread - 1) / spread),
      mean_sd = colSums(responsibilities * deviation / spread * (scaled - 3)),
      sd_sd = colSums(responsibilities * (scaled^2 - 4 * scaled + 1))))
    index <- seq_len(classes)
    full <- -crossprod(scores)
    full[cbind(index, index)] <- full[cbind(index, index)] + within_class$mean_mean
    full[cbind(index, classes + index)] <- full[cbind(index, classes + index)] +
      within_class$mean_sd
    full[cbind(classes + index, index)] <- full[cbind(classes + index, index)] +
      within_class$mean_sd
    full[cbind(classes + index, classes + index)] <-
      full[cbind(classes + index, classes + index)] + within_class$sd_sd
    crossprod(chart, full %*% chart)
  }
  units <- function(theta) three_step_ml_units(theta, z, log_weights)
  list(units = units, total = function(theta) sum(units(theta)),
       score = function(theta) colSums(unit_scores(theta)),
       unit_scores = unit_scores, hessian = hessian, kind = "closed form",
       responsibilities = function(theta) pieces(theta)$responsibilities)
}

#' The ML step-three reference
#'
#' depmixS4 fits the mixture `list(W ~ 1, z ~ 1)` (multinomial assigned
#' class, Gaussian outcome) with the class sizes fixed at the average
#' posteriors and the response probabilities of W fixed at D, from the
#' proportional and the modal moments, keeping the higher maximum; under
#' `variance = "equal"` the class standard deviations are constrained equal
#' through `fit(equal = )`. Its estimate is then carried to the fixed point by
#' Newton steps on the likelihood written here ([three_step_ml_model()],
#' closed-form score and Hessian), in a chart centred at the outcome mean.
#' The standard errors are the CR1 sandwich of the same closed-form unit
#' scores and Hessian.
#' @param z The outcome.
#' @param posteriors Units-by-classes step-one posteriors.
#' @param cluster Cluster of each unit.
#' @param variance `"varying"` or `"equal"`.
#' @return A list with `means`, `variances`, `covariance` (robust, of the
#'   means), `effective_n` (column sums of the step-three responsibilities),
#'   `log_likelihood`, `depmix_log_likelihood` and `own_at_depmix` (the
#'   likelihood written here at depmixS4's estimates), `derivative_check`
#'   (from [three_step_derivative_check()]) and `newton`.
three_step_ml_reference <- function(z, posteriors, cluster, variance = c("varying", "equal")) {
  variance <- match.arg(variance)
  classes <- ncol(posteriors)
  modal <- three_step_modal(posteriors)
  errors <- three_step_error_matrix(posteriors, modal)
  sizes <- colMeans(posteriors)
  centre <- mean(z)
  centred <- z - centre
  log_weights <- log(t(errors[, modal, drop = FALSE])) +
    matrix(log(sizes), length(z), classes, byrow = TRUE)
  frame <- data.frame(W = factor(modal, levels = seq_len(classes)), z = centred)
  moments <- function(weights) {
    means <- as.vector(crossprod(weights, centred)) / colSums(weights)
    spread <- colSums(weights * outer(centred, means, "-")^2) / colSums(weights)
    sds <- if (identical(variance, "equal")) rep(sqrt(mean(spread)), classes) else sqrt(spread)
    list(means = means, sds = sds)
  }
  indicator <- vapply(seq_len(classes), function(class) as.numeric(modal == class),
                      numeric(length(z)))
  starts <- list(moments(posteriors), moments(indicator))
  # Parameter layout of depmixS4's mix(): the prior, then per class the
  # multinomial probabilities of W and the Gaussian mean and sd.
  per_class <- classes + 2L
  free <- c(rep(0L, classes), rep(c(rep(0L, classes), 1L, 1L), classes))
  if (identical(variance, "equal")) {
    sd_positions <- classes + per_class * seq_len(classes)
    free[sd_positions] <- 2L
  }
  fits <- lapply(starts, function(start) {
    model <- depmixS4::mix(list(W ~ 1, z ~ 1), data = frame, nstates = classes,
                           family = list(depmixS4::multinomial("identity"),
                                         stats::gaussian()))
    parameters <- c(sizes, unlist(lapply(seq_len(classes), function(class) {
      c(errors[class, ], start$means[[class]], start$sds[[class]])
    })))
    model <- depmixS4::setpars(model, parameters)
    utils::capture.output(fitted <- depmixS4::fit(
      model, equal = free, method = "rsolnp", verbose = FALSE,
      solnpcntrl = list(rho = 1, outer.iter = 400, inner.iter = 800, delta = 1e-9,
                        tol = 1e-12)))
    fitted
  })
  likelihoods <- vapply(fits, function(fitted) as.numeric(depmixS4::logLik(fitted)), numeric(1))
  best <- fits[[which.max(likelihoods)]]
  parameters <- depmixS4::getpars(best)
  means <- parameters[classes + per_class * (seq_len(classes) - 1L) + classes + 1L]
  sds <- parameters[classes + per_class * seq_len(classes)]
  theta <- c(means, log(if (identical(variance, "equal")) sds[[1L]] else sds))
  model <- three_step_ml_model(centred, log_weights, variance)
  own_at_depmix <- model$total(theta)
  newton <- three_step_newton(model, theta)
  covariances <- three_step_covariances(model, newton$theta, cluster)
  list(means = centre + newton$theta[seq_len(classes)],
       variances = rep_len(exp(2 * newton$theta[-seq_len(classes)]), classes),
       covariance = covariances$robust[seq_len(classes), seq_len(classes), drop = FALSE],
       effective_n = colSums(model$responsibilities(newton$theta)),
       log_likelihood = newton$log_likelihood,
       depmix_log_likelihood = max(likelihoods), own_at_depmix = own_at_depmix,
       derivative_check = three_step_derivative_check(model, theta), newton = newton)
}

# ---------------------------------------------------------------------------
# R3STEP: membership regressed on covariates, D fixed
# ---------------------------------------------------------------------------

#' Class probabilities of a multinomial logit against the last class
#' @param design Units-by-terms design matrix.
#' @param coefficients Terms-by-(K - 1) coefficients.
#' @return Units-by-K log probabilities.
three_step_log_prior <- function(design, coefficients) {
  linear <- cbind(design %*% coefficients, 0)
  peak <- apply(linear, 1L, max)
  linear - peak - log(rowSums(exp(linear - peak)))
}

#' Per-unit log likelihood of R3STEP
#' @param theta Coefficients, terms varying fastest, in the centred chart.
#' @param design Centred design (intercept first).
#' @param log_emission Units-by-classes log D[t, W_i].
#' @return One log likelihood per unit.
three_step_r3_units <- function(theta, design, log_emission) {
  joint <- three_step_log_prior(design, matrix(theta, ncol(design))) + log_emission
  peak <- apply(joint, 1L, max)
  peak + log(rowSums(exp(joint - peak)))
}

#' The R3STEP likelihood with its score and Hessian in closed form
#'
#' l_i = log sum_t exp(eta_it) D[t, W_i] - log sum_t exp(eta_it), eta_i =
#' x_i' B with the last class's column zero. With r_i the responsibilities
#' (proportional to pi_it D[t, W_i]) and pi_i the class probabilities, the
#' score of column k of B is x_i (r_ik - pi_ik), and the Hessian block (k, l)
#' is sum_i x_i x_i' [(r_ik 1(k = l) - r_ik r_il) - (pi_ik 1(k = l) - pi_ik
#' pi_il)], the second derivative of the first log-sum-exp less that of the
#' second.
#' @param design Units-by-terms design (intercept first).
#' @param log_emission Units-by-classes log D[t, W_i].
#' @return A model, shaped like [three_step_numerical_model()]'s, with `kind =
#'   "closed form"`.
three_step_r3_model <- function(design, log_emission) {
  n_free <- ncol(log_emission) - 1L
  pieces <- function(theta) {
    log_prior <- three_step_log_prior(design, matrix(theta, ncol(design)))
    joint <- log_prior + log_emission
    peak <- apply(joint, 1L, max)
    units <- peak + log(rowSums(exp(joint - peak)))
    list(units = units, responsibilities = exp(joint - units), prior = exp(log_prior))
  }
  # Coefficients run terms fastest within a logit column, as theta does.
  unit_scores <- function(theta) {
    parts <- pieces(theta)
    do.call(cbind, lapply(seq_len(n_free), function(k) {
      (parts$responsibilities[, k] - parts$prior[, k]) * design
    }))
  }
  hessian <- function(theta) {
    parts <- pieces(theta)
    block <- function(k, l) {
      curvature <- function(p) p[, k] * (identical(k, l) - p[, l])
      crossprod(design, (curvature(parts$responsibilities) - curvature(parts$prior)) *
                  design)
    }
    do.call(rbind, lapply(seq_len(n_free), function(k) {
      do.call(cbind, lapply(seq_len(n_free), function(l) block(k, l)))
    }))
  }
  units <- function(theta) three_step_r3_units(theta, design, log_emission)
  list(units = units, total = function(theta) sum(units(theta)),
       score = function(theta) colSums(unit_scores(theta)),
       unit_scores = unit_scores, hessian = hessian, kind = "closed form")
}

#' Move multinomial-logit coefficients between a centred and a natural chart
#'
#' With x = c + x_c, the intercept in the natural chart is the centred
#' intercept minus sum_j b_j c_j; slopes are unchanged.
#' @param centres The covariate centres.
#' @param n_free Number of logit columns.
#' @return The matrix taking centred coefficients (terms varying fastest) to
#'   natural ones.
three_step_uncentre <- function(centres, n_free) {
  block <- diag(length(centres) + 1L)
  block[1L, -1L] <- -centres
  kronecker(diag(n_free), block)
}

#' The R3STEP reference
#'
#' depmixS4 fits the mixture `list(W ~ 1)` with `prior = ~ covariates` and the
#' response probabilities of W fixed at D, from zero and from the
#' proportional-assignment logit (the posterior-weighted multinomial logit,
#' concave), keeping the higher maximum. depmixS4's prior is a multinomial
#' logit against the first class; its coefficients are re-expressed against
#' the last, as the package reports them, and then carried to the fixed point
#' by Newton steps on the likelihood written here ([three_step_r3_model()],
#' closed-form score and Hessian), in the centred chart.
#' @param covariates Units-by-covariates numeric matrix with column names.
#' @param posteriors Units-by-classes step-one posteriors.
#' @param cluster Cluster of each unit.
#' @return A list with `estimates` (natural chart, terms fastest),
#'   `observed`, `robust`, `terms`, `log_likelihood`, `depmix_log_likelihood`,
#'   `own_at_depmix`, `derivative_check` and `newton`.
three_step_r3_reference <- function(covariates, posteriors, cluster) {
  classes <- ncol(posteriors)
  n_free <- classes - 1L
  modal <- three_step_modal(posteriors)
  errors <- three_step_error_matrix(posteriors, modal)
  log_emission <- log(t(errors[, modal, drop = FALSE]))
  centres <- colMeans(covariates)
  centred <- sweep(covariates, 2L, centres)
  design <- cbind(`(Intercept)` = 1, centred)
  terms <- c("(Intercept)", colnames(covariates))
  n_terms <- length(terms)
  frame <- data.frame(W = factor(modal, levels = seq_len(classes)), centred)
  prior_formula <- stats::as.formula(paste("~", paste(colnames(covariates), collapse = " + ")))
  proportional <- stats::optim(
    rep(0, n_terms * n_free),
    function(theta) -sum(posteriors * three_step_log_prior(design, matrix(theta, n_terms))),
    method = "BFGS", control = list(maxit = 1000L, reltol = 1e-14))$par
  # Against the last class (this file) <-> against the first (depmixS4):
  # depmixS4's column k is b_k - b_1 with b_K = 0 here.
  to_depmix <- function(theta) {
    against_last <- cbind(matrix(theta, n_terms), 0)
    against_last - against_last[, 1L]
  }
  from_depmix <- function(against_first) {
    as.vector(against_first[, seq_len(n_free), drop = FALSE] - against_first[, classes])
  }
  fits <- lapply(list(rep(0, n_terms * n_free), proportional), function(start) {
    model <- depmixS4::mix(list(W ~ 1), data = frame, nstates = classes,
                           family = list(depmixS4::multinomial("identity")),
                           prior = prior_formula, initdata = frame)
    # getpars() lists the prior's terms-by-classes matrix row by row.
    parameters <- c(as.vector(t(to_depmix(start))), as.vector(t(errors)))
    model <- depmixS4::setpars(model, parameters)
    fixed <- c(rep(c(TRUE, rep(FALSE, n_free)), n_terms), rep(TRUE, classes^2))
    utils::capture.output(fitted <- depmixS4::fit(
      model, fixed = fixed, method = "rsolnp", verbose = FALSE,
      solnpcntrl = list(rho = 1, outer.iter = 400, inner.iter = 800, delta = 1e-9,
                        tol = 1e-12)))
    fitted
  })
  likelihoods <- vapply(fits, function(fitted) as.numeric(depmixS4::logLik(fitted)), numeric(1))
  best <- fits[[which.max(likelihoods)]]
  against_first <- matrix(depmixS4::getpars(best)[seq_len(n_terms * classes)], n_terms,
                          classes, byrow = TRUE)
  theta <- from_depmix(against_first)
  model <- three_step_r3_model(design, log_emission)
  own_at_depmix <- model$total(theta)
  newton <- three_step_newton(model, theta)
  covariances <- three_step_covariances(model, newton$theta, cluster)
  natural <- three_step_uncentre(centres, n_free)
  list(estimates = as.vector(natural %*% newton$theta),
       observed = natural %*% covariances$observed %*% t(natural),
       robust = natural %*% covariances$robust %*% t(natural),
       terms = terms, log_likelihood = newton$log_likelihood,
       depmix_log_likelihood = max(likelihoods), own_at_depmix = own_at_depmix,
       derivative_check = three_step_derivative_check(model, theta), newton = newton)
}

#' The table r3step() reports
#'
#' Wald statistics, two-sided normal p-values, the correction across the
#' covariate terms (never the intercepts) and normal intervals.
#' @param labels A `data.frame` of the key columns, one row per coefficient.
#' @param estimates,covariance The coefficients and their covariance.
#' @param tested Which rows enter the multiplicity correction.
#' @param adjust A p.adjust() method.
#' @param ci_level Confidence level.
#' @return `labels` with the value columns appended.
three_step_coefficient_table <- function(labels, estimates, covariance, tested,
                                         adjust = "BH", ci_level = 0.95) {
  errors <- sqrt(diag(covariance))
  statistic <- estimates / errors
  p_value <- 2 * stats::pnorm(-abs(statistic))
  adjusted <- rep(NA_real_, length(p_value))
  adjusted[tested] <- stats::p.adjust(p_value[tested], method = adjust)
  quantile <- stats::qnorm(1 - (1 - ci_level) / 2)
  data.frame(labels, estimate = estimates, standard_error = errors, statistic = statistic,
             p_value = p_value, p_adjusted = adjusted,
             conf_low = estimates - quantile * errors,
             conf_high = estimates + quantile * errors, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# R3STEP with group-class intercepts
# ---------------------------------------------------------------------------

#' Per-group log likelihood of R3STEP with group-class intercepts
#'
#' log L_g = log sum_h eta_h prod_{i in g} sum_k pi_hk(x_i) D[k, W_i], with
#' pi_h a multinomial logit against the last profile whose intercepts depend
#' on the group class and whose slopes do not, and eta a logit against the
#' last group class.
#' @param theta Intercepts (group classes fastest, then profiles), slopes
#'   (covariates fastest), then the group-class logits; centred chart.
#' @param covariates Centred units-by-covariates matrix.
#' @param log_emission Units-by-profiles log D[k, W_i].
#' @param group Dense group index of each unit (1..G).
#' @param n_group_classes Number of group classes.
#' @return One log likelihood per group.
three_step_group_units <- function(theta, covariates, log_emission, group, n_group_classes) {
  n_free <- ncol(log_emission) - 1L
  n_alpha <- n_group_classes * n_free
  n_beta <- ncol(covariates) * n_free
  alpha <- matrix(theta[seq_len(n_alpha)], n_group_classes, n_free)
  beta <- matrix(theta[n_alpha + seq_len(n_beta)], ncol(covariates), n_free)
  log_eta <- c(theta[n_alpha + n_beta + seq_len(n_group_classes - 1L)], 0)
  log_eta <- log_eta - max(log_eta) - log(sum(exp(log_eta - max(log_eta))))
  by_class <- vapply(seq_len(n_group_classes), function(h) {
    linear <- cbind(sweep(covariates %*% beta, 2L, alpha[h, ], "+"), 0)
    peak <- apply(linear, 1L, max)
    joint <- linear - peak - log(rowSums(exp(linear - peak))) + log_emission
    top <- apply(joint, 1L, max)
    as.vector(rowsum(top + log(rowSums(exp(joint - top))), group, reorder = FALSE)) +
      log_eta[[h]]
  }, numeric(max(group)))
  by_class <- matrix(by_class, max(group))
  peak <- apply(by_class, 1L, max)
  peak + log(rowSums(exp(by_class - peak)))
}

#' The reference for R3STEP with group-class intercepts
#'
#' No external program fits this model, so the maximum is found here: BFGS
#' from the step-one solution (each group class's profile mix as intercepts,
#' no slopes, the group-class shares), then Newton steps to the fixed point
#' with numDeriv's derivatives ([three_step_numerical_model()]).
#' Starting from step one keeps the group classes in the step-one order.
#' @param covariates Units-by-covariates numeric matrix with column names.
#' @param posteriors Units-by-profiles step-one (marginal) posteriors.
#' @param group Group label of each unit.
#' @param profile_mix Group-classes-by-profiles step-one P(profile | group class).
#' @param group_shares Step-one group-class probabilities.
#' @return A list with `estimates` (natural chart), `observed`, `robust`,
#'   `log_likelihood` and `newton`.
three_step_group_reference <- function(covariates, posteriors, group, profile_mix,
                                       group_shares) {
  classes <- ncol(posteriors)
  n_free <- classes - 1L
  n_group_classes <- nrow(profile_mix)
  modal <- three_step_modal(posteriors)
  errors <- three_step_error_matrix(posteriors, modal)
  log_emission <- log(t(errors[, modal, drop = FALSE]))
  dense <- match(group, unique(group))
  centres <- colMeans(covariates)
  centred <- sweep(covariates, 2L, centres)
  n_alpha <- n_group_classes * n_free
  n_beta <- ncol(covariates) * n_free
  start <- c(as.vector(log(profile_mix[, seq_len(n_free), drop = FALSE] /
                             profile_mix[, classes])),
             rep(0, n_beta),
             log(group_shares[seq_len(n_group_classes - 1L)] / group_shares[[n_group_classes]]))
  model <- three_step_numerical_model(function(point) {
    three_step_group_units(point, centred, log_emission, dense, n_group_classes)
  })
  searched <- stats::optim(start, function(point) -model$total(point), method = "BFGS",
                           control = list(maxit = 5000L, reltol = 1e-14))
  newton <- three_step_newton(model, searched$par)
  covariances <- three_step_covariances(model, newton$theta, seq_len(max(dense)))
  # The centring moves every group class's intercepts by the same slope term.
  natural <- diag(length(newton$theta))
  natural[seq_len(n_alpha), n_alpha + seq_len(n_beta)] <-
    kronecker(diag(n_free), -matrix(1, n_group_classes, 1L) %*% t(centres))
  list(estimates = as.vector(natural %*% newton$theta),
       observed = natural %*% covariances$observed %*% t(natural),
       robust = natural %*% covariances$robust %*% t(natural),
       log_likelihood = newton$log_likelihood, newton = newton)
}

# ---------------------------------------------------------------------------
# The three-step latent Markov model
# ---------------------------------------------------------------------------

#' Per-person log likelihood of the step-three chain
#'
#' Forward recursion over each person's occasions, scaled at every step:
#' alpha_1 = pi(x_0) * e(W_1), alpha_t = (alpha_{t-1} A_t) * e(W_t), with e(W)
#' = D[, W] and A_t a logit against staying whose covariates are read at the
#' destination occasion t.
#' @param theta Initial coefficients (terms fastest, profiles 1..K-1), then
#'   transition coefficients (terms fastest, then destination in profile order
#'   skipping the origin, then origin); centred chart.
#' @param chain A list with `emission` (persons-by-occasions-by-profiles log
#'   D[k, W_pt] array), `initial` (persons-by-terms design) and `transition`
#'   (persons-by-occasions-by-terms design; occasion 1 unused).
#' @return One log likelihood per person.
three_step_chain_units <- function(theta, chain) {
  k <- dim(chain$emission)[3L]
  p0 <- ncol(chain$initial)
  p <- dim(chain$transition)[3L]
  n_initial <- p0 * (k - 1L)
  initial <- matrix(theta[seq_len(n_initial)], p0, k - 1L)
  moves <- array(theta[-seq_len(n_initial)], c(p, k - 1L, k))
  log_initial <- three_step_log_prior(chain$initial, initial)
  start <- log_initial + chain$emission[, 1L, ]
  forward <- function(state, occasion) {
    logs <- three_step_chain_log_transition(three_step_chain_design(chain, occasion), moves)
    peak <- apply(state, 1L, max)
    carried <- exp(state - peak)
    to <- vapply(seq_len(k), function(l) {
      rowSums(carried * exp(logs[, l, ]))
    }, numeric(nrow(state)))
    peak + log(matrix(to, nrow(state))) + chain$emission[, occasion, ]
  }
  # The recursion needs the previous occasion's alpha; Reduce() carries it.
  last <- Reduce(forward, seq_len(dim(chain$emission)[2L])[-1L], init = start)
  peak <- apply(last, 1L, max)
  peak + log(rowSums(exp(last - peak)))
}

#' Log transition probabilities of the step-three chain at one occasion
#'
#' log P(to l | from j, x): a multinomial logit against staying (j), so the
#' origin's column of the linear predictor is zero.
#' @param design Persons-by-terms transition design at the occasion.
#' @param moves Terms-by-(K - 1)-by-K coefficients (destinations other than
#'   the origin in profile order, then the origin).
#' @return A persons-by-to-by-from array of log probabilities.
three_step_chain_log_transition <- function(design, moves) {
  k <- dim(moves)[3L]
  rows <- lapply(seq_len(k), function(from) {
    linear <- matrix(0, nrow(design), k)
    linear[, -from] <- design %*% matrix(moves[, , from], dim(moves)[1L])
    peak <- apply(linear, 1L, max)
    linear - peak - log(rowSums(exp(linear - peak)))
  })
  array(unlist(rows), c(nrow(design), k, k))
}

#' The transition design of one occasion, persons by terms
#' @param chain As for [three_step_chain_units()].
#' @param occasion The occasion.
#' @return A matrix.
three_step_chain_design <- function(chain, occasion) {
  matrix(chain$transition[, occasion, ], dim(chain$transition)[1L])
}

#' The step-three chain with its score in closed form
#'
#' The score of a hidden Markov likelihood is the expected complete-data
#' score under the posterior (Fisher's identity), from a scaled
#' forward-backward pass: for the initial logit of profile m, x_0 (gamma_1m -
#' pi_m(x_0)); for the transition logit from j to l != j, sum_t x_t
#' (xi_t(j, l) - gamma_{t-1}(j) A_t(j, l)), with gamma the occasion
#' posteriors and xi the two-slice posteriors. The Hessian is numDeriv's
#' Jacobian of this score (a first derivative of a smooth function), made
#' symmetric. A numDeriv gradient of the likelihood carries errors of about
#' 1e-7 here, which is where Newton steps on it stalled.
#' @param chain As for [three_step_chain_units()].
#' @return A model, shaped like [three_step_numerical_model()]'s, with `kind =
#'   "closed form"`.
three_step_chain_model <- function(chain) {
  k <- dim(chain$emission)[3L]
  p0 <- ncol(chain$initial)
  p <- dim(chain$transition)[3L]
  n_occasions <- dim(chain$emission)[2L]
  n_initial <- p0 * (k - 1L)
  emission <- exp(chain$emission)
  unit_scores <- function(theta) {
    moves <- array(theta[-seq_len(n_initial)], c(p, k - 1L, k))
    prior <- exp(three_step_log_prior(chain$initial,
                                      matrix(theta[seq_len(n_initial)], p0, k - 1L)))
    later <- seq_len(n_occasions)[-1L]
    designs <- lapply(seq_len(n_occasions), function(t) three_step_chain_design(chain, t))
    transitions <- lapply(seq_len(n_occasions), function(t) {
      if (t == 1L) NULL else exp(three_step_chain_log_transition(designs[[t]], moves))
    })
    # Scaled forward pass: alpha_t sums to one in each row, c_t the scale.
    first <- prior * emission[, 1L, ]
    forward <- Reduce(function(state, t) {
      to <- vapply(seq_len(k), function(l) rowSums(state$alpha * transitions[[t]][, l, ]),
                   numeric(nrow(first))) |> matrix(nrow(first))
      unscaled <- to * emission[, t, ]
      list(alpha = unscaled / rowSums(unscaled), scale = rowSums(unscaled))
    }, later, init = list(alpha = first / rowSums(first), scale = rowSums(first)),
    accumulate = TRUE)
    # Backward pass, beta_T = 1, scaled by the same c_t.
    backward <- Reduce(function(t, beta) {
      carried <- emission[, t, ] * beta / forward[[t]]$scale
      vapply(seq_len(k), function(j) rowSums(transitions[[t]][, , j] * carried),
             numeric(nrow(first))) |> matrix(nrow(first))
    }, later, accumulate = TRUE, right = TRUE,
    init = matrix(1, nrow(first), k))
    gamma <- Map(function(f, b) f$alpha * b, forward, backward)
    initial_scores <- do.call(cbind, lapply(seq_len(k - 1L), function(m) {
      (gamma[[1L]][, m] - prior[, m]) * chain$initial
    }))
    transition_scores <- do.call(cbind, lapply(seq_len(k), function(j) {
      do.call(cbind, lapply(seq_len(k)[-j], function(l) {
        Reduce(`+`, lapply(later, function(t) {
          xi <- forward[[t - 1L]]$alpha[, j] * transitions[[t]][, l, j] *
            emission[, t, l] * backward[[t]][, l] / forward[[t]]$scale
          (xi - gamma[[t - 1L]][, j] * transitions[[t]][, l, j]) * designs[[t]]
        }))
      }))
    }))
    cbind(initial_scores, transition_scores)
  }
  score <- function(theta) colSums(unit_scores(theta))
  units <- function(theta) three_step_chain_units(theta, chain)
  list(units = units, total = function(theta) sum(units(theta)), score = score,
       unit_scores = unit_scores,
       hessian = function(theta) {
         hessian <- numDeriv::jacobian(score, theta)
         (hessian + t(hessian)) / 2
       },
       kind = "closed form")
}

#' The three-step latent Markov reference
#'
#' Step two classifies each occasion from its own indicators under the
#' step-one measurement model (diagonal Gaussian here) and the profile shares
#' pooled over occasions, p(k | y_it) proportional to w_k f_k(y_it); D and the
#' modal profiles follow. Step three is fitted by depmixS4's hidden Markov
#' model `depmix(list(W ~ 1), prior = ~ initial, transition = ~ transition)`
#' with the response probabilities of W fixed at D. depmixS4 reads a
#' transition covariate at the origin row, so its transition covariate column
#' carries the destination occasion's value (the next row's). Its logits
#' (against the first profile in each row) are re-expressed against the last
#' profile (initial) and against staying (transitions), then carried to the
#' fixed point by Newton steps on the recursion written here, with its
#' forward-backward score ([three_step_chain_model()]).
#' @param data The step-one data, one row per person and occasion, sorted by
#'   person and occasion, every person with the same occasions.
#' @param id,vars The person column and the indicators.
#' @param means,variances Profiles-by-indicators step-one estimates.
#' @param shares Profile shares pooled over occasions.
#' @param initial_covariates,transition_covariates Covariate names.
#' @return A list with `estimates` (natural chart), `observed`, `robust`,
#'   `labels`, `tested`, `errors` (step two's D, rows the true profile),
#'   `log_likelihood`, `depmix_log_likelihood`, `own_at_depmix`,
#'   `derivative_check` and `newton`.
three_step_chain_reference <- function(data, id, vars, means, variances, shares,
                                       initial_covariates, transition_covariates) {
  k <- nrow(means)
  log_density <- vapply(seq_len(k), function(profile) {
    log(shares[[profile]]) + rowSums(vapply(vars, function(v) {
      stats::dnorm(data[[v]], means[profile, v], sqrt(variances[profile, v]), log = TRUE)
    }, numeric(nrow(data))))
  }, numeric(nrow(data)))
  occasion_posteriors <- exp(log_density - apply(log_density, 1L, max))
  occasion_posteriors <- occasion_posteriors / rowSums(occasion_posteriors)
  modal <- three_step_modal(occasion_posteriors)
  errors <- three_step_error_matrix(occasion_posteriors, modal)
  persons <- unique(data[[id]])
  n_persons <- length(persons)
  n_occasions <- nrow(data) / n_persons
  first_rows <- !duplicated(data[[id]])
  initial_values <- as.matrix(data[first_rows, initial_covariates, drop = FALSE])
  transition_values <- as.matrix(data[transition_covariates])
  later <- stats::ave(seq_len(nrow(data)), data[[id]], FUN = seq_along) > 1L
  initial_centres <- colMeans(initial_values)
  transition_centres <- colMeans(transition_values[later, , drop = FALSE])
  initial_design <- cbind(1, sweep(initial_values, 2L, initial_centres))
  transition_rows <- cbind(1, sweep(transition_values, 2L, transition_centres))
  emission <- array(log(t(errors)[modal, , drop = FALSE]), c(n_occasions, n_persons, k)) |>
    aperm(c(2L, 1L, 3L))
  chain <- list(emission = emission, initial = initial_design,
                transition = aperm(array(transition_rows,
                                         c(n_occasions, n_persons, ncol(transition_rows))),
                                   c(2L, 1L, 3L)))
  p0 <- ncol(initial_design)
  p <- ncol(transition_rows)
  # depmixS4: transition covariates at the origin row carry the destination's
  # value; the last occasion's row is never used as an origin.
  next_value <- function(column) {
    shifted <- c(column[-1L], column[[length(column)]])
    replace(shifted, !c(later[-1L], FALSE), 0)
  }
  frame <- data.frame(W = factor(modal, levels = seq_len(k)),
                      sweep(initial_values[rep(seq_len(n_persons), each = n_occasions), ,
                                           drop = FALSE], 2L, initial_centres),
                      vapply(seq_along(transition_covariates), function(j) {
                        next_value(transition_rows[, j + 1L])
                      }, numeric(nrow(data))) |>
                        matrix(nrow(data), dimnames = list(NULL, paste0(".next_",
                                                                       transition_covariates))))
  model <- depmixS4::depmix(
    list(W ~ 1), data = frame, nstates = k,
    family = list(depmixS4::multinomial("identity")),
    prior = stats::as.formula(paste("~", paste(initial_covariates, collapse = " + "))),
    transition = stats::as.formula(paste("~", paste0(".next_", transition_covariates,
                                                     collapse = " + "))),
    ntimes = rep(n_occasions, n_persons),
    initdata = frame[first_rows, , drop = FALSE])
  n_prior <- p0 * k
  n_transition <- p * k * k
  parameters <- depmixS4::getpars(model)
  parameters[n_prior + n_transition + seq_len(k * k)] <- as.vector(t(errors))
  model <- depmixS4::setpars(model, parameters)
  # The first profile's column of every logit is depmixS4's base, fixed at 0.
  base <- function(n_terms) rep(c(TRUE, rep(FALSE, k - 1L)), n_terms)
  fixed <- c(base(p0), rep(base(p), k), rep(TRUE, k * k))
  utils::capture.output(fitted <- depmixS4::fit(
    model, fixed = fixed, method = "rsolnp", verbose = FALSE,
    solnpcntrl = list(rho = 1, outer.iter = 400, inner.iter = 800, delta = 1e-9,
                      tol = 1e-12)))
  estimates <- depmixS4::getpars(fitted)
  prior_first <- matrix(estimates[seq_len(n_prior)], p0, k, byrow = TRUE)
  initial <- prior_first[, seq_len(k - 1L), drop = FALSE] - prior_first[, k]
  moves <- unlist(lapply(seq_len(k), function(from) {
    row <- matrix(estimates[n_prior + (from - 1L) * p * k + seq_len(p * k)], p, k,
                  byrow = TRUE)
    as.vector(row[, -from, drop = FALSE] - row[, from])
  }))
  theta <- c(as.vector(initial), moves)
  model <- three_step_chain_model(chain)
  own_at_depmix <- model$total(theta)
  newton <- three_step_newton(model, theta)
  covariances <- three_step_covariances(model, newton$theta, seq_len(n_persons))
  natural <- as.matrix(three_step_block_diagonal(list(
    three_step_uncentre(initial_centres, k - 1L),
    three_step_uncentre(transition_centres, k * (k - 1L)))))
  profiles <- paste0("profile_", seq_len(k))
  initial_terms <- c("(Intercept)", initial_covariates)
  transition_terms <- c("(Intercept)", transition_covariates)
  labels <- rbind(
    data.frame(block = "initial", from = NA_character_,
               to = rep(profiles[-k], each = p0), term = rep(initial_terms, k - 1L),
               stringsAsFactors = FALSE),
    do.call(rbind, lapply(seq_len(k), function(from) {
      data.frame(block = "transition", from = profiles[[from]],
                 to = rep(profiles[-from], each = p), term = rep(transition_terms, k - 1L),
                 stringsAsFactors = FALSE)
    })))
  list(estimates = as.vector(natural %*% newton$theta),
       observed = natural %*% covariances$observed %*% t(natural),
       robust = natural %*% covariances$robust %*% t(natural),
       labels = labels, tested = labels$term != "(Intercept)", errors = errors,
       log_likelihood = newton$log_likelihood,
       depmix_log_likelihood = as.numeric(depmixS4::logLik(fitted)),
       own_at_depmix = own_at_depmix,
       derivative_check = three_step_derivative_check(model, theta), newton = newton)
}

#' Block-diagonal matrix from a list of square blocks, in base R
#' @param blocks A list of matrices.
#' @return The block-diagonal matrix.
three_step_block_diagonal <- function(blocks) {
  sizes <- vapply(blocks, nrow, integer(1))
  ends <- cumsum(sizes)
  result <- matrix(0, sum(sizes), sum(sizes))
  Reduce(function(out, b) {
    rows <- (ends[[b]] - sizes[[b]]) + seq_len(sizes[[b]])
    out[rows, rows] <- blocks[[b]]
    out
  }, seq_along(blocks), result)
}
