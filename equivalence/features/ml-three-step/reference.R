# An independent implementation of the maximum-likelihood step three for a
# continuous distal outcome (Vermunt 2010; Bakk, Tekle and Vermunt 2013):
#
#   L = sum_i log sum_t pi_t D[t, W_i] phi(z_i; mu_t, sigma_t),
#
# with W_i the modal class, D[t, s] the probability that a unit of class t is
# assigned to class s (from the step-one posteriors), and pi_t either fixed at
# the average posterior (Latent GOLD's Step3-Dependent, latents) or estimated
# (Mplus DU3STEP/DE3STEP). Written from the formulas, without calling any
# latents function: optim (BFGS) from several starts, Newton polish on a
# numDeriv Jacobian of the gradient, numDeriv::hessian for the information.
#
# Coordinates: class means, log standard deviations (one per class, or one
# shared), then for free class sizes the logits of classes 2..K against 1.
# Fitting runs on the standardised outcome; means, variances and covariances
# are carried back to the outcome's units.

ref_modal <- function(posteriors) max.col(posteriors, ties.method = "first")

ref_error_matrix <- function(posteriors, modal = ref_modal(posteriors)) {
  classes <- seq_len(ncol(posteriors))
  assigned <- outer(modal, classes, "==") * 1
  crossprod(posteriors, assigned) / colSums(posteriors)
}

ref_unpack <- function(theta, n_classes, variance, sizes) {
  means <- theta[seq_len(n_classes)]
  n_sd <- if (identical(variance, "equal")) 1L else n_classes
  log_sd <- theta[n_classes + seq_len(n_sd)]
  sds <- exp(if (identical(variance, "equal")) rep(log_sd, n_classes) else log_sd)
  logits <- if (identical(sizes, "free")) {
    c(0, theta[n_classes + n_sd + seq_len(n_classes - 1L)])
  } else NULL
  list(means = means, sds = sds, logits = logits, n_sd = n_sd)
}

# Joint log terms log pi_t + log D[t, W_i] + log phi(z_i), units by classes.
ref_joint <- function(theta, problem) {
  p <- ref_unpack(theta, problem$n_classes, problem$variance, problem$sizes)
  log_pi <- if (is.null(p$logits)) log(problem$prior) else {
    p$logits - max(p$logits) - log(sum(exp(p$logits - max(p$logits))))
  }
  n <- length(problem$z)
  log_d <- log(t(problem$errors)[problem$modal, , drop = FALSE])
  density <- vapply(seq_len(problem$n_classes),
                    \(t) stats::dnorm(problem$z, p$means[t], p$sds[t], log = TRUE),
                    numeric(n))
  list(joint = log_d + matrix(log_pi, n, problem$n_classes, byrow = TRUE) +
         matrix(density, n), unpacked = p, log_pi = log_pi)
}

ref_log_sum_exp <- function(joint) {
  top <- apply(joint, 1L, max)
  top + log(rowSums(exp(joint - top)))
}

ref_loglik <- function(theta, problem) sum(ref_log_sum_exp(ref_joint(theta, problem)$joint))

# Per-unit scores, units by coordinates.
ref_scores <- function(theta, problem) {
  state <- ref_joint(theta, problem)
  p <- state$unpacked
  r <- exp(state$joint - ref_log_sum_exp(state$joint))
  n <- length(problem$z)
  residual <- outer(problem$z, p$means, "-")
  sd_matrix <- matrix(p$sds, n, problem$n_classes, byrow = TRUE)
  d_mean <- r * residual / sd_matrix^2
  d_log_sd <- r * (residual^2 / sd_matrix^2 - 1)
  if (identical(problem$variance, "equal")) d_log_sd <- matrix(rowSums(d_log_sd), ncol = 1L)
  d_logit <- if (identical(problem$sizes, "free")) {
    pi <- exp(state$log_pi)
    (r - matrix(pi, n, problem$n_classes, byrow = TRUE))[, -1L, drop = FALSE]
  } else NULL
  cbind(d_mean, d_log_sd, d_logit)
}

ref_starts <- function(problem) {
  z <- problem$z
  k <- problem$n_classes
  moments <- function(weights) {
    w <- weights / matrix(colSums(weights), nrow(weights), k, byrow = TRUE)
    means <- as.vector(crossprod(w, z))
    variances <- colSums(w * outer(z, means, "-")^2)
    variances <- pmax(variances, 1e-4)
    log_sd <- if (identical(problem$variance, "equal")) {
      log(sqrt(sum(colSums(weights) * variances) / sum(weights)))
    } else log(sqrt(variances))
    logits <- if (identical(problem$sizes, "free")) {
      share <- colSums(weights) / sum(weights)
      log(share[-1L] / share[1L])
    } else NULL
    c(means, log_sd, logits)
  }
  modal_weights <- outer(problem$modal, seq_len(k), "==") * 1
  starts <- list(moments(problem$posteriors))
  if (all(colSums(modal_weights) > 1)) starts <- c(starts, list(moments(modal_weights)))
  # Deterministic perturbations of the proportional start: the means pushed
  # apart, pulled together and shifted by a fixed pattern.
  base <- starts[[1L]]
  pattern <- cos(seq_len(k) * 2.1)
  extra <- lapply(c(-0.75, 0.5, 1.5), \(size) {
    moved <- base
    moved[seq_len(k)] <- base[seq_len(k)] + size * pattern
    moved
  })
  # Twenty random starts (seeded locally, the caller's stream restored):
  # means at random outcome quantiles, common log standard deviation.
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) saved <- get(".Random.seed", envir = globalenv())
  on.exit(if (had_seed) assign(".Random.seed", saved, envir = globalenv()) else
    rm(".Random.seed", envir = globalenv()), add = TRUE)
  set.seed(20261020L)
  random <- lapply(seq_len(20L), \(r) {
    moved <- base
    moved[seq_len(k)] <- stats::quantile(z, stats::runif(k), names = FALSE)
    n_sd <- length(base) - k - (if (identical(problem$sizes, "free")) k - 1L else 0L)
    moved[k + seq_len(n_sd)] <- log(stats::runif(n_sd, 0.3, 1))
    moved
  })
  c(starts, extra, random)
}

ref_newton <- function(theta, problem, steps = 20L) {
  gradient <- function(x) colSums(ref_scores(x, problem))
  step_once <- function(state, index) {
    if (isTRUE(state$done)) return(state)
    g <- gradient(state$theta)
    if (max(abs(g)) < 1e-11 * length(problem$z)) return(c(state["theta"], list(done = TRUE)))
    h <- numDeriv::jacobian(gradient, state$theta)
    h <- (h + t(h)) / 2
    direction <- tryCatch(solve(h, g), error = \(e) NULL)
    if (is.null(direction)) return(c(state["theta"], list(done = TRUE)))
    current <- ref_loglik(state$theta, problem)
    sizes <- 2^-(0:30)
    values <- vapply(sizes, \(s) ref_loglik(state$theta - s * direction, problem), numeric(1))
    # Near the optimum the change is below the rounding of a sum over n
    # units; a step that loses no more than that still lowers the score.
    accepted <- which(values >= current - 1e-13 * abs(current))[1L]
    if (is.na(accepted)) return(c(state["theta"], list(done = TRUE)))
    list(theta = state$theta - sizes[accepted] * direction, done = FALSE)
  }
  Reduce(step_once, seq_len(steps), init = list(theta = theta, done = FALSE))$theta
}

#' Fit the step-three model.
#' @param posteriors Units-by-classes step-one posteriors.
#' @param z The outcome.
#' @param variance "varying" or "equal".
#' @param sizes "fixed" (average posterior) or "free".
#' @param clusters Optional cluster id per unit for the sandwich (default: units).
#' @param ci_level Confidence level for the intervals.
#' @param modal Optional modal assignment; max posterior by default.
#' @return A list of tables and scalars in the outcome's units.
ref_step3 <- function(posteriors, z, variance = c("varying", "equal"),
                      sizes = c("fixed", "free"), clusters = NULL, ci_level = 0.95,
                      modal = ref_modal(posteriors), adjust = "BH") {
  variance <- match.arg(variance)
  sizes <- match.arg(sizes)
  stopifnot("posteriors must be a matrix" = is.matrix(posteriors),
            "z must match the posteriors" = length(z) == nrow(posteriors),
            "z must vary" = stats::var(z) > 0)
  center <- mean(z)
  scale <- stats::sd(z)
  n_classes <- ncol(posteriors)
  problem <- list(z = (z - center) / scale, posteriors = posteriors, modal = modal,
                  errors = ref_error_matrix(posteriors, modal),
                  prior = colMeans(posteriors), n_classes = n_classes,
                  variance = variance, sizes = sizes)
  # A start whose variance collapses (an unbounded likelihood) fails inside
  # optim or Newton; it is recorded as failed with its message, not dropped.
  runs <- lapply(ref_starts(problem), \(start) {
    tryCatch({
      found <- stats::optim(start, ref_loglik, \(x, problem) colSums(ref_scores(x, problem)),
                            problem = problem, method = "BFGS",
                            control = list(fnscale = -1, maxit = 5000, reltol = 1e-15))
      theta <- ref_newton(found$par, problem)
      # A class standard deviation below 1% of the outcome's is a spike on
      # one or a few units: the varying-variance likelihood is unbounded
      # there, so such a maximum is degenerate and not a candidate.
      sds <- ref_unpack(theta, n_classes, variance, sizes)$sds
      list(theta = theta, loglik = ref_loglik(theta, problem), failure = NA_character_,
           degenerate = min(sds) < 0.01)
    }, error = \(e) list(theta = start, loglik = -Inf, failure = conditionMessage(e),
                         degenerate = TRUE))
  })
  logliks <- vapply(runs, `[[`, numeric(1), "loglik")
  degenerate <- vapply(runs, `[[`, logical(1), "degenerate")
  best_degenerate <- max(c(-Inf, logliks[degenerate & is.finite(logliks)]))
  logliks[!is.finite(logliks) | degenerate] <- -Inf
  if (!any(is.finite(logliks))) {
    stop(errorCondition(sprintf("every start failed: %s", runs[[1L]]$failure),
                        class = "reference_no_maximum"))
  }
  best <- runs[[which.max(logliks)]]
  theta <- best$theta
  # The information from Richardson differences of the analytic gradient; the
  # Hessian of the log likelihood itself (numDeriv::hessian) is kept as a
  # check. On a nearly unidentified model (classes with one outcome
  # distribution) the second is only accurate to about 1e-4.
  jacobian <- numDeriv::jacobian(\(x) colSums(ref_scores(x, problem)), theta)
  information <- -(jacobian + t(jacobian)) / 2
  hessian <- numDeriv::hessian(ref_loglik, theta, problem = problem)
  hessian_gap <- max(abs(-hessian - information)) / max(abs(information))
  inverse <- solve(information)
  scores <- ref_scores(theta, problem)
  units <- clusters %||% seq_along(z)
  summed <- rowsum(scores, units, reorder = FALSE)
  g <- nrow(summed)
  meat <- crossprod(summed) * g / (g - 1)
  robust <- inverse %*% meat %*% inverse
  p <- ref_unpack(theta, n_classes, variance, sizes)
  mean_block <- seq_len(n_classes)
  robust_means <- robust[mean_block, mean_block] * scale^2
  model_means <- inverse[mean_block, mean_block] * scale^2
  sd_block <- n_classes + seq_len(p$n_sd)
  variances <- (scale * p$sds)^2
  # SE of a variance from that of its log standard deviation: dv/dlog sd = 2 v.
  shared <- if (identical(variance, "equal")) variances[1L] else variances
  variance_se_model <- 2 * shared * sqrt(diag(inverse)[sd_block])
  variance_se_robust <- 2 * shared * sqrt(diag(robust)[sd_block])
  means <- center + scale * p$means
  quantile <- stats::qnorm(1 - (1 - ci_level) / 2)
  r <- exp(ref_joint(theta, problem)$joint - ref_log_sum_exp(ref_joint(theta, problem)$joint))
  se <- sqrt(diag(robust_means))
  pairs <- utils::combn(n_classes, 2L)
  contrast <- function(cov, a, b) cov[a, a] + cov[b, b] - 2 * cov[a, b]
  difference <- means[pairs[2L, ]] - means[pairs[1L, ]]
  pair_se <- sqrt(vapply(seq_len(ncol(pairs)), \(j) contrast(robust_means, pairs[1L, j],
                                                             pairs[2L, j]), numeric(1)))
  pair_se_model <- sqrt(vapply(seq_len(ncol(pairs)), \(j) contrast(model_means, pairs[1L, j],
                                                                   pairs[2L, j]), numeric(1)))
  statistic <- difference / pair_se
  p_value <- 2 * stats::pnorm(-abs(statistic))
  list(
    classes = data.frame(class = seq_len(n_classes), estimate = means,
                         standard_error = se, conf_low = means - quantile * se,
                         conf_high = means + quantile * se,
                         standard_error_model = sqrt(diag(model_means)),
                         variance = variances, effective_n = colSums(r),
                         size = exp(ref_joint(theta, problem)$log_pi)),
    pairs = data.frame(class = pairs[2L, ], reference_class = pairs[1L, ],
                       estimate = difference, standard_error = pair_se,
                       statistic = statistic, p_value = p_value,
                       p_adjusted = stats::p.adjust(p_value, adjust),
                       conf_low = difference - quantile * pair_se,
                       conf_high = difference + quantile * pair_se,
                       wald_model = (difference / pair_se_model)^2),
    variance_se_model = variance_se_model, variance_se_robust = variance_se_robust,
    log_likelihood = best$loglik - length(z) * log(scale),
    all_log_likelihoods = logliks - length(z) * log(scale),
    best_degenerate_log_likelihood = best_degenerate - length(z) * log(scale),
    max_score = max(abs(colSums(scores))), hessian_gap = hessian_gap,
    information_condition = kappa(information, exact = TRUE),
    errors = problem$errors, prior = problem$prior, modal = modal,
    robust_covariance = robust_means, model_covariance = model_means)
}
