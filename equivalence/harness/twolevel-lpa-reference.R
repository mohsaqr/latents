# The two-level latent profile model, written from its definition in base R.
#
# Sourced by the twolevel-lpa suite. Nothing here calls the package. The model
# is the nonparametric two-level mixture of Vermunt (2003) with Gaussian
# indicators:
#
#   f(y_g) = sum_w delta_w prod_{i in g} sum_c gamma_wc f_c(y_i),
#
# where f_c is the normal density of the indicators observed in row i, with
# profile c's mean and covariance restricted to those indicators (full-
# information maximum likelihood: a missing indicator is integrated out, which
# for a normal density is the marginal of the observed sub-vector, and a row
# with nothing observed contributes a factor of one). The covariance of a
# profile is diagonal with its own variances (`"varying"`), diagonal with one
# variance per indicator shared by the profiles (`"equal"`), or a full matrix
# of its own (`"full"`).
#
# Parameters are kept in their natural free coordinates: the means, the free
# variances (or the lower triangle of each covariance, by column), the first
# K - 1 probabilities of each row of gamma and the first G - 1 of delta. At a
# maximum the inverse observed information transforms with the Jacobian of
# any reparameterisation, so standard errors taken here are those of the
# natural parameters whatever scale the package estimates on.
#
# Vermunt, J. K. (2003). Multilevel latent class models. Sociological
# Methodology, 33(1), 213-239.

#' The shape of a two-level Gaussian model
#' @param indicators Indicator names.
#' @param n_profiles,n_group_classes Numbers of profiles (K) and group
#'   classes (G).
#' @param structure `"varying"`, `"equal"` or `"full"`.
#' @return A list with the sizes, the structure and the positions of each
#'   block in the free parameter vector.
twolevel_lpa_layout <- function(indicators, n_profiles, n_group_classes,
                                structure = c("varying", "equal", "full")) {
  structure <- match.arg(structure)
  j <- length(indicators)
  # The free covariance elements of one block, as (row, column) pairs: the
  # diagonal for a diagonal model, the lower triangle by column for a full one.
  cells <- if (identical(structure, "full")) {
    which(lower.tri(diag(j), diag = TRUE), arr.ind = TRUE)
  } else cbind(seq_len(j), seq_len(j))
  cells <- cells[order(cells[, 2L], cells[, 1L]), , drop = FALSE]
  n_blocks <- if (identical(structure, "equal")) 1L else n_profiles
  sizes <- c(means = n_profiles * j, covariance = nrow(cells) * n_blocks,
             gamma = n_group_classes * (n_profiles - 1L), delta = n_group_classes - 1L)
  ends <- cumsum(sizes)
  list(indicators = indicators, j = j, k = n_profiles, g = n_group_classes,
       structure = structure, cells = unname(cells), n_blocks = n_blocks,
       at = lapply(stats::setNames(names(sizes), names(sizes)), function(block) {
         ends[[block]] - sizes[[block]] + seq_len(sizes[[block]])
       }), q = sum(sizes))
}

#' Natural parameters to the free vector
#' @param parameters A list with `delta` (G), `gamma` (G x K), `means`
#'   (K x J) and `covariances` (J x J x K; for `"equal"` every slice is the
#'   same).
#' @param layout From [twolevel_lpa_layout()].
#' @return The free parameter vector.
twolevel_lpa_pack <- function(parameters, layout) {
  covariance <- vapply(seq_len(layout$n_blocks), function(k) {
    matrix(parameters$covariances[, , k], layout$j)[layout$cells]
  }, numeric(nrow(layout$cells)))
  c(as.vector(parameters$means), as.vector(covariance),
    as.vector(t(parameters$gamma[, -layout$k, drop = FALSE])),
    parameters$delta[-layout$g])
}

#' The free vector to natural parameters
#' @param theta The free parameter vector.
#' @param layout From [twolevel_lpa_layout()].
#' @return A list shaped like the `parameters` of [twolevel_lpa_pack()].
twolevel_lpa_unpack <- function(theta, layout) {
  k <- layout$k
  j <- layout$j
  values <- matrix(theta[layout$at$covariance], nrow(layout$cells), layout$n_blocks)
  covariances <- vapply(seq_len(k), function(profile) {
    m <- matrix(0, j, j)
    column <- values[, if (layout$n_blocks == 1L) 1L else profile]
    m[layout$cells] <- column
    m[layout$cells[, 2:1, drop = FALSE]] <- column
    m
  }, matrix(0, j, j))
  covariances <- array(covariances, c(j, j, k),
                       dimnames = list(layout$indicators, layout$indicators, NULL))
  first <- matrix(theta[layout$at$gamma], layout$g, k - 1L, byrow = TRUE)
  delta_first <- theta[layout$at$delta]
  list(delta = c(delta_first, 1 - sum(delta_first)),
       gamma = cbind(first, 1 - rowSums(first)),
       means = matrix(theta[layout$at$means], k, j, dimnames = list(NULL, layout$indicators)),
       covariances = covariances)
}

#' Whether natural parameters are inside the parameter space
#' @param parameters From [twolevel_lpa_unpack()].
#' @return `TRUE` when every probability is positive and every covariance
#'   positive definite.
.twolevel_lpa_admissible <- function(parameters) {
  all(parameters$delta > 0) && all(parameters$gamma > 0) &&
    all(vapply(seq_len(dim(parameters$covariances)[3L]), function(k) {
      min(eigen(parameters$covariances[, , k], symmetric = TRUE,
                only.values = TRUE)$values) > 0
    }, logical(1)))
}

#' Log density of each row's observed indicators under each profile
#'
#' Rows are grouped by their pattern of observed indicators, and each pattern
#' uses the Cholesky factor of the observed block of the covariance.
#'
#' @param parameters From [twolevel_lpa_unpack()].
#' @param x Numeric rows-by-indicators matrix, `NA` where missing.
#' @return A rows-by-profiles matrix of `log f_c(y_i)`; zero for a row with
#'   nothing observed.
twolevel_lpa_log_densities <- function(parameters, x) {
  observed <- !is.na(x)
  pattern <- apply(observed, 1L, function(row) paste(as.integer(row), collapse = ""))
  out <- matrix(0, nrow(x), nrow(parameters$means))
  by_pattern <- lapply(unique(pattern), function(key) {
    rows <- which(pattern == key)
    seen <- observed[rows[[1L]], ]
    if (!any(seen)) return(NULL)
    values <- x[rows, seen, drop = FALSE]
    densities <- vapply(seq_len(nrow(parameters$means)), function(k) {
      sigma <- matrix(parameters$covariances[seen, seen, k], sum(seen))
      root <- chol(sigma)
      standardised <- backsolve(root, t(values) - parameters$means[k, seen], transpose = TRUE)
      -0.5 * colSums(matrix(standardised, sum(seen))^2) - sum(log(diag(root))) -
        0.5 * sum(seen) * log(2 * pi)
    }, numeric(length(rows)))
    list(rows = rows, densities = matrix(densities, length(rows)))
  })
  by_pattern <- Filter(Negate(is.null), by_pattern)
  out_rows <- unlist(lapply(by_pattern, `[[`, "rows"))
  out[out_rows, ] <- do.call(rbind, lapply(by_pattern, `[[`, "densities"))
  out
}

#' Evaluate the two-level Gaussian model at given parameters
#'
#' @param parameters From [twolevel_lpa_unpack()].
#' @param x Numeric rows-by-indicators matrix.
#' @param groups Group identifier of each row.
#' @return A list: `log_likelihood`; `group_log_likelihood`, named by group in
#'   first-occurrence order; `group_posteriors`, groups by group classes with
#'   those row names; `posteriors`, each row's marginal profile posterior
#'   `P(c | y_g) = sum_w P(w | y_g) P(c | y_i, w)`; `conditional`, a list over
#'   group classes of the rows' `P(c | y_i, w)`; and `row_conditional_entropy`
#'   (see `two_level_row_conditional_entropy()`).
twolevel_lpa_evaluate <- function(parameters, x, groups) {
  group_ids <- unique(as.character(groups))
  group_index <- match(as.character(groups), group_ids)
  log_conditional <- twolevel_lpa_log_densities(parameters, x)
  peak <- apply(log_conditional, 1L, max)
  scaled <- exp(log_conditional - peak)
  mixed <- scaled %*% t(parameters$gamma)
  per_group <- rowsum(log(mixed) + peak, group_index, reorder = TRUE)
  joint <- sweep(per_group, 2L, log(parameters$delta), "+")
  top <- apply(joint, 1L, max)
  group_log_likelihood <- top + log(rowSums(exp(joint - top)))
  group_posteriors <- exp(joint - group_log_likelihood)
  conditional <- lapply(seq_along(parameters$delta), function(w) {
    within <- sweep(scaled, 2L, parameters$gamma[w, ], "*")
    within / rowSums(within)
  })
  posteriors <- Reduce(`+`, lapply(seq_along(parameters$delta), function(w) {
    conditional[[w]] * group_posteriors[group_index, w]
  }))
  row_conditional_entropy <- two_level_row_conditional_entropy(
    scaled, parameters$gamma, group_posteriors[group_index, , drop = FALSE])
  dimnames(group_posteriors) <- list(group_ids, NULL)
  list(log_likelihood = sum(group_log_likelihood),
       group_log_likelihood = stats::setNames(as.vector(group_log_likelihood), group_ids),
       group_posteriors = group_posteriors, posteriors = posteriors,
       conditional = conditional, row_conditional_entropy = row_conditional_entropy)
}

#' The log-likelihood as a function of the free vector
#' @inheritParams twolevel_lpa_evaluate
#' @param layout From [twolevel_lpa_layout()].
#' @param what `"total"` or `"groups"` (the per-group contributions).
#' @return A function of `theta`; outside the parameter space the total is
#'   `-Inf` (a refining step there is refused) and the per-group vector `NA`.
twolevel_lpa_objective <- function(layout, x, groups, what = c("total", "groups")) {
  what <- match.arg(what)
  function(theta) {
    parameters <- twolevel_lpa_unpack(theta, layout)
    if (!.twolevel_lpa_admissible(parameters)) {
      return(if (identical(what, "total")) -Inf else
        rep(NA_real_, length(unique(groups))))
    }
    evaluated <- twolevel_lpa_evaluate(parameters, x, groups)
    if (identical(what, "total")) evaluated$log_likelihood else
      evaluated$group_log_likelihood
  }
}

#' Carry estimates to the maximum by Newton steps with a fixed Hessian
#'
#' The gradient is the Richardson-extrapolated numerical gradient (numDeriv)
#' and the Hessian is taken once, at the start. The iteration therefore
#' converges linearly; each step must not lose likelihood, so the climb stays
#' on the starting mode.
#'
#' @param theta Starting free vector (Latent GOLD's printed estimates).
#' @param layout,x,groups As for [twolevel_lpa_objective()].
#' @param tolerance Largest absolute gradient accepted as a maximum.
#' @param max_steps Bound on the number of steps.
#' The numerical gradient bottoms out near 1e-7 here (a tolerance of 1e-8 is
#' never reached), so the result is resolved only to about 1e-8 in the
#' parameters; [twolevel_lpa_fixed_point()] finishes the climb.
#'
#' @return A list with `theta`, `log_likelihood`, `steps`, `gradient` (the
#'   largest absolute gradient at the result), and at the starting point
#'   `start_log_likelihood`, `start_gradient` and `start_hessian`.
twolevel_lpa_refine <- function(theta, layout, x, groups, tolerance = 1e-6,
                                max_steps = 50L) {
  objective <- twolevel_lpa_objective(layout, x, groups)
  hessian <- numDeriv::hessian(objective, theta)
  climb <- function(current, log_likelihood, steps) {
    gradient <- numDeriv::grad(objective, current)
    if (max(abs(gradient)) < tolerance) {
      return(list(theta = current, log_likelihood = log_likelihood, steps = steps,
                  gradient = max(abs(gradient))))
    }
    if (steps >= max_steps) {
      stop(sprintf("Newton refinement did not reach a gradient below %g in %d steps",
                   tolerance, max_steps), call. = FALSE)
    }
    proposal <- current - solve(hessian, gradient)
    proposed <- objective(proposal)
    stopifnot("a refining step must not lose likelihood" =
                is.finite(proposed) && proposed >= log_likelihood - 1e-9)
    climb(proposal, proposed, steps + 1L)
  }
  start <- objective(theta)
  out <- climb(theta, start, 0L)
  out$start_log_likelihood <- start
  out$start_gradient <- numDeriv::grad(objective, theta)
  out$start_hessian <- hessian
  out
}

#' One EM step of the two-level Gaussian model, from its definition
#'
#' E step: the joint posterior of row i's group class w and profile c is
#' `P(w | y_g) P(c | y_i, w)`. M step: `delta_w` is the mean group posterior;
#' `gamma_wc` the posterior mass of (w, c) over that of w; each mean and
#' variance the posterior-weighted moment over the rows where the indicator is
#' observed (exact for a diagonal covariance under FIML, where a missing value
#' is independent of the observed ones within a profile); an equal variance
#' pools the profiles; a full covariance is the weighted cross-product, which
#' is used here only on complete data.
#'
#' @param theta The free vector.
#' @param layout,x,groups As for [twolevel_lpa_objective()].
#' @return The free vector after one step.
twolevel_lpa_em_step <- function(theta, layout, x, groups) {
  parameters <- twolevel_lpa_unpack(theta, layout)
  evaluated <- twolevel_lpa_evaluate(parameters, x, groups)
  group_index <- match(as.character(groups), unique(as.character(groups)))
  row_group <- evaluated$group_posteriors[group_index, , drop = FALSE]
  joint <- lapply(seq_len(layout$g), function(w) evaluated$conditional[[w]] * row_group[, w])
  gamma <- t(vapply(seq_len(layout$g), function(w) colSums(joint[[w]]) / sum(row_group[, w]),
                    numeric(layout$k)))
  delta <- colMeans(evaluated$group_posteriors)
  posteriors <- Reduce(`+`, joint)
  observed <- !is.na(x)
  filled <- ifelse(observed, x, 0)
  mass <- crossprod(posteriors, observed)
  means <- crossprod(posteriors, filled) / mass
  covariances <- if (identical(layout$structure, "full")) {
    stopifnot("the full-covariance M step here assumes complete data" = all(observed))
    vapply(seq_len(layout$k), function(k) {
      centred <- sweep(x, 2L, means[k, ])
      crossprod(centred * posteriors[, k], centred) / sum(posteriors[, k])
    }, matrix(0, layout$j, layout$j))
  } else {
    squares <- t(vapply(seq_len(layout$k), function(k) {
      colSums(posteriors[, k] * ifelse(observed, sweep(filled, 2L, means[k, ])^2, 0))
    }, numeric(layout$j)))
    variances <- if (identical(layout$structure, "equal")) {
      matrix(colSums(squares) / colSums(observed), layout$k, layout$j, byrow = TRUE)
    } else squares / mass
    vapply(seq_len(layout$k), function(k) diag(variances[k, ], layout$j),
           matrix(0, layout$j, layout$j))
  }
  twolevel_lpa_pack(list(delta = delta, gamma = matrix(gamma, layout$g),
                         means = matrix(means, layout$k),
                         covariances = array(covariances, c(layout$j, layout$j, layout$k))),
                    layout)
}

#' Run EM from a point to its fixed point
#'
#' Newton steps on a numerical gradient resolve the maximum only to about
#' 1e-8 in the parameters (see [twolevel_lpa_refine()]); EM needs no
#' derivative, so continued to its fixed point it resolves the maximum to
#' rounding error. Steps are taken in blocks of `block` (an iteration count
#' not known in advance, written as a bounded recursion over blocks rather
#' than a loop) until one step moves no parameter by more than `tolerance`.
#'
#' @param theta Starting free vector, near the maximum.
#' @param layout,x,groups As for [twolevel_lpa_objective()].
#' @param tolerance Largest parameter change of the last step accepted.
#' @param max_blocks Bound on the number of blocks.
#' @param block Steps per block.
#' @return A list with `theta`, `log_likelihood`, `steps` and `change` (the
#'   largest parameter change of the last step).
twolevel_lpa_fixed_point <- function(theta, layout, x, groups, tolerance = 1e-13,
                                     max_blocks = 100L, block = 100L) {
  objective <- twolevel_lpa_objective(layout, x, groups)
  run <- function(current, blocks) {
    before <- objective(current)
    advanced <- Reduce(function(point, step) twolevel_lpa_em_step(point, layout, x, groups),
                       seq_len(block - 1L), current)
    last <- twolevel_lpa_em_step(advanced, layout, x, groups)
    reached <- objective(last)
    stopifnot("EM must not lose likelihood" = reached >= before - 1e-9)
    change <- max(abs(last - advanced))
    if (change < tolerance) {
      return(list(theta = last, log_likelihood = reached, steps = (blocks + 1L) * block,
                  change = change))
    }
    if (blocks + 1L >= max_blocks) {
      stop(sprintf("EM did not reach a step below %g in %d steps", tolerance,
                   max_blocks * block), call. = FALSE)
    }
    run(last, blocks + 1L)
  }
  run(theta, 0L)
}

#' Every natural parameter, labelled as the package labels it
#'
#' The free parameters plus the dependent ones (each row's last probability,
#' and for `"full"` nothing more: the upper triangle repeats the lower), with
#' the Jacobian of the map from the free vector. Labels follow
#' `parameter_inference()`: a mean or a variance has `outcome` `profile_k`
#' (`shared` for a variance of the `"equal"` model) and `term` the indicator;
#' a covariance element has `term` `"row:column"`; `P(profile k | group class
#' w)` is level `profile`, outcome `profile_k`, term `group_class_w`; `P(group
#' class w)` is level `group`, outcome `group_class_w`, term `""`. Classes are
#' numbered in the reference's own order.
#'
#' @param theta The free vector.
#' @param layout From [twolevel_lpa_layout()].
#' @return A list with `labels` (a `data.frame` with `level`, `outcome`,
#'   `term`, `parameter`, `estimate`) and `jacobian` (natural by free).
twolevel_lpa_natural <- function(theta, layout) {
  k <- layout$k
  g <- layout$g
  indicators <- layout$indicators
  cells <- layout$cells
  means <- data.frame(level = "measurement",
                      outcome = sprintf("profile_%d", rep(seq_len(k), times = layout$j)),
                      term = rep(indicators, each = k), parameter = "mean",
                      stringsAsFactors = FALSE)
  shared <- identical(layout$structure, "equal")
  full <- identical(layout$structure, "full")
  covariance <- data.frame(
    level = "measurement",
    outcome = if (shared) "shared" else
      sprintf("profile_%d", rep(seq_len(layout$n_blocks), each = nrow(cells))),
    term = if (full) paste(indicators[cells[, 1L]], indicators[cells[, 2L]], sep = ":") else
      indicators[cells[, 1L]],
    parameter = if (full) "covariance" else "variance", stringsAsFactors = FALSE)
  free_gamma <- seq_len(k - 1L)
  gamma <- data.frame(level = "profile", outcome = sprintf("profile_%d", rep(seq_len(k), times = g)),
                      term = sprintf("group_class_%d", rep(seq_len(g), each = k)),
                      parameter = "probability", stringsAsFactors = FALSE)
  delta <- data.frame(level = "group", outcome = sprintf("group_class_%d", seq_len(g)),
                      term = "", parameter = "probability", stringsAsFactors = FALSE)
  labels <- rbind(means, covariance, gamma, delta)
  # The free coordinates map to themselves; each row's last probability is one
  # minus the others.
  n_measurement <- length(layout$at$means) + length(layout$at$covariance)
  jacobian <- matrix(0, nrow(labels), layout$q)
  jacobian[cbind(seq_len(n_measurement), seq_len(n_measurement))] <- 1
  # Row of P(profile c | group class w): one at its own free coordinate for
  # c < K; minus one at each of the group class's free coordinates for c = K.
  gamma_cells <- expand.grid(profile = seq_len(k), group_class = seq_len(g))
  gamma_cells$row <- n_measurement + seq_len(k * g)
  gamma_entries <- do.call(rbind, lapply(seq_len(nrow(gamma_cells)), function(index) {
    profile <- gamma_cells$profile[[index]]
    columns <- layout$at$gamma[(gamma_cells$group_class[[index]] - 1L) * (k - 1L) + free_gamma]
    if (profile < k) cbind(gamma_cells$row[[index]], columns[[profile]], 1) else
      cbind(gamma_cells$row[[index]], columns, -1)
  }))
  jacobian[gamma_entries[, 1:2, drop = FALSE]] <- gamma_entries[, 3L]
  delta_rows <- n_measurement + k * g + seq_len(g)
  delta_free <- layout$at$delta
  jacobian[cbind(delta_rows[-g], delta_free)] <- 1
  jacobian[delta_rows[[g]], delta_free] <- -1
  parameters <- twolevel_lpa_unpack(theta, layout)
  labels$estimate <- c(as.vector(parameters$means),
                       theta[layout$at$covariance],
                       as.vector(t(parameters$gamma)), parameters$delta)
  list(labels = labels, jacobian = jacobian)
}

#' The covariances of the natural parameters, from scratch
#'
#' `observed` is the inverse of minus the Richardson-extrapolated Hessian
#' (numDeriv) of the log-likelihood in the free coordinates; `robust` is the
#' sandwich `A B A` with `B` the cross-product of the per-group scores (the
#' numDeriv Jacobian of the groups' log-likelihood contributions), times the
#' CR1 factor G / (G - 1) for G groups, which `?parameter_inference` documents
#' for `vcov_type = "robust"`. Both are carried to every natural parameter by
#' the linear map of [twolevel_lpa_natural()].
#'
#' @param theta The free vector at a maximum.
#' @param layout,x,groups As for [twolevel_lpa_objective()].
#' @return A list with `labels` (from [twolevel_lpa_natural()]), `observed`
#'   and `robust` (natural-scale covariance matrices).
twolevel_lpa_covariances <- function(theta, layout, x, groups) {
  hessian <- numDeriv::hessian(twolevel_lpa_objective(layout, x, groups), theta)
  scores <- numDeriv::jacobian(twolevel_lpa_objective(layout, x, groups, "groups"), theta)
  bread <- solve(-hessian)
  n_groups <- nrow(scores)
  natural <- twolevel_lpa_natural(theta, layout)
  carry <- function(covariance) natural$jacobian %*% covariance %*% t(natural$jacobian)
  list(labels = natural$labels, observed = carry(bread),
       robust = carry(bread %*% crossprod(scores) %*% bread * n_groups / (n_groups - 1)))
}

#' The from-scratch solution in the harness's reference shape
#'
#' @param theta The free vector.
#' @param layout,x,groups As for [twolevel_lpa_objective()].
#' @param data The indicator data frame, in the fit's row order.
#' @param covariances From [twolevel_lpa_covariances()], or `NULL` for no
#'   standard errors.
#' @return A `harness_reference`, in the reference's own class order.
twolevel_lpa_reference <- function(theta, layout, x, groups, data, covariances = NULL) {
  parameters <- twolevel_lpa_unpack(theta, layout)
  evaluated <- twolevel_lpa_evaluate(parameters, x, groups)
  k <- layout$k
  indicators <- layout$indicators
  variances <- t(vapply(seq_len(k), function(profile) diag(parameters$covariances[, , profile]),
                        numeric(layout$j)))
  variances <- matrix(variances, k, layout$j, dimnames = list(NULL, indicators))
  full <- identical(layout$structure, "full")
  standard_errors <- if (!is.null(covariances)) {
    labels <- covariances$labels
    errors <- sqrt(diag(covariances$observed))
    measurement <- labels$level == "measurement"
    grid <- expand.grid(profile = seq_len(k), indicator = indicators, stringsAsFactors = FALSE)
    variance_label <- if (full) paste(grid$indicator, grid$indicator, sep = ":") else grid$indicator
    variance_outcome <- if (identical(layout$structure, "equal")) "shared" else
      sprintf("profile_%d", grid$profile)
    key <- paste(labels$outcome, labels$term, labels$parameter)
    lookup <- function(outcome, term, parameter) {
      errors[match(paste(outcome, term, parameter), key)]
    }
    list(mean = matrix(lookup(sprintf("profile_%d", grid$profile), grid$indicator, "mean"),
                       k, layout$j, dimnames = list(NULL, indicators)),
         variance = matrix(lookup(variance_outcome, variance_label,
                                  if (full) "covariance" else "variance"),
                           k, layout$j, dimnames = list(NULL, indicators)))
  }
  harness_reference(
    log_likelihood = evaluated$log_likelihood, n_parameters = layout$q,
    posteriors = evaluated$posteriors, groups = groups,
    means = parameters$means,
    variances = if (!full) variances,
    covariances = if (full) parameters$covariances,
    standard_errors = standard_errors,
    profile_probabilities = parameters$gamma, group_probabilities = parameters$delta,
    group_posteriors = evaluated$group_posteriors,
    group_log_likelihood = evaluated$group_log_likelihood,
    data = data, row_conditional_entropy = evaluated$row_conditional_entropy)
}

#' The `parameter_inference()` table a reference implies
#'
#' Written from `?parameter_inference`: the statistic is the Wald ratio
#' `estimate / standard_error`, except for a parameter whose null value lies on
#' the boundary of its space (a variance, a covariance diagonal, a
#' probability), which has none; `p_value` is two-sided normal and, with the
#' default `adjust = "none"`, `p_adjusted` equals it. Intervals for
#' probabilities are formed on the logit scale and for variances and
#' covariance diagonals on the log scale, by the delta method, and
#' back-transformed; every other interval is `estimate +/- z se`.
#'
#' @param labels Labels and estimates, from [twolevel_lpa_natural()], already
#'   in the fit's class numbering.
#' @param errors Standard errors, one per row of `labels`.
#' @param level Confidence level.
#' @return A `data.frame` with the columns of `parameter_inference()`.
twolevel_lpa_inference_table <- function(labels, errors, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  estimate <- labels$estimate
  diagonal <- labels$parameter == "covariance" &
    vapply(strsplit(labels$term, ":", fixed = TRUE), function(pair) {
      length(pair) == 2L && identical(pair[[1L]], pair[[2L]])
    }, logical(1))
  logged <- labels$parameter == "variance" | diagonal
  logit <- labels$parameter == "probability"
  tested <- !(logged | logit)
  statistic <- ifelse(tested, estimate / errors, NA_real_)
  p_value <- 2 * stats::pnorm(-abs(statistic))
  # Each interval on its own scale; a transform is applied only to the
  # parameters whose range it maps (the logit of a mean is undefined).
  low <- estimate - z * errors
  high <- estimate + z * errors
  half_log <- z * errors[logged] / estimate[logged]
  low[logged] <- estimate[logged] * exp(-half_log)
  high[logged] <- estimate[logged] * exp(half_log)
  centre <- stats::qlogis(estimate[logit])
  half_logit <- z * errors[logit] / (estimate[logit] * (1 - estimate[logit]))
  low[logit] <- stats::plogis(centre - half_logit)
  high[logit] <- stats::plogis(centre + half_logit)
  data.frame(labels[c("level", "outcome", "term", "parameter")], estimate = estimate,
             standard_error = errors, statistic = statistic, p_value = p_value,
             p_adjusted = p_value, conf_low = low, conf_high = high,
             stringsAsFactors = FALSE)
}

#' Renumber the classes named in parameter labels
#'
#' @param labels A `data.frame` with `outcome` and `term`, in the reference's
#'   numbering.
#' @param profiles,group_classes The fit's class for each reference class.
#' @return `labels` in the fit's numbering.
twolevel_lpa_relabel <- function(labels, profiles, group_classes) {
  renumber <- function(text, prefix, map) {
    pattern <- sprintf("^%s_([0-9]+)$", prefix)
    hit <- grepl(pattern, text)
    text[hit] <- sprintf("%s_%d", prefix, map[as.integer(sub(pattern, "\\1", text[hit]))])
    text
  }
  labels$outcome <- renumber(renumber(labels$outcome, "profile", profiles),
                             "group_class", group_classes)
  labels$term <- renumber(labels$term, "group_class", group_classes)
  labels
}
