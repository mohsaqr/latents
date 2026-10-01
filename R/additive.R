# Internal engine of the group-class families of multilpa(family = ):
# additive, dispersion and additive-dispersion (Houle, Morin & Harvey, 2026).
# Groups belong to latent group classes and carry a Gaussian intercept per
# indicator; the families differ in which blocks vary across classes. See
# validation/ADDITIVE_MODEL_DESIGN.md for the specification.

#' Reduce raw ratings to the additive model's sufficient statistics
#'
#' Validates the raw-rating contract of the first increment (complete numeric
#' indicators, a known group identifier, replication within at least one
#' group) and returns, per group and indicator, the mean and the scatter about
#' it. The exact likelihood depends on the data only through these, the group
#' sizes and nothing else.
#'
#' @param data A data frame, one row per observation.
#' @param vars Names of the continuous indicator columns.
#' @param id Name of the group identifier column.
#' @return A list of class `latents_additive_data` with `ids` (group labels in
#'   first-occurrence order), `index` (group of each row), `sizes` (n_j),
#'   `n_groups` (J), `n_obs` (N), `vars`, `averages` and `scatter` (J x d
#'   matrices, rows named by group) and `pooled_scatter` (sum over groups of
#'   the scatter, one value per indicator).
#' @section Conditions:
#' `latents_bad_data` for a malformed frame, a non-numeric, missing,
#' non-finite or constant indicator, or a missing group identifier.
#' `latents_unidentified` when every group is a singleton, because the within
#' and between variances are then observed only through their sum.
#' @noRd
.additive_prepare <- function(data, vars, id, weights = NULL) {
  if (!is.data.frame(data) || !is.character(vars) || !is.character(id) ||
      nrow(data) < 2L || length(vars) < 1L || anyNA(vars) ||
      anyDuplicated(vars) || anyDuplicated(names(data)) ||
      length(id) != 1L || is.na(id) || id %in% vars ||
      !all(c(vars, id) %in% names(data))) {
    stop(errorCondition(
      "Supply a data frame with at least two rows, unique existing indicators, and one distinct group column.",
      class = "latents_bad_data", call = NULL))
  }
  frame <- data[, vars, drop = FALSE]
  if (!all(vapply(frame, is.numeric, logical(1))) ||
      any(vapply(frame, function(value) !is.null(dim(value)), logical(1)))) {
    stop(errorCondition(
      "Every indicator must be a numeric vector; the additive model takes continuous ratings.",
      class = "latents_bad_data", call = NULL))
  }
  x <- matrix(as.numeric(as.matrix(frame)), nrow = nrow(data),
              ncol = length(vars), dimnames = list(NULL, vars))
  if (!all(is.finite(x))) {
    stop(errorCondition(
      "Indicators contain missing or non-finite values; the additive model requires complete ratings.",
      class = "latents_bad_data", call = NULL))
  }
  if (any(vapply(seq_along(vars), function(r) {
    length(unique(x[, r])) < 2L
  }, logical(1)))) {
    stop(errorCondition("Constant indicators cannot identify group intercepts.",
                        class = "latents_bad_data", call = NULL))
  }
  groups <- .multilpa_prepare_groups(data[[id]])
  if (all(groups$sizes == 1L)) {
    stop(errorCondition(paste(
      "Every group has one observation, so within- and between-group variances",
      "are observed only through their sum. Supply groups with repeated ratings."),
      class = "latents_unidentified", call = NULL))
  }
  sizes <- groups$sizes
  averages <- rowsum(x, groups$index, reorder = TRUE) / sizes
  scatter <- rowsum((x - averages[groups$index, , drop = FALSE])^2,
                    groups$index, reorder = TRUE)
  dimnames(averages) <- dimnames(scatter) <- list(groups$ids, vars)
  stopifnot(
    "one sufficient-statistic row per group" =
      nrow(averages) == groups$n && nrow(scatter) == groups$n,
    "scatter must be nonnegative and finite" =
      all(is.finite(scatter)) && all(scatter >= 0),
    "group sizes must sum to the number of rows" = sum(sizes) == nrow(data)
  )
  # Sampling weights, one per group, travel with the statistics: every
  # posterior-weighted sum the M-step and the scores form is weighted by them.
  sampling_weights <- .latents_sampling_weights(data, weights, groups$index, groups$n)
  structure(list(
    ids = groups$ids, index = groups$index, sizes = sizes,
    n_groups = groups$n, n_obs = nrow(data), vars = vars,
    averages = averages, scatter = scatter,
    pooled_scatter = colSums(scatter),
    sampling_weights = sampling_weights,
    weighted_n_obs = if (is.null(sampling_weights)) NULL else
      sum(sampling_weights * sizes)
  ), class = "latents_additive_data")
}

#' Which parameter blocks are class-specific in each group-class family
#'
#' Every group-class family shares one engine; they differ only in which
#' blocks vary across group classes (Houle, Morin & Harvey, 2026):
#' additive (class means, shared within variance), dispersion (class within
#' variances, shared means and between variances) and additive-dispersion
#' (class means and class within variances).
#' @param family `"additive"`, `"dispersion"` or `"additive_dispersion"`.
#' @param between_variance `"varying"` or `"equal"`; fixed at `"equal"` for
#'   the dispersion family.
#' @return A list with `means`, `between` and `within`, each `"varying"` or
#'   `"equal"`.
#' @noRd
.additive_structure <- function(family, between_variance = "varying") {
  switch(family,
    additive = list(means = "varying", between = between_variance,
                    within = "equal"),
    dispersion = list(means = "equal", between = "equal", within = "varying"),
    additive_dispersion = list(means = "varying", between = between_variance,
                               within = "varying"),
    stop("unknown group-class family"))
}

#' Check a set of group-class parameters and bring `within` to H x d
#' @param parameters List with `means` and `between` (H x d matrices),
#'   `within` (length d, shared by every class, or an H x d matrix) and
#'   `weights` (length H).
#' @param n_indicators The number of indicators d.
#' @return `parameters` with `within` as an H x d matrix; raises
#'   `latents_bad_start` otherwise.
#' @noRd
.additive_check_parameters <- function(parameters, n_indicators) {
  valid <- is.list(parameters) &&
    all(c("means", "between", "within", "weights") %in% names(parameters)) &&
    is.matrix(parameters$means) && is.numeric(parameters$means) &&
    is.matrix(parameters$between) && is.numeric(parameters$between) &&
    identical(dim(parameters$means), dim(parameters$between)) &&
    ncol(parameters$means) == n_indicators && nrow(parameters$means) >= 1L &&
    is.numeric(parameters$within) &&
    (identical(dim(parameters$within), dim(parameters$means)) ||
       (is.null(dim(parameters$within)) &&
          length(parameters$within) == n_indicators)) &&
    is.numeric(parameters$weights) &&
    length(parameters$weights) == nrow(parameters$means) &&
    all(is.finite(parameters$means)) && all(is.finite(parameters$between)) &&
    all(parameters$between >= 0) && all(is.finite(parameters$within)) &&
    all(parameters$within > 0) && all(is.finite(parameters$weights)) &&
    all(parameters$weights > 0) &&
    abs(sum(parameters$weights) - 1) < sqrt(.Machine$double.eps)
  if (!isTRUE(valid)) {
    stop(errorCondition(paste(
      "Group-class parameters need H x d `means` and nonnegative `between`,",
      "positive `within` (length d or H x d), and positive `weights` summing",
      "to one."), class = "latents_bad_start", call = NULL))
  }
  if (is.null(dim(parameters$within))) {
    parameters$within <- matrix(parameters$within, nrow(parameters$means),
                                n_indicators, byrow = TRUE)
  }
  parameters
}

#' Row h of an H x d parameter matrix repeated over the J groups
#' @noRd
.additive_row <- function(values, h, n_groups) {
  matrix(values[h, ], n_groups, ncol(values), byrow = TRUE)
}

#' Component log densities of every group under every group class
#'
#' The raw observations of group j on indicator r have covariance
#' `W_hr I + T_hr 11'`, whose determinant is `W_hr^(n-1) (W_hr + n T_hr)` and
#' whose quadratic form splits into the within scatter over `W_hr` and the mean
#' deviation over `(W_hr + n T_hr) / n`. Every normalizing constant of the raw
#' observations is kept, so the value is comparable with any raw-observation
#' likelihood.
#'
#' @param stats Result of `.additive_prepare()`.
#' @param parameters Group-class parameters (see `.additive_check_parameters()`).
#' @return A J x H matrix of log densities, rows named by group.
#' @noRd
.additive_log_density <- function(stats, parameters) {
  stopifnot("`stats` must come from .additive_prepare()" =
              inherits(stats, "latents_additive_data"))
  parameters <- .additive_check_parameters(parameters, length(stats$vars))
  sizes <- stats$sizes
  j_count <- stats$n_groups
  densities <- vapply(seq_len(nrow(parameters$means)), function(h) {
    within <- .additive_row(parameters$within, h, j_count)
    denominator <- within + sizes * .additive_row(parameters$between, h, j_count)
    -0.5 * rowSums(sizes * log(2 * pi) + (sizes - 1) * log(within) +
                     stats$scatter / within + log(denominator) +
                     sizes * (stats$averages -
                                .additive_row(parameters$means, h, j_count))^2 /
                       denominator)
  }, numeric(j_count))
  matrix(densities, nrow = j_count,
         dimnames = list(stats$ids, rownames(parameters$means)))
}

#' Group posteriors and group-intercept moments
#'
#' The E-step. For each group class the posterior of a group's intercept is
#' Gaussian with mean `mu + n T / (W + n T) (a - mu)` and variance
#' `T W / (W + n T)`; the marginal moments mix these over the group's class
#' posterior, so they include uncertainty about class membership. All are
#' conditional on the supplied parameters.
#'
#' @param stats Result of `.additive_prepare()`.
#' @param parameters Group-class parameters.
#' @return A list with `log_likelihood` (scalar), `group_log_likelihood`
#'   (length J), `log_density` (J x H), `posterior` (J x H, rows sum to one),
#'   `conditional_mean` and `conditional_variance` (J x d x H arrays),
#'   `intercept_mean` and `intercept_variance` (J x d matrices).
#' @noRd
.additive_expectation <- function(stats, parameters) {
  parameters <- .additive_check_parameters(parameters, length(stats$vars))
  log_density <- .additive_log_density(stats, parameters)
  log_joint <- sweep(log_density, 2L, log(parameters$weights), "+")
  group_log_likelihood <- .multilpa_log_sum_exp(log_joint)
  posterior <- exp(log_joint - group_log_likelihood)
  n_classes <- nrow(parameters$means)
  j_count <- stats$n_groups
  dims <- c(j_count, length(stats$vars), n_classes)
  moments <- lapply(seq_len(n_classes), function(h) {
    within <- .additive_row(parameters$within, h, j_count)
    between <- .additive_row(parameters$between, h, j_count)
    means <- .additive_row(parameters$means, h, j_count)
    denominator <- within + stats$sizes * between
    list(mean = means + stats$sizes * between / denominator *
           (stats$averages - means),
         variance = between * within / denominator)
  })
  labels <- list(stats$ids, stats$vars, rownames(parameters$means))
  conditional_mean <- array(unlist(lapply(moments, `[[`, "mean")), dims,
                            dimnames = labels)
  conditional_variance <- array(unlist(lapply(moments, `[[`, "variance")), dims,
                                dimnames = labels)
  # weight[j, r, h] is the posterior of class h for group j, for every r.
  weight <- array(posterior[, rep(seq_len(n_classes),
                                  each = length(stats$vars))], dims)
  intercept_mean <- rowSums(weight * conditional_mean, dims = 2L)
  intercept_variance <- rowSums(weight * (conditional_variance +
    (conditional_mean - as.vector(intercept_mean))^2), dims = 2L)
  dimnames(intercept_mean) <- dimnames(intercept_variance) <- labels[1:2]
  stopifnot(
    "posterior rows must sum to one" =
      all(abs(rowSums(posterior) - 1) < sqrt(.Machine$double.eps)),
    "intercept variances must be nonnegative" =
      all(intercept_variance >= -sqrt(.Machine$double.eps))
  )
  list(log_likelihood = sum((stats$sampling_weights %||% 1) * group_log_likelihood),
       group_log_likelihood = stats::setNames(group_log_likelihood, stats$ids),
       log_density = log_density, posterior = posterior,
       conditional_mean = conditional_mean,
       conditional_variance = conditional_variance,
       intercept_mean = intercept_mean, intercept_variance = intercept_variance)
}

#' One group class's J x d slice of a J x d x H array, keeping both dimensions
#' @noRd
.additive_slice <- function(values, h) {
  matrix(values[, , h], nrow = dim(values)[1L], ncol = dim(values)[2L])
}

#' Derivative of the log likelihood with respect to each between variance
#'
#' By Fisher's identity the derivative of the marginal log likelihood is the
#' posterior-weighted derivative of each component density:
#' `0.5 n (n (a - mu)^2 / D^2 - 1 / D)` with `D = W + n T`. At `T = 0` its sign
#' is the Karush-Kuhn-Tucker test for a boundary maximum.
#'
#' @param stats Result of `.additive_prepare()`.
#' @param parameters Group-class parameters.
#' @param posterior J x H group posteriors at `parameters`.
#' @param between_variance `"varying"` gives one derivative per class and
#'   indicator; `"equal"` sums them over classes (the shared coordinate).
#' @return An H x d matrix (`"varying"`) or a length-d vector (`"equal"`).
#' @noRd
.additive_between_score <- function(stats, parameters, posterior,
                                    between_variance = c("varying", "equal")) {
  between_variance <- match.arg(between_variance)
  parameters <- .additive_check_parameters(parameters, length(stats$vars))
  posterior <- posterior * (stats$sampling_weights %||% 1)
  n_classes <- nrow(parameters$means)
  j_count <- stats$n_groups
  score <- matrix(vapply(seq_len(n_classes), function(h) {
    denominator <- .additive_row(parameters$within, h, j_count) +
      stats$sizes * .additive_row(parameters$between, h, j_count)
    deviation <- stats$averages - .additive_row(parameters$means, h, j_count)
    colSums(posterior[, h] * 0.5 * stats$sizes *
              (stats$sizes * deviation^2 / denominator^2 - 1 / denominator))
  }, numeric(length(stats$vars))), nrow = n_classes, byrow = TRUE)
  if (identical(between_variance, "equal")) colSums(score) else score
}

#' Pool an H x d matrix of class values into one shared row when a block is
#' shared, with class weights `weight` (H x d)
#' @noRd
.additive_pool <- function(values, weight, mode) {
  if (identical(mode, "varying")) return(values)
  pooled <- colSums(values * weight) / colSums(weight)
  matrix(pooled, nrow(values), ncol(values), byrow = TRUE)
}

#' Starting values from a (soft) partition of the groups
#'
#' Start 1 partitions the groups by Ward's method on the features in which
#' the family's classes differ (standardized group means when class means
#' vary, standardized log within-group variances when within variances vary),
#' which is deterministic; later starts assign every group to the nearest of H
#' randomly chosen groups on the same features. The partition is smoothed so no
#' class starts empty. Between variances start strictly positive: EM cannot
#' leave a zero between variance.
#'
#' @param structure Result of `.additive_structure()`.
#' @return Group-class parameters with H x d `within`.
#' @noRd
.additive_initialize <- function(stats, n_classes, structure, min_variance,
                                 start_index) {
  standardize <- function(values) {
    spread <- apply(values, 2L, stats::sd)
    spread[!is.finite(spread) | spread <= 0] <- 1
    sweep(sweep(values, 2L, colMeans(values)), 2L, spread, "/")
  }
  # Log within-group variance per group; a singleton group, which has none,
  # takes the indicator's average over the replicated groups.
  replicated <- stats$sizes > 1L
  log_spread <- log(pmax(stats$scatter / pmax(stats$sizes - 1, 1),
                         min_variance))
  if (any(!replicated)) {
    fill <- colMeans(log_spread[replicated, , drop = FALSE])
    log_spread[!replicated, ] <- matrix(fill, sum(!replicated), ncol(log_spread),
                                        byrow = TRUE)
  }
  features <- cbind(
    if (identical(structure$means, "varying")) standardize(stats$averages),
    if (identical(structure$within, "varying")) standardize(log_spread))
  hard <- if (n_classes == 1L) {
    rep(1L, stats$n_groups)
  } else if (start_index == 1L) {
    .multilpa_ward_assignments(features, n_classes)
  } else {
    centres <- features[sample.int(stats$n_groups, n_classes), , drop = FALSE]
    distance <- vapply(seq_len(n_classes), function(h) {
      colSums((t(features) - centres[h, ])^2)
    }, numeric(stats$n_groups))
    max.col(-matrix(distance, stats$n_groups), ties.method = "first")
  }
  posterior <- 0.9 * outer(hard, seq_len(n_classes), "==") + 0.1 / n_classes
  counts <- colSums(posterior)
  n_indicators <- length(stats$vars)
  class_weight <- matrix(counts, n_classes, n_indicators)
  means <- .additive_pool(crossprod(posterior, stats$averages) / counts,
                          class_weight, structure$means)
  pooled_within <- pmax(stats$pooled_scatter / (stats$n_obs - stats$n_groups),
                        min_variance)
  within <- if (identical(structure$within, "varying")) {
    degrees <- colSums(posterior * pmax(stats$sizes - 1, 0))
    class_within <- crossprod(posterior, stats$scatter) / pmax(degrees, 1)
    usable <- degrees > 0
    if (any(!usable)) {
      class_within[!usable, ] <- matrix(pooled_within, sum(!usable),
                                        n_indicators, byrow = TRUE)
    }
    pmax(class_within, min_variance)
  } else {
    matrix(pooled_within, n_classes, n_indicators, byrow = TRUE)
  }
  floor <- pmax(0.01 * apply(stats$averages, 2L, stats::var), min_variance)
  floor[!is.finite(floor)] <- min_variance
  raw <- matrix(vapply(seq_len(n_classes), function(h) {
    mean_size <- sum(posterior[, h] * stats$sizes) / counts[h]
    colSums(posterior[, h] * sweep(stats$averages, 2L, means[h, ])^2) /
      counts[h] - within[h, ] / mean_size
  }, numeric(n_indicators)), nrow = n_classes, byrow = TRUE)
  between <- pmax(.additive_pool(raw, class_weight, structure$between),
                  matrix(floor, n_classes, n_indicators, byrow = TRUE))
  list(means = unname(means), between = unname(between),
       within = unname(within), weights = unname(counts / stats$n_groups))
}

#' The M-step of the group-class families
#'
#' Coordinates with a positive between variance use the latent-intercept
#' augmentation; coordinates held at zero (`zero`) have no latent intercept,
#' so their complete data are the ratings themselves. Class means: the
#' posterior average of the conditional intercept means (latent coordinates)
#' or the size-weighted average of group means (held coordinates). When the
#' means are shared across classes, the class contributions are pooled with
#' their precisions (`1 / T` for latent, `n / W` for held coordinates), using
#' the current within variances. Then between variances (spread of the
#' intercepts about the new means, per class or pooled over groups) and within
#' variances (expected residual sum of squares, per class over that class's
#' observations or pooled over all observations). Each update maximizes the
#' expected complete-data log likelihood given the ones before it, so the step
#' never lowers the likelihood (an ECM step when shared means meet
#' class-specific within variances, an exact M-step otherwise).
#'
#' @param stats Result of `.additive_prepare()`.
#' @param expectation Result of `.additive_expectation()` at `parameters`.
#' @param parameters The current group-class parameters.
#' @param structure Result of `.additive_structure()`.
#' @param min_variance Lower bound on the within variances.
#' @param zero H x d logical matrix of between variances held at zero.
#' @return Updated group-class parameters.
#' @noRd
.additive_maximization <- function(stats, expectation, parameters, structure,
                                   min_variance, zero) {
  parameters <- .additive_check_parameters(parameters, length(stats$vars))
  # Weighted posteriors are weighted counts; the weights sum to the groups.
  posterior <- expectation$posterior * (stats$sampling_weights %||% 1)
  n_classes <- ncol(posterior)
  n_indicators <- length(stats$vars)
  counts <- colSums(posterior)
  # Per class: the numerator and the precision of the class-mean update.
  pieces <- lapply(seq_len(n_classes), function(h) {
    latent_sum <- colSums(posterior[, h] *
                            .additive_slice(expectation$conditional_mean, h))
    direct_weight <- sum(posterior[, h] * stats$sizes)
    direct_sum <- colSums(posterior[, h] * stats$sizes * stats$averages)
    held <- zero[h, ]
    list(value = ifelse(held, direct_sum / direct_weight, latent_sum / counts[h]),
         precision = ifelse(held, direct_weight / parameters$within[h, ],
                            counts[h] / pmax(parameters$between[h, ],
                                             .Machine$double.xmin)))
  })
  class_means <- matrix(vapply(pieces, `[[`, numeric(n_indicators), "value"),
                        nrow = n_classes, byrow = TRUE)
  precision <- matrix(vapply(pieces, `[[`, numeric(n_indicators), "precision"),
                      nrow = n_classes, byrow = TRUE)
  means <- .additive_pool(class_means, precision, structure$means)
  spread <- matrix(vapply(seq_len(n_classes), function(h) {
    conditional <- .additive_slice(expectation$conditional_mean, h)
    value <- colSums(posterior[, h] * (
      .additive_slice(expectation$conditional_variance, h) +
        sweep(conditional, 2L, means[h, ])^2))
    ifelse(zero[h, ], 0, value)
  }, numeric(n_indicators)), nrow = n_classes, byrow = TRUE)
  between <- if (identical(structure$between, "equal")) {
    matrix(colSums(spread) / stats$n_groups, n_classes, n_indicators,
           byrow = TRUE)
  } else spread / counts
  residual <- matrix(vapply(seq_len(n_classes), function(h) {
    conditional <- .additive_slice(expectation$conditional_mean, h)
    latent <- stats$scatter + stats$sizes *
      ((stats$averages - conditional)^2 +
         .additive_slice(expectation$conditional_variance, h))
    direct <- stats$scatter + stats$sizes *
      sweep(stats$averages, 2L, means[h, ])^2
    held <- matrix(zero[h, ], stats$n_groups, n_indicators, byrow = TRUE)
    colSums(posterior[, h] * ifelse(held, direct, latent))
  }, numeric(n_indicators)), nrow = n_classes, byrow = TRUE)
  observations <- colSums(posterior * stats$sizes)
  within <- if (identical(structure$within, "varying")) {
    residual / observations
  } else {
    matrix(colSums(residual) / (stats$weighted_n_obs %||% stats$n_obs),
           n_classes, n_indicators, byrow = TRUE)
  }
  list(means = means, between = between, within = pmax(within, min_variance),
       weights = counts / stats$n_groups)
}

#' Plain EM for a group-class family from one starting point
#'
#' @param zero H x d logical matrix of between variances held at zero; they
#'   must be zero in `parameters`.
#' @return A list with `parameters`, `expectation`, `converged`, `iterations`
#'   and `history` (log likelihood after each step, the start first).
#' @noRd
.additive_em <- function(stats, parameters, structure, min_variance,
                         max_iter, tol, zero) {
  parameters <- .additive_check_parameters(parameters, length(stats$vars))
  stopifnot("held between variances must be zero" =
              all(parameters$between[zero] == 0))
  expectation <- .additive_expectation(stats, parameters)
  history <- expectation$log_likelihood
  converged <- FALSE
  iteration <- 0L
  # Sequential by nature: each step starts from the last point.
  while (iteration < max_iter && !converged) {
    updated <- .additive_maximization(stats, expectation, parameters, structure,
                                      min_variance, zero)
    updated_expectation <- .additive_expectation(stats, updated)
    gain <- updated_expectation$log_likelihood - expectation$log_likelihood
    if (gain < -1e-10 * (1 + abs(expectation$log_likelihood))) {
      stop("EM likelihood decreased beyond numerical roundoff.")
    }
    converged <- abs(gain) <= tol * (1 + abs(expectation$log_likelihood))
    iteration <- iteration + 1L
    history <- c(history, updated_expectation$log_likelihood)
    parameters <- updated
    expectation <- updated_expectation
  }
  list(parameters = parameters, expectation = expectation,
       converged = converged, iterations = iteration, history = history)
}

#' Check and, where needed, reach a zero-between-variance maximum
#'
#' EM approaches a zero between variance only sublinearly and cannot reach it.
#' Each coordinate whose likelihood derivative at zero (others held at the EM
#' solution) is not positive is a boundary candidate. The candidates are held
#' at zero and the remaining parameters re-maximized; any held coordinate whose
#' derivative at the new point turns positive is released and the fit repeated.
#' The boundary solution replaces the interior one only when its likelihood is
#' at least as high (to the convergence tolerance) and the Karush-Kuhn-Tucker
#' condition holds at every held coordinate.
#'
#' @param fit Result of `.additive_em()` with no coordinate held.
#' @return `fit` or the boundary fit, with `zero` (H x d logical), `kkt`
#'   (logical: every held derivative is nonpositive) and total `iterations`.
#' @noRd
.additive_polish <- function(stats, fit, structure, min_variance, max_iter, tol) {
  n_classes <- nrow(fit$parameters$means)
  n_indicators <- length(stats$vars)
  none <- matrix(FALSE, n_classes, n_indicators)
  interior <- c(fit, list(zero = none, kkt = TRUE))
  # One coordinate per class and indicator ("varying") or per indicator
  # ("equal"): zero it alone and read the derivative there.
  coordinates <- if (identical(structure$between, "equal")) {
    lapply(seq_len(n_indicators), function(r) {
      mask <- none
      mask[, r] <- TRUE
      mask
    })
  } else {
    lapply(seq_len(n_classes * n_indicators), function(k) {
      mask <- none
      mask[k] <- TRUE
      mask
    })
  }
  score_at_zero <- vapply(coordinates, function(mask) {
    point <- fit$parameters
    point$between[mask] <- 0
    posterior <- .additive_expectation(stats, point)$posterior
    sum(.additive_between_score(stats, point, posterior, "varying")[mask])
  }, numeric(1))
  candidates <- Reduce(`|`, coordinates[score_at_zero <= 0], none)
  if (!any(candidates)) return(interior)
  slack <- sqrt(.Machine$double.eps) * (1 + abs(fit$expectation$log_likelihood))
  refit <- function(mask, point) {
    point$between[mask] <- 0
    c(.additive_em(stats, point, structure, min_variance, max_iter, tol, mask),
      list(zero = mask))
  }
  boundary <- refit(candidates, fit$parameters)
  rounds <- 1L
  # Releasing coordinates shrinks the held set, so this ends within
  # length(coordinates) rounds.
  repeat {
    score <- .additive_between_score(stats, boundary$parameters,
                                     boundary$expectation$posterior, "varying")
    held_scores <- vapply(coordinates, function(mask) {
      if (all(boundary$zero[mask])) sum(score[mask]) else -Inf
    }, numeric(1))
    release <- held_scores > slack
    if (!any(release) || rounds > length(coordinates)) break
    mask <- boundary$zero & !Reduce(`|`, coordinates[release], none)
    point <- boundary$parameters
    released <- boundary$zero & !mask
    point$between[released] <- pmax(0.01 * point$within, min_variance)[released]
    boundary <- refit(mask, point)
    rounds <- rounds + 1L
  }
  boundary$kkt <- !any(release)
  boundary$iterations <- fit$iterations + boundary$iterations
  boundary$history <- c(fit$history, boundary$history)
  better <- boundary$expectation$log_likelihood >=
    fit$expectation$log_likelihood - tol * (1 + abs(fit$expectation$log_likelihood))
  if (isTRUE(boundary$kkt) && better && any(boundary$zero)) boundary else interior
}

#' Arguments `multilpa()` accepts for a group-class family
#' @noRd
.additive_arguments <- function() {
  c("data", "vars", "id", "n_group_classes", "family", "between_variance",
    "n_starts", "max_iter", "tol", "min_variance", "seed", "weights")
}

#' Number of free rows (1 shared, or H) in a parameter block
#' @noRd
.additive_rows <- function(mode, n_classes) {
  if (identical(mode, "equal")) 1L else n_classes
}

#' Nominal free parameters of a group-class family
#'
#' Means, between variances and within variances, each one set per class or
#' one shared set, and H - 1 class weights. A variance estimated at zero stays
#' counted: it was estimated, not fixed before fitting.
#' @noRd
.additive_count_parameters <- function(n_classes, n_indicators, structure) {
  as.integer(n_indicators * (.additive_rows(structure$means, n_classes) +
                               .additive_rows(structure$between, n_classes) +
                               .additive_rows(structure$within, n_classes)) +
               n_classes - 1L)
}

#' Readable name of a group-class family
#' @noRd
.additive_family_label <- function(family) {
  switch(family, additive = "Additive", dispersion = "Dispersion",
         additive_dispersion = "Additive-dispersion")
}

#' Fit a group-class family
#'
#' Called by `multilpa(family = )` after its own argument checks. Runs
#' `n_starts` EM starts, polishes each to a Karush-Kuhn-Tucker boundary
#' maximum where a between variance goes to zero, keeps every start's record
#' (failures included), finishes an interior best start by Newton steps, and
#' returns the highest likelihood.
#'
#' @param family `"additive"`, `"dispersion"` or `"additive_dispersion"`.
#' @param between_variance `"varying"`, `"equal"`, or `NULL` for the family's
#'   default.
#' @return An object of class `multilpa_additive`.
#' @noRd
.additive_fit <- function(data, vars, id, n_group_classes, between_variance,
                          n_starts, max_iter, tol, min_variance, seed, call,
                          family = "additive", weights = NULL) {
  if (identical(family, "dispersion") && identical(between_variance, "varying")) {
    stop(errorCondition(paste(
      "The dispersion family holds between-group means and variances equal",
      "across group classes; `between_variance = \"varying\"` belongs to",
      "`family = \"additive\"` or `\"additive_dispersion\"`."),
      class = "latents_bad_argument", call = NULL))
  }
  between_variance <- between_variance %||%
    if (identical(family, "dispersion")) "equal" else "varying"
  structure <- .additive_structure(family, between_variance)
  counts <- list(n_group_classes = n_group_classes, n_starts = n_starts)
  invisible(lapply(names(counts), function(field) {
    value <- counts[[field]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value != floor(value) || value > .Machine$integer.max) {
      stop(errorCondition(sprintf("%s must be a positive integer.", field),
                          class = "latents_bad_argument", call = NULL))
    }
  }))
  if (!is.numeric(max_iter) || length(max_iter) != 1L || !is.finite(max_iter) ||
      max_iter < 0 || max_iter != floor(max_iter) ||
      max_iter > .Machine$integer.max) {
    stop(errorCondition("max_iter must be a nonnegative integer.",
                        class = "latents_bad_argument", call = NULL))
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0 ||
      !is.numeric(min_variance) || length(min_variance) != 1L ||
      !is.finite(min_variance) || min_variance <= 0) {
    stop(errorCondition("tol and min_variance must be finite positive numbers.",
                        class = "latents_bad_argument", call = NULL))
  }
  .multilpa_check_seed(seed)
  n_group_classes <- as.integer(n_group_classes)
  n_starts <- as.integer(n_starts)
  max_iter <- as.integer(max_iter)
  stats <- .additive_prepare(data, vars, id, weights)
  if (n_group_classes > stats$n_groups) {
    stop(errorCondition(sprintf(
      "%d group classes cannot be estimated from %d groups.",
      n_group_classes, stats$n_groups),
      class = "latents_unidentified", call = NULL))
  }
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)  # nolint: object_name_linter. R's name for the RNG state.
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  none <- matrix(FALSE, n_group_classes, length(vars))
  attempts <- lapply(seq_len(n_starts), function(start_index) {
    tryCatch({
      initial <- .additive_initialize(stats, n_group_classes, structure,
                                      min_variance, start_index)
      fit <- .additive_em(stats, initial, structure, min_variance, max_iter,
                          tol, none)
      if (max_iter > 0L) {
        .additive_polish(stats, fit, structure, min_variance, max_iter, tol)
      } else c(fit, list(zero = none, kkt = NA))
    }, error = function(error) list(error = conditionMessage(error)))
  })
  valid <- vapply(attempts, function(attempt) is.null(attempt$error), logical(1))
  if (!any(valid)) {
    stop(errorCondition(sprintf(
      "All %d starts failed: %s", n_starts,
      paste(unique(vapply(attempts, `[[`, character(1), "error")), collapse = "; ")),
      class = "latents_all_starts_failed", call = NULL))
  }
  scores <- vapply(attempts, function(attempt) {
    if (is.null(attempt$error)) attempt$expectation$log_likelihood else -Inf
  }, numeric(1))
  converged <- vapply(attempts, function(attempt) isTRUE(attempt$converged),
                      logical(1))
  best_start <- .multilpa_select_start(scores, converged)
  best <- attempts[[best_start]]
  # An interior solution is finished by Newton steps (see .additive_newton());
  # a boundary solution is left as the KKT pass returned it.
  newton <- list(steps = 0L, decrement = NA_real_)
  if (max_iter > 0L && isTRUE(best$converged) && !any(best$zero) &&
      all(best$parameters$within > min_variance * (1 + sqrt(.Machine$double.eps)))) {
    newton <- .additive_newton(stats, best$parameters, structure, min_variance)
    best$parameters <- newton$parameters
    best$expectation <- .additive_expectation(stats, best$parameters)
    scores[best_start] <- best$expectation$log_likelihood
  }
  # The reported likelihood must be the likelihood at the reported parameters.
  check <- .additive_expectation(stats, best$parameters)
  stopifnot("selected likelihood must reproduce at the returned parameters" =
              abs(check$log_likelihood - scores[best_start]) <=
              1e-10 * (1 + abs(check$log_likelihood)))
  class_names <- paste0("group_class_", seq_len(n_group_classes))
  # Labels are arbitrary; report the largest class first (ties by the first
  # indicator's mean, then its within variance) so the labelling does not
  # depend on which start won. Every class-indexed block moves together.
  canonical <- order(-best$parameters$weights, best$parameters$means[, 1L],
                     best$parameters$within[, 1L], seq_len(n_group_classes))
  parameters <- best$parameters
  parameters$means <- parameters$means[canonical, , drop = FALSE]
  parameters$between <- parameters$between[canonical, , drop = FALSE]
  parameters$within <- parameters$within[canonical, , drop = FALSE]
  parameters$weights <- parameters$weights[canonical]
  best$zero <- best$zero[canonical, , drop = FALSE]
  dimnames(parameters$means) <- dimnames(parameters$between) <-
    dimnames(parameters$within) <- list(class_names, vars)
  names(parameters$weights) <- class_names
  expectation <- .additive_expectation(stats, parameters)
  dimnames(best$zero) <- list(class_names, vars)
  within_floor <- parameters$within <= min_variance * (1 + sqrt(.Machine$double.eps))
  n_parameters <- .additive_count_parameters(n_group_classes, length(vars),
                                             structure)
  log_likelihood <- expectation$log_likelihood
  tolerance <- 1e-6 * (1 + abs(log_likelihood))
  starts <- data.frame(
    start = seq_len(n_starts),
    log_likelihood = ifelse(valid, scores, NA_real_),
    converged = converged,
    iterations = vapply(attempts, function(attempt) {
      if (is.null(attempt$error)) as.integer(attempt$iterations) else NA_integer_
    }, integer(1)),
    between_zero = vapply(attempts, function(attempt) {
      if (is.null(attempt$error)) sum(attempt$zero) else NA_integer_
    }, integer(1)),
    error = vapply(attempts, function(attempt) attempt$error %||% NA_character_,
                   character(1)),
    selected = seq_len(n_starts) == best_start)
  small_classes <- any(colSums(expectation$posterior) < 1)
  result <- list(
    call = call, family = family, structure = structure,
    between_variance = between_variance, vars = vars, id = id,
    means = parameters$means, between_variances = parameters$between,
    within_variances = parameters$within,
    group_probabilities = parameters$weights,
    log_likelihood = log_likelihood,
    group_log_likelihood = expectation$group_log_likelihood,
    group_posteriors = expectation$posterior,
    group_assignments = stats::setNames(
      max.col(expectation$posterior, ties.method = "first"), stats$ids),
    conditional_intercept_means = expectation$conditional_mean,
    conditional_intercept_variances = expectation$conditional_variance,
    intercept_means = expectation$intercept_mean,
    intercept_variances = expectation$intercept_variance,
    n_groups = stats$n_groups, n_obs = stats$n_obs,
    group_sizes = stats::setNames(stats$sizes, stats$ids),
    n_parameters = n_parameters,
    aic = -2 * log_likelihood + 2 * n_parameters,
    bic = -2 * log_likelihood + n_parameters * log(stats$n_groups),
    bic_individual = -2 * log_likelihood + n_parameters * log(stats$n_obs),
    converged = isTRUE(best$converged), iterations = best$iterations,
    newton_steps = newton$steps, newton_decrement = newton$decrement,
    between_zero = best$zero, kkt = best$kkt, within_floor = within_floor,
    zero_within_scatter = stats$pooled_scatter <= 0,
    small_classes = small_classes,
    starts = starts, best_start = best_start, n_failed_starts = sum(!valid),
    n_best_replicated = sum(valid & abs(scores - log_likelihood) <= tolerance),
    replication_tolerance = tolerance,
    log_likelihood_history = best$history,
    min_variance = min_variance, tol = tol, max_iter = max_iter, seed = seed,
    weights = weights,
    sampling_weights = if (is.null(stats$sampling_weights)) NULL else
      stats::setNames(stats$sampling_weights, stats$ids),
    sufficient_statistics = stats)
  class(result) <- "multilpa_additive"
  if (any(!valid)) {
    warning(warningCondition(sprintf(
      "%d of %d starts failed; the start records keep their errors.",
      sum(!valid), n_starts), class = "latents_failed_starts", call = NULL))
  }
  if (!result$converged && max_iter > 0L) {
    warning(warningCondition(
      "The best start did not converge; increase max_iter.",
      class = "latents_unconverged", call = NULL))
  }
  if (any(best$zero) || any(within_floor)) {
    warning(warningCondition(paste(
      "A between-group variance is estimated at zero or a within-group variance",
      "reached min_variance; this is a boundary fit and Wald inference does not apply."),
      class = "latents_boundary", call = NULL))
  }
  if (small_classes) {
    warning(warningCondition("A group class has effective membership below one.",
                             class = "latents_small_classes", call = NULL))
  }
  result
}

#' Print a group-class fit
#'
#' A short summary of a `multilpa(family = )` group-class fit (additive,
#' dispersion or additive-dispersion): the family and between-variance
#' restriction, sizes, likelihood, parameter count, BIC over groups,
#' convergence and whether the fit is on a boundary.
#'
#' @param x A `multilpa_additive` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.multilpa_additive <- function(x, ...) {
  cat(sprintf(paste0(
    "%s group-class model (%s between variances)\n",
    "%d group classes, %d groups, %d observations, %d indicators\n",
    "log-likelihood %.3f, %d parameters, BIC (groups) %.3f\n",
    "converged: %s; best start replicated %d of %d\n"),
    .additive_family_label(x$family), x$between_variance,
    length(x$group_probabilities), x$n_groups, x$n_obs,
    length(x$vars), x$log_likelihood, x$n_parameters, x$bic,
    if (x$converged) "yes" else "no", x$n_best_replicated, nrow(x$starts)))
  .latents_print_weights(x)
  if (any(x$between_zero) || any(x$within_floor)) {
    cat("Boundary fit: a variance is at zero or at min_variance.\n")
  }
  cat(paste("Tables: get_results(x, what = ), e.g. \"parameters\",",
            "\"group_classes\", \"groups\", \"intercepts\"; summary(x); plot(x).\n"))
  invisible(x)
}
