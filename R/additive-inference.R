# Wald inference for the group-class families (additive, dispersion,
# additive-dispersion). The estimation scale is
# class means, log within variances, log between variances and baseline class
# logits (group class 1 is the baseline); natural-scale tables come from it by
# the delta method.

#' Coefficient names of one parameter block
#'
#' A block shared by every class has one coefficient per indicator
#' (`prefix.indicator`); a class-specific block has one per class and
#' indicator (`prefix.group_class_h.indicator`), class-major.
#' @noRd
.additive_block_names <- function(prefix, mode, class_names, vars) {
  if (identical(mode, "equal")) return(paste(prefix, vars, sep = "."))
  paste(prefix, rep(class_names, each = length(vars)), vars, sep = ".")
}

#' The free values of an H x d block: its first row when shared
#' @noRd
.additive_block_pack <- function(values, mode) {
  if (identical(mode, "equal")) values[1L, ] else as.vector(t(values))
}

#' Expand a block's free values back to an H x d matrix
#' @noRd
.additive_block_unpack <- function(values, mode, n_classes, n_indicators) {
  rows <- .additive_rows(mode, n_classes)
  block <- matrix(values, rows, n_indicators, byrow = TRUE)
  block[rep(seq_len(rows), length.out = n_classes), , drop = FALSE]
}

#' Combine per-class J x d score pieces: summed when the block is shared
#' @noRd
.additive_block_combine <- function(pieces, mode) {
  if (identical(mode, "equal")) Reduce(`+`, pieces) else do.call(cbind, pieces)
}

#' Names of the free (estimation-scale) and natural coefficients
#'
#' Blocks in order: means, within variances, between variances, then class
#' logits (estimation scale; group class 1 is the baseline) or weights
#' (natural scale).
#' @param class_names,vars Labels of the group classes and indicators.
#' @param structure Result of `.additive_structure()`.
#' @param scale `"unconstrained"` or `"natural"`.
#' @return A character vector in the packing order of `.additive_pack()`.
#' @noRd
.additive_names <- function(class_names, vars, structure,
                            scale = c("unconstrained", "natural")) {
  scale <- match.arg(scale)
  log_prefix <- if (identical(scale, "unconstrained")) "log_" else ""
  weights <- if (identical(scale, "unconstrained")) {
    sprintf("logit.%s", class_names[-1L])
  } else paste0("weight.", class_names)
  c(.additive_block_names("mean", structure$means, class_names, vars),
    .additive_block_names(paste0(log_prefix, "within_variance"), structure$within,
                          class_names, vars),
    .additive_block_names(paste0(log_prefix, "between_variance"),
                          structure$between, class_names, vars),
    weights)
}

#' Pack group-class parameters onto the estimation scale
#' @return A named numeric vector.
#' @noRd
.additive_pack <- function(parameters, structure, class_names, vars) {
  parameters <- .additive_check_parameters(parameters, length(vars))
  stats::setNames(c(
    .additive_block_pack(parameters$means, structure$means),
    log(.additive_block_pack(parameters$within, structure$within)),
    log(.additive_block_pack(parameters$between, structure$between)),
    log(parameters$weights[-1L] / parameters$weights[1L])),
    .additive_names(class_names, vars, structure))
}

#' Unpack an estimation-scale vector into group-class parameters
#' @return Group-class parameters (unlabelled H x d matrices and weights).
#' @noRd
.additive_unpack <- function(theta, n_classes, n_indicators, structure) {
  theta <- unname(theta)
  sizes <- n_indicators * c(.additive_rows(structure$means, n_classes),
                            .additive_rows(structure$within, n_classes),
                            .additive_rows(structure$between, n_classes))
  at <- cumsum(c(0L, sizes))
  block <- function(k) theta[seq.int(at[k] + 1L, length.out = sizes[k])]
  logits <- c(0, theta[-seq_len(at[4L])])
  list(means = .additive_block_unpack(block(1L), structure$means, n_classes,
                                      n_indicators),
       within = exp(.additive_block_unpack(block(2L), structure$within,
                                           n_classes, n_indicators)),
       between = exp(.additive_block_unpack(block(3L), structure$between,
                                            n_classes, n_indicators)),
       weights = exp(logits - max(logits)) / sum(exp(logits - max(logits))))
}

#' Per-group scores on the estimation scale
#'
#' Fisher's identity: each group's score is its posterior-weighted component
#' score, plus the class-logit score `R_jk - omega_k`. With
#' `D = W + n T` and `e = a - mu`, the component derivatives are
#' `n e / D` (mean), `-0.5 W [(n - 1)/W + 1/D - S/W^2 - n e^2/D^2]` (log W)
#' and `0.5 T n [n e^2/D^2 - 1/D]` (log T). A shared block's score is the sum
#' of the class contributions.
#'
#' @return A J x p matrix, one row per group, columns named as the
#'   estimation-scale coefficients.
#' @noRd
.additive_group_scores <- function(stats, parameters, structure, class_names,
                                   vars) {
  parameters <- .additive_check_parameters(parameters, length(vars))
  expectation <- .additive_expectation(stats, parameters)
  # A weighted group's score is its weight times its own score.
  posterior <- expectation$posterior * (stats$sampling_weights %||% 1)
  n_classes <- nrow(parameters$means)
  j_count <- stats$n_groups
  pieces <- lapply(seq_len(n_classes), function(h) {
    within <- .additive_row(parameters$within, h, j_count)
    between <- .additive_row(parameters$between, h, j_count)
    deviation <- stats$averages - .additive_row(parameters$means, h, j_count)
    denominator <- within + stats$sizes * between
    weight <- posterior[, h]
    list(
      mean = weight * stats$sizes * deviation / denominator,
      log_within = weight * -0.5 * within *
        ((stats$sizes - 1) / within + 1 / denominator - stats$scatter / within^2 -
           stats$sizes * deviation^2 / denominator^2),
      log_between = weight * 0.5 * between * stats$sizes *
        (stats$sizes * deviation^2 / denominator^2 - 1 / denominator))
  })
  part <- function(name, mode) {
    .additive_block_combine(lapply(pieces, `[[`, name), mode)
  }
  logits <- posterior[, -1L, drop = FALSE] -
    outer(rowSums(posterior), parameters$weights[-1L])
  scores <- cbind(part("mean", structure$means),
                  part("log_within", structure$within),
                  part("log_between", structure$between), logits)
  dimnames(scores) <- list(stats$ids,
                           .additive_names(class_names, vars, structure))
  scores
}

#' Observed information of an additive fit
#'
#' The Jacobian of the analytic total score by a fourth-order central
#' difference (step scaled to each coordinate), symmetrized.
#' @return A p x p matrix: minus the Hessian of the log likelihood.
#' @noRd
.additive_information <- function(stats, theta, structure, class_names,
                                  vars, step) {
  n_classes <- length(class_names)
  total <- function(point) {
    colSums(.additive_group_scores(
      stats, .additive_unpack(point, n_classes, length(vars), structure),
      structure, class_names, vars))
  }
  # The shared five-point score differences with h = step * max(1, |theta|):
  # -(c + t(c)) / 2 equals the former -((c + t(c)) / 2) exactly, since
  # negation and halving are exact.
  .inference_score_information(total, theta, step, scheme = "five_point",
                               h_rule = "max1")
}

#' Estimation-scale covariance of an additive fit, or a classed refusal
#'
#' @param x A `multilpa_additive` fit.
#' @param vcov_type `"observed"` (inverse observed information), `"robust"`
#'   (the CR0 group-score sandwich) or `"opg"` (inverse outer product of the
#'   group scores).
#' @param step Relative finite-difference step for the information.
#' @return A list with `theta`, `vcov`, `information`, `gradient`,
#'   `group_scores`, `scaled_score`, `condition_ratio` and `vcov_type`.
#' @noRd
.additive_inference <- function(x, vcov_type = c("observed", "robust", "opg"),
                                step = 1e-4, warn_weak = TRUE) {
  vcov_type <- match.arg(vcov_type)
  stopifnot("`step` must be finite and positive" =
              is.numeric(step) && length(step) == 1L && is.finite(step) && step > 0)
  if (!isTRUE(x$converged)) {
    stop(errorCondition("Inference requires a converged fit.",
                        class = "latents_no_converge", call = NULL))
  }
  if (any(x$between_zero) || any(x$within_floor)) {
    stop(errorCondition(paste(
      "A between-group variance is estimated at zero or a within-group variance",
      "sits at min_variance; the parameter is on the edge of its space and Wald",
      "inference does not apply."), class = "latents_boundary_fit", call = NULL))
  }
  if (any(x$group_probabilities <= 1e-6) || isTRUE(x$small_classes)) {
    stop(errorCondition(
      "A group class has (nearly) no membership, so its parameters are not identified.",
      class = "latents_boundary_fit", call = NULL))
  }
  stats <- x$sufficient_statistics
  class_names <- names(x$group_probabilities)
  parameters <- list(means = unname(x$means),
                     between = unname(x$between_variances),
                     within = unname(x$within_variances),
                     weights = unname(x$group_probabilities))
  theta <- .additive_pack(parameters, x$structure, class_names, x$vars)
  scores <- .additive_group_scores(stats, parameters, x$structure,
                                   class_names, x$vars)
  information <- .additive_information(stats, theta, x$structure,
                                       class_names, x$vars, step)
  # The ratio rule (and a non-positive smallest eigenvalue) refuses with the
  # ratio in the message; the ratio is also reported with the inference.
  inverse <- .inference_invert(information, "ratio_le_eps")
  eigenvalues <- eigen(information, symmetric = TRUE, only.values = TRUE)$values
  condition_ratio <- min(eigenvalues) / max(eigenvalues)
  dimnames(inverse) <- dimnames(information)
  gradient <- colSums(scores)
  scaled_score <- sqrt(max(0, drop(crossprod(gradient, inverse %*% gradient))))
  if (scaled_score > 1e-2) {
    warning(warningCondition(sprintf(paste(
      "The score at the estimate is not near zero (scaled %.3g); the fit may",
      "not be at a maximum. Refit with a smaller `tol`."), scaled_score),
      class = "latents_unconverged", call = NULL))
  }
  if (vcov_type != "observed") {
    .inference_require_rank(scores, paste(
      "The %d group scores span %d of %d dimensions, too few to form a",
      "%s covariance."), formatted = TRUE, extra = list(vcov_type))
  }
  meat <- crossprod(scores)
  covariance <- switch(vcov_type,
    observed = inverse,
    robust = .inference_sandwich(inverse, meat),
    opg = chol2inv(chol(meat)))
  dimnames(covariance) <- dimnames(information)
  inference <- list(theta = theta, vcov = covariance, information = information,
                    gradient = gradient, group_scores = scores,
                    scaled_score = scaled_score, condition_ratio = condition_ratio,
                    vcov_type = vcov_type)
  if (isTRUE(warn_weak)) .additive_warn_weak(x, inference)
  inference
}

#' Effective groups below which Wald intervals are flagged
#'
#' Chosen by the rule declared in validation/additive-diagnostic-protocol.md
#' (largest grid value with at most 1% flags in the well-separated reference
#' simulation cells, replicates 1-500) and evaluated on held-out replicates:
#' 0.8% flags in the reference cells, 100% in the rare, weakly separated cell,
#' class-parameter coverage 0.948 unflagged against 0.800 flagged. See
#' validation/ADDITIVE_SIMULATION.md.
#' @noRd
.additive_weak_class_threshold <- 50

#' Warn when a group class rests on few effective groups
#' @return `NULL`, invisibly; raises a `latents_weak_class` warning.
#' @noRd
.additive_warn_weak <- function(x, inference) {
  effective <- .additive_effective_groups(x, inference)
  weak <- is.finite(effective$effective_groups) &
    effective$effective_groups < .additive_weak_class_threshold
  if (any(weak)) {
    # effective = expected class size x information ratio, so the two factors
    # say whether the class is small or poorly separated.
    detail <- sprintf("%s: %.1f effective of %.1f expected groups, %.0f%% of the information kept",
                      effective$group_class[weak], effective$effective_groups[weak],
                      x$n_groups * effective$weight[weak],
                      100 * effective$information_ratio[weak])
    warning(warningCondition(sprintf(paste(
      "Group classes supported by fewer than %d effective groups (%s). A small",
      "information share means poor separation; a small expected count means",
      "few groups. Wald intervals for class means and weights can be",
      "miscalibrated; see get_results(x, \"group_classes\")."),
      .additive_weak_class_threshold, paste(detail, collapse = "; ")),
      class = "latents_weak_class", call = NULL))
  }
  invisible(NULL)
}

#' Natural coefficients and their Jacobian with respect to the free ones
#' @return A list with `estimate` (named natural vector) and `jacobian`
#'   (natural x free).
#' @noRd
.additive_natural <- function(theta, n_classes, vars, structure, class_names) {
  parameters <- .additive_unpack(theta, n_classes, length(vars), structure)
  weights <- parameters$weights
  means <- .additive_block_pack(parameters$means, structure$means)
  variances <- c(.additive_block_pack(parameters$within, structure$within),
                 .additive_block_pack(parameters$between, structure$between))
  estimate <- c(means, variances, weights)
  n_head <- length(means)
  n_var <- length(variances)
  # The first class is the logit reference, so its column is dropped.
  softmax <- .inference_simplex_jacobian(weights, "first")
  jacobian <- matrix(0, length(estimate), length(theta))
  jacobian[seq_len(n_head), seq_len(n_head)] <- diag(n_head)
  positions <- n_head + seq_len(n_var)
  jacobian[positions, positions] <- diag(variances, n_var)
  jacobian[n_head + n_var + seq_len(n_classes),
           n_head + n_var + seq_len(n_classes - 1L)] <- softmax
  names(estimate) <- .additive_names(class_names, vars, structure, "natural")
  dimnames(jacobian) <- list(names(estimate), names(theta))
  list(estimate = estimate, jacobian = jacobian)
}

#' Labels of one block's rows in the parameter table
#' @noRd
.additive_block_labels <- function(mode, class_names, vars) {
  if (identical(mode, "equal")) {
    return(list(group_class = rep("shared", length(vars)), indicator = vars))
  }
  list(group_class = rep(class_names, each = length(vars)),
       indicator = rep(vars, length(class_names)))
}

#' The tidy parameter table of a group-class fit
#'
#' @param x A `multilpa_additive` fit.
#' @param inference Result of `.additive_inference()`, or `NULL` for a table
#'   without standard errors.
#' @param level Confidence level.
#' @param adjust p-value adjustment across the class means (the only rows
#'   with a Wald test).
#' @return A data frame, one row per natural parameter.
#' @noRd
.additive_parameter_table <- function(x, inference, level = 0.95,
                                      adjust = "none") {
  stopifnot("`level` must be a single finite number in (0, 1)" =
              is.numeric(level) && length(level) == 1L && is.finite(level) &&
              level > 0 && level < 1)
  class_names <- names(x$group_probabilities)
  n_classes <- length(class_names)
  vars <- x$vars
  structure <- x$structure
  blocks <- list(
    means = c(.additive_block_labels(structure$means, class_names, vars),
              list(level = "between", parameter = "mean")),
    within = c(.additive_block_labels(structure$within, class_names, vars),
               list(level = "within", parameter = "variance")),
    between = c(.additive_block_labels(structure$between, class_names, vars),
                list(level = "between", parameter = "variance")))
  table <- rbind(
    do.call(rbind, lapply(blocks, function(block) {
      data.frame(level = block$level, group_class = block$group_class,
                 indicator = block$indicator, parameter = block$parameter,
                 stringsAsFactors = FALSE)
    })),
    data.frame(level = "group_class", group_class = class_names,
               indicator = NA_character_, parameter = "weight",
               stringsAsFactors = FALSE))
  theta <- .additive_pack(
    list(means = unname(x$means), between = unname(x$between_variances),
         within = unname(x$within_variances),
         weights = unname(x$group_probabilities)),
    structure, class_names, vars)
  natural <- .additive_natural(theta, n_classes, vars, structure, class_names)
  stopifnot("one table row per natural coefficient" =
              nrow(table) == length(natural$estimate))
  table$estimate <- unname(natural$estimate)
  standard_error <- if (is.null(inference)) {
    rep(NA_real_, nrow(table))
  } else {
    .inference_delta_se(natural$jacobian, inference$vcov, "quadratic")
  }
  table$standard_error <- unname(standard_error)
  tested <- table$parameter == "mean"
  table$statistic <- ifelse(tested, table$estimate / table$standard_error,
                            NA_real_)
  table$p_value <- 2 * stats::pnorm(-abs(table$statistic))
  table$p_adjusted <- stats::p.adjust(table$p_value, method = adjust)
  kind <- ifelse(table$parameter == "weight", "probability",
                 ifelse(table$parameter == "variance", "positive", "real"))
  bounds <- .multilpa_wald_bounds(table$estimate, table$standard_error, kind,
                                  .inference_critical(level, "one_minus"))
  table$conf_low <- bounds$low
  table$conf_high <- bounds$high
  rownames(table) <- NULL
  table
}

#' Finish an interior EM solution with safeguarded Newton steps
#'
#' EM stops on the likelihood's relative change, which leaves a slowly moving
#' fit (weakly separated or rare classes) measurably short of the maximum.
#' Newton steps on the estimation scale, with the analytic score and the
#' observed information, converge quadratically from there. Each step is
#' halved until the likelihood does not fall and every within variance stays
#' at or above `min_variance`; a non-positive-definite information stops the
#' polishing and returns the EM point unchanged.
#'
#' @return A list with `parameters`, `steps` (Newton steps taken) and
#'   `decrement` (the final Newton decrement, `NA` if not computed).
#' @noRd
.additive_newton <- function(stats, parameters, structure, min_variance,
                             max_steps = 25L, step = 1e-4) {
  n_classes <- nrow(parameters$means)
  vars <- stats$vars
  class_names <- paste0("group_class_", seq_len(n_classes))
  theta <- .additive_pack(parameters, structure, class_names, vars)
  unpack <- function(point) {
    .additive_unpack(point, n_classes, length(vars), structure)
  }
  log_likelihood <- function(point) .additive_expectation(stats, unpack(point))$log_likelihood
  current <- log_likelihood(theta)
  steps <- 0L
  decrement <- NA_real_
  # Sequential by nature: each Newton step starts from the last accepted point.
  repeat {
    point <- unpack(theta)
    gradient <- colSums(.additive_group_scores(stats, point, structure,
                                               class_names, vars))
    information <- .additive_information(stats, theta, structure,
                                         class_names, vars, step)
    factor <- tryCatch(chol(information), error = function(error) NULL)
    if (is.null(factor)) break
    direction <- backsolve(factor, forwardsolve(t(factor), gradient))
    decrement <- sqrt(max(0, sum(gradient * direction)))
    if (decrement < 1e-8 || steps >= max_steps) break
    scale <- 1
    accepted <- FALSE
    while (scale >= 1 / 1024) {
      candidate <- theta + scale * direction
      if (all(unpack(candidate)$within >= min_variance)) {
        value <- log_likelihood(candidate)
        if (is.finite(value) && value >= current - 1e-12 * (1 + abs(current))) {
          accepted <- TRUE
          break
        }
      }
      scale <- scale / 2
    }
    if (!accepted) break
    theta <- candidate
    current <- value
    steps <- steps + 1L
  }
  list(parameters = unpack(theta), steps = steps, decrement = decrement)
}

#' Effective number of groups behind each group class
#'
#' With every group's class known, a weight estimate has variance
#' `w (1 - w) / J`, so `w^2 (1 - w) / Var(w)` equals the expected number of
#' groups in the class. Uncertain classification inflates `Var(w)`, so this
#' count falls when a class is small, poorly separated, or both; `Var(w)` is
#' the delta-method variance from the requested covariance.
#'
#' @param x A `multilpa_additive` fit.
#' @param inference Result of `.additive_inference()`.
#' @return A data frame, one row per group class: `group_class`, `weight`,
#'   `effective_groups` and `information_ratio` (complete-data over observed
#'   variance of the weight; one minus the fraction of missing information).
#'   Both are `NA` for a one-class model.
#' @noRd
.additive_effective_groups <- function(x, inference) {
  class_names <- names(x$group_probabilities)
  natural <- .additive_natural(inference$theta, length(class_names), x$vars,
                               x$structure, class_names)
  rows <- paste0("weight.", class_names)
  covariance <- .inference_delta_covariance(natural$jacobian, inference$vcov)
  weight <- unname(natural$estimate[rows])
  variance <- unname(diag(covariance)[rows])
  defined <- length(class_names) > 1L & variance > 0
  data.frame(
    group_class = class_names, weight = weight,
    effective_groups = ifelse(defined, weight^2 * (1 - weight) / variance, NA_real_),
    information_ratio = ifelse(defined,
                               weight * (1 - weight) / (x$n_groups * variance),
                               NA_real_),
    stringsAsFactors = FALSE)
}
