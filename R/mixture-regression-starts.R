# Starting values, start selection and class ordering for mixture_regression().

#' Neutral parameters to seed the first M-step
#' @param spec The model specification.
#' @return A parameter list.
#' @noRd
.mixture_null_params <- function(spec) {
  n_classes <- spec$n_classes
  class_names <- paste0("class_", seq_len(n_classes))
  params <- list(
    beta = matrix(0, ncol(spec$x), n_classes,
                  dimnames = list(colnames(spec$x), class_names)),
    common = stats::setNames(numeric(ncol(spec$z)), colnames(spec$z)),
    sigma2 = rep(if (identical(spec$family, "gaussian"))
      max(stats::var(spec$y), spec$min_variance) else 1, n_classes))
  if (identical(spec$family, "ordinal")) {
    params$thresholds <- .ordinal_start_thresholds(spec)
  }
  if (identical(spec$nesting, "two-level")) {
    params$delta <- matrix(0, ncol(spec$v), spec$n_group_classes)
    params$class_logits <- matrix(0, spec$n_group_classes + ncol(spec$w),
                                  n_classes)
  } else {
    params$gamma <- matrix(0, ncol(spec$w), n_classes)
  }
  params
}

#' Soften a hard partition into posterior-like weights
#' @param labels Integer labels.
#' @param n_categories Number of categories.
#' @return A matrix whose rows sum to one.
#' @noRd
.mixture_soften <- function(labels, n_categories) {
  weights <- matrix(0.2 / n_categories, length(labels), n_categories)
  weights[cbind(seq_along(labels), labels)] <-
    weights[cbind(seq_along(labels), labels)] + 0.8
  weights
}

#' Labels from the residuals of one pooled regression
#'
#' Classes that differ in intercept separate along the residuals of a pooled
#' fit, so cutting them into quantile bins is a deterministic start that
#' finds the common case directly.
#'
#' @param spec The model specification.
#' @return Integer labels per top-level unit.
#' @noRd
.mixture_residual_labels <- function(spec) {
  pooled <- spec
  pooled$n_classes <- 1L
  pooled$nesting <- "observation"
  params <- .mixture_null_params(pooled)
  params$sigma2 <- params$sigma2[1L]
  params <- .mixture_update_regression(pooled, params,
                                      matrix(1, spec$n, 1L))
  residual <- (spec$y - spec$trials *
                 switch(spec$family,
                        gaussian = drop(.mixture_linear_predictors(pooled, params)),
                        binomial = stats::plogis(
                          drop(.mixture_linear_predictors(pooled, params))),
                        poisson = ,
                        negative_binomial = exp(
                          drop(.mixture_linear_predictors(pooled, params))),
                        ordinal = drop(.mixture_class_means(pooled, params))))
  if (identical(spec$nesting, "group")) {
    residual <- drop(rowsum(residual, spec$group_index, reorder = TRUE)) /
      tabulate(spec$group_index)
  }
  ranks <- rank(residual, ties.method = "first")
  as.integer(ceiling(ranks * spec$n_classes / length(ranks)))
}

#' Pseudo E-step from labels, used to launch EM
#' @param spec The model specification.
#' @param labels Labels per top-level unit.
#' @param group_labels Group-class labels (two-level model only).
#' @return A list shaped like `.mixture_expectation()`'s result.
#' @noRd
.mixture_seed_expectation <- function(spec, labels, group_labels = NULL) {
  switch(spec$nesting,
    observation = list(tau = .mixture_soften(labels, spec$n_classes)),
    group = {
      group_tau <- .mixture_soften(labels, spec$n_classes)
      list(group_tau = group_tau,
           tau = group_tau[spec$group_index, , drop = FALSE])
    },
    "two-level" = {
      rho <- .mixture_soften(group_labels, spec$n_group_classes)
      tau <- .mixture_soften(labels, spec$n_classes)
      # Group classes start with different class mixes, or the second level
      # would begin at a symmetric point EM cannot leave.
      tau_by_group_class <- lapply(seq_len(spec$n_group_classes), function(h) {
        tilt <- matrix(1, spec$n, spec$n_classes)
        tilt[, 1L + (h - 1L) %% spec$n_classes] <- 2
        tilted <- tau * tilt
        rho[spec$group_index, h] * tilted / rowSums(tilted)
      })
      list(rho = rho, tau = Reduce(`+`, tau_by_group_class),
           tau_by_group_class = tau_by_group_class)
    })
}

#' Run EM from every start
#' @param spec The model specification.
#' @param n_starts Number of random starts.
#' @param max_iter Maximum EM iterations.
#' @param tol Convergence tolerance.
#' @return A list of EM results, the first from the residual start.
#' @noRd
.mixture_run_starts <- function(spec, n_starts, max_iter, tol) {
  units <- if (identical(spec$nesting, "observation")) spec$n else spec$n_groups
  random_labels <- function() sample(rep_len(seq_len(spec$n_classes), units))
  random_group_labels <- function() {
    if (!identical(spec$nesting, "two-level")) return(NULL)
    sample(rep_len(seq_len(spec$n_group_classes), spec$n_groups))
  }
  row_labels <- function(labels) {
    if (identical(spec$nesting, "two-level") && length(labels) != spec$n) {
      labels[spec$group_index]
    } else labels
  }
  launches <- c(
    list(list(labels = .mixture_residual_labels(spec),
              group_labels = random_group_labels())),
    lapply(seq_len(n_starts), function(start) {
      list(labels = random_labels(), group_labels = random_group_labels())
    }))
  # Two stages, as Mplus's STARTS = a b: every start is run for a short
  # screening phase, and only the more promising half is carried on to
  # convergence. A start stalled near a saddle point would otherwise spend
  # the whole iteration budget going nowhere.
  screening <- min(max_iter, 50L)
  screened <- lapply(launches, function(launch) {
    labels <- launch$labels
    if (identical(spec$nesting, "two-level")) labels <- row_labels(labels)
    seed <- .mixture_seed_expectation(spec, labels, launch$group_labels)
    params <- .mixture_maximization(spec, .mixture_null_params(spec), seed)
    result <- .mixture_em(spec, params, screening, tol)
    result$stage <- "screened"
    result
  })
  log_likelihood <- vapply(screened, function(r) r$expectation$log_likelihood,
                           numeric(1))
  n_final <- max(2L, ceiling(length(screened) / 2))
  finalists <- utils::head(order(-log_likelihood), n_final)
  lapply(seq_along(screened), function(start) {
    result <- screened[[start]]
    if (!start %in% finalists || result$converged || result$degenerate ||
        max_iter <= screening) {
      if (start %in% finalists || result$converged) result$stage <- "completed"
      return(result)
    }
    continued <- .mixture_em(spec, result$params, max_iter - screening, tol)
    continued$iterations <- continued$iterations + result$iterations
    continued$history <- c(result$history, continued$history[-1L])
    continued$stage <- "completed"
    continued
  })
}

#' Choose among starts, order the classes, and assemble the fit
#' @param spec The model specification.
#' @param results EM results from `.mixture_run_starts()`.
#' @param select_start Selection rule.
#' @return An object of class `latents_mixture_regression` without inference.
#' @noRd
.mixture_select <- function(spec, results, select_start) {
  log_likelihood <- vapply(results, function(r) r$expectation$log_likelihood,
                           numeric(1))
  converged <- vapply(results, `[[`, logical(1), "converged")
  degenerate <- vapply(results, `[[`, logical(1), "degenerate")
  usable <- is.finite(log_likelihood) & !degenerate
  if (!any(usable)) {
    stop(errorCondition(sprintf(paste(
      "Every one of the %d starts degenerated: a class emptied or collapsed",
      "onto too few rows to estimate its regression. Fewer classes, more",
      "starts, or a larger `min_variance` may help."), length(results)),
      class = "latents_no_valid_start", call = NULL))
  }
  if (any(degenerate)) {
    warning(warningCondition(sprintf(paste(
      "%d of %d starts degenerated (an emptied class or a collapsed",
      "variance) and were set aside."), sum(degenerate), length(results)),
      class = "latents_degenerate_start", call = NULL))
  }
  candidates <- which(usable)
  best <- candidates[.multilpa_select_start(log_likelihood[candidates],
                                            converged[candidates],
                                            select_start)]
  chosen <- results[[best]]
  if (!chosen$converged) {
    warning(warningCondition(sprintf(paste(
      "The selected start did not converge in %d iterations; the estimates",
      "are not a maximum. Raise `max_iter` or use",
      "`select_start = \"converged\"`."), chosen$iterations),
      class = "latents_unconverged", call = NULL))
  }
  params <- .mixture_order_classes(spec, chosen$params, chosen$expectation)
  # Reported posteriors are each unit's own; the likelihood is the weighted one.
  expectation <- .mixture_expectation(spec, params, weighted = FALSE)
  if (isTRUE(params$separation) || isTRUE(params$membership_separation) ||
      isTRUE(params$group_separation)) {
    warning(warningCondition(paste(
      "A coefficient exceeded 30 on the logit or log scale: a class",
      "separates perfectly on some predictor, and that coefficient and its",
      "standard error are not reliable."),
      class = "latents_separation", call = NULL))
  }
  if (identical(spec$family, "negative_binomial") &&
      any(params$sigma2 <= .latents_min_dispersion * (1 + 1e-6))) {
    warning(warningCondition(paste(
      "A negative-binomial dispersion is estimated at zero: that class's",
      "counts show no overdispersion, which is the Poisson limit. This is a",
      "boundary fit; its dispersion has no standard error, and",
      "`family = \"poisson\"` may suit those data."),
      class = "latents_boundary", call = NULL))
  }
  best_value <- log_likelihood[best]
  completed <- vapply(results, function(r) identical(r$stage, "completed"),
                      logical(1))
  replicated <- sum(usable & completed &
                      abs(log_likelihood - best_value) <= 1e-6 * (1 + abs(best_value)))
  starts <- data.frame(
    start = seq_along(results),
    kind = c("residuals", rep("random", length(results) - 1L)),
    log_likelihood = log_likelihood, converged = converged,
    iterations = vapply(results, `[[`, integer(1), "iterations"),
    stage = vapply(results, `[[`, character(1), "stage"),
    degenerate = degenerate, selected = seq_along(results) == best)
  structure(list(
    spec = spec, params = params, expectation = expectation,
    log_likelihood = .mixture_expectation(spec, params)$log_likelihood,
    n_parameters = .mixture_count_parameters(spec),
    converged = chosen$converged, iterations = chosen$iterations,
    history = chosen$history, starts = starts,
    n_best_replicated = replicated),
    class = "latents_mixture_regression")
}

#' Number of free parameters
#' @param spec The model specification.
#' @return A single integer.
#' @noRd
.mixture_count_parameters <- function(spec) {
  k <- spec$n_classes
  regression <- ncol(spec$x) * k + ncol(spec$z)
  dispersion <- if (spec$family %in% c("gaussian", "negative_binomial")) {
    if (identical(spec$variance, "equal")) 1L else k
  } else if (identical(spec$family, "ordinal")) {
    (spec$n_categories - 1L) * k
  } else 0L
  mixing <- if (identical(spec$nesting, "two-level")) {
    ncol(spec$v) * (spec$n_group_classes - 1L) +
      (spec$n_group_classes + ncol(spec$w)) * (k - 1L)
  } else ncol(spec$w) * (k - 1L)
  as.integer(regression + dispersion + mixing)
}

#' Put classes in a deterministic order
#'
#' Classes are ordered by decreasing size, then by their first coefficient.
#' Group classes are ordered by decreasing probability of the (largest) first
#' class among their rows: group classes exist to differ in that composition,
#' whereas their sizes are often nearly equal, and ordering on a near-tie
#' would let the labels swap between otherwise identical fits. Label
#' switching would otherwise make the numbering depend on the start that won.
#'
#' @param spec The model specification.
#' @param params The parameter list.
#' @param expectation The E-step at `params`.
#' @return The reordered parameter list.
#' @noRd
.mixture_order_classes <- function(spec, params, expectation) {
  orders <- .mixture_class_orders(spec, params, expectation)
  .mixture_permute(spec, params, orders$class_order, orders$group_order)
}

#' The canonical order of the classes and group classes
#'
#' Classes by decreasing size, ties by their first coefficient; group classes
#' by how much of the largest class they hold, then by size.
#'
#' @param spec The structure.
#' @param params The parameter list (for `beta`).
#' @param expectation The structure's E-step.
#' @return A list of `class_order` and `group_order` (`NULL` for one level).
#' @noRd
.mixture_class_orders <- function(spec, params, expectation) {
  sizes <- colSums(expectation$tau)
  # An ordinal class may have no class-specific slope; its first threshold
  # breaks the tie instead.
  tie_break <- if (nrow(params$beta) > 0L) params$beta[1L, ] else params$thresholds[1L, ]
  class_order <- order(-round(sizes, 8), tie_break)
  group_order <- NULL
  if (identical(spec$nesting, "two-level")) {
    first_class <- class_order[1L]
    composition <- vapply(seq_len(spec$n_group_classes), function(h) {
      weight <- expectation$rho[spec$group_index, h]
      sum(expectation$tau_by_group_class[[h]][, first_class]) /
        max(sum(weight), .Machine$double.xmin)
    }, numeric(1))
    group_order <- order(-round(composition, 8),
                         -colSums(expectation$rho))
  }
  list(class_order = class_order, group_order = group_order)
}

#' Relabel classes and group classes
#'
#' Every logit block is re-expressed with the new first class (and first
#' group class) as its reference, which leaves the likelihood unchanged.
#'
#' @param spec The model specification.
#' @param params The parameter list.
#' @param class_order New order of the classes (old indices).
#' @param group_order New order of the group classes, or `NULL`.
#' @return The relabelled parameter list.
#' @noRd
.mixture_permute <- function(spec, params, class_order, group_order = NULL) {
  stopifnot(
    "`class_order` must be a permutation of the classes" =
      setequal(class_order, seq_len(spec$n_classes)) &&
      length(class_order) == spec$n_classes)
  params$beta <- params$beta[, class_order, drop = FALSE]
  colnames(params$beta) <- paste0("class_", seq_len(spec$n_classes))
  params$sigma2 <- params$sigma2[class_order]
  if (!is.null(params$thresholds)) {
    params$thresholds <- params$thresholds[, class_order, drop = FALSE]
    colnames(params$thresholds) <- colnames(params$beta)
  }
  if (identical(spec$nesting, "two-level")) {
    group_order <- group_order %||% seq_len(spec$n_group_classes)
    params$delta <- .latents_rebase_logits(params$delta, group_order, "first")
    logits <- .latents_rebase_logits(params$class_logits, class_order, "first")
    intercepts <- seq_len(spec$n_group_classes)
    logits[intercepts, ] <- logits[intercepts[group_order], , drop = FALSE]
    params$class_logits <- logits
  } else {
    params$gamma <- .latents_rebase_logits(params$gamma, class_order, "first")
  }
  params
}
