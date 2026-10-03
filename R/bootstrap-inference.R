# Cluster bootstrap inference, for every covariance structure.
#
# The Wald path differentiates one log variance per profile and indicator. A
# structure that constrains the volume, the shape or the orientation has fewer
# free parameters than that chart has coordinates, so its information matrix is
# singular by construction and the Wald path refuses it. The constrained model
# is still a model, and its natural parameters -- the means, the variances or
# covariances, the mixing probabilities -- are still estimable. The bootstrap
# reaches them without needing a chart at all: resample, refit inside the same
# family, and read the spread of the estimates.
#
# Groups are the resampling unit, not rows. Observations inside a group are not
# independent, and resampling rows would treat a repeated measure as fresh
# information and return intervals that are too narrow (Field & Welsh, 2007,
# Journal of the Royal Statistical Society B, 69, 369-390; Davison & Hinkley,
# 1997, Bootstrap Methods and their Application, chapter 3).

#' The `multilpa()` arguments that refit a model as it was fitted
#'
#' A refit that reads `variance_model` and `covariance_model` alone silently
#' widens any structure that constrains the volume, the shape or the
#' orientation: a `VEI` fit comes back as `VVI`, a different model with two more
#' parameters. The three pieces of the structure carry all of it, and naming
#' them alongside `covariance_model` is an error, so a fit that records a
#' structure is refitted from the structure and a fit that does not is refitted
#' from the two switches it does record.
#'
#' @param x A fitted `multilpa` model.
#' @return A named list of arguments for [multilpa()].
#' @noRd
.multilpa_refit_arguments <- function(x) {
  stopifnot(inherits(x, "multilpa"))
  shared <- list(vars = x$vars, id = x$id, n_profiles = x$n_profiles,
                 n_group_classes = x$n_group_classes,
                 categorical = x$categorical %||% character(),
                 min_variance = x$min_variance,
                 min_probability = x$min_probability %||% 1e-10,
                 missing = x$missing %||% "error",
                 time = x$time,
                 centering = x$centering %||% "none",
                 prior = x$prior,
                 noise = isTRUE(x$noise),
                 acceleration = x$acceleration %||% "none",
                 weights = x$weights, ordinal = x$ordinal %||% character(),
                 count = x$count %||% character(),
                 count_model = x$extra_data$count_model %||% "poisson",
                 count_dispersion = x$extra_data$count_dispersion %||% "varying")
  structure <- x$covariance_structure
  if (is.null(structure) || is.na(structure) || !structure %in% .multilpa_structures()) {
    return(c(shared, list(variance_model = x$variance_model,
                          covariance_model = x$covariance_model %||% "diagonal")))
  }
  c(shared, .multilpa_structure_arguments(structure))
}

#' The `fit_staged()` arguments that refit a staged model as it was fitted
#'
#' A staged fit is refitted by staging again, so that each resample estimates
#' its own measurement before its own group-class structure. That is what
#' carries the first stage's sampling variability into the second stage's
#' replicates; refitting with the measurement held would carry none of it.
#'
#' @param x A `multilpa` fit from [fit_staged()].
#' @return A named list of arguments for [fit_staged()], without `data` and the
#'   start and convergence controls.
#' @noRd
.multilpa_staged_refit_arguments <- function(x) {
  stopifnot(inherits(x, "multilpa"), isTRUE(x$staged))
  list(vars = x$vars, id = x$id, n_profiles = x$n_profiles,
       n_group_classes = x$n_group_classes,
       variance_model = x$variance_model,
       covariance_model = x$covariance_model %||% "diagonal",
       categorical = x$categorical %||% character(),
       min_variance = x$min_variance,
       min_probability = x$min_probability %||% 1e-10,
       missing = x$missing %||% "error", time = x$time)
}

#' Bootstrap replicates of the natural coefficients
#'
#' A staged fit is refitted through [fit_staged()], both stages on every
#' resample; any other fit through [multilpa()].
#'
#' @param x A fitted `multilpa` model.
#' @param data The fitting data.
#' @param iter How many resamples to draw.
#' @param n_starts,max_iter,tol Passed to each refit.
#' @return A list with `estimates`, an `iter` by parameters matrix whose failed
#'   rows are `NA`, and `messages`, the reason each failure gave.
#' @noRd
.multilpa_bootstrap_replicates <- function(x, data, iter, n_starts, max_iter, tol) {
  staged <- isTRUE(x$staged)
  estimator <- if (staged) fit_staged else multilpa
  arguments <- if (staged) .multilpa_staged_refit_arguments(x) else
    .multilpa_refit_arguments(x)
  .latents_cluster_bootstrap(
    data, x$group_index, x$n_groups, x$id, iter,
    refit = function(resampled) {
      do.call(estimator, c(list(data = resampled), arguments,
                           list(n_starts = n_starts, max_iter = max_iter, tol = tol)))
    },
    extract = function(fit) {
      unname(.multilpa_coefficients(.multilpa_align_labels(fit, x), "natural"))
    },
    reference_names = names(.multilpa_coefficients(x, "natural")))
}

#' Percentile bootstrap inference for a fitted model
#'
#' Called by [parameter_inference()] when `method = "bootstrap"`. The interval
#' is the percentile interval of the replicates and the standard error is their
#' standard deviation, so both are read off the same replicates and neither
#' assumes the estimate is normal around its truth.
#'
#' A [fit_staged()] result is bootstrapped by staging every resample again, so
#' the replicates vary the measurement as well as the group-class structure.
#' Its table therefore has a row for every parameter, the held measurement
#' included, and the second-stage errors carry the first stage's uncertainty.
#'
#' @param x A fitted `multilpa` model.
#' @param data The fitting data.
#' @param level Interval level.
#' @param iter How many resamples to draw.
#' @param n_starts,max_iter,tol Passed to each refit.
#' @param adjust Kept for a frame of the same shape as the Wald path's.
#' @return A `data.frame`, one row per natural coefficient.
#' @noRd
.multilpa_bootstrap_inference <- function(x, data, level, iter, n_starts,
                                          max_iter, tol, adjust) {
  staged <- isTRUE(x$staged)
  if (length(x$fixed %||% character()) > 0L && !staged) {
    stop(errorCondition(paste(
      "A fit that holds a block fixed cannot be bootstrapped here: the held",
      "values came from another fit, and resampling these data does not",
      "resample them. Bootstrap the fit the measurement came from, or fit",
      "both stages with fit_staged(), whose bootstrap refits each stage."),
      class = "latents_unsupported_inference", call = NULL))
  }
  if (!isTRUE(x$converged)) {
    stop(errorCondition(
      "The original fit did not converge; bootstrap it only once it has.",
      class = "latents_no_converge", call = NULL))
  }
  ## Every coefficient of a staged fit is resampled, the held measurement
  ## included, so every one of them gets a row.
  free_natural <- if (staged) seq_along(.multilpa_coefficients(x, "natural")) else
    .multilpa_free_index(x, "natural")
  drawn <- .multilpa_bootstrap_replicates(x, data, iter, n_starts, max_iter, tol)
  estimates <- drawn$estimates[, free_natural, drop = FALSE]
  valid <- stats::complete.cases(estimates)
  if (sum(valid) < 2L) {
    reasons <- drawn$messages[!is.na(drawn$messages)]
    stop(errorCondition(sprintf(
      "Only %d of %d resamples produced a usable fit; the first reason was: %s",
      sum(valid), iter,
      if (length(reasons) > 0L) reasons[1L] else "no reason was recorded"),
      class = "latents_bootstrap_failed", call = NULL))
  }
  kept <- estimates[valid, , drop = FALSE]
  if (sum(valid) < iter) {
    warning(warningCondition(sprintf(
      "%d of %d resamples did not produce a usable fit and were dropped.",
      iter - sum(valid), iter), class = "latents_bootstrap_dropped", call = NULL))
  }
  bounds <- c((1 - level) / 2, (1 + level) / 2)
  quantiles <- apply(kept, 2L, stats::quantile, probs = bounds, names = FALSE)
  point <- .multilpa_coefficients(x, "natural")[free_natural]
  result <- data.frame(
    .multilpa_coefficient_labels(x, "natural")[free_natural, , drop = FALSE],
    estimate = unname(point),
    standard_error = apply(kept, 2L, stats::sd),
    ## A Wald statistic would put a normal approximation back on top of the
    ## replicates the interval was read from, which is the assumption this path
    ## exists to avoid. The interval is the inference.
    statistic = NA_real_, p_value = NA_real_,
    conf_low = unname(quantiles[1L, ]), conf_high = unname(quantiles[2L, ]),
    row.names = NULL, stringsAsFactors = FALSE)
  result <- .multilpa_adjust_p(result, adjust)
  covariance <- stats::cov(kept)
  dimnames(covariance) <- list(names(point), names(point))
  attributes(result) <- c(attributes(result), list(
    covariance = covariance, covariance_unconstrained = NULL,
    level = level, method = "bootstrap", iter = iter, n_valid = sum(valid),
    replicates = kept, messages = drawn$messages,
    structure = x$covariance_structure,
    fixed = x$fixed %||% character(),
    stages_resampled = if (staged) 2L else 1L))
  result
}
