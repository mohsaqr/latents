# Maximum a posteriori estimation under mclust's default conjugate prior.
#
# Fraley, C., & Raftery, A. E. (2007). Bayesian regularization for normal
# mixture estimation and model-based clustering. Journal of Classification,
# 24, 155-181. doi:10.1007/s00357-007-0004-5
#
# The hyperparameters are those of `mclust::defaultPrior()` and each M-step is
# the one mclust's EM (`me()`, hence `Mclust()`) performs --- the Fortran
# `meeiip`, `meviip`, ..., `mevvvp` of mclust 6.1 --- reproduced term by term,
# so that a fit with `prior = prior_control()` is the fit
# `mclust::me(..., prior = priorControl())` gives from the same start. Several
# of those M-steps depart from a textbook conjugate update, and the departures
# are reproduced rather than corrected, because agreeing with mclust is the
# point of the option:
#
# * EII leaves the prior scale out of the variance and always counts the mean
#   prior in the denominator, `nu + (n + G) d + 2`;
# * EEI does not shrink the means, and EEI, EVI and EEV count the mean prior
#   once (`+ 1`) rather than once per profile;
# * VEV adds the prior scale and the mean term to the scatter but keeps the
#   maximum-likelihood denominators, ignoring the degrees of freedom.
#
# mclust's standalone `mstep()` differs from its own `me()` for EII, EEI, EVI
# and EEV (different denominators; EEV even drops the shrinkage); the EM, and
# therefore every fitted model, uses `me()`'s, which are the ones here.
#
# mclust has no prior for VEE, EVE, VVE or EVV, and neither has this package.

#' Conjugate prior for Gaussian mixture estimation
#'
#' Requests maximum a posteriori (MAP) estimation of the Gaussian profile means
#' and covariances under the conjugate prior of Fraley and Raftery (2007),
#' exactly as `mclust::priorControl()` does for `mclust::Mclust()`. Pass the
#' result as `multilpa(prior = )`. The prior keeps covariances away from
#' singularity, which is what lets a many-profile or many-indicator model be
#' fitted where maximum likelihood degenerates.
#'
#' Every argument left `NULL` takes `mclust::defaultPrior()`'s value, computed
#' from the indicators being fitted: the prior mean is their column means, the
#' degrees of freedom are `d + 2` for `d` indicators, and the scale is
#' `(1/G)^(2/d)` times the sample covariance (ellipsoidal structures) or times
#' the average sample variance (axis-parallel and spherical structures), for
#' `G` profiles.
#'
#' @param shrinkage Non-negative prior precision factor for the means, `kappa`;
#'   mclust's default `0.01`. Zero leaves the means unshrunk.
#' @param mean Optional prior mean, one value per continuous indicator, in the
#'   indicators' own units. Naming one needs a positive `shrinkage`.
#' @param dof Optional prior degrees of freedom, a single number.
#' @param scale Optional prior scale: a single positive number for an
#'   axis-parallel or spherical structure, or a positive-definite `d` by `d`
#'   matrix for an ellipsoidal one.
#' @return An object of class `latents_prior`: a list with `shrinkage`, `mean`,
#'   `dof` and `scale`, the last three `NULL` where the default is to be
#'   computed from the data.
#' @section Conditions:
#'   Raises `latents_bad_argument` for a negative or non-finite `shrinkage`, a
#'   `mean` with a zero `shrinkage`, or a non-positive `dof` or `scale`.
#' @references Fraley, C., & Raftery, A. E. (2007). Bayesian regularization for
#'   normal mixture estimation and model-based clustering. *Journal of
#'   Classification*, 24, 155--181. \doi{10.1007/s00357-007-0004-5}
#' @seealso [multilpa()], whose `prior` argument takes this.
#' @examples
#' prior_control()
#' prior_control(shrinkage = 0)
#' @export
prior_control <- function(shrinkage = 0.01, mean = NULL, dof = NULL,
                          scale = NULL) {
  bad <- function(message) {
    stop(errorCondition(message, class = "latents_bad_argument", call = NULL))
  }
  if (!is.numeric(shrinkage) || length(shrinkage) != 1L ||
      !is.finite(shrinkage) || shrinkage < 0) {
    bad("`shrinkage` must be a single non-negative number.")
  }
  if (!is.null(mean) && (!is.numeric(mean) || anyNA(mean) ||
                         any(!is.finite(mean)))) {
    bad("`mean` must be NULL or a finite numeric vector.")
  }
  if (!is.null(mean) && shrinkage == 0) {
    bad("A prior `mean` needs a positive `shrinkage`, or it has no effect.")
  }
  if (!is.null(dof) && (!is.numeric(dof) || length(dof) != 1L ||
                        !is.finite(dof) || dof <= 0)) {
    bad("`dof` must be NULL or a single positive number.")
  }
  if (!is.null(scale) && (!is.numeric(scale) || any(!is.finite(scale)) ||
                          any(diag(as.matrix(scale)) <= 0))) {
    bad("`scale` must be NULL, a positive number or a positive-definite matrix.")
  }
  structure(list(shrinkage = shrinkage, mean = mean, dof = dof, scale = scale),
            class = "latents_prior")
}

#' @rdname prior_control
#' @param x A `latents_prior` object.
#' @param ... Ignored.
#' @return `print()` returns `x` invisibly.
#' @export
print.latents_prior <- function(x, ...) {
  describe <- function(value) {
    if (is.null(value)) "mclust default (from the data)" else
      if (length(value) == 1L) format(value) else
        sprintf("supplied (%s)", paste(dim(as.matrix(value)), collapse = " x "))
  }
  cat("Conjugate prior (Fraley & Raftery, 2007)\n")
  cat(sprintf("  shrinkage: %s\n  mean:      %s\n  dof:       %s\n  scale:     %s\n",
              format(x$shrinkage), describe(x$mean), describe(x$dof),
              describe(x$scale)))
  invisible(x)
}

#' The structures mclust defines a conjugate prior for
#' @return A character vector of mclust model codes.
#' @noRd
.multilpa_prior_structures <- function() {
  c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", "VEV", "VVV")
}

#' Resolve a prior request into hyperparameters
#'
#' `mclust::defaultPrior()`, computed from the indicators the EM sees.
#'
#' @param prior A `latents_prior`.
#' @param x The indicator matrix, complete, in the units the EM uses.
#' @param n_profiles `G`.
#' @param structure The mclust model code.
#' @param offset What was subtracted from each column of `x` before fitting,
#'   so that a supplied prior mean in the caller's units can be moved with it.
#' @return A list with `shrinkage`, `mean` (in the units of `x`), `dof` and
#'   `scale` (a scalar or a matrix).
#' @noRd
.multilpa_prior_parameters <- function(prior, x, n_profiles, structure,
                                       offset = numeric(ncol(x))) {
  stopifnot("`prior` must be a latents_prior" = inherits(prior, "latents_prior"),
            "`x` must be a complete numeric matrix" =
              is.matrix(x) && is.numeric(x) && !anyNA(x))
  d <- ncol(x)
  n <- nrow(x)
  ellipsoidal <- .multilpa_is_ellipsoidal(structure)
  factor <- (1 / n_profiles)^(2 / d)
  covariance <- stats::var(x)
  default_scale <- if (ellipsoidal) {
    if (n > d) factor * covariance else factor * diag(diag(covariance), d)
  } else factor * sum(diag(covariance)) / d
  scale <- prior$scale %||% default_scale
  if (ellipsoidal) {
    scale <- as.matrix(scale)
    if (!identical(dim(scale), c(d, d)) ||
        max(abs(scale - t(scale))) > 1e-10 * max(abs(scale)) ||
        min(eigen(scale, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
      stop(errorCondition(sprintf(
        "The prior `scale` for the %s structure must be a positive-definite %d x %d matrix.",
        structure, d, d), class = "latents_bad_argument", call = NULL))
    }
  } else if (length(scale) != 1L) {
    stop(errorCondition(sprintf(
      "The prior `scale` for the %s structure must be a single positive number.",
      structure), class = "latents_bad_argument", call = NULL))
  }
  mean <- if (is.null(prior$mean)) colMeans(x) else {
    if (length(prior$mean) != d) {
      stop(errorCondition(sprintf(
        "The prior `mean` must have one value per continuous indicator (%d).", d),
        class = "latents_bad_argument", call = NULL))
    }
    as.numeric(prior$mean) - offset
  }
  list(shrinkage = prior$shrinkage, mean = unname(mean),
       dof = prior$dof %||% (d + 2), scale = unname(scale))
}

#' One M-step of the Gaussian parameters under the conjugate prior
#'
#' @param x The complete indicator matrix.
#' @param posteriors Observation-by-profile posterior probabilities.
#' @param structure The mclust model code, one of `.multilpa_prior_structures()`.
#' @param hyper Hyperparameters from `.multilpa_prior_parameters()`.
#' @param min_variance The variance or eigenvalue lower bound.
#' @param max_iter,tol The inner iteration of VEI and VEV.
#' @return A list with `means`, `variances` and, for an ellipsoidal structure,
#'   `covariances`.
#' @noRd
.multilpa_prior_maximize <- function(x, posteriors, structure, hyper,
                                     min_variance, max_iter = 1000L,
                                     tol = 1e-12) {
  stopifnot("`structure` must have a conjugate prior" =
              structure %in% .multilpa_prior_structures())
  d <- ncol(x)
  k <- ncol(posteriors)
  weights <- colSums(posteriors)
  sample_means <- sweep(crossprod(posteriors, x), 1L, weights, "/")
  scatter <- lapply(seq_len(k), function(profile) {
    residuals <- sweep(x, 2L, sample_means[profile, ], "-")
    crossprod(residuals, residuals * posteriors[, profile])
  })
  kappa <- hyper$shrinkage
  bump <- as.numeric(kappa > 0)
  differences <- sweep(sample_means, 2L, hyper$mean, "-")
  contraction <- kappa * weights / (weights + kappa)
  shrunk <- sweep(sweep(sample_means, 1L, weights, "*"), 2L,
                  kappa * hyper$mean, "+") / (weights + kappa)
  mean_terms <- lapply(seq_len(k), function(profile) {
    contraction[profile] * tcrossprod(differences[profile, ])
  })
  dof <- hyper$dof
  n <- sum(weights)
  scale <- hyper$scale
  diagonals <- .multilpa_scatter_diagonals(scatter)
  mean_diagonals <- .multilpa_scatter_diagonals(mean_terms)
  means <- if (identical(structure, "EEI")) sample_means else shrunk
  normalise <- function(values) values / exp(mean(log(values)))
  spread <- switch(structure,
    EII = {
      # Without the prior scale, and counting the mean prior whatever the
      # shrinkage: meeiip.
      value <- (sum(diagonals) + sum(mean_diagonals)) / (dof + (n + k) * d + 2)
      matrix(value, k, d)
    },
    VII = {
      value <- (scale + rowSums(diagonals) + rowSums(mean_diagonals)) /
        (dof + weights * d + 2 + bump * d)
      matrix(value, k, d)
    },
    EEI = {
      value <- (scale + colSums(diagonals) + colSums(mean_diagonals)) /
        (dof + n + 1 + bump)
      matrix(value, k, d, byrow = TRUE)
    },
    VVI = (scale + diagonals + mean_diagonals) / (dof + weights + 2 + bump),
    EVI = {
      raw <- scale + diagonals + mean_diagonals
      geometric <- exp(rowMeans(log(raw)))
      volume <- sum(geometric) / (dof + n + 1 + bump)
      volume * raw / geometric
    },
    VEI = .multilpa_prior_vei(scale + diagonals + mean_diagonals,
                              dof + weights * d + 2 + bump, max_iter, tol),
    EEE = {
      pooled <- (scale + Reduce(`+`, scatter) + Reduce(`+`, mean_terms)) /
        (dof + n + d + 1 + bump * k)
      rep(list(pooled), k)
    },
    VVV = lapply(seq_len(k), function(profile) {
      (scale + scatter[[profile]] + mean_terms[[profile]]) /
        (dof + weights[profile] + d + 1 + bump)
    }),
    EEV = {
      decompositions <- lapply(seq_len(k), function(profile) {
        eigen(scale + scatter[[profile]] + mean_terms[[profile]], symmetric = TRUE)
      })
      shape <- Reduce(`+`, lapply(decompositions, `[[`, "values"))
      volume <- exp(mean(log(shape))) / (dof + d + 1 + bump + n)
      lapply(decompositions, function(decomposition) {
        volume * decomposition$vectors %*%
          (normalise(shape) * t(decomposition$vectors))
      })
    },
    VEV = .multilpa_prior_vev(lapply(seq_len(k), function(profile) {
      scale + scatter[[profile]] + mean_terms[[profile]]
    }), weights, d, max_iter, tol))
  if (.multilpa_is_ellipsoidal(structure)) {
    covariances <- lapply(spread, function(block) {
      .multilpa_bound_covariance((block + t(block)) / 2, min_variance)
    })
    return(list(means = means,
                variances = matrix(vapply(covariances, diag, numeric(d)), k, d,
                                   byrow = TRUE),
                covariances = array(unlist(covariances, use.names = FALSE),
                                    c(d, d, k))))
  }
  list(means = means, variances = pmax(spread, min_variance))
}

#' The VEI fixed point with the prior's denominators (mclust's `meveip`)
#' @param raw Profiles-by-indicators prior-augmented scatter diagonals.
#' @param denominators One volume denominator per profile.
#' @param max_iter,tol The iteration cap and relative tolerance.
#' @return A matrix of variances.
#' @noRd
.multilpa_prior_vei <- function(raw, denominators, max_iter, tol) {
  d <- ncol(raw)
  shape <- rep(1, d)
  volumes <- rep(1, nrow(raw))
  converged <- FALSE
  # A fixed point: each half maximises given the other and each step needs the
  # one before it, which no apply function expresses. The cap bounds it.
  for (iteration in seq_len(max_iter)) {
    previous_volumes <- volumes
    previous_shape <- shape
    volumes <- rowSums(sweep(raw, 2L, shape, "/")) / denominators
    updated <- colSums(raw / volumes)
    shape <- updated / exp(mean(log(updated)))
    change <- max(abs(previous_volumes - volumes) / (1 + volumes),
                  abs(previous_shape - shape) / (1 + shape))
    if (change <= tol) {
      converged <- TRUE
      break
    }
  }
  if (!converged) .multilpa_warn_structure("VEI", max_iter)
  outer(volumes, shape)
}

#' The VEV fixed point with the prior-augmented scatter (mclust's `mevevp`)
#' @param blocks Prior-augmented scatter matrices, one per profile.
#' @param weights The effective membership of each profile.
#' @param d The dimension.
#' @param max_iter,tol The iteration cap and relative tolerance.
#' @return A list of covariance matrices.
#' @noRd
.multilpa_prior_vev <- function(blocks, weights, d, max_iter, tol) {
  decompositions <- lapply(blocks, function(block) eigen(block, symmetric = TRUE))
  eigenvalues <- lapply(decompositions, `[[`, "values")
  total <- Reduce(`+`, eigenvalues)
  geometric <- exp(mean(log(total)))
  shape <- total / geometric
  volumes <- geometric / weights
  converged <- FALSE
  # The same fixed point as the maximum-likelihood VEV, from mclust's start.
  for (iteration in seq_len(max_iter)) {
    previous_volumes <- volumes
    previous_shape <- shape
    volumes <- vapply(seq_along(weights), function(profile) {
      sum(eigenvalues[[profile]] / previous_shape) / (weights[profile] * d)
    }, numeric(1))
    updated <- Reduce(`+`, lapply(seq_along(weights), function(profile) {
      eigenvalues[[profile]] / volumes[profile]
    }))
    shape <- updated / exp(mean(log(updated)))
    change <- max(abs(previous_shape - shape) / (1 + shape),
                  abs(volumes - previous_volumes) / (1 + volumes))
    if (change <= tol) {
      converged <- TRUE
      break
    }
  }
  if (!converged) .multilpa_warn_structure("VEV", max_iter)
  lapply(seq_along(weights), function(profile) {
    vectors <- decompositions[[profile]]$vectors
    volumes[profile] * vectors %*% (shape * t(vectors))
  })
}

#' Refuse a prior this package, like mclust, does not define
#'
#' @param prior `NULL` or a `latents_prior`.
#' @param structure The resolved mclust model code.
#' @param categorical,fixed The corresponding [multilpa()] arguments.
#' @param n_covariates How many membership covariates were named.
#' @return `NULL`, invisibly; raises `latents_bad_argument` for a `prior` that
#'   is not a `latents_prior` and `latents_unsupported_prior` for a combination
#'   the prior is not defined for.
#' @noRd
.multilpa_check_prior <- function(prior, structure, categorical, fixed,
                                  n_covariates) {
  if (is.null(prior)) return(invisible(NULL))
  if (!inherits(prior, "latents_prior")) {
    stop(errorCondition(
      "`prior` must be NULL or the result of prior_control().",
      class = "latents_bad_argument", call = NULL))
  }
  refuse <- function(message) {
    stop(errorCondition(message, class = "latents_unsupported_prior", call = NULL))
  }
  if (!structure %in% .multilpa_prior_structures()) {
    refuse(sprintf(paste(
      "mclust defines no conjugate prior for the %s covariance structure, and",
      "neither does this package. A prior is available for %s."),
      structure, paste(.multilpa_prior_structures(), collapse = ", ")))
  }
  if (length(categorical) > 0L) {
    refuse("`prior` is a prior on Gaussian profiles; it cannot be combined with `categorical` indicators.")
  }
  if (length(fixed) > 0L) {
    refuse("`prior` cannot be combined with `fixed`: a held block is not estimated, so it has no posterior.")
  }
  if (n_covariates > 0L) {
    refuse("`prior` cannot be combined with `profile_covariates` or `group_covariates`.")
  }
  invisible(NULL)
}
