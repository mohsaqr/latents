# The constrained diagonal covariance family.
#
# Each profile's covariance decomposes as Sigma_k = lambda_k * A_k, where
# lambda_k = |Sigma_k|^(1/d) is the volume and A_k is the shape, a diagonal
# matrix with determinant one. Constraining the two separately gives the family
# mclust names with three letters -- volume, shape, orientation -- of which the
# axis-parallel ones end in "I".
#
# This package's `variance_model` ties the two together: "equal" constrains both
# (EEI) and "varying" frees both (VVI). Between them sit the two models that
# constrain one and not the other, and below them the two spherical ones.
#
# Celeux, G., & Govaert, G. (1995). Gaussian parsimonious clustering models.
# Pattern Recognition, 28, 781-793, gives the maximum-likelihood estimates used
# here, and Table 3 of that paper is where the iterative case comes from.

#' The covariance structures this package can fit
#'
#' @return A character vector of mclust model codes.
#' @noRd
.multilpa_structures <- function() {
  c("EII", "VII", "EEI", "VEI", "EVI", "VVI",
    "EEE", "VEE", "EVE", "VVE", "EEV", "VEV", "EVV", "VVV")
}

#' Is a structure ellipsoidal rather than axis-parallel?
#'
#' The third letter is the orientation: `I` for axis-parallel, which is a
#' diagonal covariance, and `E` or `V` for an orientation shared or free across
#' profiles, which is not.
#'
#' @param structure An mclust model code.
#' @return A single logical.
#' @noRd
.multilpa_is_ellipsoidal <- function(structure) {
  !is.null(structure) && !endsWith(structure, "I")
}

#' The structures whose free parameters the inference code covers
#'
#' All fourteen. The four separable structures are charted by log variances or
#' log-Cholesky coordinates; the ten that tie profiles together by log volumes,
#' determinant-one log shapes and Cayley orientations (R/structure-inference.R).
#'
#' @return A character vector of mclust model codes.
#' @noRd
.multilpa_inferable_structures <- function() .multilpa_structures()

#' Resolve the requested covariance structure
#'
#' `variance_model` and `covariance_model` name the four structures this
#' package has always fitted. `volume` and `shape` refine the diagonal case;
#' naming either one takes precedence, and the pair that results must be one of
#' the six axis-parallel models.
#'
#' @param variance_model `"varying"` or `"equal"`.
#' @param covariance_model `"diagonal"` or `"full"`.
#' @param volume `"varying"`, `"equal"` or `NULL` to follow `variance_model`.
#' @param shape `"varying"`, `"equal"`, `"spherical"` or `NULL`.
#' @return A single mclust model code.
#' @noRd
.multilpa_resolve_structure <- function(variance_model, covariance_model,
                                        volume = NULL, shape = NULL,
                                        orientation = NULL) {
  stopifnot(
    "`variance_model` must be \"varying\" or \"equal\"" =
      length(variance_model) == 1L && variance_model %in% c("varying", "equal"),
    "`covariance_model` must be \"diagonal\" or \"full\"" =
      length(covariance_model) == 1L &&
      covariance_model %in% c("diagonal", "full"),
    "`volume` must be \"varying\", \"equal\" or NULL" =
      is.null(volume) || (length(volume) == 1L &&
                            volume %in% c("varying", "equal")),
    "`shape` must be \"varying\", \"equal\", \"spherical\" or NULL" =
      is.null(shape) || (length(shape) == 1L &&
                           shape %in% c("varying", "equal", "spherical")),
    "`orientation` must be \"varying\", \"equal\", \"axis\" or NULL" =
      is.null(orientation) || (length(orientation) == 1L &&
                                 orientation %in% c("varying", "equal", "axis")))
  refined <- !is.null(volume) || !is.null(shape) || !is.null(orientation)
  if (!refined) {
    return(switch(covariance_model,
                  diagonal = if (variance_model == "equal") "EEI" else "VVI",
                  full = if (variance_model == "equal") "EEE" else "VVV"))
  }
  volume <- volume %||% variance_model
  shape <- shape %||% variance_model
  orientation <- orientation %||%
    if (identical(covariance_model, "full")) variance_model else "axis"
  # A spherical shape is a multiple of the identity, which has no orientation
  # to constrain, so those two models carry no third letter of their own.
  if (identical(shape, "spherical")) {
    if (!identical(orientation, "axis")) {
      stop(errorCondition(paste(
        "A spherical shape is a multiple of the identity, so it has no",
        "orientation to constrain. Drop `orientation`, or name a shape that",
        "is not spherical."),
        class = "latents_bad_argument", call = NULL))
    }
    return(if (identical(volume, "equal")) "EII" else "VII")
  }
  if (identical(orientation, "axis") && identical(covariance_model, "full")) {
    stop(errorCondition(paste(
      "`orientation = \"axis\"` is a diagonal covariance and",
      "`covariance_model = \"full\"` is not. Name one of them."),
      class = "latents_bad_argument", call = NULL))
  }
  letter <- function(value) if (identical(value, "equal")) "E" else "V"
  paste0(letter(volume), letter(shape),
         if (identical(orientation, "axis")) "I" else letter(orientation))
}

#' How many free parameters a covariance structure has
#'
#' The shape matrices carry a determinant-one constraint, so a free shape costs
#' `d - 1` and not `d`; the volume it was divided out of carries the missing
#' one. Getting this wrong would move every information criterion.
#'
#' @param structure An mclust model code.
#' @param n_profiles How many profiles.
#' @param n_continuous How many continuous indicators.
#' @return A single integer.
#' @noRd
.multilpa_structure_parameters <- function(structure, n_profiles, n_continuous) {
  stopifnot("`structure` must be a supported covariance structure" =
              length(structure) == 1L && structure %in% .multilpa_structures())
  k <- n_profiles
  d <- n_continuous
  if (d == 0L) return(0L)
  # Volume, then shape, then orientation. A determinant-one shape costs
  # `d - 1` and not `d`; an orthogonal orientation costs `d(d - 1)/2`.
  volume <- if (startsWith(structure, "E")) 1L else k
  shape <- if (substr(structure, 2L, 2L) == "E") (d - 1L) else k * (d - 1L)
  # `switch()` is avoided here: a branch named `E` partially matches its own
  # `EXPR` argument, which `R CMD check` reports and a reader has to think
  # about.
  third <- substr(structure, 3L, 3L)
  orientation <- if (third == "I") 0L else if (third == "E") {
    d * (d - 1L) / 2L
  } else k * d * (d - 1L) / 2L
  if (structure %in% c("EII", "VII")) {
    return(as.integer(if (structure == "EII") 1L else k))
  }
  as.integer(volume + shape + orientation)
}

#' Maximize the diagonal covariances under a constrained structure
#'
#' @param diagonals One row per profile and column per indicator of the
#'   posterior-weighted squared deviations the M-step accumulates.
#' @param weights The effective membership of each profile.
#' @param structure An axis-parallel mclust model code.
#' @param min_variance The variance lower bound.
#' @param max_iter Iterations allowed for the one structure that needs them.
#' @param tol Relative convergence tolerance for that iteration.
#' @return A matrix of variances, one row per profile and column per indicator.
#' @noRd
.multilpa_structure_variances <- function(diagonals, weights, structure,
                                          min_variance, max_iter = 1000L,
                                          tol = sqrt(.Machine$double.eps)) {
  stopifnot(
    "`diagonals` must be one row per profile" =
      is.matrix(diagonals) && nrow(diagonals) == length(weights),
    "`structure` must be a supported covariance structure" =
      length(structure) == 1L && structure %in% .multilpa_structures())
  d <- ncol(diagonals)
  n <- sum(weights)
  tied <- structure %in% c("EVI", "VEI")
  if (tied && any(diagonals / weights <= min_variance)) {
    return(.multilpa_bounded_variances(diagonals, weights, structure,
                                      min_variance))
  }
  variances <- switch(structure,
    VVI = diagonals / weights,
    EEI = matrix(colSums(diagonals) / n, length(weights), d, byrow = TRUE),
    EII = matrix(sum(diagonals) / (n * d), length(weights), d),
    VII = matrix(rowSums(diagonals) / (weights * d), length(weights), d),
    EVI = .multilpa_evi_variances(diagonals, n, min_variance),
    VEI = .multilpa_vei_variances(diagonals, weights, min_variance, max_iter, tol))
  if (tied && any(variances < min_variance)) {
    return(.multilpa_bounded_variances(diagonals, weights, structure,
                                      min_variance))
  }
  pmax(variances, min_variance)
}

#' Equal volume, varying shape: a closed form
#'
#' Sigma_k = lambda * A_k with |A_k| = 1, so A_k is the profile's own scatter
#' rescaled to determinant one and lambda is the average of the determinants
#' those rescalings divided out. Celeux and Govaert (1995), model `[lambda B_k]`.
#'
#' @param diagonals One row of scatter diagonals per profile.
#' @param n The total effective membership.
#' @param min_variance The variance lower bound. It is not applied to the
#'   scatter: the caller sends only scatter whose every `diagonals / weights`
#'   exceeds it, so the diagonals are positive, and a profile with effective
#'   membership below one can have a diagonal below the bound while its
#'   variance is above it. Flooring that diagonal moved the maximizer.
#' @return A matrix of variances.
#' @noRd
.multilpa_evi_variances <- function(diagonals, n, min_variance) {
  d <- ncol(diagonals)
  bounded <- pmax(diagonals, .Machine$double.xmin)
  # In logs: a product of d terms underflows long before any one of them does.
  log_determinant <- rowSums(log(bounded)) / d
  shapes <- bounded / exp(log_determinant)
  volume <- sum(exp(log_determinant)) / n
  shapes * volume
}

#' Varying volume, equal shape: a bounded fixed point
#'
#' Sigma_k = lambda_k * A with |A| = 1. Each half has a closed form given the
#' other and neither has one alone, so the two are alternated. Celeux and
#' Govaert (1995), model `[lambda_k B]`.
#'
#' @param diagonals One row of scatter diagonals per profile.
#' @param weights The effective membership of each profile.
#' @param min_variance The variance lower bound. The iteration does not apply
#'   it: the caller checks the result against it and solves the bounded
#'   problem when it binds.
#' @param max_iter The iteration cap.
#' @param tol Relative convergence tolerance on the shape.
#' @return A matrix of variances.
#' @noRd
.multilpa_vei_variances <- function(diagonals, weights, min_variance, max_iter,
                                    tol) {
  d <- ncol(diagonals)
  # Positive by the caller's guard; see `.multilpa_evi_variances()`.
  bounded <- pmax(diagonals, .Machine$double.xmin)
  normalise <- function(shape) shape / exp(sum(log(shape)) / d)
  # The shape sums `colSums(W_k / lambda_k)` have no units (a count times a
  # determinant-one shape), so the variance bound, in squared data units,
  # cannot be compared with them: with indicators whose scales differ by 1e9
  # it bound here while every variance was far above it, and moved the fixed
  # point (review 2026-10-26 round 2, K2-1). They only need to stay positive.
  shape <- normalise(pmax(colSums(bounded), .Machine$double.xmin))
  converged <- FALSE
  # A fixed point: each step is a maximiser given the other half, so the
  # objective cannot decrease, and each step needs the one before it. That is
  # what a loop is, and no apply function expresses it; the cap is here because
  # "cannot decrease" is not "arrives in finite time".
  for (iteration in seq_len(max_iter)) {
    volumes <- rowSums(sweep(bounded, 2L, shape, "/")) / (weights * d)
    updated <- normalise(pmax(colSums(bounded / volumes), .Machine$double.xmin))
    converged <- max(abs(updated - shape)) <= tol * max(1, max(abs(shape)))
    shape <- updated
    if (converged) break
  }
  # Counted and reported once per fit (`latents_no_converge`), as every other
  # structure's capped M-step is; a warning here fired once per EM iteration.
  if (!converged) .multilpa_warn_structure("VEI", max_iter)
  volumes <- rowSums(sweep(bounded, 2L, shape, "/")) / (weights * d)
  outer(volumes, shape)
}

#' The scatter diagonals of a list of per-profile scatter matrices
#'
#' The observed-data M-step accumulates full matrices even for a diagonal
#' model, because a missing pair contributes an off-diagonal adjustment. Only
#' the diagonal identifies an axis-parallel structure.
#'
#' @param scatter A list of per-profile scatter matrices.
#' @return One row per profile and column per indicator.
#' @noRd
.multilpa_scatter_diagonals <- function(scatter) {
  d <- nrow(scatter[[1L]])
  # `vapply()` gives a d-by-k matrix, or a bare vector when d is one, so the
  # shape is stated rather than inferred: one row per profile either way.
  matrix(vapply(scatter, diag, numeric(d)), nrow = length(scatter), ncol = d,
         byrow = TRUE)
}

#' Refuse a standard error the coordinates cannot express
#'
#' Every structure this package fits has a chart, so this refuses only a fit
#' recording a structure code it does not know, such as one made by a newer
#' version.
#'
#' @param x A fitted model of this package.
#' @return `NULL`, invisibly, when the structure is one that is covered.
#' @noRd
.multilpa_check_structure_inference <- function(x) {
  structure <- x$covariance_structure
  if (is.null(structure) || structure %in% .multilpa_inferable_structures()) {
    return(invisible(NULL))
  }
  stop(errorCondition(sprintf(paste(
    "Wald standard errors are not available for the %s covariance structure,",
    "which this version of the package does not know. Use",
    "`method = \"bootstrap\"`, which refits inside the recorded family."),
    structure),
    class = "latents_unsupported_inference", call = NULL))
}

#' Put starting variances inside the family they will be maximized in
#'
#' EM only promises a likelihood that does not decrease *within* the family it
#' is maximizing over. A start taken from a wider structure -- a free diagonal
#' handed to an equal-volume model, which is what `start =` does when the two
#' fits differ -- sits outside it, and the first M-step projects onto the
#' family, which can lower the likelihood and trip the monotonicity guard.
#' Projecting the start first means every iteration, including the first, is a
#' step inside the family.
#'
#' @param parameters A starting parameter list.
#' @param structure An mclust model code, or `NULL` to leave it alone.
#' @param n_observations How many observations the weights should sum to.
#' @param min_variance The variance lower bound.
#' @return The parameter list, with `variances` inside the structure.
#' @noRd
.multilpa_project_start <- function(parameters, structure, n_observations,
                                    min_variance) {
  if (is.null(structure) || structure %in% .multilpa_separable_structures()) {
    return(parameters)
  }
  variances <- parameters$variances
  if (is.null(variances) || ncol(variances) == 0L) return(parameters)
  # The prevalences say how much each profile is expected to weigh, which is
  # what turns a variance back into the scatter the estimator is written for.
  prevalence <- if (is.null(parameters$profile_probabilities)) {
    rep(1 / nrow(variances), nrow(variances))
  } else colMeans(as.matrix(parameters$profile_probabilities))
  weights <- pmax(prevalence * n_observations, .Machine$double.eps)
  if (.multilpa_is_ellipsoidal(structure)) {
    blocks <- lapply(seq_along(weights), function(profile) {
      block <- if (is.null(parameters$covariances)) {
        diag(variances[profile, ], ncol(variances))
      } else {
        matrix(parameters$covariances[, , profile], ncol(variances))
      }
      block * weights[profile]
    })
    projected <- .multilpa_structure_covariances(blocks, weights, structure,
                                                 min_variance)
    parameters$covariances <- projected
    parameters$variances <- t(vapply(seq_along(weights), function(profile) {
      diag(matrix(projected[, , profile], ncol(variances)))
    }, numeric(ncol(variances))))
    return(parameters)
  }
  parameters$variances <- .multilpa_structure_variances(
    variances * weights, weights, structure, min_variance)
  parameters
}

# ---------------------------------------------------------------------------
# The ellipsoidal family.
#
# Sigma_k = lambda_k * D_k * A_k * D_k', a volume, an orientation (orthogonal
# eigenvectors) and a shape (diagonal, determinant one). Constraining the three
# separately gives the eight models whose third letter is E or V.

#' The determinant of a scatter matrix, raised to the power 1/d
#'
#' Computed in logs: a determinant of a d-by-d scatter underflows long before
#' any of its entries does.
#'
#' @param scatter A symmetric matrix.
#' @param d Its dimension.
#' @return A positive scalar.
#' @noRd
.multilpa_scaled_determinant <- function(scatter, d) {
  value <- determinant(scatter, logarithm = TRUE)
  exp(as.numeric(value$modulus) * as.numeric(value$sign) / d)
}

#' Maximize the full covariances under a constrained ellipsoidal structure
#'
#' @param scatter A list of per-profile scatter matrices, posterior weighted.
#' @param weights The effective membership of each profile.
#' @param structure An mclust model code.
#' @param min_variance The eigenvalue lower bound.
#' @param max_iter Iterations allowed for the structures that need them.
#' @param tol Relative convergence tolerance for those.
#' @param inner_steps `NULL` to solve EVE and VVE's orientation to `tol`, or a
#'   number of minorize-maximize steps to take from the warm start. The EM
#'   loop passes a small number: a generalized EM step that improves the
#'   M-step objective without solving it, which keeps every EM iteration
#'   monotone and makes the EM map a fixed function of its input.
#' @return A `d` by `d` by `k` array of covariances.
#' @noRd
.multilpa_structure_covariances <- function(scatter, weights, structure,
                                            min_variance, max_iter = 1000L,
                                            tol = sqrt(.Machine$double.eps),
                                            start = NULL, inner_steps = NULL) {
  stopifnot(
    "`structure` must be a supported covariance structure" =
      length(structure) == 1L && structure %in% .multilpa_structures(),
    "`scatter` must be one matrix per profile" =
      is.list(scatter) && length(scatter) == length(weights))
  d <- nrow(scatter[[1L]])
  k <- length(weights)
  n <- sum(weights)
  symmetric <- lapply(scatter, function(block) (block + t(block)) / 2)
  tied <- structure %in% c("EEV", "VEV", "EVV", "VEE", "EVE", "VVE")
  if (tied && any(vapply(seq_len(k), function(profile) {
    min(eigen(symmetric[[profile]], symmetric = TRUE, only.values = TRUE)$values) /
      weights[profile] <= min_variance
  }, logical(1)))) {
    return(.multilpa_bounded_covariances(symmetric, weights, structure,
      min_variance, start, inner_steps %||% max_iter, tol))
  }
  covariances <- switch(structure,
    EEE = rep(list(Reduce(`+`, symmetric) / n), k),
    VVV = lapply(seq_len(k), function(profile) symmetric[[profile]] / weights[profile]),
    EEV = .multilpa_eev_covariances(symmetric, n, d),
    VEV = .multilpa_vev_covariances(symmetric, weights, d, max_iter, tol),
    EVV = .multilpa_evv_covariances(symmetric, n, d),
    VEE = .multilpa_vee_covariances(symmetric, weights, d, max_iter, tol,
                                    .multilpa_seed_block(start, d)),
    EVE = .multilpa_eve_covariances(symmetric, weights, d,
                                    inner_steps %||% max_iter, tol, TRUE,
                                    .multilpa_seed_block(start, d),
                                    capped_by_design = !is.null(inner_steps)),
    VVE = .multilpa_eve_covariances(symmetric, weights, d,
                                    inner_steps %||% max_iter, tol, FALSE,
                                    .multilpa_seed_block(start, d),
                                    capped_by_design = !is.null(inner_steps)),
    # The axis-parallel models are diagonal, and are solved as variances.
    lapply(seq_len(k), function(profile) {
      diag(.multilpa_structure_variances(
        .multilpa_scatter_diagonals(symmetric), weights, structure,
        min_variance, max_iter, tol)[profile, ], d)
    }))
  if (tied && any(vapply(covariances, function(block) {
    min(eigen(block, symmetric = TRUE, only.values = TRUE)$values) < min_variance
  }, logical(1)))) {
    return(.multilpa_bounded_covariances(symmetric, weights, structure,
      min_variance, start, inner_steps %||% max_iter, tol))
  }
  array(unlist(lapply(covariances, function(block) {
    .multilpa_bound_covariance(block, min_variance)
  }), use.names = FALSE), c(d, d, k))
}

#' Equal volume and shape, varying orientation
#'
#' Each profile keeps its own eigenvectors; the shared shape is the normalized
#' sum of the eigenvalues and the volume is what that normalization divided out.
#' Celeux and Govaert (1995), model `[lambda D_k A D_k']`.
#'
#' @param scatter Symmetrized per-profile scatter matrices.
#' @param n The total effective membership.
#' @param d The dimension.
#' @return A list of covariance matrices.
#' @noRd
.multilpa_eev_covariances <- function(scatter, n, d) {
  decompositions <- lapply(scatter, function(block) eigen(block, symmetric = TRUE))
  eigenvalues <- Reduce(`+`, lapply(decompositions, `[[`, "values"))
  shape <- eigenvalues / exp(sum(log(eigenvalues)) / d)
  volume <- exp(sum(log(eigenvalues)) / d) / n
  lapply(decompositions, function(decomposition) {
    volume * decomposition$vectors %*% diag(shape, d) %*% t(decomposition$vectors)
  })
}

#' Varying volume, equal shape, varying orientation
#'
#' The orientations are each profile's own, so only the volume and the shared
#' shape alternate --- the same fixed point as the diagonal equal-shape model,
#' one dimension up. Celeux and Govaert (1995), model `[lambda_k D_k A D_k']`.
#'
#' @param scatter Symmetrized per-profile scatter matrices.
#' @param weights The effective membership of each profile.
#' @param d The dimension.
#' @param max_iter The iteration cap.
#' @param tol Relative convergence tolerance on the shape.
#' @return A list of covariance matrices.
#' @noRd
.multilpa_vev_covariances <- function(scatter, weights, d, max_iter, tol) {
  decompositions <- lapply(scatter, function(block) eigen(block, symmetric = TRUE))
  eigenvalues <- lapply(decompositions, `[[`, "values")
  normalise <- function(values) values / exp(sum(log(values)) / d)
  shape <- normalise(Reduce(`+`, eigenvalues))
  converged <- FALSE
  # Each half maximises given the other, and each step needs the one before it:
  # a fixed point, which no apply function expresses.
  for (iteration in seq_len(max_iter)) {
    volumes <- vapply(seq_along(weights), function(profile) {
      sum(eigenvalues[[profile]] / shape) / (weights[profile] * d)
    }, numeric(1))
    updated <- normalise(Reduce(`+`, lapply(seq_along(weights), function(profile) {
      eigenvalues[[profile]] / volumes[profile]
    })))
    converged <- max(abs(updated - shape)) <= tol * max(1, max(abs(shape)))
    shape <- updated
    if (converged) break
  }
  if (!converged) .multilpa_warn_structure("VEV", max_iter)
  volumes <- vapply(seq_along(weights), function(profile) {
    sum(eigenvalues[[profile]] / shape) / (weights[profile] * d)
  }, numeric(1))
  lapply(seq_along(weights), function(profile) {
    vectors <- decompositions[[profile]]$vectors
    volumes[profile] * vectors %*% diag(shape, d) %*% t(vectors)
  })
}

#' Equal volume, varying shape and orientation
#'
#' Each profile keeps its whole scatter, rescaled to determinant one; the
#' shared volume averages what those rescalings divided out. No eigenvectors
#' are needed. Celeux and Govaert (1995), model `[lambda D_k A_k D_k']`.
#'
#' @param scatter Symmetrized per-profile scatter matrices.
#' @param n The total effective membership.
#' @param d The dimension.
#' @return A list of covariance matrices.
#' @noRd
.multilpa_evv_covariances <- function(scatter, n, d) {
  determinants <- vapply(scatter, .multilpa_scaled_determinant, numeric(1), d = d)
  volume <- sum(determinants) / n
  lapply(seq_along(scatter), function(profile) {
    volume * scatter[[profile]] / determinants[profile]
  })
}

#' Varying volume, equal shape and orientation
#'
#' `Sigma_k = lambda_k * C` with `|C| = 1`. The shared matrix and the volumes
#' alternate, exactly as the diagonal equal-shape model does, except that `C`
#' is a full matrix and needs no orthogonal optimization.
#' Celeux and Govaert (1995), model `[lambda_k C]`.
#'
#' @param scatter Symmetrized per-profile scatter matrices.
#' @param weights The effective membership of each profile.
#' @param d The dimension.
#' @param max_iter The iteration cap.
#' @param tol Relative convergence tolerance on the shared matrix.
#' @return A list of covariance matrices.
#' @noRd
.multilpa_vee_covariances <- function(scatter, weights, d, max_iter, tol,
                                      start = NULL) {
  normalise <- function(block) block / .multilpa_scaled_determinant(block, d)
  common <- normalise(start %||% Reduce(`+`, scatter))
  converged <- FALSE
  # A fixed point again: the volumes and the shared matrix each maximise given
  # the other.
  for (iteration in seq_len(max_iter)) {
    inverse <- chol2inv(chol(common))
    volumes <- vapply(seq_along(weights), function(profile) {
      sum(inverse * scatter[[profile]]) / (weights[profile] * d)
    }, numeric(1))
    updated <- normalise(Reduce(`+`, lapply(seq_along(weights), function(profile) {
      scatter[[profile]] / volumes[profile]
    })))
    converged <- max(abs(updated - common)) <= tol * max(1, max(abs(common)))
    common <- updated
    if (converged) break
  }
  if (!converged) .multilpa_warn_structure("VEE", max_iter)
  inverse <- chol2inv(chol(common))
  volumes <- vapply(seq_along(weights), function(profile) {
    sum(inverse * scatter[[profile]]) / (weights[profile] * d)
  }, numeric(1))
  lapply(seq_along(weights), function(profile) volumes[profile] * common)
}

#' Common orientation with varying shape, at equal or varying volume
#'
#' `Sigma_k = lambda_k * D * A_k * D'` with one orthogonal `D` shared by every
#' profile. Unlike every other model here, `D` has no closed form: it minimizes
#' `sum_k tr(D A_k^-1 D' W_k) / lambda_k` over the orthogonal matrices. The
#' minorize-maximize step of Browne and McNicholas (2014) is used, which
#' majorizes each term by a linear one and solves the resulting orthogonal
#' Procrustes problem with a singular value decomposition.
#'
#' @param scatter Symmetrized per-profile scatter matrices.
#' @param weights The effective membership of each profile.
#' @param d The dimension.
#' @param max_iter The iteration cap.
#' @param tol Relative convergence tolerance on the orientation.
#' @param equal_volume `TRUE` for EVE, `FALSE` for VVE.
#' @return A list of covariance matrices.
#' @noRd
.multilpa_eve_covariances <- function(scatter, weights, d, max_iter, tol,
                                      equal_volume, start = NULL,
                                      capped_by_design = FALSE) {
  k <- length(weights)
  n <- sum(weights)
  largest <- vapply(scatter, function(block) {
    max(eigen(block, symmetric = TRUE, only.values = TRUE)$values)
  }, numeric(1))
  # Warm start: every profile shares the orientation, so the eigenvectors of
  # the previous iteration's first covariance are it. Restarting from the
  # pooled scatter each time makes the minorize-maximize step redo from far
  # away the work the last iteration already did.
  orientation <- eigen(start %||% Reduce(`+`, scatter), symmetric = TRUE)$vectors
  shapes <- lapply(seq_len(k), function(profile) rep(1, d))
  volumes <- rep(1, k)
  converged <- FALSE
  # The shapes and volumes that maximize the objective at a given orientation.
  shape_and_volume <- function(orientation) {
    rotated <- lapply(scatter, function(block) {
      diag(crossprod(orientation, block %*% orientation))
    })
    shapes <- lapply(seq_len(k), function(profile) {
      # Positive only: an absolute `.Machine$double.eps` is in squared data
      # units times counts, and bound for data in units of 1e-9 (review
      # 2026-10-26 round 2, K2-4).
      values <- pmax(rotated[[profile]], .Machine$double.xmin)
      values / exp(sum(log(values)) / d)
    })
    volumes <- if (equal_volume) {
      rep(sum(vapply(seq_len(k), function(profile) {
        .multilpa_scaled_determinant(diag(rotated[[profile]], d), d)
      }, numeric(1))) / n, k)
    } else {
      vapply(seq_len(k), function(profile) {
        sum(rotated[[profile]] / shapes[[profile]]) / (weights[profile] * d)
      }, numeric(1))
    }
    list(shapes = shapes, volumes = volumes)
  }
  # Two nested fixed points: the shapes and volumes given the orientation, and
  # the orientation given them. Each step needs the previous one.
  for (iteration in seq_len(max_iter)) {
    conditional <- shape_and_volume(orientation)
    shapes <- conditional$shapes
    volumes <- conditional$volumes
    # One MM step: majorize each trace term linearly at the current
    # orientation, then maximize the linear surrogate over the orthogonal
    # matrices, which is an orthogonal Procrustes problem.
    surrogate <- Reduce(`+`, lapply(seq_len(k), function(profile) {
      (largest[profile] * orientation - scatter[[profile]] %*% orientation) %*%
        diag(1 / (volumes[profile] * shapes[[profile]]), d)
    }))
    decomposition <- svd(surrogate)
    updated <- decomposition$u %*% t(decomposition$v)
    converged <- max(abs(abs(updated) - abs(orientation))) <= tol
    orientation <- updated
    if (converged) break
  }
  if (!converged && !capped_by_design) {
    .multilpa_warn_structure(if (equal_volume) "EVE" else "VVE", max_iter)
  }
  # The shapes and volumes that go with the final orientation: an exact
  # conditional maximization, so it can only raise the objective.
  conditional <- shape_and_volume(orientation)
  shapes <- conditional$shapes
  volumes <- conditional$volumes
  lapply(seq_len(k), function(profile) {
    volumes[profile] * orientation %*% diag(shapes[[profile]], d) %*%
      t(orientation)
  })
}

# An M-step runs once per EM iteration, hundreds of times in a fit, so a
# warning raised from inside one would fire hundreds of times for a single
# problem. The occurrences are counted here and reported once, with their
# number, when the fit that produced them finishes.
.multilpa_structure_log <- new.env(parent = emptyenv())

#' Start counting M-steps that stopped at their cap
#'
#' Called once per fit, so a count never carries from one fit into the next.
#'
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_reset_structure_log <- function() {
  .multilpa_structure_log$capped <- 0L
  .multilpa_structure_log$structure <- NA_character_
  .multilpa_structure_log$max_iter <- NA_integer_
  invisible(NULL)
}

#' Record that a structure's iteration stopped at its cap
#' @param structure The model code.
#' @param max_iter The cap it reached.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_warn_structure <- function(structure, max_iter) {
  .multilpa_structure_log$capped <- (.multilpa_structure_log$capped %||% 0L) + 1L
  .multilpa_structure_log$structure <- structure
  .multilpa_structure_log$max_iter <- max_iter
  invisible(NULL)
}

#' Report the M-steps that stopped at their cap, once
#'
#' The expectation-maximization loop keeps going after a partial M-step, and
#' usually converges anyway, so this qualifies the fit rather than rejecting
#' it. The count is the part worth reading: one M-step in three hundred is a
#' different thing from all of them.
#'
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_report_structure_log <- function() {
  capped <- .multilpa_structure_log$capped %||% 0L
  if (capped == 0L) return(invisible(NULL))
  warning(warningCondition(sprintf(paste(
    "%d M-step(s) of the %s covariance stopped at %d iterations without",
    "settling, so those steps returned the last value they reached rather",
    "than the maximizer. The fit continued from them; its own convergence is",
    "reported separately."),
    capped, .multilpa_structure_log$structure,
    .multilpa_structure_log$max_iter),
    class = "latents_no_converge"))
  invisible(NULL)
}

#' The `multilpa()` arguments that ask for a named structure
#'
#' `enumerate_classes()` grids over model codes, and `multilpa()` takes the
#' three pieces a code is made of. This turns one into the other, so the grid
#' and the fitter cannot disagree about what a code means.
#'
#' @param structure An mclust model code, or `NA_character_` to leave the
#'   caller's own arguments alone.
#' @return A named list of arguments, empty when nothing was asked for.
#' @noRd
.multilpa_structure_arguments <- function(structure) {
  if (is.na(structure)) return(list())
  stopifnot("`structure` must be a supported covariance structure" =
              structure %in% .multilpa_structures())
  spherical <- structure %in% c("EII", "VII")
  letter <- function(position) {
    if (substr(structure, position, position) == "E") "equal" else "varying"
  }
  arguments <- list(volume = letter(1L),
                    shape = if (spherical) "spherical" else letter(2L))
  if (spherical) return(arguments)
  third <- substr(structure, 3L, 3L)
  arguments$orientation <- if (third == "I") "axis" else
    if (third == "E") "equal" else "varying"
  arguments
}

#' One covariance slice of a previous iteration, as a seed
#'
#' The structures that iterate share one matrix across profiles, so any slice
#' of the previous answer carries it. `NULL` when there is nothing to start
#' from, which is the first iteration.
#'
#' @param start A `d` by `d` by `k` array, or `NULL`.
#' @param d The dimension.
#' @return A `d` by `d` matrix, or `NULL`.
#' @noRd
.multilpa_seed_block <- function(start, d) {
  if (is.null(start) || length(dim(start)) != 3L || dim(start)[1L] != d) {
    return(NULL)
  }
  block <- matrix(start[, , 1L], d, d)
  if (any(!is.finite(block))) return(NULL)
  (block + t(block)) / 2
}

#' Orientation steps per EM iteration for EVE and VVE
#'
#' Their shared orientation has no closed form. Solving it to convergence in
#' every M-step spent a median of 220 minorize-maximize steps per EM
#' iteration on work the next iteration discards; a few warm-started steps
#' per iteration reach the same maximum (a generalized EM step). Over 18
#' fits (three datasets, EVE and VVE, three seeds), 3, 5 and 10 steps landed
#' within 1e-5 to 3e-3 of the full solve's likelihood whenever they reached
#' the same maximum, and took 3.5 to 8 times less time; which maximum a fit
#' reached varied without pattern across settings, the full solve included.
#'
#' @return A single positive integer.
#' @noRd
.multilpa_orientation_steps <- function() 5L

#' Log-eigenvalue design for volume and shape constraints
#' @param k,d Numbers of profiles and dimensions.
#' @param structure Covariance structure code.
#' @return A design whose rows are profile-major log eigenvalues.
#' @noRd
.multilpa_eigenvalue_design <- function(k, d, structure) {
  stopifnot(k >= 1L, d >= 1L, structure %in% .multilpa_structures())
  equal_volume <- startsWith(structure, "E")
  equal_shape <- substr(structure, 2L, 2L) == "E"
  volume <- if (equal_volume) matrix(1, k * d, 1L) else
    diag(k)[rep(seq_len(k), each = d), , drop = FALSE]
  if (d == 1L) return(volume)
  shape <- rbind(diag(d - 1L), rep(-1, d - 1L))
  if (equal_shape) return(cbind(volume, shape[rep(seq_len(d), k), , drop = FALSE]))
  shape_design <- matrix(0, k * d, k * (d - 1L))
  invisible(lapply(seq_len(k), function(profile) {
    shape_design[(profile - 1L) * d + seq_len(d),
                 (profile - 1L) * (d - 1L) + seq_len(d - 1L)] <<- shape
  }))
  cbind(volume, shape_design)
}

#' Optimize floored variances inside a volume/shape family
#'
#' The Gaussian M-step objective is convex in log eigenvalues. The volume and
#' determinant-one shape coordinates make its equality constraints linear;
#' the eigenvalue floor is a linear inequality in those same coordinates.
#' Thus a bound never changes the requested family, unlike clipping individual
#' eigenvalues after an unrestricted update.
#' @param scatter Posterior-weighted scatter diagonals, profiles by dimensions.
#' @param weights Effective membership per profile.
#' @param structure Axis-parallel structure code.
#' @param min_variance Eigenvalue floor.
#' @return Feasible profile-by-dimension variances.
#' @noRd
.multilpa_bounded_variances <- function(scatter, weights, structure, min_variance) {
  stopifnot(is.matrix(scatter), nrow(scatter) == length(weights),
            all(weights > 0), all(is.finite(scatter)), min_variance > 0)
  k <- nrow(scatter)
  d <- ncol(scatter)
  scatter <- pmax(scatter, 0)
  if (k == 1L || structure == "VVI") return(pmax(scatter / weights, min_variance))
  if (structure == "EEI") return(matrix(pmax(colSums(scatter) / sum(weights),
    min_variance), k, d, byrow = TRUE))
  if (structure == "EII") return(matrix(max(sum(scatter) / (sum(weights) * d),
    min_variance), k, d))
  if (structure == "VII") return(matrix(pmax(rowSums(scatter) / (weights * d),
    min_variance), k, d))
  if (d == 1L) return(.multilpa_bounded_variances(scatter, weights,
    if (startsWith(structure, "E")) "EEI" else "VVI", min_variance))
  design <- .multilpa_eigenvalue_design(k, d, structure)
  lower <- log(min_variance)
  # A feasible spherical start lies inside every inequality. An active-set
  # Newton solve reaches the floor directly: an interior barrier can round a
  # positive slack to zero and then evaluate log(0) on its next outer step.
  n_volume <- if (startsWith(structure, "E")) 1L else k
  initial_volume <- if (n_volume == 1L) sum(scatter) / (sum(weights) * d) else
    rowSums(scatter) / (weights * d)
  point <- c(log(pmax(initial_volume, min_variance)) + 0.1,
             numeric(ncol(design) - n_volume))
  counts <- rep(weights, each = d)
  values <- as.vector(t(scatter))
  total <- sum(weights)
  dimension <- ncol(design)
  active <- integer()
  iteration <- 0L
  settled <- FALSE
  while (iteration < 500L && !settled) {
    iteration <- iteration + 1L
    logs <- as.vector(design %*% point)
    terms <- values * exp(-logs)
    slack <- logs - lower
    binding_tolerance <- 64 * .Machine$double.eps * dimension *
      max(1, abs(lower), max(abs(logs)))
    # Retaining a nonbinding row would solve the wrong equality-constrained
    # problem. Check its actual slack, rather than a difference of step lengths.
    if (length(active) > 0L && any(abs(slack[active]) > binding_tolerance)) {
      active <- active[abs(slack[active]) <= binding_tolerance]
    }
    gradient <- as.vector(crossprod(design, counts - terms)) / total
    hessian <- crossprod(design, design * terms) / total
    # Search along the null space of the currently active floor constraints.
    basis <- if (length(active) > 0L) {
      factor <- qr(t(design[active, , drop = FALSE]))
      qr.Q(factor, complete = TRUE)[, seq_len(dimension) > length(active), drop = FALSE]
    } else diag(dimension)
    projected <- as.vector(crossprod(basis, gradient))
    if (length(projected) == 0L ||
        max(abs(projected)) < 1e-11 * (1 + max(abs(gradient)))) {
      if (length(active) == 0L) {
        settled <- TRUE
        break
      }
      constraints <- design[active, , drop = FALSE]
      multipliers <- as.vector(solve(tcrossprod(constraints), constraints %*% gradient))
      # Certification includes primal feasibility and complementary slackness,
      # as well as stationarity and the multipliers' signs. Positive multipliers
      # alone cannot certify a constraint incorrectly frozen above its floor.
      kkt_tolerance <- 1e-11 * (1 + max(abs(gradient)))
      if (min(multipliers) >= -kkt_tolerance &&
          all(slack >= -binding_tolerance) &&
          max(abs(multipliers * slack[active])) <= kkt_tolerance) {
        settled <- TRUE
        break
      }
      release <- if (min(multipliers) < -kkt_tolerance) which.min(multipliers) else
        which.max(abs(multipliers * slack[active]))
      active <- active[-release]
      next
    }
    reduced <- crossprod(basis, hessian %*% basis)
    delta <- tryCatch({
      factor <- chol(reduced)
      -backsolve(factor, forwardsolve(t(factor), projected))
    }, error = function(condition) -projected)
    direction <- as.vector(basis %*% delta)
    velocity <- as.vector(design %*% direction)
    # A QR null-space product can leave rounding-sized velocities in rows
    # dependent on the active constraints. They are algebraically zero and
    # must not impose a zero step without adding an independent constraint.
    velocity[abs(velocity) <= 64 * .Machine$double.eps * dimension *
               max(abs(design)) * max(abs(direction))] <- 0
    limited <- setdiff(which(velocity < 0), active)
    bounds <- slack[limited] / (-velocity[limited])
    alpha <- if (length(bounds) > 0L) min(1, max(0, min(bounds))) else 1
    # Only the first limiting constraint is activated. An absolute tolerance
    # on step lengths can falsely activate a distant floor when alpha is tiny.
    # Simultaneous independent limits are added on subsequent zero-length steps.
    at_boundary <- if (length(bounds) > 0L && min(bounds) <= 1)
      limited[which.min(bounds)] else integer()
    # Evaluate the objective difference without subtracting two almost equal
    # totals. The linear term uses the gradient; expm1 retains the quadratic
    # remainder when a final Newton move is smaller than the objective's ulp.
    change <- function(alpha) {
      movement <- alpha * velocity
      sum(gradient * (alpha * direction)) +
        sum(terms * (expm1(-movement) + movement)) / total
    }
    candidate <- point + alpha * direction
    candidate_change <- change(alpha)
    halvings <- 0L
    # Only an actual decrease of the original objective is accepted; no
    # barrier objective or relaxed EM monotonicity threshold is involved.
    while ((!is.finite(candidate_change) || candidate_change > 0) && halvings < 60L) {
      alpha <- alpha / 2
      candidate <- point + alpha * direction
      candidate_change <- change(alpha)
      halvings <- halvings + 1L
      at_boundary <- integer()
    }
    point <- candidate
    if (length(at_boundary) > 0L) {
      possible <- unique(c(active, at_boundary))
      factor <- qr(t(design[possible, , drop = FALSE]))
      active <- possible[factor$pivot[seq_len(factor$rank)]]
    }
  }
  if (!settled) {
    stop(errorCondition(
      "The bounded covariance optimization did not satisfy its KKT conditions.",
      class = "latents_no_converge", call = NULL))
  }
  logs <- as.vector(design %*% point)
  stopifnot("the bounded chart must remain feasible" =
    all(is.finite(logs)), all(logs >= lower - 1e-12))
  matrix(exp(logs), k, d, byrow = TRUE)
}

#' Rotate a common orientation while retaining its eigenvalues
#' @param orientation Shared orthogonal eigenvectors.
#' @param scatter Posterior-weighted scatter matrices.
#' @param variances Profile-by-dimension eigenvalues.
#' @return Orientation after one exact sweep of pairwise rotations.
#' @noRd
.multilpa_bounded_orientation <- function(orientation, scatter, variances) {
  stopifnot(is.matrix(orientation), is.list(scatter), is.matrix(variances),
            length(scatter) == nrow(variances), all(variances > 0))
  d <- nrow(orientation)
  if (d == 1L) return(orientation)
  pairs <- which(upper.tri(diag(d)), arr.ind = TRUE)
  Reduce(function(current, pair) {
    a <- pairs[pair, 1L]
    b <- pairs[pair, 2L]
    pieces <- lapply(scatter, function(block) crossprod(current, block %*% current))
    precision_difference <- 1 / variances[, a] - 1 / variances[, b]
    cosine <- sum(precision_difference * vapply(pieces,
      function(block) (block[a, a] - block[b, b]) / 2, numeric(1)))
    sine <- sum(precision_difference * vapply(pieces,
      function(block) block[a, b], numeric(1)))
    if (cosine == 0 && sine == 0) return(current)
    angle <- atan2(-sine, -cosine) / 2
    rotation <- diag(d)
    rotation[c(a, b), c(a, b)] <- matrix(c(cos(angle), sin(angle),
      -sin(angle), cos(angle)), 2L)
    current %*% rotation
  }, seq_len(nrow(pairs)), init = orientation)
}

#' Covariance M-step with eigenvalue floors inside the ellipsoidal family
#' @param scatter,weights Posterior-weighted scatter and membership.
#' @param structure Covariance structure code.
#' @param min_variance Eigenvalue floor.
#' @param start Previous covariance array, or NULL for initial projection.
#' @param max_iter,tol Orientation iteration controls.
#' @return A dimension-by-dimension-by-profile covariance array.
#' @noRd
.multilpa_bounded_covariances <- function(scatter, weights, structure,
                                         min_variance, start, max_iter, tol) {
  stopifnot(is.list(scatter), length(scatter) == length(weights), all(weights > 0))
  k <- length(weights)
  d <- nrow(scatter[[1L]])
  common <- endsWith(structure, "E")
  eigen_structure <- paste0(substr(structure, 1L, 2L), "I")
  orientations <- if (common) {
    seed <- .multilpa_seed_block(start, d) %||% Reduce(`+`, scatter)
    rep(list(eigen(seed, symmetric = TRUE)$vectors), k)
  } else lapply(scatter, function(block) eigen(block, symmetric = TRUE)$vectors)
  solve_values <- function(rotations) {
    diagonals <- t(vapply(seq_len(k), function(profile) {
      diag(crossprod(rotations[[profile]], scatter[[profile]] %*% rotations[[profile]]))
    }, numeric(d))) |> matrix(k, d)
    .multilpa_bounded_variances(diagonals, weights, eigen_structure, min_variance)
  }
  values <- solve_values(orientations)
  objective <- function(rotations, values) sum(vapply(seq_len(k), function(profile) {
    weights[profile] * sum(log(values[profile, ])) +
      sum(diag(crossprod(rotations[[profile]], scatter[[profile]] %*%
        rotations[[profile]])) / values[profile, ])
  }, numeric(1)))
  current <- objective(orientations, values)
  iteration <- 0L
  settled <- !common
  while (common && iteration < max_iter && !settled) {
    rotated <- .multilpa_bounded_orientation(orientations[[1L]], scatter, values)
    next_orientations <- rep(list(rotated), k)
    next_values <- solve_values(next_orientations)
    next_objective <- objective(next_orientations, next_values)
    if (next_objective > current) break
    settled <- current - next_objective <= tol * (1 + abs(current))
    orientations <- next_orientations
    values <- next_values
    current <- next_objective
    iteration <- iteration + 1L
  }
  blocks <- lapply(seq_len(k), function(profile) {
    covariance <- orientations[[profile]] %*% diag(values[profile, ], d) %*%
      t(orientations[[profile]])
    (covariance + t(covariance)) / 2
  })
  # The previous feasible covariance is a valid generalized EM update if a
  # bounded orientation search cannot improve it at the available precision.
  if (!is.null(start) && identical(as.integer(dim(start)), c(d, d, k))) {
    previous <- lapply(seq_len(k), function(profile) matrix(start[, , profile], d, d))
    feasible <- all(vapply(previous, function(block) {
      min(eigen(block, symmetric = TRUE, only.values = TRUE)$values) >=
        min_variance * (1 - 1e-10)
    }, logical(1)))
    if (feasible) {
      # A SQUAREM proposal can be positive definite while lying outside the
      # nonlinear covariance family. It can seed the orientation search, but
      # cannot be retained as a feasible fallback merely because it fits well.
      eigenvalues <- t(vapply(previous, function(block) {
        eigen(block, symmetric = TRUE, only.values = TRUE)$values
      }, numeric(d))) |> matrix(k, d)
      log_volume <- rowMeans(log(eigenvalues))
      if (startsWith(structure, "E")) {
        feasible <- max(abs(log_volume - log_volume[1L])) <= 1e-10
      }
      if (structure %in% c("EEV", "VEV")) {
        log_shape <- log(eigenvalues) - log_volume
        feasible <- feasible && max(abs(sweep(log_shape, 2L, log_shape[1L, ]))) <= 1e-10
      }
      normalized <- lapply(seq_len(k), function(profile) {
        previous[[profile]] / exp(log_volume[profile])
      })
      if (structure == "VEE") {
        feasible <- feasible && all(vapply(normalized, function(block) {
          max(abs(block - normalized[[1L]])) <= 1e-10 * max(1, max(abs(block)))
        }, logical(1)))
      }
      if (structure %in% c("EVE", "VVE")) {
        first <- previous[[1L]] / max(abs(previous[[1L]]))
        feasible <- feasible && all(vapply(previous, function(block) {
          block <- block / max(abs(block))
          max(abs(first %*% block - block %*% first)) <= 1e-10
        }, logical(1)))
      }
    }
    if (feasible) {
      previous_objective <- sum(vapply(seq_len(k), function(profile) {
        factor <- chol(previous[[profile]])
        2 * weights[profile] * sum(log(diag(factor))) +
          sum(chol2inv(factor) * scatter[[profile]])
      }, numeric(1)))
      if (previous_objective < current) blocks <- previous
    }
  }
  array(unlist(blocks, use.names = FALSE), c(d, d, k))
}
