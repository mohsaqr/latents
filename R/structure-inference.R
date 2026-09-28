# Wald inference for the constrained covariance structures.
#
# Each profile's covariance decomposes as Sigma_k = lambda_k * D_k A_k D_k'.
# The ten structures that constrain the volume, the shape or the orientation
# across profiles are estimated here in free, unconstrained coordinates that
# span exactly the constrained family, so the observed information has one row
# per parameter the structure has --- the count `.multilpa_structure_parameters()`
# gives --- and nothing more:
#
# * a log volume, one per profile or one shared;
# * a log shape with the determinant-one constraint written into it: the first
#   `d - 1` log shape values are free and the last is minus their sum;
# * an orientation, as the Cayley transform R = D0 (I - S)^-1 (I + S) of a
#   skew-symmetric S around the fitted eigenvectors D0. The chart is local ---
#   S = 0 at the estimate --- which is all the curvature of the likelihood at
#   the estimate needs;
# * for the two structures whose shape and orientation are constrained
#   together (VEE shares both, EVV frees both), the determinant-one matrix
#   C = L L' in log-Cholesky coordinates whose log diagonal sums to zero,
#   which needs no eigenvectors at all.
#
# The score in these coordinates is the chain rule through the map to the
# covariance matrices, dL/dphi = sum_k <G_k, dSigma_k/dphi>, where
# G_k = (P S P - w P) / 2 is the gradient with respect to an unconstrained
# Sigma_k. Every derivative of the map is analytic.

#' The structures estimated in their own separable coordinates
#'
#' EEI, VVI, EEE and VVV constrain nothing across profiles beyond sharing or
#' freeing the whole block, so each profile's block is maximized by itself and
#' is charted as log variances or log-Cholesky coordinates. The other ten tie
#' the profiles together.
#'
#' @return A character vector of mclust model codes.
#' @noRd
.multilpa_separable_structures <- function() c("EEI", "VVI", "EEE", "VVV")

#' Is a fit estimated in a constrained-structure chart?
#' @param object A fitted `multilpa` model.
#' @return A single logical.
#' @noRd
.multilpa_uses_chart <- function(object) {
  structure <- object$covariance_structure
  !is.null(structure) && length(structure) == 1L && !is.na(structure) &&
    structure %in% .multilpa_structures() &&
    !structure %in% .multilpa_separable_structures() &&
    length(.multilpa_continuous_names(object)) > 0L
}

#' Is the natural spread block one shared block rather than one per profile?
#'
#' A block is reported once when the structure makes every profile's block the
#' same matrix, and once per profile otherwise. For the four separable
#' structures that is exactly `variance_model = "equal"`; of the other ten only
#' EII makes every block identical.
#'
#' @param object A fitted `multilpa` model.
#' @return A single logical.
#' @noRd
.multilpa_spread_shared <- function(object) {
  if (.multilpa_uses_chart(object)) {
    return(identical(object$covariance_structure, "EII"))
  }
  identical(object$variance_model, "equal")
}

#' Describe the free coordinates of a constrained structure
#'
#' @param structure An mclust model code outside the separable four.
#' @param k Number of profiles.
#' @param d Number of continuous indicators.
#' @return A list: `form` (`"spherical"`, `"diagonal"`, `"rotation"` or
#'   `"cholesky"`), the row counts and per-profile row maps of the volume, shape
#'   and rotation blocks, their widths and offsets, and `width`, the total.
#' @noRd
.multilpa_chart <- function(structure, k, d) {
  stopifnot(
    "`structure` must be a constrained structure" =
      length(structure) == 1L && structure %in% .multilpa_structures() &&
      !structure %in% .multilpa_separable_structures(),
    "`k` must be a positive whole number" = length(k) == 1L && k >= 1,
    "`d` must be a positive whole number" = length(d) == 1L && d >= 1)
  second <- substr(structure, 2L, 2L)
  third <- substr(structure, 3L, 3L)
  form <- if (structure %in% c("EII", "VII")) "spherical" else
    if (third == "I") "diagonal" else
      if (second == third) "cholesky" else "rotation"
  volume_shared <- startsWith(structure, "E")
  shape_shared <- second == "E"
  rotation_shared <- third == "E"
  map <- function(shared) if (shared) rep(1L, k) else seq_len(k)
  shape_rows <- if (form == "spherical") 0L else if (shape_shared) 1L else k
  rotation_rows <- if (form == "rotation") {
    if (rotation_shared) 1L else k
  } else 0L
  per_shape <- switch(form, spherical = 0L, diagonal = d - 1L, rotation = d - 1L,
                      cholesky = d * (d + 1L) / 2L - 1L)
  per_rotation <- d * (d - 1L) / 2L
  widths <- c(volume = if (volume_shared) 1L else k,
              shape = shape_rows * per_shape,
              rotation = rotation_rows * per_rotation)
  list(structure = structure, form = form, k = k, d = d,
       volume_shared = volume_shared, shape_shared = shape_shared,
       rotation_shared = rotation_shared,
       volume_of = map(volume_shared), shape_of = map(shape_shared),
       rotation_of = map(rotation_shared),
       shape_rows = shape_rows, rotation_rows = rotation_rows,
       per_shape = as.integer(per_shape), per_rotation = as.integer(per_rotation),
       widths = widths,
       offsets = stats::setNames(cumsum(c(0L, widths))[seq_len(3L)], names(widths)),
       width = as.integer(sum(widths)))
}

#' The (row, column) positions a Cholesky shape block is charted by
#'
#' The lower triangle in column-major order, less the last diagonal entry,
#' whose logarithm is minus the sum of the others.
#'
#' @param d The dimension.
#' @return A two-column integer matrix of positions.
#' @noRd
.multilpa_chart_cholesky_positions <- function(d) {
  positions <- which(lower.tri(matrix(0, d, d), diag = TRUE), arr.ind = TRUE)
  positions[!(positions[, 1L] == d & positions[, 2L] == d), , drop = FALSE]
}

#' The (a, b), a < b, pairs a skew-symmetric rotation block is charted by
#' @param d The dimension.
#' @return A two-column integer matrix of positions.
#' @noRd
.multilpa_chart_rotation_pairs <- function(d) {
  which(upper.tri(matrix(0, d, d)), arr.ind = TRUE)
}

#' The fitted covariance matrices of a fit, one per profile
#' @param object A fitted `multilpa` model.
#' @return A list of `d` by `d` matrices.
#' @noRd
.multilpa_fitted_blocks <- function(object) {
  d <- length(.multilpa_continuous_names(object))
  lapply(seq_len(object$n_profiles), function(profile) {
    if (is.null(object$covariances)) {
      diag(unname(object$variances[profile, ]), d)
    } else {
      block <- unname(matrix(object$covariances[, , profile], d, d))
      (block + t(block)) / 2
    }
  })
}

#' The orientations the rotation chart is centred on
#'
#' With a varying orientation each profile's own eigenvectors. With a shared one
#' the eigenvectors of a generic positive combination of the profiles'
#' covariances: every one of them is diagonal in the shared orientation, so any
#' combination is too, and unequal weights keep its eigenvalues apart.
#'
#' @param blocks The fitted covariance matrices.
#' @param chart The chart, from `.multilpa_chart()`.
#' @return A list of orthogonal matrices, one per rotation row, or `NULL`.
#' @noRd
.multilpa_chart_anchor <- function(blocks, chart) {
  if (!identical(chart$form, "rotation")) return(NULL)
  if (chart$rotation_shared) {
    combination <- Reduce(`+`, lapply(seq_along(blocks), function(profile) {
      sqrt(profile + 1) * blocks[[profile]]
    }))
    return(list(eigen(combination, symmetric = TRUE)$vectors))
  }
  lapply(blocks, function(block) eigen(block, symmetric = TRUE)$vectors)
}

#' Encode fitted covariances in the constrained-structure chart
#'
#' @param object A fitted `multilpa` model whose structure uses a chart.
#' @return A list with `values` (the named chart coordinates), `scale` (a
#'   typical magnitude for each, used only to condition the numerical Hessian),
#'   `chart` and `anchor`.
#' @noRd
.multilpa_chart_encode <- function(object) {
  continuous <- .multilpa_continuous_names(object)
  k <- object$n_profiles
  d <- length(continuous)
  chart <- .multilpa_chart(object$covariance_structure, k, d)
  blocks <- .multilpa_fitted_blocks(object)
  anchor <- .multilpa_chart_anchor(blocks, chart)
  log_volumes <- vapply(blocks, function(block) {
    as.numeric(determinant(block, logarithm = TRUE)$modulus) / d
  }, numeric(1))
  label <- function(shared, row) if (shared) "shared" else as.character(row)
  volume_values <- if (chart$volume_shared) mean(log_volumes) else log_volumes
  volume_names <- sprintf("log_volume[%s]", vapply(seq_len(chart$widths[["volume"]]),
    function(row) label(chart$volume_shared, row), character(1)))
  # The per-profile log shapes, before any sharing: log eigenvalues along the
  # anchor minus the log volume.
  profile_log_shapes <- function() {
    lapply(seq_len(k), function(profile) {
      block <- blocks[[profile]]
      values <- switch(chart$form,
        diagonal = diag(block),
        rotation = diag(crossprod(anchor[[chart$rotation_of[profile]]],
                                  block %*% anchor[[chart$rotation_of[profile]]])))
      log(values) - log_volumes[profile]
    })
  }
  shape <- if (identical(chart$form, "spherical")) {
    list(values = numeric(0), names = character(0), scale = numeric(0))
  } else if (chart$form %in% c("diagonal", "rotation")) {
    shapes <- profile_log_shapes()
    rows <- if (chart$shape_shared) list(Reduce(`+`, shapes) / k) else shapes
    free <- seq_len(d - 1L)
    list(values = unlist(lapply(rows, function(row) row[free]), use.names = FALSE),
         names = unlist(lapply(seq_along(rows), function(row) {
           sprintf("log_shape[%s,%s]", label(chart$shape_shared, row), continuous[free])
         }), use.names = FALSE),
         scale = rep(1, length(rows) * length(free)))
  } else {
    positions <- .multilpa_chart_cholesky_positions(d)
    matrices <- if (chart$shape_shared) {
      list(Reduce(`+`, lapply(seq_len(k), function(profile) {
        blocks[[profile]] / exp(log_volumes[profile])
      })) / k)
    } else {
      lapply(seq_len(k), function(profile) {
        blocks[[profile]] / exp(volume_values[chart$volume_of[profile]])
      })
    }
    pieces <- lapply(seq_along(matrices), function(row) {
      factor <- t(chol(matrices[[row]]))
      on_diagonal <- positions[, 1L] == positions[, 2L]
      values <- factor[positions]
      values[on_diagonal] <- log(values[on_diagonal])
      kind <- ifelse(on_diagonal, "log_shape_cholesky", "shape_cholesky")
      list(values = values,
           names = sprintf("%s[%s,%s,%s]", kind, label(chart$shape_shared, row),
                           continuous[positions[, 1L]], continuous[positions[, 2L]]),
           scale = ifelse(on_diagonal, 1, sqrt(diag(matrices[[row]]))[positions[, 1L]]))
    })
    list(values = unlist(lapply(pieces, `[[`, "values"), use.names = FALSE),
         names = unlist(lapply(pieces, `[[`, "names"), use.names = FALSE),
         scale = unlist(lapply(pieces, `[[`, "scale"), use.names = FALSE))
  }
  pairs <- .multilpa_chart_rotation_pairs(d)
  rotation_names <- if (chart$rotation_rows == 0L) character(0) else {
    unlist(lapply(seq_len(chart$rotation_rows), function(row) {
      sprintf("rotation[%s,%s,%s]", label(chart$rotation_shared, row),
              continuous[pairs[, 1L]], continuous[pairs[, 2L]])
    }), use.names = FALSE)
  }
  values <- c(volume_values, shape$values, rep(0, length(rotation_names)))
  names(values) <- c(volume_names, shape$names, rotation_names)
  stopifnot("the chart must have one coordinate per structure parameter" =
              length(values) == chart$width &&
              chart$width == .multilpa_structure_parameters(chart$structure, k, d))
  list(values = values,
       scale = c(rep(1, length(volume_values)), shape$scale,
                 rep(1, length(rotation_names))),
       chart = chart, anchor = anchor)
}

#' Map chart coordinates to covariance matrices, and differentiate the map
#'
#' @param phi The chart coordinates, in `.multilpa_chart_encode()` order.
#' @param chart The chart.
#' @param anchor The orientations the rotation chart is centred on.
#' @param derivatives Whether to return the Jacobian as well.
#' @return A list with `sigma`, one covariance matrix per profile, and, when
#'   asked, `jacobian`, one `d^2` by `length(phi)` matrix per profile holding
#'   `d vec(Sigma_k) / d phi`.
#' @noRd
.multilpa_chart_map <- function(phi, chart, anchor, derivatives = FALSE) {
  stopifnot("`phi` must carry one value per chart coordinate" =
              is.numeric(phi) && length(phi) == chart$width)
  phi <- unname(phi)
  d <- chart$d
  k <- chart$k
  identity <- diag(d)
  volume_block <- phi[chart$offsets[["volume"]] + seq_len(chart$widths[["volume"]])]
  shape_index <- function(row) {
    chart$offsets[["shape"]] + (row - 1L) * chart$per_shape + seq_len(chart$per_shape)
  }
  rotation_index <- function(row) {
    chart$offsets[["rotation"]] + (row - 1L) * chart$per_rotation +
      seq_len(chart$per_rotation)
  }
  # Full log shape vectors (diagonal and rotation forms) or Cholesky factors.
  log_shape <- function(row) {
    if (identical(chart$form, "spherical") || d == 1L) return(numeric(d))
    free <- phi[shape_index(row)]
    c(free, -sum(free))
  }
  positions <- .multilpa_chart_cholesky_positions(d)
  cholesky_factor <- function(row) {
    factor <- matrix(0, d, d)
    factor[positions] <- phi[shape_index(row)]
    on_diagonal <- positions[positions[, 1L] == positions[, 2L], , drop = FALSE]
    logs <- phi[shape_index(row)][positions[, 1L] == positions[, 2L]]
    factor[on_diagonal] <- exp(logs)
    factor[d, d] <- exp(-sum(logs))
    factor
  }
  pairs <- .multilpa_chart_rotation_pairs(d)
  rotation <- function(row) {
    skew <- matrix(0, d, d)
    skew[pairs] <- phi[rotation_index(row)]
    skew[pairs[, c(2L, 1L), drop = FALSE]] <- -phi[rotation_index(row)]
    inverse <- solve(identity - skew)
    cayley <- inverse %*% (identity + skew)
    list(matrix = anchor[[row]] %*% cayley, inverse = inverse, cayley = cayley)
  }
  shapes <- if (chart$form == "cholesky") {
    lapply(seq_len(chart$shape_rows), cholesky_factor)
  } else lapply(seq_len(max(chart$shape_rows, 1L)), log_shape)
  rotations <- if (chart$rotation_rows > 0L) lapply(seq_len(chart$rotation_rows), rotation)
  profile_parts <- lapply(seq_len(k), function(profile) {
    lambda <- exp(volume_block[chart$volume_of[profile]])
    shape_row <- if (chart$shape_rows == 0L) 1L else chart$shape_of[profile]
    orientation <- if (chart$rotation_rows > 0L) {
      rotations[[chart$rotation_of[profile]]]$matrix
    } else identity
    unit_sigma <- if (chart$form == "cholesky") {
      tcrossprod(shapes[[shape_row]])
    } else {
      orientation %*% (exp(shapes[[shape_row]]) * t(orientation))
    }
    list(lambda = lambda, shape_row = shape_row, orientation = orientation,
         sigma = lambda * unit_sigma)
  })
  sigma <- lapply(profile_parts, `[[`, "sigma")
  if (!derivatives) return(list(sigma = sigma))
  width <- chart$width
  jacobian <- lapply(seq_len(k), function(profile) {
    part <- profile_parts[[profile]]
    columns <- matrix(0, d * d, width)
    columns[, chart$offsets[["volume"]] + chart$volume_of[profile]] <- as.vector(part$sigma)
    if (chart$shape_rows > 0L && chart$per_shape > 0L) {
      index <- shape_index(part$shape_row)
      if (chart$form == "cholesky") {
        factor <- shapes[[part$shape_row]]
        columns[, index] <- vapply(seq_len(nrow(positions)), function(position) {
          row <- positions[position, 1L]
          column <- positions[position, 2L]
          change <- matrix(0, d, d)
          if (row == column) {
            change[row, row] <- factor[row, row]
            change[d, d] <- -factor[d, d]
          } else change[row, column] <- 1
          as.vector(part$lambda * (change %*% t(factor) + factor %*% t(change)))
        }, numeric(d * d))
      } else {
        values <- exp(shapes[[part$shape_row]])
        columns[, index] <- vapply(seq_len(d - 1L), function(free) {
          direction <- numeric(d)
          direction[free] <- 1
          direction[d] <- -1
          as.vector(part$lambda * part$orientation %*%
                      ((values * direction) * t(part$orientation)))
        }, numeric(d * d))
      }
    }
    if (chart$rotation_rows > 0L && chart$per_rotation > 0L) {
      row <- chart$rotation_of[profile]
      piece <- rotations[[row]]
      scaled <- exp(shapes[[part$shape_row]])
      columns[, rotation_index(row)] <- vapply(seq_len(nrow(pairs)), function(pair) {
        skew <- matrix(0, d, d)
        skew[pairs[pair, 1L], pairs[pair, 2L]] <- 1
        skew[pairs[pair, 2L], pairs[pair, 1L]] <- -1
        change <- anchor[[row]] %*% piece$inverse %*% skew %*% (identity + piece$cayley)
        half <- change %*% (scaled * t(part$orientation))
        as.vector(part$lambda * (half + t(half)))
      }, numeric(d * d))
    }
    columns
  })
  list(sigma = sigma, jacobian = jacobian)
}

#' Decode the spread block of a chart-coordinate parameter vector
#'
#' @param phi The chart coordinates.
#' @param object The fit defining the chart and its anchor.
#' @return A list with `variances` (profiles by indicators) and, for an
#'   ellipsoidal structure, `covariances` (`d` by `d` by profiles).
#' @noRd
.multilpa_chart_decode <- function(phi, object) {
  d <- length(.multilpa_continuous_names(object))
  chart <- .multilpa_chart(object$covariance_structure, object$n_profiles, d)
  anchor <- .multilpa_chart_anchor(.multilpa_fitted_blocks(object), chart)
  sigma <- .multilpa_chart_map(phi, chart, anchor)$sigma
  variances <- matrix(vapply(sigma, diag, numeric(d)), object$n_profiles, d,
                      byrow = TRUE)
  if (!.multilpa_is_ellipsoidal(object$covariance_structure)) {
    return(list(variances = variances))
  }
  list(variances = variances,
       covariances = array(unlist(sigma, use.names = FALSE),
                           c(d, d, object$n_profiles)))
}

#' Gradients of the log likelihood with respect to each profile's covariance
#'
#' By Fisher's identity the observed-data score is the conditional expectation
#' of the complete-data score, so with missing indicators the scatter carries
#' the conditional covariance of what was not observed. The gradient with
#' respect to an unconstrained symmetric `Sigma_k` is `(P S P - w P) / 2`.
#'
#' @param parameters Decoded parameters.
#' @param expectation The expectation step at those parameters.
#' @param x The centred indicator matrix.
#' @param group_index Group indices, or `NULL` for pooled totals.
#' @return A list with `means` (groups, or one row, by profiles-times-indicators)
#'   and `covariance`, one groups-by-`d^2` matrix per profile.
#' @noRd
.multilpa_profile_gradients <- function(parameters, expectation, x,
                                        group_index = NULL) {
  d <- ncol(x)
  k <- nrow(parameters$means)
  index <- group_index %||% rep(1L, nrow(x))
  pairs <- expand.grid(row = seq_len(d), column = seq_len(d))
  pieces <- lapply(seq_len(k), function(profile) {
    covariance <- if (is.null(parameters$covariances)) {
      diag(parameters$variances[profile, ], d)
    } else matrix(parameters$covariances[, , profile], d, d)
    precision <- chol2inv(chol(covariance))
    weights <- expectation$subject_posteriors[, profile]
    moments <- expectation$gaussian_moments[[profile]]
    expected <- if (is.null(moments)) x else moments$expected
    residuals <- sweep(expected, 2L, parameters$means[profile, ], "-")
    weighted <- residuals * weights
    scatter <- rowsum(residuals[, pairs$row, drop = FALSE] *
                        weighted[, pairs$column, drop = FALSE],
                      index, reorder = FALSE)
    adjustments <- if (is.null(moments)) list() else moments$adjustments
    if (length(adjustments) > 0L) {
      scatter <- scatter + Reduce(`+`, lapply(adjustments, function(adjustment) {
        share <- numeric(length(weights))
        share[adjustment$rows] <- weights[adjustment$rows]
        outer(as.vector(rowsum(share, index, reorder = FALSE)),
              as.vector(adjustment$covariance))
      }))
    }
    group_weight <- as.vector(rowsum(weights, index, reorder = FALSE))
    list(means = rowsum(weighted, index, reorder = FALSE) %*% precision,
         covariance = 0.5 * (scatter %*% kronecker(precision, precision) -
                               outer(group_weight, as.vector(precision))))
  })
  list(means = do.call(cbind, lapply(pieces, `[[`, "means")),
       covariance = lapply(pieces, `[[`, "covariance"))
}

#' Measurement scores in the constrained-structure chart
#'
#' @param parameters Decoded parameters.
#' @param expectation The expectation step at those parameters.
#' @param x The centred indicator matrix.
#' @param object The fit defining the chart.
#' @param phi The chart coordinates the parameters were decoded from.
#' @param group_index Group indices, or `NULL` for pooled totals.
#' @return A list with `means` and `covariances`, each a matrix with one row per
#'   group (or one row) and one column per coordinate, of the positive log
#'   likelihood's gradient.
#' @noRd
.multilpa_chart_scores <- function(parameters, expectation, x, object, phi,
                                   group_index = NULL) {
  d <- ncol(x)
  chart <- .multilpa_chart(object$covariance_structure, object$n_profiles, d)
  anchor <- .multilpa_chart_anchor(.multilpa_fitted_blocks(object), chart)
  jacobian <- .multilpa_chart_map(phi, chart, anchor, derivatives = TRUE)$jacobian
  gradients <- .multilpa_profile_gradients(parameters, expectation, x, group_index)
  list(means = gradients$means,
       covariances = Reduce(`+`, lapply(seq_len(object$n_profiles), function(profile) {
         gradients$covariance[[profile]] %*% jacobian[[profile]]
       })))
}

#' Natural spread coefficients differentiated with respect to the chart
#'
#' @param object A fit whose structure uses a chart.
#' @return A matrix with one row per natural spread coefficient, in
#'   `.multilpa_coefficients()` order, and one column per chart coordinate.
#' @noRd
.multilpa_chart_natural_jacobian <- function(object) {
  encoded <- .multilpa_chart_encode(object)
  d <- encoded$chart$d
  jacobian <- .multilpa_chart_map(encoded$values, encoded$chart, encoded$anchor,
                                  derivatives = TRUE)$jacobian
  rows <- if (.multilpa_spread_shared(object)) 1L else seq_len(object$n_profiles)
  keep <- if (.multilpa_is_ellipsoidal(object$covariance_structure)) {
    which(lower.tri(matrix(0, d, d), diag = TRUE))
  } else which(diag(d) == 1)
  do.call(rbind, lapply(rows, function(profile) {
    jacobian[[profile]][keep, , drop = FALSE]
  }))
}
