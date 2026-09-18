#' Per-group scores of the observed-data log likelihood
#'
#' Decomposes the gradient into one row per independent group. The groups are
#' the independent sampling units of the two-level likelihood, so these rows are
#' the contributions the robust cross-product matrix accumulates.
#'
#' @param theta Unconstrained parameter vector.
#' @param x Centered numeric indicator matrix, possibly containing NA.
#' @param object Fitted model defining dimensions and groups.
#' @return A groups-by-parameters matrix whose column sums equal the gradient of
#'   the positive observed-data log likelihood, that is `-.ml_lpa_score()`.
#' @noRd
.ml_lpa_group_scores <- function(theta, x, object) {
  stopifnot("`theta` must be numeric" = is.numeric(theta),
            "`x` must be a numeric matrix" = is.matrix(x) && is.numeric(x),
            "`object` must be an `ml_lpa` fit" = inherits(object, "ml_lpa"))
  parameters <- .ml_lpa_decode(theta, object)
  group_index <- object$group_index
  expectation <- .ml_lpa_expectation(x, group_index, parameters)
  n_profiles <- object$n_profiles
  n_types <- object$n_group_classes
  measurement <- if (identical(object$covariance_model, "full")) {
    .ml_lpa_full_measurement_group_score(parameters, expectation, object)
  } else {
    observed <- !is.na(x)
    blocks <- lapply(seq_len(n_profiles), function(profile) {
      residuals <- sweep(x, 2L, parameters$means[profile, ], "-")
      residuals[!observed] <- 0
      weights <- expectation$subject_posteriors[, profile]
      list(
        means = rowsum(sweep(residuals, 2L, parameters$variances[profile, ], "/") *
                         weights, group_index, reorder = FALSE),
        variances = rowsum(0.5 * (sweep(residuals^2, 2L,
          parameters$variances[profile, ], "/") - observed) * weights,
          group_index, reorder = FALSE))
    })
    variances <- if (object$variance_model == "equal") {
      Reduce(`+`, lapply(blocks, `[[`, "variances"))
    } else do.call(cbind, lapply(blocks, `[[`, "variances"))
    list(means = do.call(cbind, lapply(blocks, `[[`, "means")),
         covariances = variances)
  }
  # Each group class contributes its own multinomial score for the profile
  # proportions, with the final profile as the omitted reference category.
  profile_scores <- do.call(cbind, lapply(seq_len(n_types), function(group_type) {
    counts <- rowsum(expectation$joint[[group_type]], group_index, reorder = FALSE)
    expected <- outer(rowSums(counts), parameters$profile_probabilities[group_type, ])
    (counts - expected)[, seq_len(n_profiles - 1L), drop = FALSE]
  }))
  group_scores <- sweep(expectation$group_posteriors, 2L,
    parameters$group_probabilities, "-")[, seq_len(n_types - 1L), drop = FALSE]
  scores <- cbind(measurement$means, measurement$covariances,
                  profile_scores, group_scores)
  dimnames(scores) <- list(object$group_ids, names(theta))
  scores
}

#' Per-group Gaussian scores in log-Cholesky coordinates
#' @param parameters Decoded Gaussian and mixing parameters.
#' @param expectation Nested expectation step including Gaussian moments.
#' @param object A fitted full-covariance model.
#' @return Groups-by-parameter mean and log-Cholesky score blocks.
#' @noRd
.ml_lpa_full_measurement_group_score <- function(parameters, expectation, object) {
  stopifnot("`object` must be a full-covariance `ml_lpa` fit" =
              inherits(object, "ml_lpa") &&
              identical(object$covariance_model, "full"))
  dimension <- length(object$indicators)
  group_index <- object$group_index
  n_groups <- object$n_groups
  lower <- lower.tri(matrix(0, dimension, dimension), diag = TRUE)
  pairs <- expand.grid(row = seq_len(dimension), column = seq_len(dimension))
  component <- lapply(seq_len(object$n_profiles), function(profile) {
    covariance <- matrix(parameters$covariances[, , profile], dimension, dimension)
    precision <- chol2inv(chol(covariance))
    moments <- expectation$gaussian_moments[[profile]]
    weights <- expectation$subject_posteriors[, profile]
    residuals <- sweep(moments$expected, 2L, parameters$means[profile, ], "-")
    weighted <- residuals * weights
    group_weight <- as.vector(rowsum(weights, group_index, reorder = FALSE))
    # Weighted residual cross-products accumulated one group at a time.
    scatter <- rowsum(residuals[, pairs$row, drop = FALSE] *
                        weighted[, pairs$column, drop = FALSE],
                      group_index, reorder = FALSE)
    # Conditional-covariance corrections for each missing-data pattern enter
    # only through the posterior weight each group contributes to that pattern.
    corrections <- lapply(moments$adjustments, function(adjustment) {
      indicator <- numeric(length(weights))
      indicator[adjustment$rows] <- weights[adjustment$rows]
      list(weight = as.vector(rowsum(indicator, group_index, reorder = FALSE)),
           covariance = as.vector(adjustment$covariance))
    })
    if (length(corrections) > 0L) {
      scatter <- scatter + Reduce(`+`, lapply(corrections, function(correction) {
        outer(correction$weight, correction$covariance)
      }))
    }
    list(precision = precision, group_weight = group_weight, scatter = scatter,
         means = rowsum(weighted, group_index, reorder = FALSE) %*% precision)
  })
  covariance_blocks <- lapply(component, function(part) {
    # 0.5 * (P S P - w P) for every group, with S held as a flattened row.
    t(vapply(seq_len(n_groups), function(group) {
      scatter <- matrix(part$scatter[group, ], dimension, dimension)
      as.vector(0.5 * (part$precision %*% scatter %*% part$precision -
                         part$group_weight[group] * part$precision))
    }, numeric(dimension * dimension)))
  })
  profiles <- if (object$variance_model == "equal") 1L else seq_len(object$n_profiles)
  cholesky_blocks <- lapply(profiles, function(profile) {
    block <- if (object$variance_model == "equal") {
      Reduce(`+`, covariance_blocks)
    } else covariance_blocks[[profile]]
    factor <- t(chol(matrix(parameters$covariances[, , profile], dimension, dimension)))
    t(vapply(seq_len(n_groups), function(group) {
      factor_score <- 2 * matrix(block[group, ], dimension, dimension) %*% factor
      diag(factor_score) <- diag(factor_score) * diag(factor)
      factor_score[lower]
    }, numeric(sum(lower))))
  })
  list(means = do.call(cbind, lapply(component, `[[`, "means")),
       covariances = do.call(cbind, cholesky_blocks))
}

#' Robust cross-product matrix from per-group scores
#' @param scores Groups-by-parameters score matrix.
#' @return The summed outer product of the group scores.
#' @noRd
.ml_lpa_cross_product <- function(scores) {
  stopifnot("`scores` must be a numeric matrix" =
              is.matrix(scores) && is.numeric(scores))
  if (any(!is.finite(scores))) {
    stop(errorCondition("Per-group scores are not finite; robust inference is unavailable.",
                        class = "mllpa_bad_scores", call = NULL))
  }
  crossprod(scores)
}
