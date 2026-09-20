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
#'   the positive observed-data log likelihood, that is `-.multilpa_score()`.
#' @noRd
.multilpa_group_scores <- function(theta, x, object, codes = NULL) {
  stopifnot("`theta` must be numeric" = is.numeric(theta),
            "`x` must be a numeric matrix" = is.matrix(x) && is.numeric(x),
            "`object` must be an `multilpa` fit" = inherits(object, "multilpa"))
  parameters <- .multilpa_decode(theta, object)
  group_index <- object$group_index
  expectation <- .multilpa_expectation(x, group_index, parameters, codes)
  n_profiles <- object$n_profiles
  n_types <- object$n_group_classes
  measurement <- if (identical(object$covariance_model, "full")) {
    .multilpa_full_measurement_group_score(parameters, expectation, object)
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
  response_scores <- .multilpa_response_scores(
    codes, expectation$subject_posteriors, parameters$response_probabilities,
    group_index)
  scores <- cbind(measurement$means, measurement$covariances, response_scores,
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
.multilpa_full_measurement_group_score <- function(parameters, expectation, object) {
  stopifnot("`object` must be a full-covariance `multilpa` fit" =
              inherits(object, "multilpa") &&
              identical(object$covariance_model, "full"))
  dimension <- length(.multilpa_continuous_names(object))
  group_index <- object$group_index
  n_groups <- object$n_groups
  if (dimension == 0L) {
    empty <- matrix(numeric(0), n_groups, 0L)
    return(list(means = empty, covariances = empty))
  }
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
    t(vapply(seq_len(n_groups), function(id) {
      scatter <- matrix(part$scatter[id, ], dimension, dimension)
      as.vector(0.5 * (part$precision %*% scatter %*% part$precision -
                         part$group_weight[id] * part$precision))
    }, numeric(dimension * dimension)))
  })
  profiles <- if (object$variance_model == "equal") 1L else seq_len(object$n_profiles)
  cholesky_blocks <- lapply(profiles, function(profile) {
    block <- if (object$variance_model == "equal") {
      Reduce(`+`, covariance_blocks)
    } else covariance_blocks[[profile]]
    factor <- t(chol(matrix(parameters$covariances[, , profile], dimension, dimension)))
    t(vapply(seq_len(n_groups), function(id) {
      factor_score <- 2 * matrix(block[id, ], dimension, dimension) %*% factor
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
.multilpa_cross_product <- function(scores) {
  stopifnot("`scores` must be a numeric matrix" =
              is.matrix(scores) && is.numeric(scores))
  if (any(!is.finite(scores))) {
    stop(errorCondition("Per-group scores are not finite; robust inference is unavailable.",
                        class = "multilpa_bad_scores", call = NULL))
  }
  ## A single row cannot support a sandwich: at the estimate the rows sum to
  ## zero, so one row is itself the zero vector up to rounding and the meat is
  ## numerical dust rather than a variance.
  .multilpa_require_clusters(nrow(scores), 1L, "A robust cross-product matrix")
  crossprod(scores)
}

#' Refuse a cluster-robust variance too few independent units can support
#'
#' The meat of a cluster-robust sandwich is the sum of one outer product per
#' independent unit, formed from that unit's summed estimating-equation
#' contributions. Those contributions sum to zero across units at the estimate,
#' because that sum *is* the stationarity condition the estimate solves. The
#' meat therefore has rank at most `n_units - 1`, whatever the data: with one
#' unit it is exactly zero, and with `n_units <= n_quantities` it is singular,
#' so some contrast among the reported quantities is handed a variance of
#' exactly zero. Neither is a small standard error; both are refused.
#'
#' The rule this enforces, `n_units > n_quantities`, is the same one
#' [parameter_inference()] applies to a robust model vcov through
#' `.multilpa_check_regularity()`, so every robust path in the package asks for
#' as many independent units as the covariance it reports has dimensions, plus
#' one for the constraint that the contributions sum to zero.
#'
#' @param n_units Number of independent units the contributions were summed
#'   within: groups for a clustered variance, observations for an unclustered
#'   one.
#' @param n_quantities Number of quantities whose covariance is being reported.
#' @param what A noun phrase naming that covariance, opening the message.
#' @param unit What one independent unit is called in the message.
#' @param alternative A sentence naming a legitimate non-clustered alternative,
#'   or `NULL` when there is none.
#' @return `invisible(TRUE)` when the rule holds. Otherwise raises
#'   `multilpa_too_few_groups`.
#' @references Cameron, A. C., & Miller, D. L. (2015). A practitioner's guide
#'   to cluster-robust inference. *Journal of Human Resources*, 50, 317--372.
#' @noRd
.multilpa_require_clusters <- function(n_units, n_quantities, what,
                                       unit = "independent groups",
                                       alternative = NULL) {
  stopifnot(
    "`n_units` must be a single non-negative whole number" =
      length(n_units) == 1L && is.numeric(n_units) && is.finite(n_units) &&
      n_units >= 0,
    "`n_quantities` must be a single positive whole number" =
      length(n_quantities) == 1L && is.numeric(n_quantities) &&
      is.finite(n_quantities) && n_quantities >= 1,
    "`what` must be a single string" = is.character(what) && length(what) == 1L,
    "`unit` must be a single string" = is.character(unit) && length(unit) == 1L,
    "`alternative` must be a single string or NULL" =
      is.null(alternative) || (is.character(alternative) &&
                                 length(alternative) == 1L))
  if (n_units > n_quantities) return(invisible(TRUE))
  plural <- if (n_quantities == 1L) "y" else "ies"
  reason <- if (n_units <= 1L) {
    paste("The contributions of a single unit are the estimating equation the",
          "estimate solves, so they sum to exactly zero and the variance would",
          "be floating-point noise rather than an estimate.")
  } else {
    sprintf(paste("The contributions sum to zero at the estimate, so %d %s leave at",
                  "most %d independent contributions for %d quantit%s: the",
                  "covariance is singular and some contrast among them would be",
                  "given a variance of exactly zero."),
            n_units, unit, n_units - 1L, n_quantities, plural)
  }
  stop(errorCondition(
    paste(c(sprintf("%s needs more %s than the %d quantit%s it covers; this fit has %d.",
                    what, unit, n_quantities, plural, n_units),
            reason, alternative), collapse = " "),
    class = "multilpa_too_few_groups", call = NULL))
}
