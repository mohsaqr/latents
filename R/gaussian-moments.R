#' Evaluate observed Gaussian densities and conditional sufficient statistics
#' @param x Numeric observation matrix, possibly containing NA.
#' @param parameters Gaussian mixture parameter list.
#' @return Log densities and conditional first moments and covariance adjustments.
#' @noRd
.multilpa_gaussian_moments <- function(x, parameters) {
  stopifnot(is.matrix(x), is.numeric(x), is.list(parameters),
            !any(is.infinite(x)))
  patterns <- split(seq_len(nrow(x)), apply(!is.na(x), 1L, paste0, collapse = ""))
  components <- lapply(seq_len(nrow(parameters$means)), function(profile) {
    covariance <- if (is.null(parameters$covariances)) {
      diag(parameters$variances[profile, ], ncol(x))
    } else matrix(parameters$covariances[, , profile], ncol(x), ncol(x))
    expected <- x
    log_density <- numeric(nrow(x))
    adjustments <- lapply(patterns, function(rows) {
      observed <- which(!is.na(x[rows[1L], ]))
      missing <- which(is.na(x[rows[1L], ]))
      conditional_covariance <- matrix(0, ncol(x), ncol(x))
      if (length(observed)) {
        residuals <- sweep(x[rows, observed, drop = FALSE], 2L,
                           parameters$means[profile, observed], "-")
        chol_observed <- chol(covariance[observed, observed, drop = FALSE])
        standardized <- residuals %*% backsolve(chol_observed, diag(length(observed)))
        log_density[rows] <<- -0.5 * (length(observed) * log(2 * pi) +
          2 * sum(log(diag(chol_observed))) + rowSums(standardized^2))
        if (length(missing)) {
          regression <- covariance[missing, observed, drop = FALSE] %*%
            chol2inv(chol_observed)
          expected[rows, missing] <<- sweep(residuals %*% t(regression), 2L,
                                            parameters$means[profile, missing], "+")
          conditional_covariance[missing, missing] <-
            covariance[missing, missing, drop = FALSE] -
            regression %*% covariance[observed, missing, drop = FALSE]
        }
      } else {
        expected[rows, ] <<- matrix(parameters$means[profile, ], length(rows),
                                    ncol(x), byrow = TRUE)
        conditional_covariance <- covariance
      }
      list(rows = rows, covariance = conditional_covariance)
    })
    list(expected = expected, adjustments = adjustments, log_density = log_density)
  })
  list(log_density = matrix(vapply(components, `[[`, numeric(nrow(x)), "log_density"),
                             nrow(x), length(components)), moments = components)
}

#' Constrain a covariance matrix by its minimum eigenvalue
#' @param covariance Symmetric covariance matrix.
#' @param min_variance Positive eigenvalue lower bound.
#' @return Symmetric positive definite covariance matrix.
#' @noRd
.multilpa_bound_covariance <- function(covariance, min_variance) {
  stopifnot(is.matrix(covariance), is.numeric(covariance), min_variance > 0,
            nrow(covariance) == ncol(covariance), all(is.finite(covariance)))
  decomposition <- eigen((covariance + t(covariance)) / 2, symmetric = TRUE)
  if (min(decomposition$values) >= min_variance) return(covariance)
  bounded <- sweep(decomposition$vectors, 2L,
                   pmax(decomposition$values, min_variance), "*") %*%
    t(decomposition$vectors)
  (bounded + t(bounded)) / 2
}

#' Maximize Gaussian parameters from conditional moments
#' @param x Numeric observation matrix.
#' @param expectation Nested posteriors containing conditional Gaussian moments.
#' @param variance_model Equal or varying covariance across profiles.
#' @param min_variance Variance or eigenvalue lower bound.
#' @param covariance_model Diagonal or full covariance.
#' @param held Parameters to hold at supplied values, or `NULL` to estimate
#'   every one of them. Only `means`, `variances` and `covariances` are read.
#' @return Means, diagonal variances, and optional covariance array. A held
#'   block is returned as supplied.
#' @noRd
.multilpa_maximize_moments <- function(x, expectation, variance_model,
                                     min_variance, covariance_model,
                                     held = NULL, structure = NULL,
                                     previous = NULL) {
  stopifnot(is.matrix(x), is.list(expectation), min_variance > 0,
            variance_model %in% c("varying", "equal"),
            covariance_model %in% c("diagonal", "full"),
            is.null(held) || is.list(held),
            is.null(structure) || structure %in% .multilpa_structures())
  weights <- colSums(expectation$subject_posteriors)
  n_indicators <- ncol(x)
  components <- lapply(seq_along(weights), function(profile) {
    moments <- expectation$gaussian_moments[[profile]]
    responsibility <- expectation$subject_posteriors[, profile]
    # A held mean is the point the spread is measured around. Estimating the
    # mean first and substituting afterwards would centre the residuals on the
    # free weighted mean and report a spread that is too small.
    means <- if (is.null(held$means)) {
      colSums(moments$expected * responsibility) / weights[profile]
    } else held$means[profile, ]
    residuals <- sweep(moments$expected, 2L, means, "-")
    covariance_sum <- crossprod(residuals, residuals * responsibility) +
      Reduce(`+`, lapply(moments$adjustments, function(adjustment) {
        adjustment$covariance * sum(responsibility[adjustment$rows])
      }))
    list(means = means, covariance_sum = covariance_sum)
  })
  means <- t(matrix(vapply(components, `[[`, numeric(n_indicators), "means"),
                     n_indicators, length(weights)))
  covariance_sums <- lapply(components, `[[`, "covariance_sum")
  shared <- if (variance_model == "equal") Reduce(`+`, covariance_sums) / sum(weights) else NULL
  # The constrained diagonals tie the profiles together -- an equal volume is
  # an average over all of them -- so they are solved once for every profile
  # rather than profile by profile.
  constrained_diagonal <- !is.null(structure) &&
    identical(covariance_model, "diagonal") && !structure %in% c("EEI", "VVI")
  constrained_full <- .multilpa_is_ellipsoidal(structure) &&
    !structure %in% c("EEE", "VVV")
  structured <- if (constrained_diagonal) {
    .multilpa_structure_variances(.multilpa_scatter_diagonals(covariance_sums),
                                  weights, structure, min_variance)
  } else if (constrained_full) {
    # The previous iteration's answer seeds the structures that iterate, so
    # each M-step continues where the last one stopped rather than restarting.
    .multilpa_structure_covariances(covariance_sums, weights, structure,
                                    min_variance, start = previous$covariances)
  } else NULL
  covariances <- lapply(seq_along(weights), function(profile) {
    held_covariance <- if (!is.null(held$covariances)) {
      matrix(held$covariances[, , profile], n_indicators, n_indicators)
    } else if (!is.null(held$variances)) {
      diag(held$variances[profile, ], n_indicators)
    } else NULL
    if (!is.null(held_covariance)) return(held_covariance)
    if (constrained_diagonal) return(diag(structured[profile, ], n_indicators))
    if (constrained_full) {
      return(matrix(structured[, , profile], n_indicators, n_indicators))
    }
    covariance <- shared %||% (covariance_sums[[profile]] / weights[profile])
    if (covariance_model == "full") .multilpa_bound_covariance(covariance, min_variance)
    else diag(pmax(diag(covariance), min_variance), n_indicators)
  })
  result <- list(means = means,
    variances = t(matrix(vapply(covariances, diag, numeric(n_indicators)),
                          n_indicators, length(weights))))
  if (covariance_model == "full") {
    result$covariances <- array(unlist(covariances, use.names = FALSE),
                                c(n_indicators, n_indicators, length(weights)))
  }
  result
}

#' Detect constrained Gaussian covariance estimates
#' @param parameters Gaussian parameter list.
#' @param min_variance Variance or eigenvalue lower bound.
#' @return Whether a variance or covariance eigenvalue is on its lower bound.
#' @noRd
.multilpa_covariance_boundary <- function(parameters, min_variance) {
  stopifnot(is.list(parameters), min_variance > 0)
  if (is.null(parameters$covariances)) {
    return(any(parameters$variances <= min_variance * (1 + 1e-8)))
  }
  any(vapply(seq_len(nrow(parameters$means)), function(profile) {
    covariance <- matrix(parameters$covariances[, , profile], ncol(parameters$means))
    min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values) <=
      min_variance * (1 + 1e-8)
  }, logical(1)))
}
