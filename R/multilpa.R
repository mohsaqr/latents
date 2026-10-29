#' Row-wise log-sum-exp
#' @param log_values Numeric matrix of log weights.
#' @return Vector of log sums.
#' @noRd
.multilpa_log_sum_exp <- function(log_values) {
  stopifnot(is.matrix(log_values), is.numeric(log_values),
            nrow(log_values) > 0L, ncol(log_values) > 0L,
            !anyNA(log_values))
  row_max <- .multilpa_row_max(log_values)
  if (any(!is.finite(row_max))) {
    stop(errorCondition("All component densities vanished, or a density overflowed.",
                        class = "latents_vanished_density", call = NULL))
  }
  row_max + log(rowSums(exp(sweep(log_values, 1L, row_max, "-"))))
}

#' Compute nested posterior probabilities
#' @param x Observation-by-indicator matrix of continuous indicators, optionally
#'   containing NA. It may have zero columns when every indicator is categorical.
#' @param group_index Integer group indices in first-occurrence order.
#' @param parameters List of model parameters.
#' @param codes Optional integer matrix of categorical indicator codes.
#' @param weights Optional sampling weight per group, already scaled.
#' @return Log likelihood, posteriors, and conditional class probabilities.
#'   With `weights`, the log likelihood is the weighted pseudo log likelihood
#'   and every posterior is multiplied by its group's weight, so the M-step
#'   and the scores read weighted counts without knowing about weights;
#'   `group_log_likelihood` stays per group and unweighted.
#' @noRd
.multilpa_expectation <- function(x, group_index, parameters, codes = NULL,
                                  weights = NULL, extra = NULL) {
  stopifnot(is.matrix(x), is.numeric(x), !any(is.infinite(x)),
            length(group_index) == nrow(x), is.list(parameters),
            is.null(weights) || length(weights) == max(group_index))
  # A noise component is one more column: a constant density over the data's
  # hypervolume, with its own mixing proportion last in the profile row.
  measurement <- .latents_measurement_log_density(
    x, parameters, codes, extra, noise_log_density = parameters$noise_log_density)
  result <- c(.multilpa_nested_structure(measurement$log_density, measurement$offset,
                                         group_index, parameters),
              list(gaussian_moments = measurement$moments))
  if (is.null(weights)) result else .multilpa_weigh(result, weights, group_index)
}

#' Maximize the expected complete-data log likelihood
#' @param x Observation-by-indicator matrix.
#' @param expectation Nested posterior calculation.
#' @param variance_model Either varying or equal across profiles.
#' @param min_variance Lower bound on every indicator variance.
#' @param covariance_model Diagonal or full residual covariance.
#' @param covariance_pattern `NULL`, or a logical indicators-by-indicators
#'   matrix marking the free covariance entries (the diagonal and the
#'   `local_dependence` pairs); every other entry is held at zero.
#' @return Updated model parameters.
#' @noRd
.multilpa_maximization <- function(x, expectation, variance_model, min_variance,
                                  covariance_model = "diagonal", codes = NULL,
                                  n_categories = NULL, min_probability = 1e-10,
                                  held = NULL, structure = NULL,
                                  previous = NULL, prior = NULL, extra = NULL,
                                  covariance_pattern = NULL) {
  stopifnot(is.matrix(x), is.list(expectation),
            variance_model %in% c("varying", "equal"), min_variance > 0,
            is.null(held) || is.list(held))
  # With direct effects this is an ECM step: the means and variances are
  # maximized for the indicators with the current effects removed, then the
  # effects for the new means and variances (R/direct-effects.R).
  unshifted <- x
  x <- .latents_direct_shift(x, extra, previous)
  weights <- colSums(expectation$subject_posteriors)
  group_weights <- colSums(expectation$group_posteriors)
  if (any(weights <= 0) || any(group_weights <= 0)) {
    stop("A profile or group class has zero effective membership.")
  }
  # A noise component (the last posterior column) has a mixing proportion but
  # no Gaussian parameters, so the Gaussian M-step sees the other columns only
  # and its shared blocks average over their total weight, not over n.
  noisy <- !is.null(previous$noise_log_density)
  gaussian_expectation <- expectation
  if (noisy) {
    gaussian_expectation$subject_posteriors <- expectation$subject_posteriors[
      , -ncol(expectation$subject_posteriors), drop = FALSE]
    weights <- colSums(gaussian_expectation$subject_posteriors)
  }
  gaussian <- if (!is.null(prior)) {
    # Maximum a posteriori under mclust's conjugate prior; complete data only.
    .multilpa_prior_maximize(x, gaussian_expectation$subject_posteriors, structure,
                             prior, min_variance)
  } else if (!is.null(expectation$gaussian_moments)) {
    .multilpa_maximize_moments(x, gaussian_expectation, variance_model, min_variance,
                               covariance_model, held, structure, previous,
                               covariance_pattern = covariance_pattern)
  } else NULL
  if (ncol(x) == 0L) {
    means <- matrix(numeric(0), length(weights), 0L)
    variances <- matrix(numeric(0), length(weights), 0L)
  } else if (is.null(gaussian)) {
    # A held mean is the point the spread is measured around, so it replaces the
    # weighted mean before the residuals are formed rather than afterwards.
    means <- held$means %||%
      sweep(crossprod(gaussian_expectation$subject_posteriors, x), 1L, weights, "/")
    variance_sums <- t(matrix(vapply(seq_along(weights), function(profile) {
      residuals <- sweep(x, 2L, means[profile, ], "-")
      colSums(residuals^2 * gaussian_expectation$subject_posteriors[, profile])
    }, numeric(ncol(x))), nrow = ncol(x), ncol = length(weights)))
    # A constrained structure ties the profiles together, so it is solved for
    # all of them at once rather than profile by profile.
    constrained <- !is.null(structure) && identical(covariance_model, "diagonal") &&
      !structure %in% c("EEI", "VVI")
    variances <- if (!is.null(held$variances)) held$variances else
      if (constrained) {
        .multilpa_structure_variances(variance_sums, weights, structure,
                                      min_variance)
      } else if (variance_model == "varying") {
        sweep(variance_sums, 1L, weights, "/")
      } else {
        matrix(colSums(variance_sums) /
                 if (noisy) sum(weights) else expectation$n_total %||% nrow(x),
               length(weights), ncol(x), byrow = TRUE)
      }
    # This is the exact M-step for the stated variance >= min_variance constraint.
    # It is not a hidden ridge, and bound-active estimates are reported to callers.
    # A held variance is returned as supplied, bound or not.
    if (is.null(held$variances)) variances <- pmax(variances, min_variance)
  } else {
    means <- gaussian$means
    variances <- gaussian$variances
  }
  n_columns <- ncol(expectation$subject_posteriors)
  profile_counts <- t(matrix(vapply(expectation$joint, colSums, numeric(n_columns)),
                             nrow = n_columns, ncol = length(group_weights)))
  profile_probabilities <- profile_counts / rowSums(profile_counts)
  parameters <- list(means = means, variances = variances,
                     profile_probabilities = profile_probabilities,
                     group_probabilities = group_weights / sum(group_weights))
  if (!is.null(gaussian$covariances)) parameters$covariances <- gaussian$covariances
  if (noisy) parameters$noise_log_density <- previous$noise_log_density
  if (!is.null(codes)) {
    parameters$response_probabilities <- held$response_probabilities %||%
      .multilpa_categorical_maximize(codes, expectation$subject_posteriors,
                                     n_categories, min_probability)
  }
  if (!is.null(extra)) {
    blocks <- .latents_extra_maximize(extra, expectation$subject_posteriors, previous)
    parameters[names(blocks)] <- blocks
  }
  if (!is.null(extra$direct)) {
    parameters$direct_effects <- .latents_direct_maximize(
      unshifted, extra, expectation$subject_posteriors, means, variances)
  }
  if (any(!is.finite(unlist(parameters, use.names = FALSE)))) {
    stop("Non-finite parameters in the M-step; check indicator scales.")
  }
  parameters
}

#' Initialize a two-level Gaussian mixture
#' @param x Observation-by-indicator matrix.
#' @param group_index Integer group indices.
#' @param n_profiles Number of individual profiles.
#' @param n_types Number of group classes.
#' @param variance_model Variance constraint.
#' @param min_variance Variance lower bound.
#' @param start_index Restart number.
#' @param covariance_model Diagonal or full residual covariance.
#' @param row_weights Optional sampling weights in observation order. Only
#'   positive-weight rows supply starting values; their relative weights do
#'   not change the existing initialization algorithm.
#' @return Strictly positive initial mixture probabilities and Gaussian parameters.
#' @noRd
.multilpa_initialize <- function(x, group_index, n_profiles, n_types,
                               variance_model, min_variance, start_index,
                               covariance_model = "diagonal", codes = NULL,
                               n_categories = NULL, min_probability = 1e-10,
                               hierarchical = FALSE, extra = NULL,
                               kind = "kmeans", row_weights = NULL) {
  stopifnot(is.matrix(x), length(group_index) == nrow(x),
            n_profiles >= 1L, n_types >= 1L, start_index >= 1L,
            variance_model %in% c("varying", "equal"), min_variance > 0,
            kind %in% c("kmeans", "partition", "centres"))
  stopifnot(is.null(row_weights) ||
              (is.numeric(row_weights) && length(row_weights) == nrow(x) &&
                 all(is.finite(row_weights)) && all(row_weights >= 0)))
  if (!is.null(row_weights) && any(row_weights == 0)) {
    # Initial states contain parameters only. Keep the full observation and
    # group mapping for EM, but let excluded units influence none of the
    # imputation, clustering, group compositions or measurement parameters.
    rows <- row_weights > 0
    x <- x[rows, , drop = FALSE]
    group_index <- group_index[rows]
    if (!is.null(codes)) codes <- codes[rows, , drop = FALSE]
    extra <- .latents_extra_rows(extra, rows)
    if (n_profiles > nrow(unique(cbind(x, codes, .latents_extra_matrix(extra))))) {
      stop(errorCondition(paste(
        "n_profiles cannot exceed the number of distinct observed indicator",
        "rows with positive sampling weight."),
        class = "latents_unidentified", call = NULL))
    }
  }
  # The random kinds of the emEM pool describe each cluster by its own
  # spread, as Ward's start does, and draw the group classes at random.
  random_kind <- !identical(kind, "kmeans")
  if (ncol(x) > 0L && anyNA(x)) {
    initial_means <- colMeans(x, na.rm = TRUE)
    if (any(!is.finite(initial_means))) {
      stop(errorCondition(paste(
        "Every continuous indicator needs an observed row with positive",
        "sampling weight for initialization."),
        class = "latents_bad_data", call = NULL))
    }
    x <- matrix(vapply(seq_len(ncol(x)), function(indicator) {
      value <- x[, indicator]
      value[is.na(value)] <- initial_means[indicator]
      value
    }, numeric(nrow(x))), nrow(x), ncol(x))
  }
  overall_means <- colMeans(x)
  overall_variances <- colMeans(sweep(x, 2L, overall_means, "-")^2)
  if (n_profiles == 1L) {
    means <- matrix(overall_means, 1L)
    assignments <- rep(1L, nrow(x))
  } else {
    # Standardizing only for initialization stops a large-scale indicator from
    # determining all initial clusters; fitting remains in the original units.
    scales <- sqrt(overall_variances)
    if (!is.null(row_weights) && any(row_weights == 0)) scales[scales == 0] <- 1
    standardized <- sweep(sweep(x, 2L, overall_means, "-"),
                          2L, scales, "/")
    # Categorical indicators join the clustering as centered dummy columns, so
    # a mixed-mode model is initialized from all of its indicators at once.
    clustering_input <- if (is.null(codes)) standardized else
      cbind(standardized, .multilpa_categorical_design(codes, n_categories))
    clustering_input <- cbind(clustering_input, .latents_extra_clustering(extra))
    if (hierarchical || random_kind) {
      assignments <- if (hierarchical) {
        .multilpa_ward_assignments(clustering_input, n_profiles)
      } else .multilpa_random_assignments(clustering_input, n_profiles, kind)
      centers <- t(vapply(seq_len(n_profiles), function(profile) {
        colMeans(clustering_input[assignments == profile, seq_len(ncol(x)),
                                  drop = FALSE])
      }, numeric(ncol(x))))
      centers <- matrix(centers, n_profiles, ncol(x))
    } else {
      initial_fit <- .multilpa_start_kmeans(clustering_input, n_profiles)
      if (!is.null(initial_fit$ifault) && initial_fit$ifault != 0L) {
        stop(paste("The k-means clustering that draws this start's initial values",
                   "did not converge in 100 iterations; the start was skipped."))
      }
      assignments <- initial_fit$cluster
      centers <- initial_fit$centers[, seq_len(ncol(x)), drop = FALSE]
    }
    means <- sweep(sweep(centers, 2L, sqrt(overall_variances), "*"),
                   2L, overall_means, "+")
  }
  variances <- matrix(pmax(overall_variances, min_variance),
                      n_profiles, ncol(x), byrow = TRUE)
  # The hierarchical start describes each cluster by its own spread, which is
  # what lets it reach solutions whose profiles differ in shape. A cluster too
  # small to estimate a spread keeps the overall one.
  cluster_spread <- function(profile) {
    rows <- assignments == profile
    if (sum(rows) <= ncol(x) + 1L) return(NULL)
    stats::cov(x[rows, , drop = FALSE]) * (sum(rows) - 1) / sum(rows)
  }
  spreads <- if ((hierarchical || random_kind) && n_profiles > 1L && ncol(x) > 0L) {
    lapply(seq_len(n_profiles), cluster_spread)
  } else NULL
  if (!is.null(spreads)) {
    variances <- t(vapply(seq_len(n_profiles), function(profile) {
      if (is.null(spreads[[profile]])) variances[profile, ] else
        pmax(diag(spreads[[profile]]), min_variance)
    }, numeric(ncol(x))))
    variances <- matrix(variances, n_profiles, ncol(x))
    if (identical(variance_model, "equal")) {
      variances <- matrix(colMeans(variances), n_profiles, ncol(x),
                          byrow = TRUE)
    }
  }
  indicator <- vapply(seq_len(n_profiles), function(profile) {
    as.numeric(assignments == profile)
  }, numeric(nrow(x)))
  group_counts <- rowsum(indicator, group_index, reorder = FALSE)
  # Pseudocounts are used only to initialize, never to estimate the final model.
  group_composition <- (group_counts + 0.5) /
    (rowSums(group_counts) + 0.5 * n_profiles)
  if (n_types == 1L) {
    profile_probabilities <- matrix(colMeans(indicator), 1L)
  } else if (hierarchical) {
    # Ward's start is deterministic at both levels: its group classes are
    # Ward's clusters of the group compositions, not k-means from random
    # centres, so start 1 is the same whatever the seed (review 2026-10-26
    # round 2, ST-1).
    profile_probabilities <- .multilpa_group_class_seed(group_composition, n_types)
  } else if (!random_kind && start_index %% 2L == 1L &&
             nrow(unique(group_composition)) >= n_types) {
    group_fit <- .multilpa_start_kmeans(group_composition, n_types)
    if (!is.null(group_fit$ifault) && group_fit$ifault != 0L) {
      stop(paste("The k-means clustering that draws this start's initial group-class",
                 "values did not converge in 100 iterations; the start was skipped."))
    }
    profile_probabilities <- group_fit$centers
  } else {
    profile_probabilities <- matrix(stats::runif(n_types * n_profiles, 0.1, 1),
                                    n_types, n_profiles)
  }
  profile_probabilities <- profile_probabilities / rowSums(profile_probabilities)
  result <- list(means = means, variances = variances,
       profile_probabilities = profile_probabilities,
       group_probabilities = rep(1 / n_types, n_types))
  if (!is.null(codes)) {
    result$response_probabilities <- .multilpa_categorical_initialize(
      codes, assignments, n_profiles, n_categories, min_probability)
  }
  if (!is.null(extra)) {
    result <- c(result, .latents_extra_initialize(extra, assignments, n_profiles))
  }
  if (ncol(x) > 0L && covariance_model == "full") {
    residuals <- sweep(x, 2L, overall_means, "-")
    covariance <- .multilpa_bound_covariance(crossprod(residuals) / nrow(x), min_variance)
    blocks <- rep(list(covariance), n_profiles)
    if (!is.null(spreads)) {
      blocks <- lapply(seq_len(n_profiles), function(profile) {
        if (is.null(spreads[[profile]])) covariance else
          .multilpa_bound_covariance(spreads[[profile]], min_variance)
      })
      if (identical(variance_model, "equal")) {
        blocks <- rep(list(Reduce(`+`, blocks) / n_profiles), n_profiles)
      }
    }
    result$covariances <- array(unlist(blocks, use.names = FALSE),
                                c(ncol(x), ncol(x), n_profiles))
    result$variances <- t(vapply(blocks, diag, numeric(ncol(x))))
    result$variances <- matrix(result$variances, n_profiles, ncol(x))
  }
  result
}

#' Deterministic group-class starting values from group compositions
#'
#' Ward's clustering of the groups' profile compositions, each group class
#' starting at its cluster's mean composition. With fewer distinct
#' compositions than group classes, each class tilts the overall composition
#' towards one profile by a different amount, so the classes still start
#' apart.
#'
#' @param group_composition Groups-by-profiles matrix of smoothed profile
#'   shares, rows summing to one.
#' @param n_types Number of group classes, at least two.
#' @return A group-classes-by-profiles matrix of positive shares.
#' @noRd
.multilpa_group_class_seed <- function(group_composition, n_types) {
  stopifnot("`group_composition` must be a numeric matrix" =
              is.matrix(group_composition) && is.numeric(group_composition),
            "`n_types` must be at least two" = n_types >= 2L)
  n_profiles <- ncol(group_composition)
  if (nrow(unique(group_composition)) >= n_types) {
    labels <- .multilpa_ward_assignments(group_composition, n_types)
    return(rowsum(group_composition, labels, reorder = TRUE) /
             as.vector(tabulate(labels, n_types)))
  }
  overall <- colMeans(group_composition)
  t(vapply(seq_len(n_types), function(group_class) {
    tilt <- group_class / (n_types + 1)
    target <- (group_class - 1L) %% n_profiles + 1L
    (1 - tilt) * overall + tilt * as.numeric(seq_len(n_profiles) == target)
  }, numeric(n_profiles)))
}

#' Internal Gaussian origin from the rows that enter the objective
#' @param x Continuous indicator matrix.
#' @param row_weights Optional nonnegative weights in observation order.
#' @return Observed column means among positive-weight rows.
#' @noRd
.multilpa_initialization_center <- function(x, row_weights = NULL) {
  stopifnot(is.matrix(x), is.null(row_weights) ||
              (is.numeric(row_weights) && length(row_weights) == nrow(x) &&
                 all(is.finite(row_weights)) && all(row_weights >= 0)))
  if (ncol(x) == 0L) return(numeric(0))
  rows <- if (is.null(row_weights)) rep(TRUE, nrow(x)) else row_weights > 0
  centers <- colMeans(x[rows, , drop = FALSE], na.rm = TRUE)
  if (any(!is.finite(centers))) {
    stop(errorCondition(paste(
      "Every continuous indicator needs an observed row with positive",
      "sampling weight for initialization."),
      class = "latents_bad_data", call = NULL))
  }
  centers
}

#' Assign rows to clusters at random, for the emEM pool
#'
#' `"partition"` deals the rows into the clusters in a random order, so every
#' cluster holds a share; `"centres"` draws distinct rows as the cluster
#' centres and assigns every row to the nearest, each drawn row keeping its
#' own cluster (so none is empty even when rows repeat).
#'
#' @param input Standardized numeric matrix, one row per observation.
#' @param n_clusters Number of clusters.
#' @param kind `"partition"` or `"centres"`.
#' @return Integer cluster labels, one per row, every cluster non-empty.
#' @noRd
.multilpa_random_assignments <- function(input, n_clusters, kind) {
  stopifnot("`input` must be a numeric matrix" = is.matrix(input),
            "`n_clusters` must not exceed the rows" = n_clusters <= nrow(input),
            "`kind` must be \"partition\" or \"centres\"" =
              kind %in% c("partition", "centres"))
  if (identical(kind, "partition")) {
    return(sample(rep_len(seq_len(n_clusters), nrow(input))))
  }
  drawn <- sample.int(nrow(input), n_clusters)
  distance <- vapply(drawn, function(row) {
    colSums((t(input) - input[row, ])^2)
  }, numeric(nrow(input)))
  assignments <- max.col(-matrix(distance, nrow(input)), ties.method = "first")
  assignments[drawn] <- seq_len(n_clusters)
  assignments
}

#' Assign rows to clusters by Ward's hierarchical clustering
#'
#' Ward's Euclidean agglomeration on standardized rows. Agglomeration needs
#' every pairwise distance, so at most 2000 evenly spaced rows are clustered
#' and every row is then assigned to the nearest cluster centre.
#'
#' @param input Standardized numeric matrix, one row per observation.
#' @param n_clusters Number of clusters.
#' @return Integer cluster labels, one per row, every cluster non-empty.
#' @noRd
.multilpa_ward_assignments <- function(input, n_clusters) {
  stopifnot("`input` must be a numeric matrix" = is.matrix(input),
            "`n_clusters` must not exceed the rows" = n_clusters <= nrow(input))
  used <- unique(round(seq(1, nrow(input),
                           length.out = min(nrow(input), 2000L))))
  tree <- stats::hclust(stats::dist(input[used, , drop = FALSE]), "ward.D2")
  labels <- stats::cutree(tree, n_clusters)
  centres <- rowsum(input[used, , drop = FALSE], labels) /
    as.vector(tabulate(labels, n_clusters))
  distance <- vapply(seq_len(n_clusters), function(cluster) {
    colSums((t(input) - centres[cluster, ])^2)
  }, numeric(nrow(input)))
  assignments <- max.col(-matrix(distance, nrow(input)), ties.method = "first")
  # Nearest-centre assignment can leave a small cluster empty on the full
  # data; its sampled members keep it.
  empty <- setdiff(seq_len(n_clusters), assignments)
  if (length(empty) > 0L) {
    keep <- labels %in% empty
    assignments[used[keep]] <- labels[keep]
  }
  assignments
}

#' Validate user-supplied initial parameters
#' @param start Initial parameter list.
#' @param n_profiles Number of profiles.
#' @param n_types Number of group classes.
#' @param n_indicators Number of indicators.
#' @param variance_model Variance constraint.
#' @param min_variance Variance lower bound.
#' @param covariance_model Diagonal or full residual covariance.
#' @param n_categories Integer vector of category counts, one per categorical
#'   indicator, or `NULL` when every indicator is continuous.
#' @param min_probability Lower bound applied to response probabilities.
#' @param categorical Indicator names, in the order this fit consumes them, or
#'   `NULL` when no alignment by label is possible.
#' @param levels Named list of category labels per categorical indicator, in
#'   code order, or `NULL`.
#' @return A validated copy without dimension names, with any labelled
#'   categorical block aligned to this fit's indicator and category order.
#' @noRd
.multilpa_validate_start <- function(start, n_profiles, n_types, n_indicators,
                                   variance_model, min_variance, covariance_model = "diagonal",
                                   n_categories = NULL, min_probability = 1e-10,
                                   categorical = NULL, levels = NULL) {
  stopifnot(is.list(start), n_profiles >= 1L, n_types >= 1L,
            n_indicators >= 0L, min_variance > 0,
            "An all-categorical model needs `n_categories`" =
              n_indicators >= 1L || !is.null(n_categories))
  covariances <- NULL
  if (covariance_model == "full" && n_indicators > 0L) {
    covariances <- start$covariances
    if (!is.numeric(covariances) ||
        !identical(as.integer(dim(covariances)), c(as.integer(n_indicators),
          as.integer(n_indicators), as.integer(n_profiles))) || any(!is.finite(covariances))) {
      stop(errorCondition(
      "start$covariances must be a finite indicators x indicators x profiles array.",
      class = "latents_bad_start", call = NULL))
    }
    invisible(lapply(seq_len(n_profiles), function(profile) {
      covariance <- matrix(covariances[, , profile], n_indicators, n_indicators)
      if (max(abs(covariance - t(covariance))) > 1e-12 ||
          min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values) < min_variance * (1 - 1e-8)) {
        stop(errorCondition(
          "Starting covariances must be symmetric with eigenvalues at least min_variance.",
          class = "latents_bad_start", call = NULL))
      }
      if (variance_model == "equal" && max(abs(covariance - covariances[, , 1L])) > 1e-12) {
        stop(errorCondition("Starting covariances must be equal across profiles.",
                            class = "latents_bad_start", call = NULL))
      }
    }))
    variances <- t(matrix(vapply(seq_len(n_profiles), function(profile) {
      diag(matrix(covariances[, , profile], n_indicators, n_indicators))
    }, numeric(n_indicators)), n_indicators, n_profiles))
    if (!is.null(start$variances) && !isTRUE(all.equal(unname(start$variances), variances, tolerance = 1e-12))) {
      stop(errorCondition("Starting variances must agree with covariance diagonals.",
                          class = "latents_bad_start", call = NULL))
    }
    start$variances <- variances
    start$covariances <- NULL
  }
  # An all-categorical model has no Gaussian block; accept its absence rather
  # than making the caller assemble empty matrices.
  if (n_indicators == 0L) {
    if (is.null(start$means)) start$means <- matrix(0, n_profiles, 0L)
    if (is.null(start$variances)) start$variances <- matrix(0, n_profiles, 0L)
  }
  required <- c("means", "variances", "profile_probabilities", "group_probabilities")
  if (!is.null(n_categories)) required <- c(required, "response_probabilities")
  if (!setequal(names(start), required) || anyDuplicated(names(start))) {
    stop(errorCondition(sprintf("start must contain exactly %s.",
                                paste(required, collapse = ", ")),
                        class = "latents_bad_start", call = NULL))
  }
  dimensions <- list(means = c(n_profiles, n_indicators),
                     variances = c(n_profiles, n_indicators),
                     profile_probabilities = c(n_types, n_profiles))
  invisible(lapply(names(dimensions), function(field) {
    value <- start[[field]]
    if (!is.matrix(value) || !is.numeric(value) ||
        !identical(as.integer(dim(value)), as.integer(dimensions[[field]])) ||
        any(!is.finite(value))) {
      stop(errorCondition(sprintf(
        "start$%s must be a finite numeric matrix with dimensions %s.",
        field, paste(dimensions[[field]], collapse = " x ")),
        class = "latents_bad_start", call = NULL))
    }
  }))
  if (any(start$variances < min_variance)) {
    stop(errorCondition("start$variances must be at least min_variance.",
                        class = "latents_bad_start", call = NULL))
  }
  if (variance_model == "equal" &&
      any(abs(sweep(start$variances, 2L, start$variances[1L, ], "-")) > 1e-12)) {
    stop(errorCondition(
      "start$variances must be equal across profiles for variance_model = 'equal'.",
      class = "latents_bad_start", call = NULL))
  }
  group_probabilities <- start$group_probabilities
  if (!is.numeric(group_probabilities) || !is.null(dim(group_probabilities)) ||
      length(group_probabilities) != n_types || any(!is.finite(group_probabilities)) ||
      any(group_probabilities <= 0) || abs(sum(group_probabilities) - 1) > 1e-8 ||
      any(start$profile_probabilities <= 0) ||
      any(abs(rowSums(start$profile_probabilities) - 1) > 1e-8)) {
    stop(errorCondition(
      "Starting probabilities must be strictly positive and sum to one (within each profile-probability row).",
      class = "latents_bad_start", call = NULL))
  }
  # Normalize only roundoff admitted by the validation tolerance.
  result <- list(means = unname(start$means), variances = unname(start$variances),
       profile_probabilities = unname(start$profile_probabilities /
                                       rowSums(start$profile_probabilities)),
       group_probabilities = unname(group_probabilities / sum(group_probabilities)))
  if (!is.null(covariances)) result$covariances <- unname(covariances)
  if (!is.null(n_categories)) {
    # Labels first, shape second: a block list carried over from another fit is
    # attached to the item and the category its labels name, not to whichever
    # position it happens to occupy here.
    result$response_probabilities <- .multilpa_validate_response_start(
      .multilpa_align_response_start(start$response_probabilities,
                                     categorical, levels),
      n_profiles, n_categories, min_probability)
  }
  result
}

#' Align a labelled item-response start with this fit's own encoding
#'
#' Response blocks are consumed by position: the first block belongs to the
#' first `categorical` indicator and its columns to that indicator's categories
#' in code order. A block list carried over from another fit, or from a first
#' stage, can therefore be attached to the wrong item, or a category to the
#' wrong column, and every later number is silently wrong rather than visibly
#' broken. Whenever the caller supplies labels -- names on the list, column
#' names on a block -- the labels are the contract: the block is reordered to
#' this fit's encoding, and a label naming an item or a category this fit does
#' not have is refused, because there is no alignment to make.
#'
#' Unlabelled input keeps the documented positional contract unchanged.
#'
#' @param response_probabilities The caller's block list, possibly labelled.
#' @param categorical Indicator names, in the order this fit consumes them.
#' @param levels Named list of category labels per indicator, in code order, or
#'   `NULL` to align the blocks only.
#' @return The block list, unnamed, in this fit's item order with each block's
#'   columns in this fit's category order. Raises `latents_bad_start` when a
#'   supplied label cannot be aligned.
#' @noRd
.multilpa_align_response_start <- function(response_probabilities, categorical,
                                           levels = NULL) {
  # A list of the wrong length, or an object that is not a list at all, is
  # reported by the shape validation that follows, in its own words.
  if (!is.list(response_probabilities) || length(categorical) == 0L ||
      length(response_probabilities) != length(categorical)) {
    return(response_probabilities)
  }
  block_names <- names(response_probabilities)
  if (!is.null(block_names) && !anyNA(block_names) && all(nzchar(block_names))) {
    if (anyDuplicated(block_names) > 0L || !setequal(block_names, categorical)) {
      stop(errorCondition(sprintf(
        "start$response_probabilities is labelled %s, but this fit's categorical indicators are %s.",
        paste(sprintf("\"%s\"", block_names), collapse = ", "),
        paste(sprintf("\"%s\"", categorical), collapse = ", ")),
        class = "latents_bad_start", call = NULL))
    }
    response_probabilities <- response_probabilities[categorical]
  }
  if (is.null(levels)) return(unname(response_probabilities))
  unname(lapply(seq_along(categorical), function(index) {
    block <- response_probabilities[[index]]
    wanted <- levels[[categorical[[index]]]]
    observed <- if (is.matrix(block)) colnames(block) else NULL
    if (is.null(observed) || is.null(wanted)) return(block)
    if (anyDuplicated(observed) > 0L || !setequal(observed, wanted)) {
      stop(errorCondition(sprintf(
        "start$response_probabilities for `%s` is labelled %s, but `%s` has categories %s here.",
        categorical[[index]],
        paste(sprintf("\"%s\"", observed), collapse = ", "),
        categorical[[index]],
        paste(sprintf("\"%s\"", wanted), collapse = ", ")),
        class = "latents_bad_start", call = NULL))
    }
    block[, match(wanted, observed), drop = FALSE]
  }))
}

#' Validate a user-supplied item-response starting block
#'
#' @param response_probabilities List of profiles-by-categories matrices, one
#'   per categorical indicator.
#' @param n_profiles Number of profiles.
#' @param n_categories Integer vector of category counts per indicator.
#' @param min_probability Lower bound applied to every cell.
#' @return The validated list, unnamed, with rows renormalized to sum to one.
#' @noRd
.multilpa_validate_response_start <- function(response_probabilities, n_profiles,
                                            n_categories, min_probability) {
  if (!is.list(response_probabilities) ||
      length(response_probabilities) != length(n_categories)) {
    stop(errorCondition(sprintf(
      "start$response_probabilities must be a list of %d matrices, one per categorical indicator.",
      length(n_categories)), class = "latents_bad_start", call = NULL))
  }
  unname(lapply(seq_along(n_categories), function(indicator) {
    block <- response_probabilities[[indicator]]
    if (!is.matrix(block) || !is.numeric(block) ||
        !identical(as.integer(dim(block)),
                   c(as.integer(n_profiles), as.integer(n_categories[[indicator]]))) ||
        any(!is.finite(block))) {
      stop(errorCondition(sprintf(
        "start$response_probabilities[[%d]] must be a finite numeric %d x %d matrix.",
        indicator, n_profiles, n_categories[[indicator]]),
        class = "latents_bad_start", call = NULL))
    }
    if (any(block < 0) || any(abs(rowSums(block) - 1) > 1e-8)) {
      stop(errorCondition(sprintf(
        "start$response_probabilities[[%d]] rows must be nonnegative and sum to one.",
        indicator), class = "latents_bad_start", call = NULL))
    }
    bounded <- t(apply(block / rowSums(block), 1L, .multilpa_bound_probabilities,
                       min_probability = min_probability))
    matrix(bounded, n_profiles, n_categories[[indicator]])
  }))
}

#' Run one nested EM optimization
#' @param x Observation-by-indicator matrix.
#' @param group_index Integer group indices.
#' @param parameters Starting model parameters.
#' @param variance_model Variance constraint.
#' @param min_variance Variance lower bound.
#' @param max_iter Maximum M-steps.
#' @param tol Relative likelihood convergence tolerance.
#' @param covariance_model Diagonal or full residual covariance.
#' @param covariance_pattern `NULL`, or a logical indicators-by-indicators
#'   matrix marking the free covariance entries (the diagonal and the
#'   `local_dependence` pairs); every other entry is held at zero.
#' @param pause,resume As for `.latents_em()`: pause after that many
#'   iterations, or resume a paused run (then `parameters` is not used).
#' @return Final parameters, posterior calculation, and convergence history;
#'   with `pause` or `resume`, also the `checkpoint`.
#' @noRd
.multilpa_em <- function(x, group_index, parameters, variance_model,
                        min_variance, max_iter, tol, covariance_model = "diagonal",
                        codes = NULL, n_categories = NULL, min_probability = 1e-10,
                        held = NULL, structure = NULL, prior = NULL,
                        accelerate_em = TRUE, sampling_weights = NULL,
                        extra = NULL, covariance_pattern = NULL, pause = Inf,
                        resume = NULL) {
  stopifnot(is.matrix(x), is.list(parameters), max_iter >= 0L, tol > 0,
            length(group_index) == nrow(x), min_variance > 0,
            variance_model %in% c("varying", "equal"))
  # Under a prior the iteration climbs the posterior, not the likelihood, so
  # the likelihood may legitimately fall; convergence is still read off the
  # likelihood's relative change, which is mclust's rule. SQUAREM needs a
  # likelihood objective, so a prior fit takes plain steps.
  accelerator <- if (isTRUE(accelerate_em) && is.null(prior)) {
    .multilpa_accelerator(x, parameters,
      row_weights = if (is.null(sampling_weights)) NULL else
        sampling_weights[group_index], extra = extra)
  }
  fit <- .latents_em(
    parameters,
    evaluate = function(point) {
      .multilpa_expectation(x, group_index, point, codes, sampling_weights, extra)
    },
    # `point` warm-starts the covariance structures that iterate.
    maximize = function(point, expectation) {
      .multilpa_maximization(x, expectation, variance_model, min_variance,
                             covariance_model, codes, n_categories,
                             min_probability, held, structure, point, prior, extra,
                             covariance_pattern = covariance_pattern)
    },
    max_iter = max_iter, tol = tol,
    settings = .latents_em_settings(
      guard_decrease = is.null(prior), accelerate = accelerator,
      # A unit-free stop: the change in log likelihood against the number of
      # observed values (review 2026-10-26 round 2, K2-6).
      stop_scale = .latents_observed_values(
        x, codes, .latents_extra_matrix(extra),
        row_weights = if (is.null(sampling_weights)) NULL else
          sampling_weights[group_index])),
    pause = pause, resume = resume)
  result <- list(parameters = fit$state, expectation = fit$expectation,
                 converged = fit$converged, iterations = fit$iterations,
                 history = fit$history)
  if (!is.null(fit$checkpoint)) result$checkpoint <- fit$checkpoint
  result
}

#' Validate and extract an optional time column
#'
#' The model never uses the ordering; it is carried so that the assignments can
#' be read back in sequence without the caller supplying the data again.
#'
#' @param data The data frame passed to the fit.
#' @param time `NULL`, or the name of a single column giving each observation's
#'   position within its group.
#' @param id Name of the group column, used for the exclusion and the
#'   uniqueness check. It is known to name a column of `data` by the time this
#'   is called, because the data contract is checked first.
#' @param vars The indicator column names, which `time` may not name.
#'   Defaults to none so that callers which have not yet been updated keep
#'   working; `multilpa()` passes them.
#' @return `NULL` when `time` is `NULL`, otherwise the column in input row order.
#' @noRd
.multilpa_time_values <- function(data, time, id, vars = character()) {
  if (is.null(time)) return(NULL)
  .latents_require(
    "`time` must be NULL or a single column name" =
      is.character(time) && length(time) == 1L && !is.na(time))
  if (!time %in% names(data)) {
    stop(errorCondition(sprintf("`time` names `%s`, which `data` does not have.%s", time,
                                .latents_did_you_mean(time, names(data))),
                        class = "latents_bad_data", call = NULL))
  }
  .latents_require(
    "`time` must not be one of `vars` or the `id` column" =
      !identical(time, id) && !time %in% vars, .class = "latents_bad_data")
  values <- data[[time]]
  if (anyNA(values)) {
    stop(errorCondition("`time` must not contain missing values.",
                        class = "latents_bad_time", call = NULL))
  }
  # split() first converts numeric keys to factor labels, which can collapse
  # distinct numeric identifiers that print identically. Match the native keys.
  group_index <- match(data[[id]], unique(data[[id]]))
  duplicated_within <- any(vapply(split(values, group_index),
                                  function(v) anyDuplicated(v) > 0L, logical(1)))
  if (duplicated_within) {
    stop(errorCondition(
      "`time` must be unique within each group; a position cannot occur twice.",
      class = "latents_bad_time", call = NULL))
  }
  # Every reader of the ordering sorts these values, so numbers stored as text
  # whose alphabetical order is wrong are refused here, once, for every model.
  .lta_check_time_order(values, time)
  values
}

#' Resolve a single-level request into the two-level machinery
#'
#' `id = NULL` says the observations are independent: each row is its own unit,
#' which is the ordinary one-level mixture this model reduces to when every
#' group holds one row and there is one group class. Nothing about the
#' likelihood changes; only the grouping does.
#'
#' The column is fabricated here rather than by the caller. A caller who writes
#' `data$unit <- seq_len(nrow(data))` has to know that the fit stores that
#' column, that `get_results(x, "data")` will hand it back, and that a second call
#' must spell it the same way -- plumbing the package should own. This
#' package's own validation harness wrote that line once per comparison before
#' this argument existed.
#'
#' @param data The caller's data frame.
#' @param n_group_classes What the caller asked for, or `NULL` if they did not.
#' @param latent_class Whether every indicator is categorical, which makes the
#'   model a latent class rather than a latent profile analysis.
#' @return A list with `data`, `id` and `n_group_classes`.
#' @noRd
.multilpa_single_level <- function(data, n_group_classes,
                                   latent_class = FALSE) {
  stopifnot("`data` must be a data frame" = is.data.frame(data))
  name <- ".observation"
  if (name %in% names(data)) {
    stop(errorCondition(sprintf(paste(
      "`id = NULL` gives every observation its own unit and names that unit",
      "`%s`, which these data already have as a column. Rename it, or pass it",
      "as `id` if it is the nesting unit."), name),
      class = "latents_bad_data", call = NULL))
  }
  # Exactly one: `as.integer()` would let 1.5 or "1" through as a request for
  # one group class, which a fit with `id` refuses as not a positive integer.
  if (!is.null(n_group_classes) &&
      !(is.numeric(n_group_classes) && length(n_group_classes) == 1L &&
        identical(as.numeric(n_group_classes), 1))) {
    stop(errorCondition(sprintf(paste(
      "`id = NULL` puts one observation in each unit, so there is no",
      "composition for a second-level class to differ in and",
      "`n_group_classes = %s` has nothing to estimate. Pass `id` to fit group",
      "classes, or leave `n_group_classes` alone for a single-level fit."),
      paste(deparse(n_group_classes), collapse = "")),
      class = "latents_bad_argument", call = NULL))
  }
  .multilpa_single_level_notice(latent_class)
  data[[name]] <- seq_len(nrow(data))
  list(data = data, id = name, n_group_classes = 1L)
}

#' Fit a two-level latent profile model
#'
#' @description
#' Multilevel latent profile analysis describes heterogeneity among
#' observations and among the groups in which they are nested. Examples
#' include repeated measurements within individuals and students within
#' schools. Observation-level profiles describe patterns in the indicators,
#' while group-level classes describe differences in the probabilities of
#' those profiles.
#'
#' In the standard model, profile definitions are shared across group classes:
#' a given profile has the same indicator means and covariance structure in
#' every group class. Group classes therefore differ in their profile
#' composition. Observations within the same group can belong to different
#' profiles, with their membership probabilities depending on the group's
#' latent class.
#'
#' `multilpa()` estimates the measurement parameters, conditional profile
#' probabilities, group-class proportions, and posterior membership
#' probabilities at both levels. The default specification uses Gaussian
#' indicators with diagonal within-profile covariance matrices. Alternative
#' specifications support other covariance structures, discrete or mixed
#' indicators, and predictors of latent membership.
#'
#' @param data A data frame containing the indicators and the group identifier.
#' @param vars Unique character vector of indicator column names. Columns not
#'   named in `categorical`, `ordinal` or `count` are continuous indicators.
#' @param id Name of the group identifier column. Character, factor and
#'   numeric identifiers are supported; missing identifiers are not. `id` has
#'   no default: omitting it raises `latents_bad_argument`.
#'
#'   `id = NULL` fits a **single-level** model: each row is its own unit, the
#'   observations are independent, and `n_group_classes` is one. The fit raises
#'   a `latents_single_level` message. Methods for extracting results,
#'   visualization and inference remain available. The internal unit column is
#'   named `.observation`; it is not
#'   returned by `get_results(x, "data")`, and a `data` that already has a
#'   column of that name raises `latents_bad_data`. Requesting more than one
#'   group class with `id = NULL` raises `latents_bad_argument`.
#' @param n_profiles Positive integer number of individual profiles.
#' @param n_group_classes Positive integer number of latent group classes
#'   (default 2). With `id = NULL` it is one, and any other value is an error.
#' @param profile_covariates,group_covariates Names of numeric columns of
#'   `data` that predict individual-profile and group-class membership through
#'   multinomial logits, with the last class as reference. Naming either fits
#'   the one-step maximum-likelihood covariate model and returns a
#'   `multilpa_covariates` object. Profile intercepts differ by group class;
#'   profile slopes are shared across group classes unless
#'   `profile_slopes = "group_class"`. A group covariate must be constant
#'   within each group. Covariates enter in their supplied units, without
#'   centring or scaling. The covariate model rejects `start` and `fixed`; it
#'   accepts `missing = "fiml"` for the indicators, but the covariates must be
#'   complete. [three_step()] and [r3step()] are the staged alternatives.
#' @param profile_slopes `"shared"` (the default) gives each profile
#'   covariate one slope per profile logit, common to every group class.
#'   `"group_class"` estimates separate slopes in each group class (a
#'   cross-level interaction), reported with terms such as `z:group_class_1`;
#'   it multiplies the number of slope parameters by `n_group_classes`, each
#'   set estimated from the groups in its class. Naming it without
#'   `profile_covariates` raises `latents_bad_argument`.
#' @param variance_model `"varying"` (the default; profile-specific variances
#'   or covariance matrices) or `"equal"` (shared across profiles).
#' @param n_starts Positive integer number of EM starts (default 10). The
#'   first start is deterministic: Ward's hierarchical clustering of the
#'   standardized indicators (on at most 2000 evenly spaced rows), each cluster
#'   initialized with its own means and variances, and, with more than one group
#'   class, Ward's clustering of the groups' profile compositions. The
#'   remaining starts are k-means from random centres. A supplied `start`
#'   replaces the first start, and the remaining starts are random. With
#'   `start` and `max_iter = 0`, exactly one start is evaluated regardless of
#'   `n_starts`. Zero-weight units do not contribute to initialization. For
#'   `search = "emEM"`, `n_starts` is the number of final optimization runs,
#'   including the protected first start (see `search`).
#' @param max_iter Nonnegative integer maximum number of EM updates per start
#'   (default 1000). `max_iter = 0` performs no update: with `start`, the model
#'   is evaluated at exactly those values, for example to score them with
#'   `logLik()`; without `start`, the random initializations are evaluated and
#'   the one with the highest log likelihood is returned.
#' @param tol Positive convergence tolerance (default `1e-8`). EM has
#'   converged when one step (or SQUAREM cycle) changes the log likelihood by
#'   no more than `tol` times the number of observed indicator values (each
#'   non-missing cell in a row of positive weight, with a minimum count of one).
#'   The tolerance is therefore expressed per observed indicator value.
#' @param min_variance Lower bound on each variance, or on each covariance
#'   eigenvalue under full covariance; the fit is then constrained maximum
#'   likelihood. `NULL` (the default) sets it to `1e-6` times the smallest
#'   sample variance of the continuous indicators (observed values in rows of
#'   positive weight, on the scale fitted, after any centring and scaling), so
#'   that the default bound adapts to the measurement scale. A numeric value
#'   specifies an absolute bound in squared indicator units. The bound used is
#'   stored as `min_variance` on the fit.
#'   Estimates at the bound are reported and raise a warning.
#' @param seed Optional random seed: any integer `set.seed()` accepts,
#'   including negative integers. When supplied, the caller's random-number state
#'   is restored on exit.
#' @param start Optional list of `means` and `variances` (profiles by
#'   indicators), `profile_probabilities` (group classes by profiles), and
#'   `group_probabilities` (vector). Starting probabilities must be positive.
#'   For full covariance, supply `covariances` (indicators by indicators by
#'   profiles); `variances` may then be omitted, or must match their diagonals.
#'   A categorical model also takes `response_probabilities`, a list of one
#'   profiles-by-categories matrix per indicator named in `categorical`. The
#'   list is matched by position unless it is labelled: when its elements are
#'   named after the indicators, or its columns after the categories, each
#'   block is matched to the indicator and category its labels name. A label
#'   naming an indicator or category the fit does not have raises
#'   `latents_bad_start`.
#' @param categorical Character vector naming indicators to treat as
#'   categorical. Each is modelled by unrestricted, profile-specific response
#'   probabilities over its observed categories (the latent class measurement
#'   model); binary, ordinal and unordered indicators are handled by the same
#'   parameterization. Numeric, integer, logical, character and factor columns
#'   are accepted. Indicators not named here stay Gaussian, so naming a subset
#'   fits a mixed-mode model.
#' @param model The covariance structure as an mclust code, one of
#'   `"EII"`, `"VII"`, `"EEI"`, `"VEI"`, `"EVI"`, `"VVI"`, `"EEE"`, `"VEE"`,
#'   `"EVE"`, `"VVE"`, `"EEV"`, `"VEV"`, `"EVV"`, `"VVV"`. The letters give the
#'   volume, shape and orientation, each equal (`E`) or varying (`V`) across
#'   profiles, with `I` for axis-aligned (no covariances). These are the codes
#'   [enumerate_classes()] takes. Use it instead of
#'   `variance_model`/`covariance_model` or `volume`/`shape`/`orientation`, not
#'   together with them.
#' @param volume,shape,orientation The covariance structure in its three
#'   components. Each profile's covariance decomposes as
#'   `Sigma_k = lambda_k * D_k * A_k * D_k'`, with *volume*
#'   `lambda_k = |Sigma_k|^(1/d)`, *orientation* `D_k` (eigenvectors) and
#'   *shape* `A_k` (diagonal, determinant one). Constraining the three across
#'   profiles gives the fourteen models mclust names with three letters, all
#'   of which are available.
#'
#'   `volume` is `"equal"` or `"varying"`. `shape` is `"equal"`, `"varying"`
#'   or `"spherical"`; `"spherical"` makes the variance equal across indicators
#'   within a profile and leaves no orientation to constrain. `orientation` is
#'   `"axis"` (a diagonal covariance), `"equal"` (one orientation shared by
#'   every profile) or `"varying"`. Each defaults to `NULL`, which follows
#'   `variance_model` and `covariance_model`.
#'
#'   | `volume` | `shape` | `orientation` | model |
#'   |---|---|---|---|
#'   | equal | spherical | --- | EII |
#'   | varying | spherical | --- | VII |
#'   | equal/varying | equal/varying | axis | EEI, VEI, EVI, VVI |
#'   | equal/varying | equal/varying | equal | EEE, VEE, EVE, VVE |
#'   | equal/varying | equal/varying | varying | EEV, VEV, EVV, VVV |
#'
#'   The estimates are those of Celeux and Govaert (1995); EVE and VVE, which
#'   combine a shared orientation with a free shape, have no closed form and
#'   use the minorize-maximize step of Browne and McNicholas (2014). Parameter
#'   counts match mclust's for all fourteen.
#'
#'   Structures other than EEI, VVI, EEE and VVV are maximized across all
#'   profiles jointly and cannot be combined with a held `variances` block.
#'   [parameter_inference()] reports Wald standard errors for all fourteen,
#'   computed in the structure's free coordinates (log volumes,
#'   determinant-one log shapes and orientations) and transformed to the
#'   reported variances and covariances by the delta method;
#'   `method = "bootstrap"` is also available for all fourteen.
#' @param select_start Rule for selecting a fitted solution. `"likelihood"`
#'   (the default) selects among starts within `tol * max(1, n_observed)` of
#'   the highest log likelihood, preferring converged starts within this set.
#'   Here `n_observed` counts observed indicator values in positive-weight rows.
#'   `"converged"` restricts selection to converged starts when any are
#'   available; otherwise it uses all starts. Remaining ties are resolved by
#'   start order. Convergence is required for
#'   Wald inference but does not establish identification or global
#'   optimality. `summary()` and `get_results(fit, "starts")` report the
#'   individual starts. The reported replication
#'   count includes final-stage starts within `1e-6 * max(1, n_observed)` of
#'   the selected log likelihood, where `n_observed` counts observed indicator
#'   values in rows of positive sampling weight, as in the EM stopping rule.
#'   This absolute likelihood-gap tolerance is independent of measurement
#'   units and of `tol`, since a start stopped by the EM rule can still lie
#'   several steps below its maximum. It is stored as
#'   `replication_tolerance` and does not select a start.
#' @param boundary Whether a start with a variance at `min_variance` (a
#'   covariance eigenvalue, under full covariance) may be reported. `"allow"`
#'   (the default) permits it, as Mplus does, with a `latents_boundary`
#'   warning. `"exclude"` sets such starts aside and applies `select_start` to
#'   the rest; for models without covariates, the `boundary` column of the
#'   starts table identifies them. When every start is at the bound, the best
#'   is reported with a `latents_all_starts_boundary` warning. A Gaussian
#'   profile can increase the likelihood without limit by shrinking around one
#'   observation, so the size of that gain is set by `min_variance` rather than
#'   by the data, and more starts make such a solution more likely to be found.
#'   Categorical indicators have no variances and are unaffected.
#' @param prior `NULL` (the default) for maximum likelihood, or a
#'   [prior_control()] object for maximum a posteriori estimation of the
#'   Gaussian means and covariances under the conjugate prior of Fraley and
#'   Raftery (2007). Hyperparameters use the data-dependent defaults of
#'   `mclust::defaultPrior()`, and updates follow the conventions described in
#'   [prior_control()]. Available for EII, VII,
#'   EEI, VEI, EVI, VVI, EEE, EEV, VEV and VVV (mclust has no prior for VEE,
#'   EVE, VVE and EVV), with complete continuous indicators and without `fixed`
#'   or membership covariates; otherwise `latents_unsupported_prior` is raised.
#'   The mixing proportions have no prior. As in mclust, the fit's
#'   `log_likelihood` is the *unpenalized* log likelihood at the posterior
#'   mode, and `aic`, `bic` and their variants are computed from it with the
#'   usual parameter count. Convergence uses the change described under `tol`;
#'   because a MAP step can lower that likelihood, the decrease check made
#'   under maximum likelihood is not applied. `parameter_inference(method =
#'   "wald")` raises `latents_unsupported_inference` for such a fit, because
#'   the score is not zero at a posterior mode; `method = "bootstrap"` refits
#'   every resample with the same prior. The fit stores the request as `prior`
#'   and the resolved hyperparameters, in the indicators' units (centred,
#'   under `centering`), as `prior_parameters`.
#' @param noise `FALSE` (the default), or `TRUE` to add a noise component: a
#'   mixture component with constant density `1 / V` over the hypervolume `V`
#'   of the data, which absorbs outliers and scatter that no Gaussian profile
#'   explains (Banfield & Raftery, 1993; mclust's
#'   `initialization = list(noise = )` and `Vinv`). `V` is mclust's
#'   `hypvol()`: the smaller of the volumes of the data's axis-aligned and
#'   principal-component-aligned bounding boxes. The component has no
#'   measurement parameters; its mixing proportion is estimated with the
#'   others, and, as in mclust, it adds two parameters (the proportion and the
#'   hypervolume), following mclust's parameter-count convention. A random
#'   start assigns probability 0.1 to noise. A supplied `start` includes
#'   `noise_probability` alongside
#'   `profile_probabilities`, the two summing to one. Available for complete
#'   continuous indicators with one group class (`id = NULL` or
#'   `n_group_classes = 1`), under any of the fourteen covariance structures,
#'   with or without `prior`; otherwise `latents_unsupported_noise` is raised.
#'   In the fit the noise component is profile `0`, as in mclust's
#'   classification: `subject_profiles` is `0` for an observation assigned to
#'   it, `get_results(fit, "assignments")` adds a `posterior_noise` column,
#'   `get_results(fit, "posteriors")` and
#'   `get_results(fit, "profile_probabilities")` have profile-`0` rows, and the
#'   classification diagnostics count it as a class. `means`, `variances`,
#'   `covariances`, `profile_probabilities` and `subject_posteriors` describe
#'   the Gaussian profiles only, so the last two sum to one minus the noise
#'   share; the fit adds `noise_probability`, `noise_posteriors`, `n_noise` and
#'   `hypervolume`. [parameter_inference()], [three_step()], [r3step()],
#'   [bootstrap_lrt()], [starting_values()], bivariate residuals,
#'   [fit_staged()] and the posterior plots raise `latents_unsupported_noise`
#'   for such a fit.
#' @param acceleration `"squarem"` (the default) accelerates EM with SQUAREM
#'   (Varadhan & Roland, 2008): each cycle extrapolates along the last two EM
#'   steps in unit-free coordinates (means divided by reference standard
#'   deviations, variances and covariances by the corresponding products;
#'   probabilities unchanged). Reference standard deviations use observed values in rows
#'   of positive sampling weight. Each cycle ends with an ordinary EM step
#'   and reverts to the EM update when
#'   the extrapolated point would lower the likelihood, so the iteration stays
#'   monotone. Convergence is assessed over a complete cycle. `"none"` runs
#'   ordinary EM. `iterations` and `max_iter`
#'   count EM steps under either setting. Fits with `prior`, and models with
#'   membership covariates, always use plain EM. For EVE and VVE, whose shared
#'   orientation has no closed form, each EM iteration takes a fixed number of
#'   warm-started orientation steps (a generalized EM step) under either
#'   setting, which keeps every iteration monotone and the EM map a fixed
#'   function, as SQUAREM requires.
#' @param centering How to centre the continuous indicators before fitting.
#'   `"none"` (the default) fits them as supplied. `"person"` subtracts each
#'   group's own mean from its rows, so values are deviations from the unit's
#'   average, the profiles describe within-unit change rather than level, and
#'   the group classes are types of change pattern rather than types of unit
#'   (the person-mean-centred within-person design; Quintana, 2021; Voelkle,
#'   Brose, Schmiedek, & Lindenberger, 2014). `"grand"` subtracts one mean per
#'   indicator, which moves the origin without changing the within-group
#'   structure. The offsets are stored on the fit, so
#'   `get_results(x, "data")` returns the supplied columns and row-alignment
#'   checks still apply. `"person"` raises `latents_bad_data` when it leaves an
#'   indicator constant, as happens when a unit has one observation of it.
#' @param scaling How to scale the continuous indicators before fitting, at
#'   the same levels as `centering`. `"none"` (the default) fits them as
#'   supplied. `"grand"` divides each indicator by its standard deviation over
#'   all observations; `"person"` divides each group's values by that group's
#'   own standard deviation. With the matching centring these are z-scores:
#'   `centering = "grand", scaling = "grand"` standardizes each indicator over
#'   the sample, and `centering = "person", scaling = "person"` standardizes it
#'   within each unit. Standard deviations use the observed values with
#'   denominator `n - 1`. The scales are stored on the fit with the centring
#'   offsets, so `get_results(x, "data")` returns the supplied columns, and
#'   `predict()` scales new units by their own standard deviation under `"person"` and by
#'   the fitted standard deviation under `"grand"`. Means, variances and covariances are
#'   reported in standard-deviation units. The log likelihood and every
#'   information criterion are those of the scaled values, so only fits with
#'   the same centring and scaling are comparable. A standard deviation that is
#'   zero or undefined (under `"person"`, a unit observed once on an indicator
#'   or always giving the same value) raises `latents_bad_data`. Under
#'   `EEI`, `VEI`, `EVI`, `VVI`, `EEE`, `VEE`, `EVV` and `VVV`, whose family
#'   is unchanged when each indicator is rescaled, grand scaling changes only
#'   the units of the estimates and adds the density Jacobian to each start's
#'   log likelihood. SQUAREM uses unit-free coordinates, so its search path
#'   is equivariant in exact arithmetic when variance bounds are inactive or
#'   transformed with the model. Rounding near convergence can still change
#'   iteration counts, and nearly tied starts can change the selected labelling.
#'   A scalar variance or eigenvalue floor is not generally preserved by
#'   unequal indicator rescaling. Under the spherical `EII` and `VII`, which
#'   impose one variance on every indicator, and under `EEV`, `VEV`, `EVE` and
#'   `VVE`, which share eigenvalues or eigenvectors across profiles (neither
#'   survives rescaling the indicators by different factors), it changes the
#'   solution. Under `centering = "person"` with `scaling = "person"`, a group
#'   with only two observations of an indicator becomes exactly
#'   `-1/sqrt(2)` and `1/sqrt(2)` on it: such groups raise
#'   `latents_atomic_groups`, and when every group is like that the fit is
#'   rejected with `latents_bad_data`.
#' @param drop_empty `TRUE` (the default) or `FALSE`: how the reported
#'   statistics treat rows whose every indicator is missing, and groups
#'   consisting only of such rows. This setting affects reporting, not
#'   estimation. With `TRUE` they are excluded, as in Mplus, from the individual and
#'   group sample sizes of the information criteria, from relative entropy,
#'   AvePP and the entropy term of ICL, AWE and CLC, from the classification
#'   tables and diagnostic plots, and from the effective class counts.
#'   Row-aligned assignment and posterior tables retain these units but report
#'   their classifications and posterior probabilities as `NA`. A partially
#'   observed row remains classified, and an empty row in an informative group
#'   can retain its group-class assignment. With
#'   `FALSE` every row and group is counted in those statistics, except in the
#'   individual sample size of a
#'   `multilpa()` fit, which excludes such rows either way. Units of sampling
#'   weight zero are excluded from every count and entropy term under either
#'   setting. [compare_models()]
#'   and [bootstrap_lrt()] do not compare fits made with different
#'   `drop_empty`. The setting is stored on the fit and passed on by
#'   enumerations and refits; it has no effect on data without such rows.
#' @param cores Number of processor cores for the random starts (default 1).
#'   The starting values are drawn in order in the current process and only
#'   the EM iterations are distributed, so the result is identical to
#'   `cores = 1`. Under `search = "emEM"` both the screening of the pool and
#'   the final runs are distributed; the screening is most of the time of an
#'   emEM fit, so a large model or data set gains most. Workers are forked
#'   processes on macOS and Linux, and a socket cluster on Windows and on
#'   macOS with Apple's Accelerate BLAS, as in [enumerate_classes()]. Lost
#'   worker results are retried once with half as many workers (rounded up).
#'   One `latents_worker_lost` warning reports losses and recovery; starts still
#'   lost are recorded as failed. Socket workers retain completed chunks, but
#'   tasks within a failed chunk may also be retried. Models with
#'   membership covariates run their starts on one core. To fit many models,
#'   give `cores` to [enumerate_classes()], which shares its worker pool
#'   across candidates and their starts for models without membership
#'   covariates or joined categorical `local_dependence` pairs. Other fitting
#'   engines use one task per candidate.
#' @param local_dependence Pairs of indicators that stay associated within
#'   each class, written `"item1:item2"` (several as a character vector).
#'   A continuous indicator pair has a covariance estimated within each class
#'   (one shared covariance with `variance_model = "equal"`), every other
#'   covariance staying zero: the covariance Latent GOLD fits with
#'   `item1 <-> item2 | Cluster`. It is not combined with `covariance_model`,
#'   `model` or the structure arguments, covariates, `time`, `prior`, `noise`
#'   or held variances. A pair of categorical indicators is modelled jointly
#'   within each class. Each such pair becomes one categorical indicator, named
#'   `"item1:item2"`, whose categories are every combination of the two items'
#'   categories (the full grid: an unobserved combination is estimated on the
#'   `min_probability` floor and counted as a parameter, so each class has
#'   `K1 K2 - 1` of them), so every class has its own joint distribution of
#'   the pair: the class-specific association Latent GOLD fits with
#'   `item1 <-> item2 | Cluster`. The other indicators stay independent given
#'   the class. `get_results(fit, "local_dependence")` reports each class's
#'   joint probabilities beside their values under independence; `predict()`,
#'   `parameter_inference(data = )` and the other functions that take `data`
#'   join the pair in the caller's data. Under `missing = "fiml"` a row
#'   observing one item of a pair and not the other contributes the observed
#'   item, the pair's joint distribution summed over the missing one, which is
#'   the observed-data likelihood under MAR. A pair that mixes a categorical
#'   and a continuous indicator, that names an indicator twice or that shares
#'   an indicator with another pair raises `latents_bad_argument`; a
#'   categorical item of a pair whose categories contain `":"` raises
#'   `latents_bad_data`. Default: none.
#' @param direct_effects Direct effects of covariates on continuous
#'   indicators, written `"indicator ~ covariate"` (several covariates joined
#'   with `+`, several indicators as a character vector). The covariate moves
#'   the indicator's mean by the same amount in every profile,
#'   `y = mean_k + b z + e` (Bandeen-Roche et al., 1997), so a difference in
#'   the indicator that the covariate explains is not read as a difference
#'   between profiles: the direct effect Latent GOLD fits with
#'   `y <- 1 + Cluster + z`. Covariates are numeric columns of `data`,
#'   observed in every row (code a categorical covariate as indicator
#'   columns), and enter with their values as given, so `means` are the
#'   profile means at a covariate value of zero; centre a covariate first to
#'   read them at its mean. The effects are estimated with the profiles (a
#'   conditional maximization step in EM; Meng & Rubin, 1993), have standard
#'   errors in [parameter_inference()] (parameter `direct_effect`, term
#'   `"indicator:covariate"`), and are reported by
#'   `get_results(fit, "direct_effects")`. [predict()] and the functions that
#'   take `data` read the covariates from it. Direct effects need a diagonal
#'   covariance (EEI or VVI) and are not combined with `local_dependence`
#'   between continuous indicators, `profile_covariates`, `group_covariates`,
#'   `time`, `prior`, `noise`, `fixed`, `start`, `centering` or `scaling`;
#'   those, a term that is not `"indicator ~ covariate"`, and an indicator
#'   that is not continuous raise `latents_bad_argument`, and a covariate
#'   that is absent, not numeric, missing or constant raises
#'   `latents_bad_data`. Default: none.
#' @param time Optional name of a column giving each observation's position
#'   within its group, such as a wave, occasion or course number. It is not
#'   used in estimation; it is stored for `get_results(x, "sequences")`,
#'   `get_results(x, "sequence_summary")` and `plot(what = "sequences")`.
#'   Values must be complete and unique within each group, and are ordered as
#'   [sort()] orders them: numbers numerically, a factor by its levels, and
#'   character strings alphabetically. Numbers stored as text whose
#'   alphabetical order is not their numeric order (`"2"`, `"10"`), and
#'   numbered labels sharing their text around the number whose alphabetical
#'   order is not the numbers' (`"T2"`, `"T10"`; `"week 2"`, `"week 10"`),
#'   raise `latents_bad_time` (also `latents_bad_data`); give the order with
#'   `factor(levels = )` or zero-pad the numbers. Zero-padded numbers and ISO
#'   dates are accepted.
#' @param fixed Character vector naming measurement blocks to hold at the
#'   values in `start` instead of estimating them: any of `"means"`,
#'   `"variances"` and `"response_probabilities"`, or `"measurement"` for every
#'   block the model has. Under `covariance_model = "full"`, `"variances"`
#'   holds the residual covariance matrices. A fixed block retains its supplied values in
#'   every start and is not counted in `n_parameters`, so the result is a
#'   different model, not a different starting point. `start` must contain the
#'   named blocks; `starting_values(fit, what = "measurement")` produces them,
#'   and [fit_staged()] performs the two-stage workflow in one call. The mixing
#'   parameters cannot be held.
#' @param min_probability Positive lower bound on every categorical response
#'   probability (default `1e-10`), constraining the maximum-likelihood problem
#'   as `min_variance` does for Gaussian indicators.
#' @param missing `"error"` (the default) rejects missing indicators; `"fiml"`
#'   maximizes the observed-data likelihood under an ignorable missingness
#'   mechanism (MAR), integrating missing indicators out rather than filling
#'   them in. This also holds with membership covariates; the covariates
#'   themselves must be complete (`latents_bad_data`).
#' @param covariance_model `"diagonal"` (the default) assumes conditional
#'   independence; `"full"` estimates within-profile residual covariances.
#' @param family The model family. `"profiles"` (the default) is the model
#'   described above: individual profiles whose prevalences differ across
#'   group classes. Three group-class-only families have no individual
#'   profiles (Houle, Morin & Harvey, 2026): each group has a Gaussian
#'   intercept per indicator, drawn from its group class's distribution, and
#'   its members vary around it.
#'   * `"additive"`: group classes differ in their means (and, optionally, in
#'     between-group variances); one within-group variance per indicator is
#'     shared by every class.
#'   * `"dispersion"`: group classes differ in their within-group variances;
#'     means and between-group variances are shared.
#'   * `"additive_dispersion"`: group classes differ in means and in
#'     within-group variances.
#'
#'   These families take only `data`, `vars`, `id`, `n_group_classes`,
#'   `between_variance`, `n_starts`, `max_iter`, `tol`, `min_variance`,
#'   `seed`, `weights`, `search` and `n_screen`; any other argument raises
#'   `latents_bad_argument`. Groups of sampling weight zero are not part of the
#'   sample: the identification checks, starting values, counts, entropy and
#'   classification tables use only the groups of positive weight, though
#'   every group still receives a posterior.
#'   Indicators must be complete and continuous, and at least one group must
#'   have more than one row (`latents_unidentified` otherwise). Between
#'   variances estimated at zero are reached by a boundary maximization with a
#'   Karush-Kuhn-Tucker check and flagged (`latents_boundary`); interior fits
#'   are finished by Newton steps. The result has class `multilpa_additive`;
#'   `bic` uses the number of groups and `bic_individual` the number of rows.
#'   It is read with [get_results()]; `summary()`, `plot()`, `coef()`,
#'   `vcov()`, `confint()` and [parameter_inference()] work on it,
#'   [enumerate_classes()] compares families through its `family` argument,
#'   and [bootstrap_lrt()] tests one number of group classes against another.
#'   These families are experimental: `parameter_inference(method =
#'   "bootstrap")` raises `latents_unsupported_inference`, and missing
#'   indicators (`latents_bad_data`), `missing`, membership covariates and
#'   ordinal or count indicators (`latents_bad_argument`) are not available.
#'   See the [group-class and cross-level models
#'   article](https://pak.dynasite.org/latents/articles/additive.html).
#'
#'   Two cross-level families estimate individual profiles and group classes
#'   from the same indicators, following the manifest-aggregation
#'   specification of Houle et al. (2026, supplement): each group's means of
#'   the indicators are between-level indicators of its group class, with
#'   class-specific means and variances.
#'   * `"restricted_cross_level"`: profile prevalences are the same in every
#'     group class, so profiles and group classes are estimated separately;
#'     the profile composition of each class is reported descriptively.
#'   * `"full_cross_level"`: profile prevalences differ across group classes
#'     (as in `"profiles"`), and the group means also inform the classes.
#'
#'   These take `data`, `vars`, `id`, `n_profiles`, `n_group_classes`,
#'   `variance_model` (individual profiles), `between_variance` (group
#'   classes), `n_starts`, `max_iter`, `tol`, `min_variance`, `seed` and
#'   `weights`. The group means are computed from the same ratings, so the
#'   likelihood is the specification's working likelihood: comparable between
#'   cross-level fits of the same data, not with the other families. Standard
#'   errors are not available (`latents_unsupported_inference`). The result
#'   has class `multilpa_cross_level`, read with [get_results()]. Its
#'   `get_results(fit, "profiles")` table ends with a `prevalence` column:
#'   each individual profile's prevalence (estimated and shared for the
#'   restricted family, the class-weighted marginal for the full family) and,
#'   on the group-class rows, the class weight.
#' @param between_variance For the additive and additive-dispersion families:
#'   `"varying"` (the default) estimates between-group variances per group
#'   class, `"equal"` one set shared by all classes. The dispersion family
#'   holds them equal and rejects `"varying"`. For the cross-level families,
#'   the variances of the group means within each group class.
#' @param weights `NULL`, or the name of a numeric column of `data` holding a
#'   sampling weight for each independent unit: each `id` group of a two-level
#'   fit, each row of a single-level one. A two-level weight must be constant
#'   within its group (`latents_bad_weights` otherwise); within-unit weights
#'   are not supported. The fit maximizes the pseudo log likelihood
#'   `sum_j w_j log L_j` (Skinner, 1989) with the weights scaled to sum to the
#'   number of units of positive weight, as Mplus scales to the sample size,
#'   so integer weights give the fit to the data with each unit repeated that
#'   many times and the criteria stay on the sample's scale. A unit of weight
#'   zero is not part of the sample: it changes no estimate, likelihood,
#'   sample size, entropy or criterion, sets neither the grand offset nor the
#'   grand scale of `centering`/`scaling = "grand"`, and is not classified.
#'   Nor does it define a category: a categorical, ordinal or paired
#'   (`local_dependence`) value seen only on units of weight zero is not a
#'   category of the fit and has no parameter, and an indicator constant among
#'   the positive-weight units is rejected regardless of zero-weight values. Under
#'   `missing = "error"`, a missing categorical, ordinal or count value on a
#'   unit of weight zero is not rejected; a missing continuous value is.
#'   Standard errors are sandwich estimates:
#'   `parameter_inference()` defaults to `vcov_type = "robust"` and rejects
#'   `"observed"` and `"opg"` (`latents_unsupported_weights`), and
#'   `method = "bootstrap"` resamples the units with their weights.
#'   Classification tables report each unit's own posterior; effective counts
#'   and proportions are weighted. [bootstrap_lrt()], [three_step()],
#'   [r3step()], `prior` and `noise` reject a weighted fit. Weights are supported with
#'   every `family`, membership covariates, categorical indicators and
#'   `missing = "fiml"`.
#' @param ordinal Names of indicators in `vars` to model as ordinal: an
#'   (ordered) factor, or integers, whose categories are ordered by factor
#'   level or by value. Each follows an adjacent-category logit with category
#'   intercepts shared by every profile and one location per profile,
#'   `log P(k | c) / P(k - 1 | c) = a_k - a_(k-1) + eta_c`, with the last
#'   profile's location fixed at zero: `(K - 1) + (C - 1)` parameters per
#'   indicator, against `C (K - 1)` for a `categorical` one. This is Latent
#'   GOLD's default ordinal model, against which it is checked. The category
#'   probabilities and locations are returned by
#'   `get_results(fit, "ordinal")`. A factor level never observed in the data
#'   raises the `latents_unused_category` warning; `droplevels()` removes it.
#' @param count Names of indicators in `vars` holding non-negative whole
#'   numbers, modelled as Poisson with one mean per profile (or negative
#'   binomial, see `count_model`); the means are returned by
#'   `get_results(fit, "count_means")`. Ordinal and count indicators support
#'   `missing = "fiml"`, `weights`, two-level fits, membership covariates,
#'   Wald and bootstrap inference, [bootstrap_lrt()], [predict()] and [lta()].
#'   `start`, `fixed`, `prior` and `noise` do not accept them
#'   (`latents_unsupported_indicator`), and the group-class and cross-level
#'   families do not support them
#'   (`latents_bad_argument`).
#' @param count_model `"poisson"` (the default) or `"negative_binomial"`: NB2,
#'   with mean `mu` and variance `mu + alpha mu^2`, for counts with greater variance
#'   than a Poisson within a profile (Latent GOLD's `poisson overdispersed`).
#'   A dispersion estimated at zero is the Poisson limit, a boundary fit
#'   (`latents_boundary` warning; Wald inference is rejected there).
#' @param count_dispersion For `count_model = "negative_binomial"`: one
#'   dispersion per profile (`"varying"`, the default) or one shared by every
#'   profile (`"equal"`).
#' @param search How the starts are found: `NULL` (the default), `"kmeans"`
#'   or `"emEM"`. `NULL` uses `"emEM"` when every profile shares one
#'   covariance matrix (EII, EEI, EEE with continuous indicators), which removes
#'   the isolated-profile covariance-collapse mechanism. Shared covariances can
#'   still hit their floor on degenerate or over-extracted data; inspect the
#'   boundary diagnostics. The default is `"kmeans"` otherwise, including
#'   categorical-only models and the group-class families. The selected search
#'   strategy is recorded on the fit. `"kmeans"` runs
#'   `n_starts` starts to convergence: Ward's clustering first, then k-means
#'   from random centres. `"emEM"` screens a pool of `n_screen` starts by a
#'   short EM run each and runs `n_starts` selected starts to convergence
#'   (Biernacki, Celeux & Govaert, 2003). The first start (Ward's clustering,
#'   or supplied starting values) is retained for final optimization if its
#'   screening succeeds; remaining starts are selected by screening likelihood.
#'   Thus `n_starts = 1` finalizes Ward's start whenever it screens
#'   successfully, even if another screened start has a higher likelihood.
#'   With a supplied start and `n_starts = 1`, screening is skipped entirely.
#'   The pool begins with the kmeans search's own
#'   starts, drawn in the same order, followed by random partitions of the
#'   rows and random rows as class centres, and `get_results(fit, "starts")`
#'   gives every start of the pool with its `kind`, its `stage` and its
#'   likelihood after screening. Either search reports the start with the
#'   highest likelihood among those run to convergence (subject to
#'   `select_start` and `boundary`). With profile-specific variances that
#'   maximum can correspond to a nearly degenerate component concentrated on
#'   a few observations (McLachlan & Peel, 2000, section 3.10). More extensive
#'   searches can identify such solutions more often. A variance at
#'   `min_variance` produces a `latents_boundary` warning; use
#'   `boundary = "exclude"` to exclude these starts when an interior solution
#'   is available. Low effective Gaussian profile membership produces the
#'   heuristic warning `latents_thin_profiles`. `get_results(fit, "gaussian_diagnostics")` reports continuous
#'   covariance and support measures without an automatic spike cutoff or
#'   changing the selected start. Neither a warning's absence nor a repeated
#'   likelihood establishes a regular or globally optimal solution.
#' @param n_screen For `search = "emEM"`, the number of starts screened, at
#'   least `n_starts`. `NULL` (the default) screens five times `n_starts`.
#'   The protected first start is included in the `n_starts` final optimization
#'   runs. Increasing `n_screen` does not increase that number. A single supplied
#'   start with `n_starts = 1` disables screening.
#'   A value with `search = "kmeans"`, or one below `n_starts`, raises
#'   `latents_bad_argument`.
#' @return A `multilpa` object containing `means`, `variances`, optional
#'   `covariances` (indicators by indicators by profiles),
#'   `profile_probabilities`, `group_probabilities`, posterior matrices,
#'   classifications, log likelihood, information criteria, restart
#'   diagnostics and convergence history. Individual rows retain their input
#'   order; groups keep first-occurrence order. `bic` and `bic_groups` use the
#'   observed group count; `bic_individual` uses `n_informative`, the number of
#'   rows with at least one observed indicator. These are alternative
#'   conventions, not interchangeable criteria. `boundary` identifies variance
#'   bounds; `small_classes` flags a profile or group class whose effective
#'   membership (`effective_profile_counts`, `effective_group_counts`, weighted
#'   for a fit with `weights`) is below one, in which case the fit warns with
#'   `latents_small_classes`. The fit carries no standard errors:
#'   [parameter_inference()] computes them, and [bootstrap_lrt()] tests nested
#'   models. Global optimality is not guaranteed for any `n_starts`.
#'   `as.data.frame()` returns the tidy tables, selected by `what`.
#' @details Infinite and constant observed indicators are rejected; each
#'   indicator must have at least two distinct observed values. No rows are
#'   silently dropped. Under FIML, fully missing individuals contribute no
#'   direct measurement likelihood but receive posterior probabilities from
#'   their group's information. The individual-level BIC excludes these rows
#'   through `n_informative`; the treatment of fully missing groups follows
#'   `drop_empty`. Missing patterns can prevent
#'   parameter identification; no general identification guarantee is made.
#'   Only initialization fills missing values with indicator means; EM uses
#'   conditional Gaussian sufficient statistics and observed marginal
#'   densities. More than one group class requires more than one profile and at
#'   least one group with more than one individual; these checks are necessary
#'   but do not establish identification. A model of complete categorical and
#'   ordinal indicators alone, with one group class, warns
#'   (`latents_unidentified`) when it has more free parameters than its
#'   contingency table has free cells (negative degrees of freedom), which no
#'   data can identify (Goodman, 1974). Start selection follows `select_start`;
#'   the default can return a solution that has not converged. Inspect the
#'   `converged` flag and the `starts` table before interpreting the model.
#'   Profile and group-class labels are arbitrary.
#' @references Vermunt, J. K. (2003). Multilevel latent class models.
#'   Sociological Methodology, 33, 213--239.
#'   doi:10.1111/j.0081-1750.2003.t01-1-00131.x.
#'
#'   Goodman, L. A. (1974). Exploratory latent structure analysis using both
#'   identifiable and unidentifiable models. Biometrika, 61(2), 215--231.
#'   doi:10.1093/biomet/61.2.215.
#'
#'   Banfield, J. D., & Raftery, A. E. (1993). Model-based Gaussian and
#'   non-Gaussian clustering. Biometrics, 49, 803--821. doi:10.2307/2532201.
#'
#'   Bandeen-Roche, K., Miglioretti, D. L., Zeger, S. L., & Rathouz, P. J.
#'   (1997). Latent variable regression for multiple discrete outcomes.
#'   Journal of the American Statistical Association, 92(440), 1375--1386.
#'
#'   Meng, X.-L., & Rubin, D. B. (1993). Maximum likelihood estimation via the
#'   ECM algorithm: A general framework. Biometrika, 80(2), 267--278.
#'   doi:10.1093/biomet/80.2.267.
#'
#'   Houle, S. A., Morin, A. J. S., & Harvey, J.-F. (2026). Multilevel latent
#'   profile analyses: A comprehensive guide. Organizational Research Methods.
#'   doi:10.1177/10944281261469432.
#'
#'   Skinner, C. J. (1989). Domain means, regression and multivariate
#'   analysis. In C. J. Skinner, D. Holt, & T. M. F. Smith (Eds.), Analysis
#'   of complex surveys (pp. 59--87). Wiley.
#'
#'   Asparouhov, T. (2005). Sampling weights in latent variable modeling.
#'   Structural Equation Modeling, 12, 411--434.
#'   doi:10.1207/s15328007sem1203_4.
#'
#'   Fraley, C., & Raftery, A. E. (2007). Bayesian regularization for normal
#'   mixture estimation and model-based clustering. Journal of
#'   Classification, 24, 155--181. doi:10.1007/s00357-007-0004-5.
#' @examples
#' # Course sessions nested in students: engagement profiles of sessions, and
#' # classes of students that differ in how often they are engaged.
#' activity <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
#' fit <- multilpa(course_engagement, activity, "student",
#'                 n_profiles = 2, n_group_classes = 2, n_starts = 4,
#'                 seed = 42)
#' summary(fit)
#' as.data.frame(fit)
#' get_results(fit, what = "profile_probabilities")
#' @seealso [multilca()] for a model in which every indicator is categorical.
#' @export
#' @importFrom stats setNames
multilpa <- function(data, vars, id, n_profiles,
                       n_group_classes = 2L,
                       profile_covariates = character(),
                       group_covariates = character(),
                       variance_model = c("varying", "equal"),
                       n_starts = 10L, max_iter = 1000L, tol = 1e-8,
                       min_variance = NULL, seed = NULL, start = NULL,
                       missing = c("error", "fiml"),
                       profile_slopes = c("shared", "group_class"),
                       covariance_model = c("diagonal", "full"),
                       categorical = character(), min_probability = 1e-10,
                       time = NULL, fixed = character(),
                       centering = c("none", "person", "grand"),
                       volume = NULL, shape = NULL, orientation = NULL,
                       model = NULL,
                       select_start = c("likelihood", "converged"),
                       prior = NULL, noise = FALSE,
                       acceleration = c("squarem", "none"),
                       family = c("profiles", "additive", "dispersion",
                                  "additive_dispersion", "restricted_cross_level",
                                  "full_cross_level"),
                       between_variance = c("varying", "equal"),
                       weights = NULL, ordinal = character(),
                       count = character(),
                       count_model = c("poisson", "negative_binomial"),
                       count_dispersion = c("varying", "equal"),
                       scaling = c("none", "person", "grand"),
                       drop_empty = TRUE, cores = 1L,
                       local_dependence = character(),
                       direct_effects = character(),
                       search = NULL, n_screen = NULL,
                       boundary = c("allow", "exclude")) {
  # Forward only supplied arguments so setup retains every missing() check.
  fitting_call <- match.call()
  # The setup runs as a call of symbols: its head is the function's name and
  # any argument that arrived as a value (through do.call()) is bound to a
  # name first, so a condition raised while validating the input never
  # carries the data in its call (review 2026-10-27).
  setup_env <- new.env(parent = parent.frame())
  setup_env$.multilpa_setup <- .multilpa_setup
  setup_call <- .latents_bind_values(fitting_call, setup_env, parent.frame())
  setup_call[[1L]] <- quote(.multilpa_setup)
  setup_call$.call <- substitute(quote(value), list(value = fitting_call))
  plan <- eval(setup_call, setup_env)
  if (!inherits(plan, "latents_multilpa_plan")) return(plan)
  .multilpa_run_plan(plan, plan$cores)
}

#' Prepare a profile-model fit without running its starts
#' @inheritParams multilpa
#' @param .call Public fitting call to retain on the result.
#' @return A prepared plan, or the existing fit for another fitting engine.
#' @noRd
.multilpa_setup <- function(data, vars, id, n_profiles,
                       n_group_classes = 2L,
                       profile_covariates = character(),
                       group_covariates = character(),
                       variance_model = c("varying", "equal"),
                       n_starts = 10L, max_iter = 1000L, tol = 1e-8,
                       min_variance = NULL, seed = NULL, start = NULL,
                       missing = c("error", "fiml"),
                       profile_slopes = c("shared", "group_class"),
                       covariance_model = c("diagonal", "full"),
                       categorical = character(), min_probability = 1e-10,
                       time = NULL, fixed = character(),
                       centering = c("none", "person", "grand"),
                       volume = NULL, shape = NULL, orientation = NULL,
                       model = NULL,
                       select_start = c("likelihood", "converged"),
                       boundary = c("allow", "exclude"),
                       prior = NULL, noise = FALSE,
                       acceleration = c("squarem", "none"),
                       family = c("profiles", "additive", "dispersion",
                                  "additive_dispersion", "restricted_cross_level",
                                  "full_cross_level"),
                       between_variance = c("varying", "equal"),
                       weights = NULL, ordinal = character(),
                       count = character(),
                       count_model = c("poisson", "negative_binomial"),
                       count_dispersion = c("varying", "equal"),
                       scaling = c("none", "person", "grand"),
                       drop_empty = TRUE, cores = 1L,
                       local_dependence = character(),
                       direct_effects = character(),
                       search = NULL, n_screen = NULL, .call = NULL) {
  .call <- .call %||% match.call()
  # Read before any argument is reassigned: missing() is FALSE afterwards.
  covariance_model_given <- !missing(covariance_model)
  # Locally dependent categorical pairs are joined into one indicator each,
  # and the model is fitted on the joint indicators with the caller's other
  # arguments; continuous pairs go on to a covariance pattern below.
  if (length(local_dependence) > 0L &&
      any(.latents_parse_pairs(local_dependence)$item_1 %in% categorical)) {
    return(.latents_fit_local_dependence(.call, environment(),
                                         local_dependence))
  }
  family <- .latents_match_formal(family)
  cores <- .latents_check_cores(cores)
  count_model <- .latents_match_formal(count_model)
  count_dispersion <- .latents_match_formal(count_dispersion)
  if (family %in% c("restricted_cross_level", "full_cross_level")) {
    supplied <- setdiff(names(.call)[-1L], .cross_level_arguments())
    if (length(supplied) > 0L || missing(id) || is.null(id) || missing(n_profiles)) {
      stop(errorCondition(sprintf(paste(
        "`family = \"%s\"` needs `id` and `n_profiles`; it takes %s.%s"), family,
        paste(sprintf("`%s`", setdiff(.cross_level_arguments(), "family")),
              collapse = ", "),
        if (length(supplied) > 0L) sprintf(" Not applicable: %s.",
          paste(sprintf("`%s`", supplied), collapse = ", ")) else ""),
        class = "latents_bad_argument", call = NULL))
    }
    return(.cross_level_fit(data, vars, id, n_profiles, n_group_classes, family,
                            .latents_match_formal(variance_model),
                            if (missing(between_variance)) "varying" else
                              .latents_match_formal(between_variance),
                            n_starts, max_iter, tol, min_variance, seed,
                            .call, weights = weights,
                            search = .latents_check_search(search) %||% "kmeans",
                            n_screen = n_screen))
  }
  if (!identical(family, "profiles")) {
    # Group-class families have no individual profiles, so most arguments
    # have no meaning for them; refuse those rather than ignore them.
    supplied <- setdiff(names(.call)[-1L], .additive_arguments())
    if (length(supplied) > 0L || missing(id) || is.null(id)) {
      stop(errorCondition(sprintf(paste(
        "`family = \"%s\"` has group classes only and needs `id`; it takes",
        "%s.%s"), family,
        paste(sprintf("`%s`", setdiff(.additive_arguments(), "family")),
              collapse = ", "),
        if (length(supplied) > 0L) sprintf(" Not applicable: %s.",
          paste(sprintf("`%s`", supplied), collapse = ", ")) else ""),
        class = "latents_bad_argument", call = NULL))
    }
    return(.additive_fit(data, vars, id, n_group_classes,
                         if (missing(between_variance)) NULL else
                           .latents_match_formal(between_variance),
                         n_starts, max_iter, tol, min_variance, seed,
                         .call, family = family, weights = weights,
                         search = .latents_check_search(search) %||% "kmeans",
                         n_screen = n_screen))
  }
  if (!missing(between_variance)) {
    stop(errorCondition(paste(
      "`between_variance` belongs to the group-class families",
      "(`family = \"additive\"`); the profile model does not have it."),
      class = "latents_bad_argument", call = NULL))
  }
  select_start <- .latents_match_formal(select_start)
  boundary <- .latents_match_formal(boundary)
  # Read once here: `boundary` is reused below for the selected fit's flag.
  exclude_boundary <- identical(boundary, "exclude")
  acceleration <- .latents_match_formal(acceleration)
  search <- .latents_check_search(search)
  ## Before `stopifnot()`, which reads `id` and would otherwise force the
  ## missing argument into R's own bare "argument \"id\" is missing" error.
  if (missing(id)) .multilpa_missing_id("multilpa", "`lpa()`")
  .latents_check_data_frame(data, substitute(data))
  .latents_require(
    "`vars` must be a character vector of column names" =
      is.character(vars),
    "`id` must be a single column name, or NULL for a single-level fit" =
      is.null(id) || (is.character(id) && length(id) == 1L),
    "`profile_covariates` must be a character vector of column names" =
      is.character(profile_covariates) && !anyNA(profile_covariates),
    "`group_covariates` must be a character vector of column names" =
      is.character(group_covariates) && !anyNA(group_covariates),
    "`categorical` must be a character vector of indicator names" =
      is.character(categorical) && !anyNA(categorical),
    "`min_probability` must be a single number in (0, 1)" =
      is.numeric(min_probability) && length(min_probability) == 1L &&
      is.finite(min_probability) && min_probability > 0 &&
      min_probability < 1)
  # Every later reader of `start` indexes it as a list, so anything else would
  # surface as an unclassed `$` error from deep inside the fit.
  if (!is.null(start) && !is.list(start)) {
    stop(errorCondition(
      "`start` must be NULL or a list of parameter blocks, such as starting_values(fit).",
      class = "latents_bad_start", call = NULL))
  }
  call <- .call
  single_level <- is.null(id)
  if (single_level) {
    resolved <- .multilpa_single_level(
      data, if (missing(n_group_classes)) NULL else n_group_classes,
      latent_class = length(vars) > 0L && all(vars %in% categorical))
    data <- resolved$data
    id <- resolved$id
    n_group_classes <- resolved$n_group_classes
  }
  # `model` names the whole covariance structure in mclust's three letters,
  # the same codes enumerate_classes() takes; it replaces the switches and the
  # pieces rather than combining with them.
  if (!is.null(model)) {
    if (!is.character(model) || length(model) != 1L ||
        !model %in% .multilpa_structures()) {
      stop(errorCondition(sprintf(
        "`model` must be one of the covariance structures %s.",
        paste(sprintf("\"%s\"", .multilpa_structures()), collapse = ", ")),
        class = "latents_bad_argument", call = NULL))
    }
    if (!missing(variance_model) || !missing(covariance_model) ||
        !is.null(volume) || !is.null(shape) || !is.null(orientation)) {
      stop(errorCondition(paste(
        "Give the covariance structure either as `model` or through",
        "`variance_model`/`covariance_model` or `volume`/`shape`/`orientation`,",
        "not both."), class = "latents_bad_argument", call = NULL))
    }
    pieces <- .multilpa_structure_arguments(model)
    volume <- pieces$volume
    shape <- pieces$shape
    orientation <- pieces$orientation
  }
  variance_model <- .latents_match_formal(variance_model)
  covariance_model <- .latents_match_formal(covariance_model)
  missing <- .latents_match_formal(missing)
  centering <- .latents_match_formal(centering)
  scaling <- .latents_match_formal(scaling)
  .multilpa_check_drop_empty(drop_empty)
  profile_slopes <- .latents_match_formal(profile_slopes)
  if (identical(profile_slopes, "group_class") && length(profile_covariates) == 0L) {
    stop(errorCondition(paste(
      "`profile_slopes = \"group_class\"` lets profile-covariate slopes differ",
      "across group classes, so it needs `profile_covariates`."),
      class = "latents_bad_argument", call = NULL))
  }
  # Continuous local dependence: a covariance for each named pair within every
  # class, every other covariance zero, fitted as a patterned full covariance.
  covariance_pattern <- NULL
  if (length(local_dependence) > 0L) {
    covariance_pattern <- .latents_covariance_pattern(
      .latents_parse_pairs(local_dependence), vars, categorical, ordinal, count,
      conflicts = c(
        "covariance_model" = covariance_model_given, "model" = !is.null(model),
        "volume, shape or orientation" =
          !is.null(volume) || !is.null(shape) || !is.null(orientation),
        "profile_covariates or group_covariates" =
          length(profile_covariates) > 0L || length(group_covariates) > 0L,
        "time" = !is.null(time), "prior" = !is.null(prior),
        "fixed variances" = any(c("variances", "measurement") %in% fixed),
        "noise" = isTRUE(noise)))
    covariance_model <- "full"
  }
  structure <- .multilpa_resolve_structure(variance_model, covariance_model,
                                           volume, shape, orientation)
  # An ellipsoidal structure has an orientation, which a diagonal parameter
  # block cannot carry, so the fit keeps a full covariance array whatever
  # `covariance_model` said.
  if (.multilpa_is_ellipsoidal(structure)) covariance_model <- "full"
  # The four structures the package has always fitted are maximized by the
  # `variance_model` branch of the M-step rather than by a constrained solver,
  # so reaching one of them through `volume`/`shape`/`orientation` has to set
  # it. Without this, asking for EEI by its pieces reported EEI and counted its
  # parameters while fitting free per-profile variances.
  if (structure %in% c("EEI", "EEE")) variance_model <- "equal"
  if (structure %in% c("VVI", "VVV")) variance_model <- "varying"
  search <- search %||%
    .latents_default_search(structure, setdiff(vars, c(categorical, ordinal, count)))
  # Direct effects are checked here, before the covariate and transition
  # engines take over: those engines do not fit them.
  direct_terms <- NULL
  if (length(direct_effects) > 0L) {
    direct_terms <- .latents_parse_direct_effects(direct_effects)
    .latents_check_direct_conflicts(c(
      "a covariance structure other than EEI or VVI" = !structure %in% c("EEI", "VVI"),
      "`local_dependence` between continuous indicators" = length(local_dependence) > 0L,
      "`profile_covariates` or `group_covariates`" =
        length(profile_covariates) > 0L || length(group_covariates) > 0L,
      "`time`" = !is.null(time), "`prior`" = !is.null(prior),
      "`noise`" = isTRUE(noise),
      "`fixed` or `start`" = length(fixed) > 0L || !is.null(start),
      "`centering` or `scaling`" = !identical(centering[[1L]], "none") ||
        !identical(scaling[[1L]], "none")))
  }
  # Counted across every M-step of every start, and reported once at the end.
  .multilpa_reset_structure_log()
  report_structure <- TRUE
  on.exit(if (report_structure) .multilpa_report_structure_log(), add = TRUE)
  if (!structure %in% .multilpa_separable_structures() &&
      any(c("variances", "measurement") %in% fixed)) {
    stop(errorCondition(paste(
      "A constrained covariance structure is maximized across every profile at",
      "once, so holding one profile's variances would not leave the others at",
      "their maximum. Use `fixed = \"means\"`, or the unconstrained structure."),
      class = "latents_bad_argument", call = NULL))
  }
  .multilpa_check_prior(prior, structure, categorical, fixed,
                        length(profile_covariates) + length(group_covariates))
  .multilpa_check_noise(noise, n_group_classes, categorical, fixed,
                        length(profile_covariates) + length(group_covariates))
  if (!is.null(weights)) {
    if (!is.null(prior)) .latents_refuse_weights("`prior`")
    if (isTRUE(noise)) .latents_refuse_weights("`noise = TRUE`")
  }
  .latents_check_extra_arguments(vars, categorical, ordinal, count, list(
    "`prior`" = !is.null(prior), "`noise = TRUE`" = isTRUE(noise),
    "`fixed`" = length(fixed) > 0L, "`start`" = !is.null(start)))
  if (length(profile_covariates) > 0L || length(group_covariates) > 0L) {
    # The covariate model is fitted through `variance_model` and
    # `covariance_model` alone, which reach these four structures only; any
    # other used to be fitted silently as the nearest of them.
    if (!structure %in% c("EEI", "VVI", "EEE", "VVV")) {
      stop(errorCondition(sprintf(paste(
        "The covariance structure \"%s\" is not available with",
        "`profile_covariates` or `group_covariates`, which fit EEI, VVI, EEE",
        "and VVV only. Fit the covariate-free model with this structure and use",
        "three_step() or r3step(), or choose one of the four."), structure),
        class = "latents_unsupported_structure", call = NULL))
    }
    result <- .multilpa_covariate_model(
      data = data, vars = vars, id = id, n_profiles = n_profiles,
      n_group_classes = n_group_classes,
      profile_covariates = profile_covariates,
      group_covariates = group_covariates, variance_model = variance_model,
      n_starts = n_starts, max_iter = max_iter, tol = tol,
      min_variance = min_variance, seed = seed, start = start,
      missing = missing, profile_slopes = profile_slopes,
      covariance_model = covariance_model,
      categorical = categorical, min_probability = min_probability,
      time = time, fixed = fixed, select_start = select_start, boundary = boundary,
      call = call,
      weights = weights, ordinal = ordinal, count = count,
      count_model = count_model, count_dispersion = count_dispersion,
      centering = centering, scaling = scaling, drop_empty = drop_empty,
      search = search, n_screen = n_screen)
    result$single_level <- single_level
    return(result)
  }
  # The data contract comes first: `time` is checked against the `id` column,
  # so an `id` that does not name a column of `data` must be reported as
  # that, not as an opaque failure inside the time check.
  .multilpa_check_arguments(data, vars, id, n_profiles, n_group_classes,
                          n_starts, max_iter, tol, min_variance, min_probability,
                          seed, categorical)
  time_values <- .multilpa_time_values(data, time, id, vars)
  measurement <- .multilpa_prepare_indicators(data, vars, categorical,
                                            missing, min_probability,
                                            other = c(ordinal, count),
                                            weights = weights)
  # Ordinal categories and the count checks are judged on the rows of
  # positive weight, as the other indicators are (review 2026-10-27, CORE-2).
  extra <- .latents_set_count_model(
    .latents_prepare_extra(data, ordinal, count, missing,
                           contributing = .multilpa_contributing_rows(data, weights)),
    count_model, count_dispersion)
  if (!is.null(direct_terms)) {
    extra <- .latents_add_direct(extra, .latents_prepare_direct(
      data, direct_terms, vars, measurement$continuous))
  }
  continuous <- measurement$continuous
  indicator_frame <- measurement$frame
  encoded <- measurement$encoded
  codes <- measurement$codes
  n_categories <- measurement$n_categories
  x <- measurement$x
  groups <- .multilpa_prepare_groups(data[[id]])
  group_values <- groups$values
  group_index <- groups$index
  group_ids <- groups$ids
  n_groups <- groups$n
  sampling_weights <- .latents_sampling_weights(data, weights, group_index, n_groups)
  # The groups that enter the likelihood. A group of weight zero is neither
  # counted nor classified, so no data check may lean on it (review
  # 2026-10-26, ML-2).
  carrying <- if (is.null(sampling_weights)) rep(TRUE, n_groups) else
    sampling_weights > 0
  # Centring happens here, between validating the indicators and fitting them,
  # so the measurement model sees deviations and every verb downstream can put
  # a supplied frame on the same scale from the offsets kept on the fit.
  centred <- .multilpa_center_indicators(x, group_index, n_groups, centering,
                                         scaling, carrying = carrying)
  if ((!identical(centering, "none") || !identical(scaling, "none")) &&
      ncol(x) > 0L) {
    .multilpa_check_centered(centred$x, centering, raw = x,
                             group_index = group_index, carrying = carrying)
    .multilpa_check_atomic_groups(x, group_index, centering, scaling,
                                  carrying = carrying)
    x <- centred$x
    indicator_frame <- as.data.frame(x)
  }
  # The default variance bound is relative to the indicators as they are
  # fitted, once, here (review 2026-10-26 round 2, K2-2).
  min_variance <- .latents_resolve_min_variance(
    min_variance, x,
    if (is.null(sampling_weights)) NULL else sampling_weights[group_index])
  group_sizes <- groups$sizes
  # Identification is judged on the carrying groups: a group of weight zero
  # can supply neither a distinct row nor a group for a group class to hold.
  carried_rows <- carrying[group_index]
  distinct_rows <- nrow(unique(cbind(x, codes, .latents_extra_matrix(extra))[
    carried_rows, , drop = FALSE]))
  if (n_profiles > distinct_rows) {
    stop(errorCondition(
      "n_profiles cannot exceed the number of distinct observed indicator rows.",
      class = "latents_unidentified", call = NULL))
  }
  if (n_group_classes > sum(carrying)) {
    stop(errorCondition(
      "n_group_classes cannot exceed the number of groups (of positive weight).",
      class = "latents_unidentified", call = NULL))
  }
  if (n_group_classes > 1L &&
      (n_profiles == 1L || all(group_sizes[carrying] == 1L))) {
    stop(errorCondition(
      "Multiple group classes are not identifiable with one profile or only singleton groups.",
      class = "latents_unidentified", call = NULL))
  }
  if (length(fixed) > 0L) {
    start <- .multilpa_complete_start(start, n_profiles, n_group_classes)
  }
  # `fixed` is checked against the start before the start's own shape is, so
  # that a request to hold a block the start does not carry is reported as
  # that, rather than as a generic complaint about the start's contents.
  fixed <- .multilpa_validate_fixed(fixed, start, covariance_model, ncol(x),
                                    n_categories)
  # Identification by counting: a complete table of discrete indicators has
  # fewer free cells than the model has parameters (review 2026-10-27, CORE-5).
  if (ncol(x) == 0L && n_group_classes == 1L && !isTRUE(noise)) {
    .multilpa_warn_negative_df(
      if (is.null(codes)) NULL else codes[carried_rows, , drop = FALSE], n_categories,
      .latents_extra_rows(extra, which(carried_rows)),
      .multilpa_count_parameters(n_profiles, n_group_classes, 0L, n_categories,
                                 variance_model, covariance_model) +
        .latents_extra_n_parameters(extra, n_profiles) -
        .multilpa_fixed_parameters(fixed, n_profiles, 0L, n_categories,
                                   variance_model, covariance_model))
  }
  # A noise start carries its proportion beside the profile probabilities;
  # it is taken out here so the ordinary validation sees the ordinary shape.
  noise_start <- NULL
  if (isTRUE(noise) && !is.null(start)) {
    split_start <- .multilpa_noise_start(start)
    start <- split_start$start
    noise_start <- split_start$noise_probability
  }
  if (!is.null(start)) {
    start <- .multilpa_validate_start(start, n_profiles, n_group_classes,
                                    ncol(x), variance_model, min_variance,
                                    covariance_model, n_categories, min_probability,
                                    categorical, encoded$levels)
  }
  .latents_local_seed(seed)
  # Center to reduce cancellation in weighted means when indicators have large
  # offsets. This is a translation, with no density Jacobian or scale change.
  centers <- .multilpa_initialization_center(x, sampling_weights[group_index])
  x <- sweep(x, 2L, centers, "-")
  # Direct effects are fitted on centred covariates too: an intercept and the
  # effect of a covariate far from zero are nearly collinear, which makes the
  # conditional M-steps crawl (4911 EM iterations for a covariate centred at
  # 41, against 21 at 1). The model is the same; the means are moved back to
  # a covariate value of zero after EM.
  fitting_extra <- .latents_center_direct(extra, x, sampling_weights[group_index])
  if (any(!is.finite(x[!is.na(x)]^2))) {
    stop(errorCondition(
      "Indicator scales overflow squared residuals; rescale the data.",
      class = "latents_bad_data", call = NULL))
  }
  if (!is.null(start)) start$means <- sweep(start$means, 2L, centers, "-")
  held <- .multilpa_held_parameters(start, fixed, covariance_model)
  # The prior is resolved on the indicators the EM sees, so its default mean is
  # their (centred) column means and a supplied mean moves with the centring.
  prior_parameters <- if (is.null(prior)) NULL else {
    if (anyNA(x)) {
      stop(errorCondition(paste(
        "`prior` is defined for complete indicators, as in mclust; these have",
        "missing values. Drop `prior`, or use complete rows."),
        class = "latents_unsupported_prior", call = NULL))
    }
    .multilpa_prior_parameters(prior, x, n_profiles, structure, offset = centers)
  }
  # The noise component's density is the reciprocal of mclust's hypervolume of
  # the data, which a translation does not change.
  noise_log_density <- if (isTRUE(noise)) {
    if (anyNA(x)) {
      stop(errorCondition(paste(
        "A noise component is a uniform density over the complete indicator",
        "space, as in mclust; these indicators have missing values."),
        class = "latents_unsupported_noise", call = NULL))
    }
    -.multilpa_log_hypervolume(x)
  } else NULL
  # `max_iter = 0` performs no update, so a restart cannot improve on anything:
  # scoring extra random initializations and keeping the highest would return a
  # parameter set nobody supplied in place of the one the caller asked to have
  # evaluated. Evaluate-only with a `start` therefore uses that start alone, and
  # `get_results(fit, "starts")` shows the single start that was run.
  if (max_iter == 0L && !is.null(start)) n_starts <- 1L
  # Under emEM a pool is screened and the best `n_starts` are run to the end;
  # the pool's first `n_starts` starts are the kmeans search's own, drawn in
  # the same order, so emEM considers every start that search would run.
  screen <- .latents_search_screen(search, n_starts, n_screen, max_iter,
                                   supplied = !is.null(start))
  start_kind <- function(start_index) .latents_start_kind(start_index, n_starts)
  # Each start runs to completion or is recorded as failed; the best is
  # chosen by the shared rule (see .latents_run_starts()).
  prepare <- function(start_index) {
    initial <- if (start_index == 1L && !is.null(start)) start else {
      # The first start is Ward's hierarchical clustering, deterministic
      # and good at profiles that differ in shape; the rest are k-means
      # from random centres, and in an emEM pool random partitions and
      # random centres after those.
      .multilpa_initialize(x, group_index, n_profiles, n_group_classes,
                         variance_model, min_variance, start_index,
                         covariance_model, codes, n_categories, min_probability,
                         hierarchical = start_index == 1L, extra = fitting_extra,
                         kind = start_kind(start_index),
                         row_weights = sampling_weights[group_index])
    }
    # Every restart begins from the held values, so a random start cannot
    # report a measurement solution that was neither estimated nor supplied.
    initial <- .multilpa_apply_held(initial, held)
    initial <- .multilpa_project_start(initial, structure, nrow(x), min_variance)
    # A random start begins with a tenth of the mass on the noise component;
    # a supplied start with the proportion it supplied.
    if (!is.null(noise_log_density)) {
      initial <- .multilpa_add_noise(
        initial,
        if (start_index == 1L && !is.null(noise_start)) noise_start else 0.1,
        noise_log_density)
    }
    initial
  }
  initial <- .latents_prepare_starts(
    if (is.null(screen)) n_starts else screen$n_pool, prepare)
  context <- as.list(environment())
  context[c("prepare", "start_kind", "initial")] <- NULL
  plan <- .lta_plan(context, initial, .multilpa_run_start, .multilpa_finalise,
    c("x", "group_index", "variance_model", "min_variance", "max_iter", "tol",
      "covariance_model", "codes", "n_categories", "min_probability", "held",
      "structure", "prior_parameters", "acceleration", "sampling_weights",
      "fitting_extra", "covariance_pattern"),
    n_observed = .latents_observed_values(
      x, codes, .latents_extra_matrix(fitting_extra),
      row_weights = if (is.null(sampling_weights)) NULL else
        sampling_weights[group_index]))
  plan$interior <- if (exclude_boundary) .lta_bind_context(
    .multilpa_start_interior, list(min_variance = min_variance))
  plan$cores <- cores
  class(plan) <- "latents_multilpa_plan"
  report_structure <- FALSE
  plan
}

#' Warn when a discrete-indicator model has negative degrees of freedom
#'
#' A latent class model of complete categorical (and ordinal) indicators is a
#' model of their contingency table, which has `prod(categories) - 1` free
#' cell probabilities. With more free parameters than that, the degrees of
#' freedom are negative and the model is not identified: different class sizes
#' and response probabilities reproduce the same table (Goodman, 1974). This
#' is a counting condition, necessary but not sufficient: a model with
#' non-negative degrees of freedom can still be locally unidentified. It is
#' judged only where the count is exact -- no continuous, count or
#' direct-effect part, one group class, every row of positive weight complete
#' -- and leaves every other model as it was.
#'
#' @param codes Categorical codes of the rows of positive weight, or `NULL`.
#' @param n_categories Categories per categorical indicator, or `NULL`.
#' @param extra Prepared ordinal and count indicators of the same rows, or
#'   `NULL`.
#' @param n_parameters The model's free parameters.
#' @return `NULL`, invisibly; raises a `latents_unidentified` warning when the
#'   degrees of freedom are negative.
#' @references Goodman, L. A. (1974). Exploratory latent structure analysis
#'   using both identifiable and unidentifiable models. Biometrika, 61(2),
#'   215-231.
#' @noRd
.multilpa_warn_negative_df <- function(codes, n_categories, extra, n_parameters) {
  ordinal <- extra$ordinal
  counted <- is.null(extra) ||
    ((ncol(extra$count) %||% 0L) == 0L && is.null(extra$direct))
  # A half-observed pair row has an NA code, so it too makes the table
  # incomplete.
  if (!counted || (is.null(codes) && is.null(ordinal)) ||
      anyNA(codes) || anyNA(ordinal)) {
    return(invisible(NULL))
  }
  # In doubles: a product of many indicators' categories overflows an integer.
  cells <- prod(as.numeric(c(n_categories, extra$ordinal_categories)))
  degrees <- cells - 1 - n_parameters
  if (degrees >= 0) return(invisible(NULL))
  warning(warningCondition(sprintf(paste(
    "The model is not identified: it has %d free parameters, but its",
    "indicators' %s response patterns leave only %s free cell probabilities",
    "(%s degrees of freedom). Its class sizes and response probabilities are",
    "not determined by the data; fit fewer classes or use more indicators."),
    as.integer(n_parameters), format(cells, big.mark = ","),
    format(cells - 1, big.mark = ","), format(degrees)),
    class = "latents_unidentified", call = NULL))
  invisible(NULL)
}

#' `match.arg()` on a formal of the calling function
#'
#' Reads the choices from the caller's formals, as `match.arg(arg)` does, and
#' matches with `.latents_match_arg()`, so an invalid value raises
#' `latents_bad_argument` naming the argument and its choices where base
#' `match.arg()` raised an unclassed error about `'arg'` (review 2026-10-27,
#' CORE-7).
#' @param arg The argument, written as the bare name of a formal of the caller.
#' @return The matched choice, a single string.
#' @noRd
.latents_match_formal <- function(arg) {
  name <- deparse(substitute(arg))
  caller <- sys.function(sys.parent())
  choices <- eval(formals(caller)[[name]], envir = environment(caller))
  .latents_match_arg(arg, choices, name)
}

#' Check that the indicators and the group column exist and are distinct
#'
#' Each broken contract has its own message naming the columns concerned,
#' with the likely spelling of a misspelt one.
#' @param data The data frame.
#' @param vars The indicator names.
#' @param id The group column (a synthetic column for single-level fits).
#' @return `NULL`, invisibly, or `latents_bad_data`.
#' @noRd
.multilpa_check_columns <- function(data, vars, id) {
  refuse <- function(message) {
    stop(errorCondition(message, class = "latents_bad_data", call = NULL))
  }
  if (nrow(data) < 2L) {
    refuse(sprintf("`data` has %d %s; a fit needs at least two.", nrow(data),
                   if (nrow(data) == 1L) "row" else "rows"))
  }
  if (length(vars) < 1L || anyNA(vars)) refuse("`vars` must name at least one column, and no NA.")
  repeated <- unique(vars[duplicated(vars)])
  if (length(repeated) > 0L) {
    refuse(sprintf("`vars` lists %s more than once.", .latents_name_list(repeated)))
  }
  doubled <- unique(names(data)[duplicated(names(data))])
  if (length(doubled) > 0L) {
    refuse(sprintf("`data` has more than one column named %s.", .latents_name_list(doubled)))
  }
  absent <- setdiff(vars, names(data))
  if (length(absent) > 0L) {
    refuse(sprintf("`vars` names %s, which `data` does not have.%s",
                   .latents_name_list(absent), .latents_did_you_mean(absent, names(data))))
  }
  if (length(id) != 1L || is.na(id)) refuse("`id` must be a single column name.")
  if (id %in% vars) {
    refuse(sprintf("`id` (`%s`) is also listed in `vars`; the group column cannot be an indicator.",
                   id))
  }
  if (!id %in% names(data)) {
    refuse(sprintf("`id` names `%s`, which `data` does not have.%s", id,
                   .latents_did_you_mean(id, names(data))))
  }
  invisible(NULL)
}

#' Validate the scalar arguments of a fit
#' @return `NULL`, invisibly; raises on the first broken contract.
#' @noRd
.multilpa_check_arguments <- function(data, vars, id, n_profiles,
                                    n_group_classes, n_starts, max_iter, tol,
                                    min_variance, min_probability, seed,
                                    categorical) {
  .multilpa_check_columns(data, vars, id)
  counts <- list(n_profiles = n_profiles, n_group_classes = n_group_classes,
                 n_starts = n_starts)
  invisible(lapply(names(counts), function(field) {
    value <- counts[[field]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value != floor(value) || value > .Machine$integer.max) {
      stop(errorCondition(sprintf("%s must be a positive integer.", field),
                          class = "latents_bad_argument", call = NULL))
    }
  }))
  # max_iter = 0 evaluates the likelihood at the supplied start without moving.
  if (!is.numeric(max_iter) || length(max_iter) != 1L || !is.finite(max_iter) ||
      max_iter < 0 || max_iter != floor(max_iter) ||
      max_iter > .Machine$integer.max) {
    stop(errorCondition("max_iter must be a nonnegative integer.",
                        class = "latents_bad_argument", call = NULL))
  }
  # `min_variance = NULL` is resolved from the indicators once they are
  # prepared (`.latents_resolve_min_variance()`).
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0 ||
      (!is.null(min_variance) && (!is.numeric(min_variance) ||
        length(min_variance) != 1L || !is.finite(min_variance) ||
        min_variance <= 0))) {
    stop(errorCondition("tol and min_variance must be finite positive numbers.",
                        class = "latents_bad_argument", call = NULL))
  }
  .multilpa_check_seed(seed)
  if (anyDuplicated(categorical) || !all(categorical %in% vars)) {
    stop(errorCondition("`categorical` must name distinct indicators listed in `vars`.",
                        class = "latents_bad_categorical", call = NULL))
  }
  invisible(NULL)
}

#' Split indicators into Gaussian and categorical blocks
#'
#' Returns the continuous indicator matrix, which may have zero columns, and the
#' encoded categorical codes, which may be `NULL`. Every measurement-side
#' contract is checked here so the fitting path can assume clean input.
#'
#' The checks that ask whether an indicator carries information -- observed
#' at all, more than one value, two or more categories -- are made on the rows
#' of positive sampling weight only. A unit of weight zero is absent from the
#' fit (it is neither counted nor classified), so it must not be able to rescue
#' an indicator that is constant among the units that are: dropping those rows
#' and keeping them give the same refusal (review 2026-10-26, ML-2). The
#' checks that protect the arithmetic (finite values, numeric type) still see
#' every row.
#'
#' @param weights `NULL`, or the name of the sampling-weight column of `data`.
#' @return A list with `continuous`, `frame`, `x`, `encoded`, `codes` and
#'   `n_categories`.
#' @noRd
.multilpa_prepare_indicators <- function(data, vars, categorical, missing,
                                       min_probability, other = character(),
                                       weights = NULL) {
  continuous <- setdiff(vars, c(categorical, other))
  contributing <- .multilpa_contributing_rows(data, weights)
  encoded <- if (length(categorical) > 0L) {
    .multilpa_encode_categorical(data[, categorical, drop = FALSE],
                                 contributing = contributing)
  } else NULL
  .latents_warn_suspect_categorical(data[, categorical, drop = FALSE], encoded)
  codes <- encoded$codes
  n_categories <- encoded$n_categories
  if (!is.null(codes) && min_probability * max(n_categories) >= 1) {
    stop(errorCondition(sprintf(
      "`min_probability` of %g leaves no room for an indicator with %d categories.",
      min_probability, max(n_categories)),
      class = "latents_bad_categorical", call = NULL))
  }
  if (!is.null(codes)) {
    .latents_refuse_missing_values(codes[contributing, , drop = FALSE], missing)
  }
  frame <- data[, continuous, drop = FALSE]
  not_numeric <- names(frame)[!vapply(frame, function(value) {
    is.numeric(value) && is.null(dim(value))
  }, logical(1))]
  if (length(not_numeric) > 0L) {
    kinds <- vapply(frame[not_numeric], function(value) class(value)[[1L]], character(1))
    stop(errorCondition(sprintf(paste(
      "%s %s not numeric (%s). A continuous indicator must be a numeric vector;",
      "declare categorical ones in `categorical`, or fit them with lca() or multilca()."),
      .latents_name_list(not_numeric), if (length(not_numeric) == 1L) "is" else "are",
      paste(sprintf("%s: %s", not_numeric, kinds), collapse = ", ")),
      class = "latents_bad_data", call = NULL))
  }
  x <- matrix(as.matrix(frame), nrow = nrow(data), ncol = length(continuous),
              dimnames = list(NULL, continuous))
  .latents_refuse_missing_values(x, missing)
  carried <- x[contributing, , drop = FALSE]
  if (any(colSums(!is.na(carried)) == 0L)) {
    stop(errorCondition(
      "Every indicator needs an observed row with positive sampling weight.",
      class = "latents_bad_data", call = NULL))
  }
  # Constancy is decided exactly, from the range of the observed values, as
  # the centring checks decide it (see .multilpa_indicator_scales()).
  constant <- vapply(seq_len(ncol(carried)), function(column) {
    values <- carried[!is.na(carried[, column]), column]
    max(values) == min(values)
  }, logical(1))
  if (any(constant)) {
    stop(errorCondition(sprintf(paste(
      "Constant indicators cannot identify Gaussian profiles: %s takes one",
      "value in every row of positive weight."),
      paste(sprintf("`%s`", continuous[constant]), collapse = ", ")),
      class = "latents_bad_data", call = NULL))
  }
  list(continuous = continuous, frame = frame, x = x, encoded = encoded,
       codes = codes, n_categories = n_categories, other = other)
}

#' Which rows carry sampling weight
#'
#' The data checks of a fit judge the rows that enter the likelihood, and a row
#' of weight zero does not. The weights themselves are validated later, by
#' `.latents_sampling_weights()`; a column that would fail there (missing,
#' not numeric, negative, all zero) leaves every row in, so the data checks
#' behave as unweighted ones and that refusal is still the one reported.
#'
#' @param data The data frame of the fit.
#' @param weights `NULL`, or the name of the sampling-weight column.
#' @return A logical vector with one element per row of `data`.
#' @noRd
.multilpa_contributing_rows <- function(data, weights = NULL) {
  stopifnot("`data` must be a data frame" = is.data.frame(data))
  every <- rep(TRUE, nrow(data))
  if (is.null(weights)) return(every)
  # Either the name of the weight column, or (the covariate engine, which has
  # already resolved them) one weight per row.
  values <- if (is.numeric(weights) && length(weights) == nrow(data)) weights else {
    if (!is.character(weights) || length(weights) != 1L || is.na(weights) ||
        !weights %in% names(data)) return(every)
    data[[weights]]
  }
  if (!is.numeric(values) || !is.null(dim(values)) || anyNA(values) ||
      any(!is.finite(values)) || any(values < 0) || !any(values > 0)) {
    return(every)
  }
  values > 0
}

#' Index the observed groups in first-occurrence order
#' @return A list with `values`, `index`, `ids`, `n` and `sizes`.
#' @noRd
.multilpa_prepare_groups <- function(raw_groups) {
  if (!(is.character(raw_groups) || is.factor(raw_groups) || is.numeric(raw_groups)) ||
      !is.null(dim(raw_groups)) || anyNA(raw_groups) ||
      (is.numeric(raw_groups) && any(!is.finite(raw_groups)))) {
    stop(errorCondition(
      "group must contain nonmissing, finite numeric, factor, or character identifiers.",
      class = "latents_bad_data", call = NULL))
  }
  values <- unique(raw_groups)
  index <- match(raw_groups, values)
  ids <- make.unique(as.character(values))
  list(values = values, index = index, ids = ids, n = length(ids),
       sizes = tabulate(index, nbins = length(ids)))
}

#' Attach readable dimension names to the fitted parameters
#' @return The parameter list, with profiles, indicators and categories named.
#' @noRd
.multilpa_label_parameters <- function(parameters, profile_names, continuous,
                                     categorical, levels) {
  dimnames(parameters$means) <- dimnames(parameters$variances) <-
    list(profile_names, continuous)
  if (!is.null(parameters$covariances)) {
    dimnames(parameters$covariances) <- list(continuous, continuous, profile_names)
  }
  if (!is.null(parameters$response_probabilities)) {
    parameters$response_probabilities <- stats::setNames(
      lapply(seq_along(categorical), function(index) {
        block <- parameters$response_probabilities[[index]]
        dimnames(block) <- list(profile_names, levels[[categorical[index]]])
        block
      }), categorical)
  }
  parameters
}

#' Count the free parameters of a fitted model
#'
#' Mixing proportions contribute `(H - 1) + H (K - 1)`. Gaussian indicators add
#' means plus variances or covariances, shared or profile-specific. Categorical
#' measurement is always profile-specific, so `variance_model` has no
#' categorical analogue and never reduces that block.
#'
#' @return A single number of free parameters.
#' @noRd
.multilpa_count_parameters <- function(n_profiles, n_group_classes, n_continuous,
                                     n_categories, variance_model,
                                     covariance_model, structure = NULL) {
  covariance_parameters <- if (covariance_model == "full") {
    n_continuous * (n_continuous + 1) / 2
  } else n_continuous
  gaussian <- if (n_continuous == 0L) 0L else if (!is.null(structure)) {
    # The constrained structures divide a determinant-one shape out of a
    # volume, so their spread block is not `n_profiles` copies of anything.
    n_profiles * n_continuous +
      .multilpa_structure_parameters(structure, n_profiles, n_continuous)
  } else {
    n_profiles * n_continuous +
      if (variance_model == "varying") n_profiles * covariance_parameters else
        covariance_parameters
  }
  categorical <- if (is.null(n_categories)) 0L else
    .multilpa_categorical_parameters(n_profiles, n_categories)
  (n_group_classes - 1L) + n_group_classes * (n_profiles - 1L) +
    gaussian + categorical
}

#' The parameters of one profile's own Gaussian density
#'
#' Its free means plus what its spread adds to the model: every variance or
#' covariance when profiles vary, nothing when they share one, and the
#' profile's share of a constrained structure (the parameters one more
#' profile adds). Covariances held at zero by a pattern are excluded, as
#' are held means; held variances cannot collapse during this fit.
#'
#' @param n_continuous Number of continuous indicators.
#' @param variance_model,covariance_model,structure As in [multilpa()].
#' @param covariance_pattern The free covariance entries, or `NULL`.
#' @param fixed The expanded names of held measurement blocks.
#' @return A single number; `0` when there are no continuous indicators or
#'   the spread is shared by every profile or held fixed, so no profile has
#'   a free variance of its own to collapse.
#' @noRd
.latents_profile_parameters <- function(n_continuous, variance_model,
                                        covariance_model, structure = NULL,
                                        covariance_pattern = NULL,
                                        fixed = character()) {
  stopifnot(is.numeric(n_continuous), length(n_continuous) == 1L,
            is.finite(n_continuous), n_continuous >= 0,
            is.character(fixed), !anyNA(fixed),
            is.null(covariance_pattern) ||
              (is.matrix(covariance_pattern) && is.logical(covariance_pattern) &&
                 identical(dim(covariance_pattern), rep(as.integer(n_continuous), 2L))))
  if (n_continuous == 0L || any(c("variances", "measurement") %in% fixed)) return(0)
  own_spread <- if (!is.null(structure)) {
    .multilpa_structure_parameters(structure, 2L, n_continuous) -
      .multilpa_structure_parameters(structure, 1L, n_continuous)
  } else if (identical(variance_model, "varying")) {
    if (identical(covariance_model, "full")) n_continuous * (n_continuous + 1) / 2 else
      n_continuous
  } else 0
  if (own_spread == 0) return(0)
  if (!is.null(covariance_pattern)) {
    own_spread <- own_spread - .latents_pattern_zeros(covariance_pattern)
  }
  own_spread + if ("means" %in% fixed) 0L else n_continuous
}

#' Warn when a profile has fewer members than its own parameters
#'
#' This is a conservative membership-versus-parameter heuristic for profiles
#' with class-specific spread. It is not a rank or identification test:
#' positive soft memberships can yield full-rank covariance at arbitrarily
#' small total mass, and a larger profile can be nearly flat. Under sampling
#' weights the membership is normalized mass, not an independent sample count.
#' The estimate is unchanged; continuous covariance diagnostics are reported
#' separately without using a rejection threshold.
#'
#' @param effective_counts Effective membership of each profile.
#' @param n_parameters The output of `.latents_profile_parameters()`.
#' @param context Optional text identifying a separately estimated measurement
#'   block, such as `"at occasion 2"`.
#' @return The positions of the thin profiles, invisibly; a
#'   `latents_thin_profiles` warning when there are any.
#' @noRd
.latents_warn_thin_profiles <- function(effective_counts, n_parameters, context = NULL) {
  stopifnot(is.numeric(effective_counts), is.numeric(n_parameters),
            length(n_parameters) == 1L, is.finite(n_parameters), n_parameters >= 0,
            is.null(context) ||
              (is.character(context) && length(context) == 1L &&
                 !is.na(context) && nzchar(context)))
  thin <- which(effective_counts < n_parameters)
  if (n_parameters == 0 || length(effective_counts) < 2L || length(thin) == 0L) {
    return(invisible(integer()))
  }
  warning(warningCondition(sprintf(paste(
    "%s %s%s: %s Gaussian membership mass for %d free means and variances%s.",
    "This is a conservative diagnostic, not an identification test;",
    "inspect the Gaussian covariance diagnostics."),
    if (length(thin) == 1L) "Profile" else "Profiles",
    paste(thin, collapse = ", "),
    if (is.null(context)) "" else paste0(" ", context),
    paste(format(round(effective_counts[thin], 1), nsmall = 1), collapse = ", "),
    as.integer(n_parameters), if (length(thin) == 1L) "" else " each"),
    class = "latents_thin_profiles", call = NULL))
  invisible(thin)
}

#' Hand a covariate request to the covariate estimator
#'
#' `multilpa()` is one verb over one model, and membership covariates are part
#' of that model rather than a different one, so they are arguments and not a
#' separate entry point. Two of `multilpa()`'s own arguments have no meaning
#' on this path: the covariate likelihood has no starting-value contract of
#' the shape the covariate-free EM uses, and no held-measurement machinery.
#' Each is refused by name instead of being accepted and ignored, which would
#' return a fit that is not the one asked for. `missing = "fiml"` is supported:
#' the measurement density is the covariate-free model's observed-data one.
#'
#' @param data,vars,id,n_profiles,n_group_classes As in [multilpa()].
#' @param profile_covariates,group_covariates The requested predictors.
#' @param variance_model,n_starts,max_iter,tol,min_variance As in [multilpa()].
#' @param seed,time,covariance_model,categorical,min_probability As in
#'   [multilpa()].
#' @param missing,profile_slopes As in [multilpa()].
#' @param start,fixed Refused when they are anything but their default.
#' @param boundary As in [multilpa()].
#' @param call The user's call, recorded on the result.
#' @return An object of class `multilpa_covariates`.
#' @noRd
.multilpa_covariate_model <- function(data, vars, id, n_profiles,
                                      n_group_classes, profile_covariates,
                                      group_covariates, variance_model,
                                      n_starts, max_iter, tol, min_variance,
                                      seed, start, missing,
                                      profile_slopes = "shared", covariance_model,
                                      categorical, min_probability, time,
                                      fixed, select_start = "likelihood",
                                      boundary = "allow",
                                      call, weights = NULL, ordinal = character(),
                                      count = character(), count_model = "poisson",
                                      count_dispersion = "varying",
                                      centering = "none", scaling = "none",
                                      drop_empty = TRUE, search = "kmeans",
                                      n_screen = NULL) {
  unsupported <- c(
    start = !is.null(start),
    fixed = length(fixed) > 0L)
  if (any(unsupported)) {
    stop(errorCondition(sprintf(paste(
      "%s cannot be combined with `profile_covariates` or",
      "`group_covariates`. Fit the covariate model without %s, or fit the",
      "covariate-free model and use three_step() or r3step()."),
      paste(sprintf("`%s`", names(unsupported)[unsupported]), collapse = ", "),
      if (sum(unsupported) > 1L) "them" else "it"),
      class = "latents_bad_argument", call = NULL))
  }
  .multilpa_fit_covariates(
    data = data, vars = vars, id = id, n_profiles = n_profiles,
    n_group_classes = n_group_classes,
    profile_covariates = profile_covariates,
    group_covariates = group_covariates, variance_model = variance_model,
    n_starts = n_starts, max_iter = max_iter, tol = tol,
    min_variance = min_variance, seed = seed, time = time,
    covariance_model = covariance_model, categorical = categorical,
    min_probability = min_probability, missing = missing,
    profile_slopes = profile_slopes, select_start = select_start,
    boundary = boundary, call = call,
    weights = weights, ordinal = ordinal, count = count,
    count_model = count_model, count_dispersion = count_dispersion,
    centering = centering, scaling = scaling, drop_empty = drop_empty,
    search = search, n_screen = n_screen)
}

#' One k-means run for a random start
#'
#' Lloyd's algorithm warns "empty cluster" when a start's random centres leave
#' a cluster without members. That start then fails or converges like any
#' other, and the fit reports failed starts itself ("k of n starts failed"),
#' so base R's warning is muffled here rather than repeated to the caller.
#' @param x Numeric matrix to cluster.
#' @param centers Number of clusters.
#' @return The `stats::kmeans()` result.
#' @noRd
.multilpa_start_kmeans <- function(x, centers) {
  withCallingHandlers(
    stats::kmeans(x, centers = centers, iter.max = 100L, algorithm = "Lloyd"),
    warning = function(condition) {
      # Lloyd's k-means warns when it stops at `iter.max` and also sets
      # `ifault` to 2, which the caller turns into a failed start that the
      # starts table records; the bare warning would only repeat it, worded as
      # though the model had not converged.
      message <- conditionMessage(condition)
      if (startsWith(message, "empty cluster") ||
          startsWith(message, "did not converge in")) {
        invokeRestart("muffleWarning")
      }
    })
}
