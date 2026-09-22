#' Row-wise log-sum-exp
#' @param log_values Numeric matrix of log weights.
#' @return Vector of log sums.
#' @noRd
.multilpa_log_sum_exp <- function(log_values) {
  stopifnot(is.matrix(log_values), is.numeric(log_values),
            nrow(log_values) > 0L, ncol(log_values) > 0L,
            !anyNA(log_values))
  row_max <- apply(log_values, 1L, max)
  if (any(!is.finite(row_max))) {
    stop("All component densities vanished, or a density overflowed.")
  }
  row_max + log(rowSums(exp(sweep(log_values, 1L, row_max, "-"))))
}

#' Compute nested posterior probabilities
#' @param x Observation-by-indicator matrix of continuous indicators, optionally
#'   containing NA. It may have zero columns when every indicator is categorical.
#' @param group_index Integer group indices in first-occurrence order.
#' @param parameters List of model parameters.
#' @param codes Optional integer matrix of categorical indicator codes.
#' @return Log likelihood, posteriors, and conditional class probabilities.
#' @noRd
.multilpa_expectation <- function(x, group_index, parameters, codes = NULL) {
  stopifnot(is.matrix(x), is.numeric(x), !any(is.infinite(x)),
            length(group_index) == nrow(x), is.list(parameters))
  n_profiles <- nrow(parameters$means)
  n_types <- length(parameters$group_probabilities)
  gaussian <- if (ncol(x) > 0L && (anyNA(x) || !is.null(parameters$covariances))) {
    .multilpa_gaussian_moments(x, parameters)
  } else NULL
  log_density <- if (!is.null(gaussian)) gaussian$log_density else vapply(seq_len(n_profiles), function(profile) {
    residuals <- sweep(x, 2L, parameters$means[profile, ], "-")
    -0.5 * rowSums(sweep(residuals^2, 2L,
                        parameters$variances[profile, ], "/") +
                    matrix(log(2 * pi) + log(parameters$variances[profile, ]),
                           nrow(x), ncol(x), byrow = TRUE))
  }, numeric(nrow(x)))
  # Remove the common measurement offset before adding log priors. Otherwise
  # a very small density (for example log(f) = -5e15) rounds away the priors,
  # and subtracting the absolute log marginal can produce posteriors above one.
  density_offset <- apply(log_density, 1L, max)
  if (any(!is.finite(density_offset))) {
    stop("All component densities vanished, or a density overflowed.")
  }
  log_density <- sweep(log_density, 1L, density_offset, "-")
  # Categorical indicators are conditionally independent of the continuous ones
  # given the profile. Add their log densities after removing the Gaussian
  # offset so that small categorical contributions retain their precision too.
  if (!is.null(codes)) {
    log_density <- log_density +
      .multilpa_categorical_log_density(codes, parameters$response_probabilities)
  }
  # Each group type supplies a different prior over the same Gaussian profiles.
  conditional <- lapply(seq_len(n_types), function(group_type) {
    log_joint <- sweep(log_density, 2L,
                       log(parameters$profile_probabilities[group_type, ]), "+")
    log_marginal <- .multilpa_log_sum_exp(log_joint)
    posterior <- exp(sweep(log_joint, 1L, apply(log_joint, 1L, max), "-"))
    list(log_marginal = log_marginal,
         posterior = posterior / rowSums(posterior))
  })
  group_scores <- matrix(vapply(seq_len(n_types), function(group_type) {
    as.numeric(rowsum(conditional[[group_type]]$log_marginal,
                      group_index, reorder = FALSE))
  }, numeric(max(group_index))), nrow = max(group_index), ncol = n_types)
  group_offset <- apply(group_scores, 1L, max)
  group_scores <- sweep(sweep(group_scores, 1L, group_offset, "-"), 2L,
                        log(parameters$group_probabilities), "+")
  group_log_likelihood <- .multilpa_log_sum_exp(group_scores) + group_offset +
    as.numeric(rowsum(density_offset, group_index, reorder = FALSE))
  group_posteriors <- exp(sweep(group_scores, 1L, apply(group_scores, 1L, max), "-"))
  group_posteriors <- group_posteriors / rowSums(group_posteriors)
  joint <- lapply(seq_len(n_types), function(group_type) {
    conditional[[group_type]]$posterior * group_posteriors[group_index, group_type]
  })
  subject_posteriors <- Reduce(`+`, joint)
  log_likelihood <- sum(group_log_likelihood)
  if (!is.finite(log_likelihood) || any(!is.finite(subject_posteriors))) {
    stop("Non-finite likelihood or posterior probabilities.")
  }
  list(log_likelihood = log_likelihood,
       group_log_likelihood = group_log_likelihood,
       group_posteriors = group_posteriors,
       subject_posteriors = subject_posteriors, joint = joint,
       gaussian_moments = gaussian$moments)
}

#' Maximize the expected complete-data log likelihood
#' @param x Observation-by-indicator matrix.
#' @param expectation Nested posterior calculation.
#' @param variance_model Either varying or equal across profiles.
#' @param min_variance Lower bound on every indicator variance.
#' @param covariance_model Diagonal or full residual covariance.
#' @return Updated model parameters.
#' @noRd
.multilpa_maximization <- function(x, expectation, variance_model, min_variance,
                                  covariance_model = "diagonal", codes = NULL,
                                  n_categories = NULL, min_probability = 1e-10,
                                  held = NULL, structure = NULL,
                                  previous = NULL) {
  stopifnot(is.matrix(x), is.list(expectation),
            variance_model %in% c("varying", "equal"), min_variance > 0,
            is.null(held) || is.list(held))
  weights <- colSums(expectation$subject_posteriors)
  group_weights <- colSums(expectation$group_posteriors)
  if (any(weights <= 0) || any(group_weights <= 0)) {
    stop("A profile or group class has zero effective membership.")
  }
  gaussian <- if (!is.null(expectation$gaussian_moments)) {
    .multilpa_maximize_moments(x, expectation, variance_model, min_variance,
                               covariance_model, held, structure, previous)
  } else NULL
  if (ncol(x) == 0L) {
    means <- matrix(numeric(0), length(weights), 0L)
    variances <- matrix(numeric(0), length(weights), 0L)
  } else if (is.null(gaussian)) {
    # A held mean is the point the spread is measured around, so it replaces the
    # weighted mean before the residuals are formed rather than afterwards.
    means <- held$means %||%
      sweep(crossprod(expectation$subject_posteriors, x), 1L, weights, "/")
    variance_sums <- t(matrix(vapply(seq_along(weights), function(profile) {
      residuals <- sweep(x, 2L, means[profile, ], "-")
      colSums(residuals^2 * expectation$subject_posteriors[, profile])
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
        matrix(colSums(variance_sums) / nrow(x),
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
  profile_counts <- t(matrix(vapply(expectation$joint, colSums, numeric(length(weights))),
                             nrow = length(weights), ncol = length(group_weights)))
  profile_probabilities <- profile_counts / rowSums(profile_counts)
  parameters <- list(means = means, variances = variances,
                     profile_probabilities = profile_probabilities,
                     group_probabilities = group_weights / sum(group_weights))
  if (!is.null(gaussian$covariances)) parameters$covariances <- gaussian$covariances
  if (!is.null(codes)) {
    parameters$response_probabilities <- held$response_probabilities %||%
      .multilpa_categorical_maximize(codes, expectation$subject_posteriors,
                                     n_categories, min_probability)
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
#' @return Strictly positive initial mixture probabilities and Gaussian parameters.
#' @noRd
.multilpa_initialize <- function(x, group_index, n_profiles, n_types,
                               variance_model, min_variance, start_index,
                               covariance_model = "diagonal", codes = NULL,
                               n_categories = NULL, min_probability = 1e-10) {
  stopifnot(is.matrix(x), length(group_index) == nrow(x),
            n_profiles >= 1L, n_types >= 1L, start_index >= 1L,
            variance_model %in% c("varying", "equal"), min_variance > 0)
  if (ncol(x) > 0L && anyNA(x)) {
    initial_means <- colMeans(x, na.rm = TRUE)
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
    standardized <- sweep(sweep(x, 2L, overall_means, "-"),
                          2L, sqrt(overall_variances), "/")
    # Categorical indicators join the clustering as centered dummy columns, so
    # a mixed-mode model is initialized from all of its indicators at once.
    clustering_input <- if (is.null(codes)) standardized else
      cbind(standardized, .multilpa_categorical_design(codes, n_categories))
    initial_fit <- stats::kmeans(clustering_input, centers = n_profiles,
                                iter.max = 100L, algorithm = "Lloyd")
    if (!is.null(initial_fit$ifault) && initial_fit$ifault != 0L) {
      stop("Initial k-means did not converge.")
    }
    assignments <- initial_fit$cluster
    centers <- initial_fit$centers[, seq_len(ncol(x)), drop = FALSE]
    means <- sweep(sweep(centers, 2L, sqrt(overall_variances), "*"),
                   2L, overall_means, "+")
  }
  variances <- matrix(pmax(overall_variances, min_variance),
                      n_profiles, ncol(x), byrow = TRUE)
  indicator <- vapply(seq_len(n_profiles), function(profile) {
    as.numeric(assignments == profile)
  }, numeric(nrow(x)))
  group_counts <- rowsum(indicator, group_index, reorder = FALSE)
  # Pseudocounts are used only to initialize, never to estimate the final model.
  group_composition <- (group_counts + 0.5) /
    (rowSums(group_counts) + 0.5 * n_profiles)
  if (n_types == 1L) {
    profile_probabilities <- matrix(colMeans(indicator), 1L)
  } else if (start_index %% 2L == 1L &&
             nrow(unique(group_composition)) >= n_types) {
    group_fit <- stats::kmeans(group_composition, centers = n_types,
                              iter.max = 100L, algorithm = "Lloyd")
    if (!is.null(group_fit$ifault) && group_fit$ifault != 0L) {
      stop("Initial group k-means did not converge.")
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
  if (ncol(x) > 0L && covariance_model == "full") {
    residuals <- sweep(x, 2L, overall_means, "-")
    covariance <- .multilpa_bound_covariance(crossprod(residuals) / nrow(x), min_variance)
    result$covariances <- array(rep(covariance, n_profiles), c(ncol(x), ncol(x), n_profiles))
    result$variances <- matrix(diag(covariance), n_profiles, ncol(x), byrow = TRUE)
  }
  result
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
      class = "multilpa_bad_start", call = NULL))
    }
    invisible(lapply(seq_len(n_profiles), function(profile) {
      covariance <- matrix(covariances[, , profile], n_indicators, n_indicators)
      if (max(abs(covariance - t(covariance))) > 1e-12 ||
          min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values) < min_variance * (1 - 1e-8)) {
        stop(errorCondition(
          "Starting covariances must be symmetric with eigenvalues at least min_variance.",
          class = "multilpa_bad_start", call = NULL))
      }
      if (variance_model == "equal" && max(abs(covariance - covariances[, , 1L])) > 1e-12) {
        stop(errorCondition("Starting covariances must be equal across profiles.",
                            class = "multilpa_bad_start", call = NULL))
      }
    }))
    variances <- t(matrix(vapply(seq_len(n_profiles), function(profile) {
      diag(matrix(covariances[, , profile], n_indicators, n_indicators))
    }, numeric(n_indicators)), n_indicators, n_profiles))
    if (!is.null(start$variances) && !isTRUE(all.equal(unname(start$variances), variances, tolerance = 1e-12))) {
      stop(errorCondition("Starting variances must agree with covariance diagonals.",
                          class = "multilpa_bad_start", call = NULL))
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
                        class = "multilpa_bad_start", call = NULL))
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
        class = "multilpa_bad_start", call = NULL))
    }
  }))
  if (any(start$variances < min_variance)) {
    stop(errorCondition("start$variances must be at least min_variance.",
                        class = "multilpa_bad_start", call = NULL))
  }
  if (variance_model == "equal" &&
      any(abs(sweep(start$variances, 2L, start$variances[1L, ], "-")) > 1e-12)) {
    stop(errorCondition(
      "start$variances must be equal across profiles for variance_model = 'equal'.",
      class = "multilpa_bad_start", call = NULL))
  }
  group_probabilities <- start$group_probabilities
  if (!is.numeric(group_probabilities) || !is.null(dim(group_probabilities)) ||
      length(group_probabilities) != n_types || any(!is.finite(group_probabilities)) ||
      any(group_probabilities <= 0) || abs(sum(group_probabilities) - 1) > 1e-8 ||
      any(start$profile_probabilities <= 0) ||
      any(abs(rowSums(start$profile_probabilities) - 1) > 1e-8)) {
    stop(errorCondition(
      "Starting probabilities must be strictly positive and sum to one (within each profile-probability row).",
      class = "multilpa_bad_start", call = NULL))
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
#'   columns in this fit's category order. Raises `multilpa_bad_start` when a
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
        class = "multilpa_bad_start", call = NULL))
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
        class = "multilpa_bad_start", call = NULL))
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
      length(n_categories)), class = "multilpa_bad_start", call = NULL))
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
        class = "multilpa_bad_start", call = NULL))
    }
    if (any(block < 0) || any(abs(rowSums(block) - 1) > 1e-8)) {
      stop(errorCondition(sprintf(
        "start$response_probabilities[[%d]] rows must be nonnegative and sum to one.",
        indicator), class = "multilpa_bad_start", call = NULL))
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
#' @return Final parameters, posterior calculation, and convergence history.
#' @noRd
.multilpa_em <- function(x, group_index, parameters, variance_model,
                        min_variance, max_iter, tol, covariance_model = "diagonal",
                        codes = NULL, n_categories = NULL, min_probability = 1e-10,
                        held = NULL, structure = NULL) {
  # `parameters` is the current point, which the M-step uses to warm start the
  # covariance structures that iterate.
  stopifnot(is.matrix(x), is.list(parameters), max_iter >= 0L, tol > 0,
            length(group_index) == nrow(x), min_variance > 0,
            variance_model %in% c("varying", "equal"))
  expectation <- .multilpa_expectation(x, group_index, parameters, codes)
  history <- expectation$log_likelihood
  converged <- FALSE
  iteration <- 0L
  while (iteration < max_iter && !converged) {
    parameters <- .multilpa_maximization(x, expectation, variance_model, min_variance,
                                       covariance_model, codes, n_categories,
                                       min_probability, held, structure,
                                       parameters)
    updated <- .multilpa_expectation(x, group_index, parameters, codes)
    improvement <- updated$log_likelihood - expectation$log_likelihood
    if (improvement < -1e-10 * (1 + abs(expectation$log_likelihood))) {
      stop("EM likelihood decreased beyond numerical roundoff.")
    }
    converged <- abs(improvement) <= tol * (1 + abs(expectation$log_likelihood))
    iteration <- iteration + 1L
    history <- c(history, updated$log_likelihood)
    expectation <- updated
  }
  list(parameters = parameters, expectation = expectation,
       converged = converged, iterations = iteration, history = history)
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
  stopifnot(
    "`time` must be NULL or a single column name" =
      is.character(time) && length(time) == 1L && !is.na(time),
    "`time` must name a column of `data`" = time %in% names(data),
    "`time` must not be one of `vars` or the `id` column" =
      !identical(time, id) && !time %in% vars
  )
  values <- data[[time]]
  if (anyNA(values)) {
    stop(errorCondition("`time` must not contain missing values.",
                        class = "multilpa_bad_time", call = NULL))
  }
  # split() first converts numeric keys to factor labels, which can collapse
  # distinct numeric identifiers that print identically. Match the native keys.
  group_index <- match(data[[id]], unique(data[[id]]))
  duplicated_within <- any(vapply(split(values, group_index),
                                  function(v) anyDuplicated(v) > 0L, logical(1)))
  if (duplicated_within) {
    stop(errorCondition(
      "`time` must be unique within each group; a position cannot occur twice.",
      class = "multilpa_bad_time", call = NULL))
  }
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
#' column, that `get_data(x, "data")` will hand it back, and that a second call
#' must spell it the same way -- plumbing the package should own. This
#' package's own validation harness wrote that line once per comparison before
#' this argument existed.
#'
#' @param data The caller's data frame.
#' @param n_group_classes What the caller asked for, or `NULL` if they did not.
#' @return A list with `data`, `id` and `n_group_classes`.
#' @noRd
.multilpa_single_level <- function(data, n_group_classes) {
  stopifnot("`data` must be a data frame" = is.data.frame(data))
  name <- ".observation"
  if (name %in% names(data)) {
    stop(errorCondition(sprintf(paste(
      "`id = NULL` gives every observation its own unit and names that unit",
      "`%s`, which these data already have as a column. Rename it, or pass it",
      "as `id` if it is the nesting unit."), name),
      class = "multilpa_bad_data", call = NULL))
  }
  if (!is.null(n_group_classes) && !isTRUE(as.integer(n_group_classes) == 1L)) {
    stop(errorCondition(sprintf(paste(
      "`id = NULL` puts one observation in each unit, so there is no",
      "composition for a second-level class to differ in and",
      "`n_group_classes = %s` has nothing to estimate. Pass `id` to fit group",
      "classes, or leave `n_group_classes` alone for a single-level fit."),
      format(n_group_classes)),
      class = "multilpa_bad_argument", call = NULL))
  }
  warning(warningCondition(paste(
    "`id = NULL` fits a single-level model: every observation is its own unit",
    "and no second level is estimated. This package exists for the two-level",
    "model; pass the nesting column as `id` to fit it."),
    class = "multilpa_single_level", call = NULL))
  data[[name]] <- seq_len(nrow(data))
  list(data = data, id = name, n_group_classes = 1L)
}

#' Fit a two-level latent profile model
#'
#' Fits individual Gaussian profiles nested within observed groups. A discrete
#' latent group class determines the profile proportions. Profile means and
#' variances are shared across group classes (measurement invariance), and
#' indicators are independent conditional on individual profile membership by
#' default. Full residual covariance and observed-data maximum likelihood for
#' missing indicators are available.
#'
#' @param data A data frame containing indicators and a group identifier.
#' @param vars Unique character vector of continuous indicator column names.
#' @param id Name of the observed group identifier column. Character, factor, or
#'   numeric identifiers are supported; missing identifiers are not.
#'
#'   `id` has no default: omitting it raises `multilpa_bad_argument`, because a
#'   forgotten grouping would otherwise be fitted as a different model without
#'   saying so. Passing `id = NULL` explicitly fits a **single-level** model:
#'   the observations are treated as independent, each row is its own unit, and
#'   `n_group_classes` becomes one. That fit raises a `multilpa_single_level`
#'   warning, since this package exists for the two-level model. This is the ordinary Gaussian or latent-class mixture that the
#'   two-level model reduces to, and every verb of this package works on it. The
#'   unit column is fabricated internally as `.observation`; it is not returned
#'   by `get_data(x, "data")`, and a `data` that already has a column of that
#'   name raises `multilpa_bad_data`. Asking for more than one group class
#'   without an `id` raises `multilpa_bad_argument`, because one observation per
#'   unit leaves no composition for a second-level class to differ in.
#' @param n_profiles Positive integer number of individual profiles.
#' @param n_group_classes Positive integer number of latent group classes. With
#'   `id = NULL` it is one, and naming anything else is an error.
#' @param profile_covariates,group_covariates Names of numeric columns of
#'   `data` predicting individual-profile and group-class membership through
#'   multinomial logits, with the final class as reference. Naming either one
#'   fits the one-step covariate model and returns a `multilpa_covariates`
#'   object: profile slopes are shared across group classes, profile
#'   intercepts differ by group class, and a `group_covariate` must be
#'   constant within each group. Covariates enter in their supplied units, so
#'   centre or scale them beforehand if that is what you want. This is one-step
#'   maximum likelihood, not a regression on assigned classes; [three_step()]
#'   and [r3step()] are the staged alternatives. The covariate path supports
#'   neither `start`, nor `missing = "fiml"`, nor `fixed`, and refuses them by
#'   name rather than ignoring them.
#'   Use one for an ordinary, pooled LPA.
#' @param variance_model Either `"varying"` (profile-specific variances or
#'   covariance matrices) or `"equal"` (shared across profiles).
#' @param n_starts Positive integer number of EM starts. When `start` is supplied,
#'   it supplies the first start; remaining starts are random initializations.
#'   It is ignored when `max_iter = 0` and `start` is supplied: that call
#'   evaluates the supplied parameters and nothing else, so exactly one start is
#'   run whatever `n_starts` says.
#' @param max_iter Nonnegative integer maximum number of EM updates per start.
#'   `max_iter = 0` performs no update. With `start`, it returns the model
#'   evaluated at exactly those values, whatever `n_starts` is, so that
#'   `logLik()` scores a parameter set supplied from elsewhere rather than the
#'   best of some random initializations that were never asked for. Without
#'   `start` there is nothing to evaluate at, so the random initializations are
#'   scored and the highest is returned.
#' @param tol Positive relative log-likelihood tolerance. Convergence requires
#'   absolute change no greater than `tol * (1 + abs(previous log likelihood))`.
#' @param min_variance Positive lower bound on each variance, or each covariance
#'   eigenvalue for full covariance, in squared input units. This defines a constrained maximum-likelihood problem. Bound-active
#'   estimates are explicitly reported and generate a warning.
#' @param seed Optional nonnegative integer random seed. With a supplied seed,
#'   the caller's random-number state is restored on exit.
#' @param start Optional list of `means` and `variances` (profiles by indicators),
#'   `profile_probabilities` (group classes by profiles), and
#'   `group_probabilities` (vector). Starting probabilities must be positive.
#'   For full covariance, supply `covariances` (indicators by indicators by
#'   profiles); `variances` may be omitted or must match their diagonals.
#'   A categorical model also takes `response_probabilities`, a list of one
#'   profiles-by-categories matrix per indicator named in `categorical`. That
#'   list is read by position unless it is labelled: name its elements after the
#'   indicators, or its columns after the categories, and each block is matched
#'   to the indicator and category its labels name, so a differently ordered
#'   `categorical` or a differently ordered set of factor levels cannot attach a
#'   distribution to the wrong item. A label naming an indicator or a category
#'   this fit does not have raises `multilpa_bad_start` rather than being
#'   aligned by position.
#' @param categorical Character vector naming indicators to treat as
#'   categorical. Each is modelled by unrestricted, profile-specific response
#'   probabilities over its observed categories, which is the latent class
#'   measurement model. Binary, ordinal and unordered indicators are all handled
#'   by the same unrestricted parameterization; numeric, integer, logical,
#'   character and factor columns are accepted. Indicators not named here stay
#'   Gaussian, so naming a subset fits a mixed-mode model.
#' @param volume,shape,orientation The covariance structure, in the three
#'   pieces it is made of. Each profile's covariance decomposes as
#'   `Sigma_k = lambda_k * D_k * A_k * D_k'`: a *volume*
#'   `lambda_k = |Sigma_k|^(1/d)`, an *orientation* `D_k` of eigenvectors, and
#'   a *shape* `A_k`, diagonal with determinant one. Constraining the three
#'   across profiles gives the fourteen models `mclust` names with three
#'   letters, and this package fits all of them.
#'
#'   `volume` is `"equal"` or `"varying"`. `shape` is `"equal"`, `"varying"`
#'   or `"spherical"`, the last making every indicator's spread equal within a
#'   profile, which leaves no orientation to constrain. `orientation` is
#'   `"axis"` (axis-parallel, a diagonal covariance), `"equal"` (one
#'   orientation shared by every profile) or `"varying"`. Each is `NULL` by
#'   default, which follows `variance_model` and `covariance_model`, so a call
#'   that names none of them fits exactly what it always did.
#'
#'   | `volume` | `shape` | `orientation` | model |
#'   |---|---|---|---|
#'   | equal | spherical | --- | EII |
#'   | varying | spherical | --- | VII |
#'   | equal/varying | equal/varying | axis | EEI, VEI, EVI, VVI |
#'   | equal/varying | equal/varying | equal | EEE, VEE, EVE, VVE |
#'   | equal/varying | equal/varying | varying | EEV, VEV, EVV, VVV |
#'
#'   The estimates are those of Celeux and Govaert (1995); the two models with
#'   a shared orientation and a free shape, EVE and VVE, have no closed form
#'   and use the minorize-maximize step of Browne and McNicholas (2014).
#'   Parameter counts match `mclust`'s own for all fourteen.
#'
#'   Anything other than EEI, VVI, EEE or VVV is maximized across every profile
#'   at once, so it cannot be combined with a held `variances` block, and
#'   `parameter_inference(method = "wald")` refuses it with
#'   `multilpa_unsupported_inference`: the free coordinates are log variances,
#'   which is the wrong chart for a constrained volume, shape or orientation.
#'   `parameter_inference(method = "bootstrap")` reports all fourteen: it
#'   resamples groups and refits inside the same family, so it needs no chart.
#' @param centering How to centre the continuous indicators before fitting.
#'   `"none"`, the default, fits them as supplied. `"person"` subtracts each
#'   group's own mean from its rows, so a value reads as a deviation from that
#'   unit's average and the profiles become profiles of *change* rather than of
#'   level: this is the within-person, person-mean-centred design (Quintana,
#'   2021; Voelkle, Brose, Schmiedek, & Lindenberger, 2014). `"grand"`
#'   subtracts one mean per indicator, which moves the origin without touching
#'   the within-group structure. The offsets are kept on the fit, so
#'   `get_data(x, "data")` still returns the columns you supplied and every
#'   verb that checks row alignment still checks it. Centring removes exactly
#'   the between-unit variation, so `"person"` refuses with
#'   `multilpa_bad_data` when it leaves an indicator constant --- which is what
#'   happens when a unit has one observation of it. With `"person"` the group
#'   classes become types of *change pattern*, not types of unit.
#' @param time Optional name of a column giving each observation's position
#'   within its group, such as a wave, occasion or course number. The model does
#'   not use it; it is stored so that `get_data(x, "sequences")`,
#'   `get_data(x, "sequence_summary")` and
#'   `plot(what = "sequences")` can read the assignments back in order. Values
#'   must be complete and unique within each group.
#' @param fixed Character vector naming measurement blocks to hold at the
#'   values `start` supplies, instead of estimating them: any of `"means"`,
#'   `"variances"` and `"response_probabilities"`, or `"measurement"` for every
#'   block the model has. Under `covariance_model = "full"`, `"variances"`
#'   holds the residual covariance matrices. A held block stays exactly as
#'   supplied, in every restart, and stops counting towards `n_parameters`, so
#'   this is a different model rather than a different starting point for the
#'   same one. `start` must carry the named blocks;
#'   `starting_values(fit, what = "measurement")` produces them, and
#'   [fit_staged()] wraps the whole two-stage workflow in one call. The mixing
#'   parameters cannot be held: they are what a fixed-measurement fit is for.
#' @param min_probability Positive lower bound on every categorical response
#'   probability, defining a constrained maximum-likelihood problem in the same
#'   way `min_variance` does for Gaussian indicators.
#' @param missing `"error"` rejects missing indicators; `"fiml"` maximizes the
#'   observed-data likelihood under an ignorable missingness mechanism (MAR).
#'   Missing indicators are integrated out, not filled in for likelihood fitting.
#' @param covariance_model `"diagonal"` assumes conditional independence;
#'   `"full"` estimates within-profile residual covariances.
#' @return An `multilpa` object containing `means`, `variances`, optional
#'   `covariances` (indicators by indicators by profiles),
#'   `profile_probabilities`, `group_probabilities`, posterior matrices,
#'   classifications, log likelihood, information criteria, restart diagnostics,
#'   and convergence history. Individual rows retain their input order; groups
#'   retain first-occurrence order. `bic` and `bic_groups` use the observed group
#'   count; `bic_individual` uses `n_informative`, the number of rows carrying at
#'   least one observed indicator. These are alternative
#'   conventions, not interchangeable criteria. `boundary` identifies variance
#'   bounds; `small_classes` flags effective memberships below one. The fit
#'   itself carries no standard errors: [parameter_inference()] computes them
#'   from the fit and the data it was fitted to, and [bootstrap_lrt()] tests
#'   nested models. No guarantee of global optimality is
#'   given, whatever `n_starts` is used.
#'   Read the tidy form with `as.data.frame()`; `what` selects which table.
#' @details Infinite and constant observed indicators are rejected. Each
#'   indicator must have at least two distinct observed values. No rows are
#'   silently dropped. In FIML mode, fully missing individuals contribute no
#'   direct measurement likelihood but receive posterior probabilities from
#'   their group's information. Observation counts retain these individuals;
#'   the individual-level BIC excludes them through `n_informative`, while the
#'   group-level BIC uses all observed groups. Missing patterns can prevent parameter
#'   identification; no general identification guarantee is made. Initialization
#'   alone uses indicator-mean filling. EM uses conditional Gaussian sufficient
#'   statistics and observed marginal densities. More than one group class requires more than one profile
#'   and at least one group with multiple individuals; these checks are necessary
#'   but do not establish identification. The highest finite likelihood across
#'   starts is returned, even if that start did not converge; inspect `converged`
#'   and `starts`. Profile and group-class labels are arbitrary.
#' @references Vermunt, J. K. (2003). Multilevel latent class models.
#'   Sociological Methodology, 33, 213--239.
#'   doi:10.1111/j.0081-1750.2003.t01-1-00131.x.
#' @examples
#' set.seed(7)
#' # Two kinds of school, differing only in how often a pupil scores highly.
#' example_data <- data.frame(
#'   school = rep(seq_len(24), each = 10),
#'   school_type = rep(c("mixed", "high"), each = 120)
#' )
#' example_data$high <- rbinom(240, 1L,
#'   ifelse(example_data$school_type == "high", 0.8, 0.2))
#' example_data$score_a <- rnorm(240, mean = 2 * example_data$high)
#' example_data$score_b <- rnorm(240, mean = 2 * example_data$high)
#'
#' fit <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                 n_profiles = 2, n_group_classes = 2, n_starts = 4,
#'                 seed = 42)
#' summary(fit)
#' as.data.frame(fit)
#' get_data(fit, what = "profile_probabilities")
#' @export
#' @importFrom stats setNames
multilpa <- function(data, vars, id, n_profiles,
                       n_group_classes = 2L,
                       profile_covariates = character(),
                       group_covariates = character(),
                       variance_model = c("varying", "equal"),
                       n_starts = 10L, max_iter = 1000L, tol = 1e-8,
                       min_variance = 1e-6, seed = NULL, start = NULL,
                       missing = c("error", "fiml"),
                       covariance_model = c("diagonal", "full"),
                       categorical = character(), min_probability = 1e-10,
                       time = NULL, fixed = character(),
                       centering = c("none", "person", "grand"),
                       volume = NULL, shape = NULL, orientation = NULL) {
  ## Before `stopifnot()`, which reads `id` and would otherwise force the
  ## missing argument into R's own bare "argument \"id\" is missing" error.
  if (missing(id)) {
    stop(errorCondition(paste(
      "`id` names the column the observations are nested in, and this model",
      "needs it. Pass it, or pass `id = NULL` to say the observations are",
      "independent and fit a single-level model."),
      class = "multilpa_bad_argument", call = NULL))
  }
  stopifnot(
    "`data` must be a data frame" = is.data.frame(data),
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
  call <- match.call()
  single_level <- is.null(id)
  if (single_level) {
    resolved <- .multilpa_single_level(
      data, if (missing(n_group_classes)) NULL else n_group_classes)
    data <- resolved$data
    id <- resolved$id
    n_group_classes <- resolved$n_group_classes
  }
  variance_model <- match.arg(variance_model)
  covariance_model <- match.arg(covariance_model)
  missing <- match.arg(missing)
  centering <- match.arg(centering)
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
  # Counted across every M-step of every start, and reported once at the end.
  .multilpa_reset_structure_log()
  on.exit(.multilpa_report_structure_log(), add = TRUE)
  if (!structure %in% .multilpa_inferable_structures() &&
      any(c("variances", "measurement") %in% fixed)) {
    stop(errorCondition(paste(
      "A constrained covariance structure is maximized across every profile at",
      "once, so holding one profile's variances would not leave the others at",
      "their maximum. Use `fixed = \"means\"`, or the unconstrained structure."),
      class = "multilpa_bad_argument", call = NULL))
  }
  if (length(profile_covariates) > 0L || length(group_covariates) > 0L) {
    return(.multilpa_covariate_model(
      data = data, vars = vars, id = id, n_profiles = n_profiles,
      n_group_classes = n_group_classes,
      profile_covariates = profile_covariates,
      group_covariates = group_covariates, variance_model = variance_model,
      n_starts = n_starts, max_iter = max_iter, tol = tol,
      min_variance = min_variance, seed = seed, start = start,
      missing = missing, covariance_model = covariance_model,
      categorical = categorical, min_probability = min_probability,
      time = time, fixed = fixed, call = call))
  }
  # The data contract comes first: `time` is checked against the `id` column,
  # so an `id` that does not name a column of `data` must be reported as
  # that, not as an opaque failure inside the time check.
  .multilpa_check_arguments(data, vars, id, n_profiles, n_group_classes,
                          n_starts, max_iter, tol, min_variance, min_probability,
                          seed, categorical)
  time_values <- .multilpa_time_values(data, time, id, vars)
  measurement <- .multilpa_prepare_indicators(data, vars, categorical,
                                            missing, min_probability)
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
  # Centring happens here, between validating the indicators and fitting them,
  # so the measurement model sees deviations and every verb downstream can put
  # a supplied frame on the same scale from the offsets kept on the fit.
  centred <- .multilpa_center_indicators(x, group_index, n_groups, centering)
  if (!identical(centering, "none") && ncol(x) > 0L) {
    x <- centred$x
    .multilpa_check_centered(x, centering)
    indicator_frame <- as.data.frame(x)
  }
  group_sizes <- groups$sizes
  distinct_rows <- if (is.null(codes)) nrow(unique(x)) else
    nrow(unique(cbind(x, codes)))
  if (n_profiles > distinct_rows) {
    stop(errorCondition(
      "n_profiles cannot exceed the number of distinct observed indicator rows.",
      class = "multilpa_unidentified", call = NULL))
  }
  if (n_group_classes > n_groups) {
    stop(errorCondition("n_group_classes cannot exceed the number of groups.",
                        class = "multilpa_unidentified", call = NULL))
  }
  if (n_group_classes > 1L && (n_profiles == 1L || all(group_sizes == 1L))) {
    stop(errorCondition(
      "Multiple group classes are not identifiable with one profile or only singleton groups.",
      class = "multilpa_unidentified", call = NULL))
  }
  if (length(fixed) > 0L) {
    start <- .multilpa_complete_start(start, n_profiles, n_group_classes)
  }
  # `fixed` is checked against the start before the start's own shape is, so
  # that a request to hold a block the start does not carry is reported as
  # that, rather than as a generic complaint about the start's contents.
  fixed <- .multilpa_validate_fixed(fixed, start, covariance_model, ncol(x),
                                    n_categories)
  if (!is.null(start)) {
    start <- .multilpa_validate_start(start, n_profiles, n_group_classes,
                                    ncol(x), variance_model, min_variance,
                                    covariance_model, n_categories, min_probability,
                                    categorical, encoded$levels)
  }
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  # Center to reduce cancellation in weighted means when indicators have large
  # offsets. This is a translation, with no density Jacobian or scale change.
  centers <- if (ncol(x) > 0L) colMeans(x, na.rm = TRUE) else numeric(0)
  x <- sweep(x, 2L, centers, "-")
  if (any(!is.finite(x[!is.na(x)]^2))) {
    stop(errorCondition(
      "Indicator scales overflow squared residuals; rescale the data.",
      class = "multilpa_bad_data", call = NULL))
  }
  if (!is.null(start)) start$means <- sweep(start$means, 2L, centers, "-")
  held <- .multilpa_held_parameters(start, fixed, covariance_model)
  # `max_iter = 0` performs no update, so a restart cannot improve on anything:
  # scoring extra random initializations and keeping the highest would return a
  # parameter set nobody supplied in place of the one the caller asked to have
  # evaluated. Evaluate-only with a `start` therefore uses that start alone, and
  # `get_data(fit, "starts")` shows the single start that was run.
  if (max_iter == 0L && !is.null(start)) n_starts <- 1L
  attempts <- lapply(seq_len(n_starts), function(start_index) {
    tryCatch({
      initial <- if (start_index == 1L && !is.null(start)) start else {
        .multilpa_initialize(x, group_index, n_profiles, n_group_classes,
                           variance_model, min_variance, start_index,
                           covariance_model, codes, n_categories, min_probability)
      }
      # Every restart begins from the held values, so a random start cannot
      # report a measurement solution that was neither estimated nor supplied.
      initial <- .multilpa_apply_held(initial, held)
      initial <- .multilpa_project_start(initial, structure, nrow(x), min_variance)
      .multilpa_em(x, group_index, initial, variance_model, min_variance, max_iter,
                 tol, covariance_model, codes, n_categories, min_probability,
                 held, structure)
    }, error = function(error) list(error = conditionMessage(error)))
  })
  valid <- vapply(attempts, function(attempt) is.null(attempt$error), logical(1))
  if (!any(valid)) {
    stop(errorCondition(sprintf(
      "All %d starts failed: %s", n_starts,
      paste(unique(vapply(attempts, `[[`, character(1), "error")), collapse = "; ")),
      class = "multilpa_all_starts_failed", call = NULL))
  }
  scores <- vapply(attempts, function(attempt) {
    if (is.null(attempt$error)) attempt$expectation$log_likelihood else -Inf
  }, numeric(1))
  best_start <- which.max(scores)
  best <- attempts[[best_start]]
  parameters <- best$parameters
  parameters$means <- sweep(parameters$means, 2L, centers, "+")
  profile_names <- paste0("profile_", seq_len(n_profiles))
  type_names <- paste0("group_class_", seq_len(n_group_classes))
  parameters <- .multilpa_label_parameters(parameters, profile_names, continuous,
                                         categorical, encoded$levels)
  dimnames(parameters$profile_probabilities) <- list(type_names, profile_names)
  names(parameters$group_probabilities) <- type_names
  subject_posteriors <- best$expectation$subject_posteriors
  group_posteriors <- best$expectation$group_posteriors
  dimnames(subject_posteriors) <- list(rownames(data), profile_names)
  dimnames(group_posteriors) <- list(group_ids, type_names)
  n_parameters_unconstrained <- .multilpa_count_parameters(structure = structure,
    n_profiles, n_group_classes, ncol(x), n_categories, variance_model,
    covariance_model)
  n_parameters <- n_parameters_unconstrained -
    .multilpa_fixed_parameters(fixed, n_profiles, ncol(x), n_categories,
                               variance_model, covariance_model)
  ## A held fit is comparable with other held fits on `n_parameters` and with a
  ## joint fit on the count that includes the measurement it was handed. Only
  ## `fit_staged()` used to record the second one, so `print()` fell through
  ## `%||%` on a `multilpa(fixed = )` fit and reported the same number twice.
  n_parameters_with_measurement <- if (length(fixed) == 0L) NULL else
    n_parameters_unconstrained
  log_likelihood <- best$expectation$log_likelihood
  starts <- do.call(rbind, lapply(seq_along(attempts), function(start_index) {
    attempt <- attempts[[start_index]]
    data.frame(start = start_index, log_likelihood = scores[start_index],
               converged = if (valid[start_index]) attempt$converged else FALSE,
               iterations = if (valid[start_index]) attempt$iterations else NA_integer_,
               error = if (valid[start_index]) NA_character_ else attempt$error,
               boundary = if (valid[start_index])
                 .multilpa_covariance_boundary(attempt$parameters, min_variance) else NA)
  }))
  # A row with no observed indicator contributes nothing to the likelihood and
  # nothing to any estimate, so it must not enlarge the sample size the
  # individual-level BIC penalizes against.
  observed_per_row <- rowSums(!is.na(x)) +
    if (is.null(codes)) 0L else rowSums(!is.na(codes))
  # Reported by print() and summary() rather than warned about: resampling verbs
  # refit the same data many times, and a warning would fire once per replicate.
  n_informative <- sum(observed_per_row > 0L)
  boundary <- .multilpa_covariance_boundary(parameters, min_variance)
  small_classes <- any(colSums(subject_posteriors) < 1) || any(colSums(group_posteriors) < 1)
  result <- c(parameters, list(
    call = call, vars = vars, continuous = continuous,
    categorical = categorical,
    categorical_levels = encoded$levels, min_probability = min_probability,
    indicator_data = as.matrix(indicator_frame),
    centering = centering, centering_offsets = centred$offsets,
    categorical_data = codes,
    id = id, group_ids = group_ids, single_level = single_level,
    time = time, time_values = time_values,
    group_values = group_values, group_index = group_index,
    group_sizes = setNames(group_sizes, group_ids),
    n_observations = nrow(x), n_informative = n_informative,
    n_groups = n_groups, n_profiles = as.integer(n_profiles),
    n_group_classes = as.integer(n_group_classes), variance_model = variance_model,
    covariance_model = covariance_model, covariance_structure = structure,
    missing = missing,
    n_observed_by_indicator = setNames(
      c(colSums(!is.na(x)), if (is.null(codes)) NULL else colSums(!is.na(codes))),
      c(continuous, categorical)),
    min_variance = min_variance, standard_deviations = sqrt(parameters$variances),
    fixed = fixed,
    measurement_model = if (is.null(codes)) "gaussian" else
      if (ncol(x) == 0L) "categorical" else "mixed",
    subject_posteriors = subject_posteriors, group_posteriors = group_posteriors,
    subject_profiles = max.col(subject_posteriors, ties.method = "first"),
    group_classes = max.col(group_posteriors, ties.method = "first"),
    log_likelihood = log_likelihood,
    group_log_likelihood = setNames(best$expectation$group_log_likelihood, group_ids),
    n_parameters = n_parameters,
    n_parameters_with_measurement = n_parameters_with_measurement,
    aic = -2 * log_likelihood + 2 * n_parameters,
    bic = -2 * log_likelihood + log(n_groups) * n_parameters,
    bic_groups = -2 * log_likelihood + log(n_groups) * n_parameters,
    bic_individual = -2 * log_likelihood + log(n_informative) * n_parameters,
    converged = best$converged, iterations = best$iterations,
    log_likelihood_history = best$history, starts = starts, best_start = best_start,
    n_failed_starts = sum(!valid), boundary = boundary, small_classes = small_classes,
    effective_profile_counts = colSums(subject_posteriors),
    effective_group_counts = colSums(group_posteriors),
    n_best_replicated = sum(valid & abs(scores - log_likelihood) <=
                             1e-6 * (1 + abs(log_likelihood))),
    replication_tolerance = 1e-6 * (1 + abs(log_likelihood))))
  class(result) <- "multilpa"
  if (any(!valid)) {
    warning(warningCondition(sprintf(
      "%d of %d starts failed; see get_data(fit, \"starts\").",
      sum(!valid), n_starts), class = "multilpa_failed_starts", call = NULL))
  }
  # max_iter = 0 is a deliberate evaluate-only call, so non-convergence is
  # expected rather than an anomaly worth reporting.
  if (!best$converged && max_iter > 0L) {
    warning(warningCondition(
      "The best start did not converge; increase max_iter and see get_data(fit, \"starts\").",
      class = "multilpa_unconverged", call = NULL))
  }
  if (boundary) {
    warning(warningCondition(if (covariance_model == "full")
      "A covariance eigenvalue reached min_variance; this is a bound-active constrained fit." else
      "A variance reached min_variance; this is a bound-active constrained fit.",
      class = "multilpa_boundary", call = NULL))
  }
  if (small_classes) {
    warning(warningCondition(
      "A profile or group class has effective membership below one.",
      class = "multilpa_small_classes", call = NULL))
  }
  result
}

#' Validate the scalar arguments of a fit
#' @return `NULL`, invisibly; raises on the first broken contract.
#' @noRd
.multilpa_check_arguments <- function(data, vars, id, n_profiles,
                                    n_group_classes, n_starts, max_iter, tol,
                                    min_variance, min_probability, seed,
                                    categorical) {
  if (nrow(data) < 2L || length(vars) < 1L || anyNA(vars) ||
      anyDuplicated(vars) || anyDuplicated(names(data)) ||
      length(id) != 1L || is.na(id) || id %in% vars ||
      !all(c(vars, id) %in% names(data))) {
    stop(errorCondition(
      "Supply at least two rows, unique existing indicators, and one distinct group column.",
      class = "multilpa_bad_data", call = NULL))
  }
  counts <- list(n_profiles = n_profiles, n_group_classes = n_group_classes,
                 n_starts = n_starts)
  invisible(lapply(names(counts), function(field) {
    value <- counts[[field]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value != floor(value) || value > .Machine$integer.max) {
      stop(errorCondition(sprintf("%s must be a positive integer.", field),
                          class = "multilpa_bad_argument", call = NULL))
    }
  }))
  # max_iter = 0 evaluates the likelihood at the supplied start without moving.
  if (!is.numeric(max_iter) || length(max_iter) != 1L || !is.finite(max_iter) ||
      max_iter < 0 || max_iter != floor(max_iter) ||
      max_iter > .Machine$integer.max) {
    stop(errorCondition("max_iter must be a nonnegative integer.",
                        class = "multilpa_bad_argument", call = NULL))
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0 ||
      !is.numeric(min_variance) || length(min_variance) != 1L ||
      !is.finite(min_variance) || min_variance <= 0) {
    stop(errorCondition("tol and min_variance must be finite positive numbers.",
                        class = "multilpa_bad_argument", call = NULL))
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
      !is.finite(seed) || seed < 0 || seed != floor(seed) ||
      seed > .Machine$integer.max)) {
    stop(errorCondition("seed must be a nonnegative integer or NULL.",
                        class = "multilpa_bad_argument", call = NULL))
  }
  if (anyDuplicated(categorical) || !all(categorical %in% vars)) {
    stop(errorCondition("`categorical` must name distinct indicators listed in `vars`.",
                        class = "multilpa_bad_categorical", call = NULL))
  }
  invisible(NULL)
}

#' Split indicators into Gaussian and categorical blocks
#'
#' Returns the continuous indicator matrix, which may have zero columns, and the
#' encoded categorical codes, which may be `NULL`. Every measurement-side
#' contract is checked here so the fitting path can assume clean input.
#'
#' @return A list with `continuous`, `frame`, `x`, `encoded`, `codes` and
#'   `n_categories`.
#' @noRd
.multilpa_prepare_indicators <- function(data, vars, categorical, missing,
                                       min_probability) {
  continuous <- setdiff(vars, categorical)
  encoded <- if (length(categorical) > 0L) {
    .multilpa_encode_categorical(data[, categorical, drop = FALSE])
  } else NULL
  codes <- encoded$codes
  n_categories <- encoded$n_categories
  if (!is.null(codes) && min_probability * max(n_categories) >= 1) {
    stop(errorCondition(sprintf(
      "`min_probability` of %g leaves no room for an indicator with %d categories.",
      min_probability, max(n_categories)),
      class = "multilpa_bad_categorical", call = NULL))
  }
  if (!is.null(codes) && missing == "error" && anyNA(codes)) {
    stop(errorCondition("Indicators contain missing or non-finite values.",
                        class = "multilpa_bad_data", call = NULL))
  }
  frame <- data[, continuous, drop = FALSE]
  if (!all(vapply(frame, is.numeric, logical(1))) ||
      any(vapply(frame, function(value) !is.null(dim(value)), logical(1)))) {
    stop(errorCondition(
      "Every indicator must be a numeric vector; factors are not continuous indicators.",
      class = "multilpa_bad_data", call = NULL))
  }
  x <- matrix(as.matrix(frame), nrow = nrow(data), ncol = length(continuous),
              dimnames = list(NULL, continuous))
  if (any(is.infinite(x)) || any(is.nan(x)) || (missing == "error" && anyNA(x))) {
    stop(errorCondition("Indicators contain missing or non-finite values.",
                        class = "multilpa_bad_data", call = NULL))
  }
  if (any(colSums(!is.na(x)) == 0L)) {
    stop(errorCondition("Every indicator must have observed values.",
                        class = "multilpa_bad_data", call = NULL))
  }
  if (any(vapply(frame, function(value) length(unique(value[!is.na(value)])) < 2L,
                 logical(1)))) {
    stop(errorCondition("Constant indicators cannot identify Gaussian profiles.",
                        class = "multilpa_bad_data", call = NULL))
  }
  list(continuous = continuous, frame = frame, x = x, encoded = encoded,
       codes = codes, n_categories = n_categories)
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
      class = "multilpa_bad_data", call = NULL))
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

#' Hand a covariate request to the covariate estimator
#'
#' `multilpa()` is one verb over one model, and membership covariates are part
#' of that model rather than a different one, so they are arguments and not a
#' separate entry point. Three of `multilpa()`'s own arguments have no meaning
#' on this path: the covariate likelihood has no observed-data form, no
#' starting-value contract of the shape the covariate-free EM uses, and no
#' held-measurement machinery. Each is refused by name instead of being
#' accepted and ignored, which would return a fit that is not the one asked
#' for.
#'
#' @param data,vars,id,n_profiles,n_group_classes As in [multilpa()].
#' @param profile_covariates,group_covariates The requested predictors.
#' @param variance_model,n_starts,max_iter,tol,min_variance As in [multilpa()].
#' @param seed,time,covariance_model,categorical,min_probability As in
#'   [multilpa()].
#' @param start,missing,fixed Refused when they are anything but their default.
#' @param call The user's call, recorded on the result.
#' @return An object of class `multilpa_covariates`.
#' @noRd
.multilpa_covariate_model <- function(data, vars, id, n_profiles,
                                      n_group_classes, profile_covariates,
                                      group_covariates, variance_model,
                                      n_starts, max_iter, tol, min_variance,
                                      seed, start, missing, covariance_model,
                                      categorical, min_probability, time,
                                      fixed, call) {
  unsupported <- c(
    start = !is.null(start),
    missing = !identical(missing, "error"),
    fixed = length(fixed) > 0L)
  if (any(unsupported)) {
    stop(errorCondition(sprintf(paste(
      "%s cannot be combined with `profile_covariates` or",
      "`group_covariates`. Fit the covariate model without %s, or fit the",
      "covariate-free model and use three_step() or r3step()."),
      paste(sprintf("`%s`", names(unsupported)[unsupported]), collapse = ", "),
      if (sum(unsupported) > 1L) "them" else "it"),
      class = "multilpa_bad_argument", call = NULL))
  }
  .multilpa_fit_covariates(
    data = data, vars = vars, id = id, n_profiles = n_profiles,
    n_group_classes = n_group_classes,
    profile_covariates = profile_covariates,
    group_covariates = group_covariates, variance_model = variance_model,
    n_starts = n_starts, max_iter = max_iter, tol = tol,
    min_variance = min_variance, seed = seed, time = time,
    covariance_model = covariance_model, categorical = categorical,
    min_probability = min_probability, call = call)
}
