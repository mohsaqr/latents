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

#' Every permutation of `seq_len(k)`
#'
#' @param k A single positive integer.
#' @return A matrix with `factorial(k)` rows, one permutation per row.
#' @noRd
.multilpa_permutations <- function(k) {
  stopifnot("`k` must be a single positive integer" =
              length(k) == 1L && is.finite(k) && k >= 1 && k == floor(k))
  k <- as.integer(k)
  if (k == 1L) return(matrix(1L, 1L, 1L))
  smaller <- .multilpa_permutations(k - 1L)
  orders <- do.call(rbind, lapply(seq_len(k), function(position) {
    remaining <- seq_len(k)[-position]
    cbind(position, matrix(remaining[smaller], nrow(smaller)))
  }))
  ## `cbind()` names the first column after the scalar it was given, and a
  ## dimname here would travel into every row `apply()` hands out.
  dimnames(orders) <- NULL
  orders
}

#' What a profile looks like, for matching one fit's profiles to another's
#'
#' Bootstrap replicates come back with their profiles in whatever order the EM
#' happened to label them, so the replicates have to be matched to the original
#' before their spread means anything. The signature is what the matching reads:
#' the profile means, each categorical indicator's response probabilities, and
#' the mixing probability, side by side.
#'
#' Every column is put in units where a difference of one is a large difference.
#' Means are divided by the indicator's own within-profile standard deviation,
#' taken from `reference` so that the two signatures are measured the same way;
#' probabilities are already on that scale and are left alone. Standardizing a
#' column by its spread *across profiles* instead is what an earlier version did,
#' and it is wrong: two profiles that happen to be mixed half and half give that
#' column a spread near zero, and dividing by it turns a difference of 0.02 in a
#' mixing proportion into the term that decides the match.
#'
#' @param x A fitted `multilpa` model.
#' @param reference The fit whose scale both signatures are measured in.
#' @return A profiles-by-features numeric matrix.
#' @noRd
.multilpa_profile_signature <- function(x, reference = x) {
  stopifnot(inherits(x, "multilpa"), inherits(reference, "multilpa"))
  blocks <- list()
  if (length(.multilpa_continuous_names(x)) > 0L) {
    spread <- sqrt(colMeans(reference$variances))
    spread[!is.finite(spread) | spread <= 0] <- 1
    blocks$means <- sweep(unname(x$means), 2L, spread, "/")
  }
  responses <- x$response_probabilities %||% list()
  if (length(responses) > 0L) {
    blocks$response <- do.call(cbind, lapply(responses, unname))
  }
  ## A weak identifier on its own, and a tie-breaker next to the measurement:
  ## kept on its natural scale so it can never outvote the means.
  blocks$mixing <- matrix(colMeans(x$profile_probabilities), ncol = 1L)
  do.call(cbind, blocks)
}

#' The order that best matches one set of profiles to another
#'
#' Exhaustive for the profile counts latent profile analysis actually uses; a
#' greedy nearest match above that, because `factorial(k)` stops being free.
#' Both signatures arrive already in comparable units, so nothing is rescaled
#' here; see [.multilpa_profile_signature()] for why that scaling is its job.
#'
#' @param reference,candidate Profiles-by-features matrices with the same shape.
#' @return An integer vector `order` with `candidate[order, ]` matched to
#'   `reference`.
#' @noRd
.multilpa_match_order <- function(reference, candidate) {
  stopifnot("the two signatures must have the same shape" =
              is.matrix(reference) && is.matrix(candidate) &&
              identical(dim(reference), dim(candidate)))
  k <- nrow(reference)
  if (k == 1L) return(1L)
  cost <- vapply(seq_len(k), function(j) {
    rowSums(sweep(reference, 2L, candidate[j, ], "-")^2)
  }, numeric(k))
  if (k <= 7L) {
    orders <- .multilpa_permutations(k)
    total <- vapply(seq_len(nrow(orders)), function(i) {
      sum(cost[cbind(seq_len(k), orders[i, ])])
    }, numeric(1))
    return(as.integer(orders[which.min(total), ]))
  }
  # Greedy: give each reference profile its closest unclaimed candidate, taking
  # the most decisive pairing first so an obvious match is not spent elsewhere.
  Reduce(function(state, i) {
    choice <- which(state$available)[which.min(cost[i, state$available])]
    state$order[i] <- choice
    state$available[choice] <- FALSE
    state
  }, order(apply(cost, 1L, min)),
  init = list(order = integer(k), available = rep(TRUE, k)))$order
}

#' Reorder a fit's profile-indexed blocks
#' @param x A fitted `multilpa` model.
#' @param order An integer permutation of the profiles.
#' @return The fit, with every profile-indexed block in the new order.
#' @noRd
.multilpa_permute_profiles <- function(x, order) {
  stopifnot(inherits(x, "multilpa"),
            "`order` must be a permutation of the profiles" =
              identical(sort(as.integer(order)), seq_len(x$n_profiles)))
  x$means <- x$means[order, , drop = FALSE]
  x$variances <- x$variances[order, , drop = FALSE]
  if (!is.null(x$covariances)) x$covariances <- x$covariances[, , order, drop = FALSE]
  x$profile_probabilities <- x$profile_probabilities[, order, drop = FALSE]
  responses <- x$response_probabilities %||% list()
  if (length(responses) > 0L) {
    x$response_probabilities <- lapply(responses, function(block) {
      block[order, , drop = FALSE]
    })
  }
  x
}

#' Reorder a fit's group-class-indexed blocks
#' @param x A fitted `multilpa` model.
#' @param order An integer permutation of the group classes.
#' @return The fit, with every group-class-indexed block in the new order.
#' @noRd
.multilpa_permute_group_classes <- function(x, order) {
  stopifnot(inherits(x, "multilpa"),
            "`order` must be a permutation of the group classes" =
              identical(sort(as.integer(order)), seq_len(x$n_group_classes)))
  x$profile_probabilities <- x$profile_probabilities[order, , drop = FALSE]
  x$group_probabilities <- x$group_probabilities[order]
  x
}

#' Put a replicate's labels back on the original's
#'
#' Profiles first, because a group class is described by the profile mixture it
#' carries and that description only means something once the profiles agree.
#'
#' @param replicate A fitted `multilpa` model from a resample.
#' @param reference The original fit.
#' @return The replicate, relabelled to the reference.
#' @noRd
.multilpa_align_labels <- function(replicate, reference) {
  replicate <- .multilpa_permute_profiles(
    replicate,
    .multilpa_match_order(.multilpa_profile_signature(reference, reference),
                          .multilpa_profile_signature(replicate, reference)))
  if (reference$n_group_classes == 1L) return(replicate)
  signature <- function(x) cbind(x$profile_probabilities, x$group_probabilities)
  .multilpa_permute_group_classes(
    replicate, .multilpa_match_order(signature(reference), signature(replicate)))
}

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
                 centering = x$centering %||% "none")
  structure <- x$covariance_structure
  if (is.null(structure) || is.na(structure) || !structure %in% .multilpa_structures()) {
    return(c(shared, list(variance_model = x$variance_model,
                          covariance_model = x$covariance_model %||% "diagonal")))
  }
  c(shared, .multilpa_structure_arguments(structure))
}

#' One resample of the groups, as a data frame
#'
#' A group drawn twice has to become two groups, or the refit would pool the two
#' copies into one unit and the resample would carry fewer groups than it drew.
#'
#' @param data The fitting data.
#' @param rows_by_group Row indices of each group, in group order.
#' @param id The group column's name.
#' @param drawn Which groups were drawn, with replacement.
#' @return A data frame with one relabelled group per draw.
#' @noRd
.multilpa_resample_groups <- function(data, rows_by_group, id, drawn) {
  stopifnot(is.data.frame(data), is.list(rows_by_group), is.character(id),
            length(id) == 1L)
  taken <- rows_by_group[drawn]
  resampled <- data[unlist(taken, use.names = FALSE), , drop = FALSE]
  resampled[[id]] <- rep(seq_along(drawn), lengths(taken))
  row.names(resampled) <- NULL
  resampled
}

#' Bootstrap replicates of the natural coefficients
#'
#' @param x A fitted `multilpa` model.
#' @param data The fitting data.
#' @param iter How many resamples to draw.
#' @param n_starts,max_iter,tol Passed to each refit.
#' @return A list with `estimates`, an `iter` by parameters matrix whose failed
#'   rows are `NA`, and `messages`, the reason each failure gave.
#' @noRd
.multilpa_bootstrap_replicates <- function(x, data, iter, n_starts, max_iter, tol) {
  rows_by_group <- split(seq_len(nrow(data)), x$group_index)
  arguments <- .multilpa_refit_arguments(x)
  reference_names <- names(.multilpa_coefficients(x, "natural"))
  replicates <- lapply(seq_len(iter), function(i) {
    drawn <- sample.int(x$n_groups, x$n_groups, replace = TRUE)
    resampled <- .multilpa_resample_groups(data, rows_by_group, x$id, drawn)
    fit <- tryCatch(
      do.call(multilpa, c(list(data = resampled), arguments,
                          list(n_starts = n_starts, max_iter = max_iter, tol = tol))),
      error = function(error) conditionMessage(error))
    if (is.character(fit)) {
      return(list(estimate = rep(NA_real_, length(reference_names)), message = fit))
    }
    if (!isTRUE(fit$converged)) {
      return(list(estimate = rep(NA_real_, length(reference_names)),
                  message = "the replicate did not converge"))
    }
    aligned <- .multilpa_align_labels(fit, x)
    list(estimate = unname(.multilpa_coefficients(aligned, "natural")),
         message = NA_character_)
  })
  estimates <- do.call(rbind, lapply(replicates, `[[`, "estimate"))
  colnames(estimates) <- reference_names
  list(estimates = estimates,
       messages = vapply(replicates, `[[`, character(1), "message"))
}

#' Percentile bootstrap inference for a fitted model
#'
#' Called by [parameter_inference()] when `method = "bootstrap"`. The interval
#' is the percentile interval of the replicates and the standard error is their
#' standard deviation, so both are read off the same replicates and neither
#' assumes the estimate is normal around its truth.
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
  if (length(x$fixed %||% character()) > 0L) {
    stop(errorCondition(paste(
      "A fit that holds a block fixed cannot be bootstrapped here: the held",
      "values came from another fit, and resampling these data does not",
      "resample them. Bootstrap the fit the measurement came from."),
      class = "multilpa_unsupported_inference", call = NULL))
  }
  if (!isTRUE(x$converged)) {
    stop(errorCondition(
      "The original fit did not converge; bootstrap it only once it has.",
      class = "multilpa_no_converge", call = NULL))
  }
  free_natural <- .multilpa_free_index(x, "natural")
  drawn <- .multilpa_bootstrap_replicates(x, data, iter, n_starts, max_iter, tol)
  estimates <- drawn$estimates[, free_natural, drop = FALSE]
  valid <- stats::complete.cases(estimates)
  if (sum(valid) < 2L) {
    reasons <- drawn$messages[!is.na(drawn$messages)]
    stop(errorCondition(sprintf(
      "Only %d of %d resamples produced a usable fit; the first reason was: %s",
      sum(valid), iter,
      if (length(reasons) > 0L) reasons[1L] else "no reason was recorded"),
      class = "multilpa_bootstrap_failed", call = NULL))
  }
  kept <- estimates[valid, , drop = FALSE]
  if (sum(valid) < iter) {
    warning(warningCondition(sprintf(
      "%d of %d resamples did not produce a usable fit and were dropped.",
      iter - sum(valid), iter), class = "multilpa_bootstrap_dropped", call = NULL))
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
    fixed = x$fixed %||% character()))
  result
}
