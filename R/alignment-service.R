# The alignment service: matching one fit's class labels to another's, and
# relabelling a fit exactly.
#
# A mixture's labels are arbitrary. Bootstrap replicates, imputations and
# refits from other seeds come back with their profiles and group classes in
# whatever order their EM happened to settle, so every verb that compares or
# averages fits has to put them on one labelling first. That is three steps,
# and each lives here once:
#
#   signature   what a class looks like, on a scale-free footing
#               (`.latents_measurement_signature()`, plus the group-class
#               signatures);
#   matching    the permutation that best maps one signature onto another
#               (`.multilpa_match_order()`), exhaustive up to seven classes;
#   permuting   relabelling every class-indexed block so the likelihood is
#               unchanged, which for baseline-category logits means rebasing
#               on the new reference class (`.latents_rebase_logits()`).
#
# The engines keep what is particular to their parameters: the transition
# models relabel their own logits (lta-bootstrap.R) and the regression and
# growth mixtures their own blocks (mixture-regression-starts.R,
# growth-mixture.R); they call the primitives here.

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

#' The order that best matches one set of profiles to another
#'
#' Exhaustive for the profile counts latent profile analysis actually uses; a
#' greedy nearest match above that, because `factorial(k)` stops being free.
#' Both signatures arrive already in comparable units, so nothing is rescaled
#' here; see `.latents_measurement_signature()` for why that scaling is its job.
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

#' What a class looks like, for matching one fit's classes to another's
#'
#' Replicates, imputations and refits come back with their classes in whatever
#' order the EM happened to label them, so they have to be matched to a
#' reference before they can be compared. The signature is what the matching
#' reads: per measurement block, the class means, each categorical indicator's
#' response probabilities, the ordinal category probabilities and the count
#' means, side by side, and optionally a mixing column.
#'
#' Every column is put in units where a difference of one is a large
#' difference. Means are divided by the indicator's own within-class standard
#' deviation, taken from `scales` so that the two signatures are measured the
#' same way; probabilities are already on that scale and are left alone; count
#' means are divided by their standard deviation under the reference's count
#' model. Standardizing a column by its spread *across classes* instead is
#' what an earlier version did, and it is wrong: two classes that happen to be
#' mixed half and half give that column a spread near zero, and dividing by it
#' turns a difference of 0.02 in a mixing proportion into the term that decides
#' the match.
#'
#' The profile, covariate and transition models measure their classes with the
#' same blocks and differed only in the mixing column, so they share this one
#' builder: the constant mixing probabilities of a profile fit, the posterior
#' share of a covariate fit (which has no constant mixing), and none for a
#' transition fit (whose class shares vary by occasion).
#'
#' @param blocks A list of measurement blocks, each a list with any of
#'   `means`, `response_probabilities`, `ordinal_intercepts`,
#'   `ordinal_locations`, `count_means` (a fitted profile model is such a
#'   block).
#' @param scales The reference's blocks, in the same order, carrying the
#'   `variances`, `count_means`, `count_dispersion` and `n_profiles` the
#'   columns are scaled by.
#' @param mixing `NULL`, or one number per class appended as the last column.
#'   A weak identifier on its own and a tie-breaker next to the measurement, so
#'   it is kept on its natural scale where it can never outvote the means.
#' @return A classes-by-features numeric matrix, or `NULL` when there is
#'   neither measurement nor mixing to match on.
#' @noRd
.latents_measurement_signature <- function(blocks, scales, mixing = NULL) {
  stopifnot("`blocks` and `scales` must be lists of the same length" =
              is.list(blocks) && is.list(scales) && length(blocks) == length(scales))
  features <- unlist(lapply(seq_along(blocks), function(b) {
    block <- blocks[[b]]
    scale <- scales[[b]]
    parts <- list()
    if (length(block$means) > 0L && ncol(block$means) > 0L) {
      spread <- sqrt(colMeans(scale$variances))
      spread[!is.finite(spread) | spread <= 0] <- 1
      parts$means <- sweep(unname(block$means), 2L, spread, "/")
    }
    if (length(block$response_probabilities) > 0L) {
      parts$response <- do.call(cbind, lapply(block$response_probabilities, unname))
    }
    c(parts, .multilpa_extra_signature(block, scale))
  }), recursive = FALSE)
  if (!is.null(mixing)) features$mixing <- matrix(mixing, ncol = 1L)
  do.call(cbind, features)
}

#' Ordinal probabilities and standardized count means for label matching
#' @param x,reference Fitted models with the same extra indicators.
#' @return A list of profiles-by-features matrices.
#' @noRd
.multilpa_extra_signature <- function(x, reference) {
  blocks <- list()
  if (length(x$ordinal_intercepts) > 0L) {
    blocks$ordinal <- do.call(cbind, lapply(seq_along(x$ordinal_intercepts), function(j) {
      exp(.latents_ordinal_log_probabilities(x$ordinal_intercepts[[j]],
                                            x$ordinal_locations[, j]))
    }))
  }
  if (length(x$count_means) > 0L) {
    variance <- reference$count_means
    if (length(reference$count_dispersion) > 0L) {
      dispersion <- reference$count_dispersion
      if (nrow(dispersion) == 1L) {
        dispersion <- dispersion[rep(1L, reference$n_profiles), , drop = FALSE]
      }
      variance <- variance + dispersion * reference$count_means^2
    }
    spread <- sqrt(colMeans(variance))
    spread[!is.finite(spread) | spread <= 0] <- 1
    blocks$count <- sweep(unname(x$count_means), 2L, spread, "/")
  }
  blocks
}


#' What a profile model's profile looks like, for matching
#'
#' The measurement signature of a profile model, with its mixing
#' probabilities as the tie-breaking column.
#'
#' @param x A fitted `multilpa` model.
#' @param reference The fit whose scale both signatures are measured in.
#' @return A profiles-by-features numeric matrix.
#' @noRd
.multilpa_profile_signature <- function(x, reference = x) {
  stopifnot(inherits(x, "multilpa"), inherits(reference, "multilpa"))
  .latents_measurement_signature(list(x), list(reference),
                                 mixing = colMeans(x$profile_probabilities))
}

#' What a covariate fit's profile looks like, for matching
#'
#' As `.multilpa_profile_signature()`, with the profile's overall share taken
#' from the posteriors, since a covariate fit has no constant mixing
#' probabilities.
#' @param x,reference `multilpa_covariates` fits.
#' @return A profiles-by-features matrix.
#' @noRd
.multilpa_cov_profile_signature <- function(x, reference = x) {
  .latents_measurement_signature(list(x), list(reference),
                                 mixing = colMeans(x$subject_posteriors))
}

#' What a profile model's group class looks like, for matching
#'
#' The profile mixture the group class carries, next to its share of the
#' groups. Read after the profiles are aligned, since the mixture is only
#' comparable once they agree.
#' @param x A fitted `multilpa` model.
#' @return A group-classes-by-features matrix.
#' @noRd
.multilpa_group_signature <- function(x) {
  cbind(x$profile_probabilities, x$group_probabilities)
}

#' What a covariate fit's group class looks like, for matching
#'
#' The profile mixture a group class implies, averaged over the observed
#' covariates, next to the class's share of the groups. Read after the profiles
#' are aligned, since the mixture is only comparable once they agree.
#' @param x A `multilpa_covariates` fit.
#' @return A group-classes-by-features matrix.
#' @noRd
.multilpa_cov_group_signature <- function(x) {
  mixtures <- t(vapply(x$profile_design, function(design) {
    colMeans(.multilpa_softmax(design, x$profile_coefficients))
  }, numeric(x$n_profiles)))
  cbind(matrix(mixtures, x$n_group_classes), colMeans(x$group_posteriors))
}

#' Relabel baseline-category logits, rebasing on the new reference class
#'
#' Membership and initial-state logits are contrasts against a reference
#' class, the last (profile and covariate models, transitions) or the first
#' (regression and growth mixtures). Relabelling can move a different class
#' into the reference position, so reordering the columns is not enough: the
#' full logits (the reference's column being zero) are reordered and the new
#' reference's column is subtracted from every column, which leaves every
#' class probability unchanged. Each engine had its own copy of this; the
#' arithmetic here is theirs, one subtraction per entry, so the relabelled
#' values are bit-identical to what they computed.
#'
#' @param full A matrix of full logits, one column per class, the reference's
#'   column included.
#' @param order Integer permutation: new class `j` is old class `order[j]`.
#' @param reference `"last"` or `"first"`: which position is the reference
#'   after relabelling.
#' @return `full` with its columns reordered and rebased, the reference
#'   column zero.
#' @noRd
.latents_rebase_logits <- function(full, order, reference = c("last", "first")) {
  reference <- match.arg(reference)
  stopifnot("`full` must be a matrix" = is.matrix(full),
            "`order` must be a permutation of the columns" =
              length(order) == ncol(full) && setequal(order, seq_len(ncol(full))))
  moved <- full[, order, drop = FALSE]
  moved - moved[, if (identical(reference, "last")) ncol(moved) else 1L]
}

#' Relabel extra-indicator parameters and restore the ordinal reference
#' @param x A fitted profile or covariate model.
#' @param order An integer permutation of its profiles.
#' @return The fit with extra parameters in the new chart.
#' @noRd
.multilpa_permute_extra <- function(x, order) {
  labels <- paste0("profile_", seq_len(x$n_profiles))
  if (length(x$ordinal_intercepts) > 0L) {
    locations <- x$ordinal_locations[order, , drop = FALSE]
    reference <- locations[nrow(locations), ]
    x$ordinal_intercepts <- lapply(seq_along(x$ordinal_intercepts), function(j) {
      intercepts <- x$ordinal_intercepts[[j]]
      intercepts + seq_along(intercepts) * reference[j]
    })
    names(x$ordinal_intercepts) <- x$ordinal
    x$ordinal_locations <- sweep(locations, 2L, reference, "-")
    rownames(x$ordinal_locations) <- labels
  }
  if (length(x$count_means) > 0L) {
    x$count_means <- x$count_means[order, , drop = FALSE]
    rownames(x$count_means) <- labels
  }
  if (length(x$count_dispersion) > 0L && nrow(x$count_dispersion) == x$n_profiles) {
    x$count_dispersion <- x$count_dispersion[order, , drop = FALSE]
    rownames(x$count_dispersion) <- labels
  }
  x
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
  x <- .multilpa_permute_extra(x, order)
  labels <- paste0("profile_", seq_len(x$n_profiles))
  x$means <- x$means[order, , drop = FALSE]
  rownames(x$means) <- labels
  x$variances <- x$variances[order, , drop = FALSE]
  rownames(x$variances) <- labels
  if (!is.null(x$standard_deviations)) {
    x$standard_deviations <- x$standard_deviations[order, , drop = FALSE]
    rownames(x$standard_deviations) <- labels
  }
  if (!is.null(x$covariances)) {
    x$covariances <- x$covariances[, , order, drop = FALSE]
    dimnames(x$covariances)[[3L]] <- labels
  }
  x$profile_probabilities <- x$profile_probabilities[, order, drop = FALSE]
  colnames(x$profile_probabilities) <- labels
  responses <- x$response_probabilities %||% list()
  if (length(responses) > 0L) {
    x$response_probabilities <- lapply(responses, function(block) {
      reordered <- block[order, , drop = FALSE]
      rownames(reordered) <- labels
      reordered
    })
  }
  # Everything indexed by profile moves, not only the parameter blocks. The
  # posteriors are columns and reorder like the rest; `subject_profiles` holds
  # labels rather than positions, so it takes the inverse permutation -- a case
  # on old profile `order[j]` is on new profile `j`. Permuting the parameters
  # while leaving the assignments behind would leave the fit self-inconsistent,
  # which is invisible to a caller that reads only parameters and wrong for one
  # that compares assignments.
  if (!is.null(x$subject_posteriors)) {
    x$subject_posteriors <- x$subject_posteriors[, order, drop = FALSE]
    colnames(x$subject_posteriors) <- labels
  }
  if (!is.null(x$subject_profiles)) {
    # A noise unit is profile 0, which no permutation of the profiles moves.
    x$subject_profiles <- match(x$subject_profiles, order, nomatch = 0L)
  }
  if (!is.null(x$effective_profile_counts)) {
    x$effective_profile_counts <- x$effective_profile_counts[order]
    names(x$effective_profile_counts) <- labels
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
  labels <- paste0("group_class_", seq_len(x$n_group_classes))
  x$profile_probabilities <- x$profile_probabilities[order, , drop = FALSE]
  rownames(x$profile_probabilities) <- labels
  x$group_probabilities <- x$group_probabilities[order]
  names(x$group_probabilities) <- labels
  if (!is.null(x$group_posteriors)) {
    x$group_posteriors <- x$group_posteriors[, order, drop = FALSE]
    colnames(x$group_posteriors) <- labels
  }
  if (!is.null(x$group_classes)) {
    x$group_classes <- match(x$group_classes, order)
  }
  if (!is.null(x$effective_group_counts)) {
    x$effective_group_counts <- x$effective_group_counts[order]
    names(x$effective_group_counts) <- labels
  }
  x
}

#' Relabel a covariate fit's profiles and group classes
#'
#' New profile `j` is old profile `profile_order[j]`, and likewise for group
#' classes. Measurement blocks and posteriors are reordered. The membership
#' logits are reparameterized, not merely reordered: they are contrasts against
#' the last profile (or group class), and relabelling may move a different one
#' into that reference position. The full logits (the reference's being zero)
#' are reordered and the new reference's are subtracted. The profile design
#' needs no change: group class `h`'s design is the same matrix whichever
#' fitted class carries that label.
#'
#' @param x A `multilpa_covariates` fit.
#' @param profile_order,group_order Integer permutations.
#' @return The relabelled fit, whose likelihood is unchanged.
#' @noRd
.multilpa_cov_permute <- function(x, profile_order, group_order) {
  stopifnot(inherits(x, "multilpa_covariates"),
            identical(sort(as.integer(profile_order)), seq_len(x$n_profiles)),
            identical(sort(as.integer(group_order)), seq_len(x$n_group_classes)))
  x <- .multilpa_permute_extra(x, profile_order)
  n_profiles <- x$n_profiles
  n_groups <- x$n_group_classes
  profiles <- paste0("profile_", seq_len(n_profiles))
  reorder_rows <- function(block) {
    if (is.null(block)) return(NULL)
    names_kept <- colnames(block)
    reordered <- block[profile_order, , drop = FALSE]
    dimnames(reordered) <- list(if (!is.null(rownames(block))) profiles, names_kept)
    reordered
  }
  x$means <- reorder_rows(x$means)
  x$variances <- reorder_rows(x$variances)
  x$standard_deviations <- reorder_rows(x$standard_deviations)
  if (!is.null(x$covariances)) {
    x$covariances <- x$covariances[, , profile_order, drop = FALSE]
  }
  if (!is.null(x$response_probabilities)) {
    x$response_probabilities <- lapply(x$response_probabilities, reorder_rows)
  }
  ## Baseline-category logits: complete them with the reference's zeros,
  ## rebase on the new last category, and drop its (zero) column.
  recontrast <- function(coefficients, order) {
    full <- .latents_rebase_logits(cbind(coefficients, 0), order)
    result <- full[, -ncol(full), drop = FALSE]
    dimnames(result) <- dimnames(coefficients)
    result
  }
  beta <- recontrast(x$profile_coefficients, profile_order)
  ## Rows: one intercept per group class, then the slopes, which come once
  ## (shared) or in one block per group class.
  n_slopes <- nrow(beta) - n_groups
  slope_rows <- if (identical(x$profile_slopes %||% "shared", "group_class")) {
    width <- n_slopes / n_groups
    unlist(lapply(group_order, function(h) n_groups + (h - 1L) * width + seq_len(width)))
  } else n_groups + seq_len(n_slopes)
  kept_names <- dimnames(beta)
  beta <- beta[c(group_order, slope_rows), , drop = FALSE]
  dimnames(beta) <- kept_names
  x$profile_coefficients <- beta
  x$group_coefficients <- recontrast(x$group_coefficients, group_order)
  ## Posteriors and assignments follow the labels.
  x$subject_posteriors <- x$subject_posteriors[, profile_order, drop = FALSE]
  x$group_posteriors <- x$group_posteriors[, group_order, drop = FALSE]
  x$subject_profiles <- match(x$subject_profiles, profile_order)
  x$group_classes <- match(x$group_classes, group_order)
  x$effective_profile_counts <- x$effective_profile_counts[profile_order]
  x$effective_group_counts <- x$effective_group_counts[group_order]
  if (!is.null(x$group_priors)) {
    x$group_priors <- x$group_priors[, group_order, drop = FALSE]
  }
  if (!is.null(x$profile_priors)) {
    x$profile_priors <- lapply(x$profile_priors[group_order], function(prior) {
      prior[, profile_order, drop = FALSE]
    })
  }
  x
}


#' Put one fit's labels on a reference fit's
#'
#' Profiles first, because a group class is described by the profile mixture
#' it carries and that description only means something once the profiles
#' agree. A covariate fit is matched on its posterior shares and relabelled
#' with its logits rebased; any other profile model on its mixing
#' probabilities. Used by the cluster bootstrap, seed sensitivity and
#' multiple-imputation pooling, which all need the same alignment.
#'
#' @param fit A fit to relabel.
#' @param reference The fit whose labels it takes.
#' @return A list with the aligned `fit` and the `profile_order` and
#'   `group_order` that aligned it (new label `j` is old label `order[j]`).
#' @noRd
.multilpa_pool_align <- function(fit, reference) {
  if (inherits(fit, "multilpa_covariates")) {
    profile_order <- .multilpa_match_order(
      .multilpa_cov_profile_signature(reference, reference),
      .multilpa_cov_profile_signature(fit, reference))
    fit <- .multilpa_cov_permute(fit, profile_order, seq_len(fit$n_group_classes))
    group_order <- .multilpa_match_order(.multilpa_cov_group_signature(reference),
                                         .multilpa_cov_group_signature(fit))
    fit <- .multilpa_cov_permute(fit, seq_len(fit$n_profiles), group_order)
    return(list(fit = fit, profile_order = profile_order, group_order = group_order))
  }
  profile_order <- .multilpa_match_order(
    .multilpa_profile_signature(reference, reference),
    .multilpa_profile_signature(fit, reference))
  fit <- .multilpa_permute_profiles(fit, profile_order)
  # One group class has nothing to match, and relabelling it to itself would
  # only rewrite its names.
  if (reference$n_group_classes == 1L) {
    return(list(fit = fit, profile_order = profile_order, group_order = 1L))
  }
  group_order <- .multilpa_match_order(.multilpa_group_signature(reference),
                                       .multilpa_group_signature(fit))
  list(fit = .multilpa_permute_group_classes(fit, group_order),
       profile_order = profile_order, group_order = group_order)
}

#' Put a replicate's labels back on the original's
#'
#' @param replicate A fitted model from a resample or another seed.
#' @param reference The original fit.
#' @return The replicate, relabelled to the reference.
#' @noRd
.multilpa_align_labels <- function(replicate, reference) {
  .multilpa_pool_align(replicate, reference)$fit
}
