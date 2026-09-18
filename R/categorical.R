#' Encode categorical indicators as consecutive integer codes
#'
#' Maps each categorical indicator's observed values onto `1..C`, keeping the
#' original values so that results can be reported on the user's own scale.
#' Numeric, integer, logical, character and factor indicators are all accepted;
#' factor levels are preserved in their declared order, and other types are
#' ordered by `sort()`, so the coding is deterministic and does not depend on
#' row order.
#'
#' @param frame A data frame of categorical indicator columns.
#' @return A list with `codes` (an integer matrix of the same shape, `NA`
#'   preserved), `levels` (a list of the observed values per indicator, in code
#'   order), and `n_categories` (an integer vector of category counts).
#' @noRd
.ml_lpa_encode_categorical <- function(frame) {
  stopifnot("`frame` must be a data frame" = is.data.frame(frame),
            "`frame` must have at least one column" = ncol(frame) >= 1L)
  encoded <- lapply(names(frame), function(indicator) {
    value <- frame[[indicator]]
    if (!is.null(dim(value))) {
      stop(errorCondition(sprintf("Categorical indicator `%s` must be a vector.", indicator),
                          class = "mllpa_bad_categorical", call = NULL))
    }
    if (is.numeric(value) && any(!is.finite(value[!is.na(value)]))) {
      stop(errorCondition(sprintf("Categorical indicator `%s` contains non-finite values.", indicator),
                          class = "mllpa_bad_categorical", call = NULL))
    }
    levels_observed <- if (is.factor(value)) {
      levels(droplevels(value))
    } else sort(unique(value[!is.na(value)]))
    if (length(levels_observed) < 2L) {
      stop(errorCondition(sprintf("Categorical indicator `%s` has fewer than two observed categories.", indicator),
                          class = "mllpa_bad_categorical", call = NULL))
    }
    code <- if (is.factor(value)) {
      match(as.character(value), levels_observed)
    } else match(value, levels_observed)
    list(code = as.integer(code), levels = as.character(levels_observed))
  })
  names(encoded) <- names(frame)
  codes <- matrix(unlist(lapply(encoded, `[[`, "code"), use.names = FALSE),
                  nrow = nrow(frame), ncol = ncol(frame),
                  dimnames = list(NULL, names(frame)))
  list(codes = codes,
       levels = lapply(encoded, `[[`, "levels"),
       n_categories = vapply(encoded, function(part) length(part$levels), integer(1)))
}

#' Conditional log density of categorical indicators
#'
#' Sums the log response probabilities of the observed categories. Missing
#' cells contribute nothing, which is exactly the observed-data likelihood under
#' an ignorable missingness mechanism, so no imputation is involved.
#'
#' @param codes Integer matrix of category codes, `NA` where unobserved.
#' @param response_probabilities List of profiles-by-categories matrices, one
#'   per categorical indicator.
#' @return An observations-by-profiles matrix of log densities.
#' @noRd
.ml_lpa_categorical_log_density <- function(codes, response_probabilities) {
  stopifnot("`codes` must be an integer matrix" = is.matrix(codes),
            "`response_probabilities` must be a list" =
              is.list(response_probabilities),
            "one response matrix per categorical indicator is required" =
              length(response_probabilities) == ncol(codes))
  n_profiles <- nrow(response_probabilities[[1L]])
  contributions <- lapply(seq_len(ncol(codes)), function(indicator) {
    log_probabilities <- t(log(response_probabilities[[indicator]]))
    code <- codes[, indicator]
    observed <- !is.na(code)
    contribution <- matrix(0, nrow(codes), n_profiles)
    contribution[observed, ] <- log_probabilities[code[observed], , drop = FALSE]
    contribution
  })
  Reduce(`+`, contributions)
}

#' Maximize the expected complete-data log likelihood for categorical indicators
#'
#' The unrestricted multinomial M-step: each profile's response distribution is
#' the posterior-weighted share of observed categories. This is closed form,
#' unlike the Gaussian variance step, and needs no iteration.
#'
#' @param codes Integer matrix of category codes.
#' @param posteriors Observations-by-profiles posterior probability matrix.
#' @param n_categories Integer vector of category counts per indicator.
#' @param min_probability Lower bound on every response probability.
#' @return A list of profiles-by-categories matrices with rows summing to one.
#' @noRd
.ml_lpa_categorical_maximize <- function(codes, posteriors, n_categories,
                                          min_probability) {
  stopifnot("`codes` must be a matrix" = is.matrix(codes),
            "`posteriors` must have one row per observation" =
              nrow(posteriors) == nrow(codes),
            "`min_probability` must be in (0, 1)" =
              is.numeric(min_probability) && length(min_probability) == 1L &&
              is.finite(min_probability) && min_probability > 0 &&
              min_probability < 1)
  lapply(seq_len(ncol(codes)), function(indicator) {
    code <- codes[, indicator]
    observed <- !is.na(code)
    categories <- seq_len(n_categories[[indicator]])
    membership <- outer(code[observed], categories, "==")
    counts <- crossprod(posteriors[observed, , drop = FALSE], membership)
    totals <- rowSums(counts)
    if (any(totals <= 0)) {
      stop(errorCondition("A profile has no observed responses for an indicator.",
                          class = "mllpa_empty_profile", call = NULL))
    }
    # Normalizing after a simple pmax would push a bounded share back below the
    # bound, so the exact constrained solution is used instead.
    t(apply(counts, 1L, .ml_lpa_bound_probabilities,
            min_probability = min_probability))
  })
}

#' Exact bounded multinomial maximum-likelihood shares
#'
#' Solves `maximize sum(counts * log(p))` subject to `sum(p) == 1` and
#' `p >= min_probability`. The Karush-Kuhn-Tucker solution holds every category
#' whose unconstrained share falls below the bound exactly at it, and leaves the
#' remaining categories in their unconstrained proportions, rescaled to the
#' probability that is left over.
#'
#' @param counts Nonnegative weighted counts for one profile and indicator.
#' @param min_probability Lower bound on every returned probability.
#' @return A probability vector summing to one with no entry below the bound.
#' @noRd
.ml_lpa_bound_probabilities <- function(counts, min_probability) {
  stopifnot("`counts` must be nonnegative" = all(counts >= 0),
            "`min_probability` must leave room for every category" =
              min_probability * length(counts) < 1)
  active <- rep(FALSE, length(counts))
  # The active set is not known in advance and can only grow, so this loop runs
  # at most once per category. Vectorizing it would mean enumerating every
  # possible active set, which is exponential in the number of categories.
  repeat {
    free_total <- sum(counts[!active])
    remaining <- 1 - sum(active) * min_probability
    share <- numeric(length(counts))
    share[active] <- min_probability
    share[!active] <- if (free_total > 0) {
      remaining * counts[!active] / free_total
    } else remaining / max(sum(!active), 1L)
    newly_bound <- !active & share < min_probability
    if (!any(newly_bound)) return(share)
    active <- active | newly_bound
  }
}

#' Initialize categorical response probabilities from a hard assignment
#' @param codes Integer matrix of category codes.
#' @param assignments Integer profile assignment per observation.
#' @param n_profiles Number of profiles.
#' @param n_categories Integer vector of category counts per indicator.
#' @param min_probability Lower bound on every response probability.
#' @return A list of profiles-by-categories matrices.
#' @noRd
.ml_lpa_categorical_initialize <- function(codes, assignments, n_profiles,
                                           n_categories, min_probability) {
  stopifnot("`assignments` must have one value per observation" =
              length(assignments) == nrow(codes))
  posteriors <- vapply(seq_len(n_profiles), function(profile) {
    as.numeric(assignments == profile)
  }, numeric(nrow(codes)))
  # Pseudocounts are used only to initialize, never to estimate the final model.
  smoothed <- lapply(seq_len(ncol(codes)), function(indicator) {
    code <- codes[, indicator]
    observed <- !is.na(code)
    categories <- seq_len(n_categories[[indicator]])
    membership <- outer(code[observed], categories, "==")
    counts <- crossprod(posteriors[observed, , drop = FALSE], membership) + 0.5
    counts / rowSums(counts)
  })
  lapply(smoothed, function(block) {
    t(apply(block, 1L, .ml_lpa_bound_probabilities,
            min_probability = min_probability))
  })
}

#' One-hot representation used to initialize categorical clustering
#'
#' Expands categorical codes into centered indicator columns so that the shared
#' k-means initialization can see them alongside continuous indicators.
#'
#' @param codes Integer matrix of category codes.
#' @param n_categories Integer vector of category counts per indicator.
#' @return A numeric matrix with `sum(n_categories - 1)` columns.
#' @noRd
.ml_lpa_categorical_design <- function(codes, n_categories) {
  stopifnot("`codes` must be a matrix" = is.matrix(codes))
  blocks <- lapply(seq_len(ncol(codes)), function(indicator) {
    code <- codes[, indicator]
    # The final category is the redundant reference, dropped to keep the design
    # full rank; the missing cells take the indicator's observed shares.
    categories <- seq_len(n_categories[[indicator]] - 1L)
    membership <- outer(code, categories, "==") * 1
    shares <- colMeans(membership[!is.na(code), , drop = FALSE])
    membership[is.na(code), ] <- rep(shares, each = sum(is.na(code)))
    membership
  })
  do.call(cbind, blocks)
}

#' Count free categorical measurement parameters
#' @param n_profiles Number of profiles.
#' @param n_categories Integer vector of category counts per indicator.
#' @return Number of free response parameters.
#' @noRd
.ml_lpa_categorical_parameters <- function(n_profiles, n_categories) {
  as.integer(n_profiles * sum(n_categories - 1L))
}

#' Response probabilities as logit thresholds
#'
#' Converts unrestricted response probabilities to the cumulative-logit
#' thresholds that mixture software reports for categorical indicators, so that
#' native estimates can be compared with external output directly.
#'
#' @param probabilities A profiles-by-categories matrix with rows summing to one.
#' @return A profiles-by-(categories minus one) matrix of thresholds. Entry
#'   `[k, c]` is `qlogis(P(y <= c | profile k))`.
#' @noRd
.ml_lpa_categorical_thresholds <- function(probabilities) {
  stopifnot("`probabilities` must be a numeric matrix" =
              is.matrix(probabilities) && is.numeric(probabilities))
  cumulative <- t(apply(probabilities, 1L, cumsum))
  keep <- seq_len(ncol(probabilities) - 1L)
  thresholds <- stats::qlogis(cumulative[, keep, drop = FALSE])
  dimnames(thresholds) <- NULL
  thresholds
}
