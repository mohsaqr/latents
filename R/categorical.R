#' Validate native categorical values and obtain their ordered labels
#'
#' Validation precedes conversion to text so distinct numeric values cannot
#' silently become one category. Prediction subsets may have no variation.
#' @param value A categorical column in its original type.
#' @param indicator Its column name, for errors.
#' @param min_categories Minimum number of observed categories; fitting uses
#'   two, while preparing prediction subsets uses zero.
#' @param contributing `NULL`, or one logical per element of `value`: `FALSE`
#'   for a row of sampling weight zero. The categories are still read off
#'   every row, so each row has a code, but `min_categories` is counted on the
#'   contributing rows only: a row that enters no likelihood cannot supply the
#'   second category of an indicator (review 2026-10-26, ML-2).
#' @return The observed categories, in their native type and declared order;
#'   joined local-dependence items retain their full category grid.
#' @noRd
.multilpa_categorical_levels <- function(value, indicator, min_categories = 2L,
                                         contributing = NULL) {
  if (!(is.factor(value) || is.numeric(value) || is.logical(value) ||
        is.character(value))) {
    stop(errorCondition(sprintf(
      "Categorical indicator `%s` must be numeric, logical, character or a factor.", indicator),
      class = "latents_bad_categorical", call = NULL))
  }
  if (!is.null(dim(value))) {
    stop(errorCondition(sprintf("Categorical indicator `%s` must be a vector.", indicator),
                        class = "latents_bad_categorical", call = NULL))
  }
  if (is.numeric(value) && any(!is.finite(value[!is.na(value)]))) {
    stop(errorCondition(sprintf("Categorical indicator `%s` contains non-finite values.", indicator),
                        class = "latents_bad_categorical", call = NULL))
  }
  levels_observed <- if (is.factor(value)) {
    # The joined item's declared levels include unobserved combinations,
    # each of which has an estimated probability at the model's floor. Its
    # markers of half-observed rows ("a::", ":::b") are not categories: such
    # a row is a set of joint categories (`.latents_partial_map()`).
    if (length(attr(value, "latents_pair", exact = TRUE)) == 2L) {
      declared <- levels(value)
      declared[!.latents_partial_levels(declared)]
    } else levels(droplevels(value))
  } else sort(unique(value[!is.na(value)]))
  if (length(levels_observed) < min_categories) {
    stop(errorCondition(sprintf("Categorical indicator `%s` has fewer than two observed categories.", indicator),
                        class = "latents_bad_categorical", call = NULL))
  }
  # A joined item keeps its full grid whatever was observed, so only a plain
  # indicator is counted again on the rows that carry weight.
  if (!is.null(contributing) && !all(contributing) &&
      length(attr(value, "latents_pair", exact = TRUE)) != 2L) {
    stopifnot("`contributing` must hold one logical per value" =
                is.logical(contributing) && !anyNA(contributing) &&
                length(contributing) == length(value))
    carried <- value[contributing & !is.na(value)]
    n_carried <- if (is.factor(carried)) nlevels(droplevels(carried)) else
      length(unique(carried))
    if (n_carried < min_categories) {
      stop(errorCondition(sprintf(paste(
        "Categorical indicator `%s` has fewer than two observed categories in",
        "the rows of positive weight."), indicator),
        class = "latents_bad_categorical", call = NULL))
    }
  }
  # Categories are reported, aligned and mapped back to values by their
  # labels, so two distinct values with one label (0.3 and 0.1 + 0.2 both
  # print as "0.3") would be two categories nobody can tell apart.
  labels <- as.character(levels_observed)
  if (anyDuplicated(labels) > 0L) {
    stop(errorCondition(sprintf(paste(
      "Categorical indicator `%s` has distinct values that print alike (%s);",
      "round the column so that each category has its own label."),
      indicator, paste(sprintf("\"%s\"", unique(labels[duplicated(labels)])),
                       collapse = ", ")),
      class = "latents_bad_categorical", call = NULL))
  }
  levels_observed
}

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
#' @param contributing `NULL`, or one logical per row, `FALSE` for a row of
#'   sampling weight zero; passed to `.multilpa_categorical_levels()`.
#' @return A list with `codes` (an integer matrix of the same shape, `NA`
#'   preserved), `levels` (a list of the observed values per indicator, in code
#'   order), and `n_categories` (an integer vector of category counts). When a
#'   joined `local_dependence` item has rows observing one of its two items,
#'   those rows are `NA` in `codes`, which carries the attribute
#'   `"latents_partial"` (see `.multilpa_attach_partial()`).
#' @noRd
.multilpa_encode_categorical <- function(frame, contributing = NULL) {
  stopifnot("`frame` must be a data frame" = is.data.frame(frame),
            "`frame` must have at least one column" = ncol(frame) >= 1L)
  encoded <- lapply(names(frame), function(indicator) {
    value <- frame[[indicator]]
    levels_observed <- .multilpa_categorical_levels(value, indicator,
                                                    contributing = contributing)
    code <- if (is.factor(value)) {
      match(as.character(value), levels_observed)
    } else match(value, levels_observed)
    list(code = as.integer(code), levels = as.character(levels_observed))
  })
  names(encoded) <- names(frame)
  codes <- matrix(unlist(lapply(encoded, `[[`, "code"), use.names = FALSE),
                  nrow = nrow(frame), ncol = ncol(frame),
                  dimnames = list(NULL, names(frame)))
  codes <- .multilpa_attach_partial(codes, frame,
                                    lapply(encoded, `[[`, "levels"))
  list(codes = codes,
       levels = lapply(encoded, `[[`, "levels"),
       n_categories = vapply(encoded, function(part) length(part$levels), integer(1)))
}

#' Record the half-observed rows of joined pair items on a code matrix
#'
#' A row observing one item of a `local_dependence` pair is not missing: it
#' observes the set of joint categories that agree with that item, and its
#' likelihood is the sum of their probabilities. Its code stays `NA`, so every
#' consumer that reads codes as categories sees no category, and the sets
#' travel as the attribute `"latents_partial"`: a list with one entry per
#' column, `NULL` or a `.latents_partial_map()` result. The attribute is set
#' only when some row is half observed, so other code matrices are unchanged.
#'
#' @param codes The integer code matrix.
#' @param frame The columns the codes were read from (a data frame or list).
#' @param levels The categories of each column, in code order.
#' @return `codes`, with the attribute when a joined item has such rows.
#' @noRd
.multilpa_attach_partial <- function(codes, frame, levels) {
  maps <- lapply(seq_len(ncol(codes)), function(column) {
    value <- frame[[column]]
    if (!is.factor(value) ||
        length(attr(value, "latents_pair", exact = TRUE)) != 2L) return(NULL)
    .latents_partial_map(value, as.character(levels[[column]]))
  })
  if (all(vapply(maps, is.null, logical(1)))) return(codes)
  names(maps) <- colnames(codes)
  attr(codes, "latents_partial") <- maps
  codes
}

#' The half-observed set map of one column of a code matrix
#' @param codes A code matrix, possibly carrying `"latents_partial"`.
#' @param column The column position.
#' @return `NULL`, or the column's `.latents_partial_map()` result.
#' @noRd
.multilpa_partial_of <- function(codes, column) {
  maps <- attr(codes, "latents_partial", exact = TRUE)
  if (is.null(maps)) NULL else maps[[column]]
}

#' Rows-by-sets membership of a half-observed set map
#' @param map A `.latents_partial_map()` result.
#' @param n_categories The joint item's number of categories.
#' @return A sets-by-categories 0/1 matrix.
#' @noRd
.multilpa_partial_sets <- function(map, n_categories) {
  t(vapply(map$sets, function(set) {
    as.numeric(seq_len(n_categories) %in% set)
  }, numeric(n_categories)))
}

#' Warn when numeric categorical indicators look continuous
#'
#' Every distinct value of a categorical indicator is fitted as a category,
#' so a continuous column declared categorical (directly, through
#' `type = "lca"` or through [multilca()]) becomes many categories without
#' any error. Numeric indicators with non-integer values or more than
#' `max_categories` categories are reported in one classed warning.
#' @param frame The categorical indicator columns.
#' @param encoded The result of `.multilpa_encode_categorical()`, or `NULL`.
#' @param max_categories Largest category count not reported.
#' @return NULL, invisibly.
#' @noRd
.latents_warn_suspect_categorical <- function(frame, encoded, max_categories = 10L) {
  if (is.null(encoded) || ncol(frame) == 0L) return(invisible(NULL))
  findings <- vapply(names(frame), function(indicator) {
    value <- frame[[indicator]]
    if (!is.numeric(value)) return(NA_character_)
    observed <- value[!is.na(value)]
    n_categories <- encoded$n_categories[[indicator]]
    if (any(abs(observed - round(observed)) > sqrt(.Machine$double.eps))) {
      sprintf("`%s` has non-integer values (%d categories)", indicator, n_categories)
    } else if (n_categories > max_categories) {
      sprintf("`%s` has %d categories (%s to %s)", indicator, n_categories,
              format(min(observed)), format(max(observed)))
    } else NA_character_
  }, character(1))
  findings <- findings[!is.na(findings)]
  if (length(findings) == 0L) return(invisible(NULL))
  message <- if (length(findings) == 1L) {
    sprintf(paste("Categorical indicator %s; each distinct value is fitted as a",
                  "category. If it is continuous, a model with continuous (Gaussian)",
                  "indicators, such as a latent profile model, may be more appropriate."),
            findings)
  } else {
    sprintf(paste("%d categorical indicators may be continuous: %s. Each distinct",
                  "value is fitted as a category. If they are continuous, a model with",
                  "continuous (Gaussian) indicators, such as a latent profile model, may",
                  "be more appropriate."),
            length(findings), paste(findings, collapse = "; "))
  }
  warning(warningCondition(message, class = "latents_suspect_categorical",
                           indicators = names(findings), call = NULL))
  invisible(NULL)
}

#' Let a data warning through once while an enumeration fits its candidates
#'
#' Every candidate prepares the same indicators, so a warning about the data
#' would otherwise be repeated once per candidate.
#' @param expr The fitting call.
#' @param class The class of the warning to report once.
#' @return The value of `expr`.
#' @noRd
.latents_report_once <- function(expr, class = "latents_suspect_categorical") {
  seen <- FALSE
  withCallingHandlers(expr, warning = function(condition) {
    if (!inherits(condition, class)) return(invisible(NULL))
    if (seen) invokeRestart("muffleWarning")
    seen <<- TRUE
  })
}

#' One original, typed value for each encoded category
#'
#' Category labels alone do not reconstruct a numeric or ordered-factor
#' indicator: converting them to text changes the data a refit sees, and
#' sorting the text can change the codes of a factor with declared levels.
#' Keeping one original value per code preserves the input type and factor
#' attributes without copying every observation into the fit a second time.
#'
#' @param data The original fitting data frame.
#' @param categorical Names of its categorical indicators.
#' @param encoded The result of `.multilpa_encode_categorical()`, or `NULL`.
#' @return A named list of typed vectors, one value per observed category.
#' @noRd
.multilpa_categorical_values <- function(data, categorical, encoded) {
  if (length(categorical) == 0L) return(NULL)
  stats::setNames(lapply(categorical, function(name) {
    column <- data[[name]]
    values <- column[match(encoded$levels[[name]], as.character(column))]
    # Factor subsetting drops custom attributes. Keep the explicit pair
    # provenance: labels alone can also be ordinary user categories.
    attr(values, "latents_pair") <- attr(column, "latents_pair", exact = TRUE)
    values
  }), categorical)
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
.multilpa_categorical_log_density <- function(codes, response_probabilities) {
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
    # A half-observed pair row contributes the probability of its observed
    # item: the sum of the joint probabilities that agree with it.
    partial <- .multilpa_partial_of(codes, indicator)
    if (!is.null(partial)) {
      sets <- .multilpa_partial_sets(partial, ncol(response_probabilities[[indicator]]))
      set_log <- log(tcrossprod(sets, response_probabilities[[indicator]]))
      half <- !is.na(partial$index)
      contribution[half, ] <- set_log[partial$index[half], , drop = FALSE]
    }
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
.multilpa_categorical_maximize <- function(codes, posteriors, n_categories,
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
    partial <- .multilpa_partial_of(codes, indicator)
    if (!is.null(partial)) {
      return(.multilpa_partial_maximize(counts, partial, posteriors,
                                        length(categories), min_probability))
    }
    totals <- rowSums(counts)
    if (any(totals <= 0)) {
      stop(errorCondition("A profile has no observed responses for an indicator.",
                          class = "latents_empty_profile", call = NULL))
    }
    # Normalizing after a simple pmax would push a bounded share back below the
    # bound, so the exact constrained solution is used instead.
    t(apply(counts, 1L, .multilpa_bound_probabilities,
            min_probability = min_probability))
  })
}

#' M-step of a joined pair item with half-observed rows
#'
#' Each profile maximizes `sum_j n_j log p_j + sum_s m_s log(sum_{j in S_s}
#' p_j)` over the bounded simplex, where `n` are the posterior-weighted counts
#' of fully observed rows and `m_s` those of the rows observing set `S_s` (one
#' item of the pair). The objective is concave -- a sum of logarithms of
#' linear functions of `p` -- and the bounded simplex is convex, so the
#' maximum is global. It has no closed form when both items have
#' half-observed rows, so it is reached by an inner EM over the missing item:
#' each step spreads `m_s` over `S_s` in proportion to the current `p` and
#' solves the bounded multinomial exactly, which raises the objective
#' (a minorize-maximize step). It stops when no probability moves by more than
#' `tolerance` (1e-13, below the rounding of probabilities of order one), and
#' after `max_steps` steps at most, with a `latents_no_converge` warning.
#'
#' @param counts Profiles-by-categories counts of the fully observed rows.
#' @param partial The column's `.latents_partial_map()` result.
#' @param posteriors Observations-by-profiles posterior weights.
#' @param n_categories The number of joint categories.
#' @param min_probability Lower bound on every probability.
#' @param tolerance,max_steps Stopping rule of the inner EM.
#' @return A profiles-by-categories probability matrix.
#' @noRd
.multilpa_partial_maximize <- function(counts, partial, posteriors, n_categories,
                                       min_probability, tolerance = 1e-13,
                                       max_steps = 10000L) {
  half <- !is.na(partial$index)
  sets <- .multilpa_partial_sets(partial, n_categories)
  set_membership <- outer(partial$index[half], seq_len(nrow(sets)), "==") * 1
  # Profiles-by-sets weights of the half-observed rows.
  set_counts <- crossprod(posteriors[half, , drop = FALSE], set_membership)
  if (any(rowSums(counts) + rowSums(set_counts) <= 0)) {
    stop(errorCondition("A profile has no observed responses for an indicator.",
                        class = "latents_empty_profile", call = NULL))
  }
  bound <- function(expected) {
    t(apply(expected, 1L, .multilpa_bound_probabilities,
            min_probability = min_probability))
  }
  spread <- function(probability) {
    counts + (set_counts / tcrossprod(probability, sets)) %*% sets * probability
  }
  # Start inside the simplex: each set's weight spread evenly over its cells.
  probability <- bound(counts + sweep(set_counts, 2L, rowSums(sets), "/") %*% sets)
  step <- 0L
  moved <- Inf
  # Each inner EM step depends on the previous one; nothing to vectorize.
  while (moved > tolerance && step < max_steps) {
    updated <- bound(spread(probability))
    moved <- max(abs(updated - probability))
    probability <- updated
    step <- step + 1L
  }
  if (moved > tolerance) {
    warning(warningCondition(sprintf(paste(
      "The joint distribution of a local-dependence pair with half-observed",
      "rows did not converge in %d inner EM steps (largest change %.2g)."),
      max_steps, moved), class = "latents_no_converge", call = NULL))
  }
  probability
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
.multilpa_bound_probabilities <- function(counts, min_probability) {
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
.multilpa_categorical_initialize <- function(codes, assignments, n_profiles,
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
    t(apply(block, 1L, .multilpa_bound_probabilities,
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
.multilpa_categorical_design <- function(codes, n_categories) {
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
.multilpa_categorical_parameters <- function(n_profiles, n_categories) {
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
.multilpa_categorical_thresholds <- function(probabilities) {
  stopifnot("`probabilities` must be a numeric matrix" =
              is.matrix(probabilities) && is.numeric(probabilities))
  cumulative <- t(apply(probabilities, 1L, cumsum))
  keep <- seq_len(ncol(probabilities) - 1L)
  thresholds <- stats::qlogis(cumulative[, keep, drop = FALSE])
  dimnames(thresholds) <- NULL
  thresholds
}
