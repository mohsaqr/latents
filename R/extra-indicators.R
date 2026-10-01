# Ordinal and count indicators.
#
# Both are conditionally independent of the other indicators given the profile,
# like the categorical ones, and enter the E-step as one more n x K block of
# log densities. They are kept in one "extra" object so that every engine sees
# a single data argument and a single set of parameters:
#
#   ordinal   adjacent-category logit with category intercepts shared by every
#             profile and one location per profile (Latent GOLD's default
#             ordinal model): log P(k) / P(k - 1) = a_k - a_(k-1) + eta_c, so
#             log P(k | c) = a_k + eta_c (k - 1) - normaliser, with a_1 = 0
#             and the last profile's location fixed at 0. (K - 1) + (C - 1)
#             parameters per indicator, against C (K - 1) for an unrestricted
#             categorical indicator.
#   count     Poisson, one mean per profile.
#
# Missing values are integrated out (they contribute nothing), as for the
# other indicators under `missing = "fiml"`.

#' Prepare ordinal and count indicators
#'
#' @param data The data frame of the fit.
#' @param ordinal,count Names of the ordinal and count indicators.
#' @param missing `"error"` or `"fiml"`.
#' @return `NULL` when neither is given; otherwise a list with `ordinal` (an
#'   integer code matrix, codes 1..K), `ordinal_levels` (the category labels
#'   per indicator), `ordinal_categories` (K per indicator) and `count` (a
#'   numeric matrix of non-negative integers).
#' @noRd
.latents_prepare_extra <- function(data, ordinal = character(),
                                   count = character(), missing = "error") {
  if (length(ordinal) == 0L && length(count) == 0L) return(NULL)
  levels_of <- function(value) {
    if (is.factor(value)) levels(value) else as.character(sort(unique(value)))
  }
  ordinal_levels <- stats::setNames(lapply(ordinal, function(name) {
    value <- data[[name]]
    if (!(is.factor(value) || is.numeric(value)) ||
        (is.numeric(value) && any(value[!is.na(value)] != round(value[!is.na(value)])))) {
      stop(errorCondition(sprintf(paste(
        "Ordinal indicator `%s` must be an (ordered) factor or whole numbers;",
        "its categories are ordered by factor level or by value."), name),
        class = "latents_bad_data", call = NULL))
    }
    levels_of(value)
  }), ordinal)
  ordinal_codes <- if (length(ordinal) == 0L) NULL else
    matrix(vapply(ordinal, function(name) {
      as.integer(match(as.character(data[[name]]), ordinal_levels[[name]]))
    }, integer(nrow(data))), nrow(data), length(ordinal),
    dimnames = list(NULL, ordinal))
  categories <- lengths(ordinal_levels)
  if (any(categories < 2L)) {
    stop(errorCondition(sprintf(
      "Ordinal indicator %s has fewer than two observed categories.",
      paste(sprintf("`%s`", ordinal[categories < 2L]), collapse = ", ")),
      class = "latents_bad_data", call = NULL))
  }
  counts <- if (length(count) == 0L) NULL else
    matrix(vapply(count, function(name) {
      value <- data[[name]]
      observed <- value[!is.na(value)]
      if (!is.numeric(value) || any(!is.finite(observed)) || any(observed < 0) ||
          any(observed != round(observed))) {
        stop(errorCondition(sprintf(
          "Count indicator `%s` must hold non-negative whole numbers.", name),
          class = "latents_bad_data", call = NULL))
      }
      if (length(observed) > 0L && all(observed == 0)) {
        stop(errorCondition(sprintf(
          "Count indicator `%s` is zero throughout; a Poisson mean is not identified.",
          name), class = "latents_bad_data", call = NULL))
      }
      as.numeric(value)
    }, numeric(nrow(data))), nrow(data), length(count),
    dimnames = list(NULL, count))
  if (identical(missing, "error") && (anyNA(ordinal_codes) || anyNA(counts))) {
    stop(errorCondition("Indicators contain missing or non-finite values.",
                        class = "latents_bad_data", call = NULL))
  }
  observed <- c(if (is.null(ordinal_codes)) NULL else colSums(!is.na(ordinal_codes)),
                if (is.null(counts)) NULL else colSums(!is.na(counts)))
  if (any(observed == 0L)) {
    stop(errorCondition("Every indicator must have observed values.",
                        class = "latents_bad_data", call = NULL))
  }
  list(ordinal = ordinal_codes, ordinal_levels = ordinal_levels,
       ordinal_categories = as.integer(categories), count = counts)
}

#' The names of a fit's ordinal and count indicators
#' @noRd
.latents_extra_names <- function(extra) {
  list(ordinal = colnames(extra$ordinal) %||% character(),
       count = colnames(extra$count) %||% character())
}

#' Log category probabilities of one ordinal indicator, profiles x categories
#'
#' @param intercepts The K - 1 intercepts a_2..a_K (a_1 = 0).
#' @param locations The C profile locations (the last is 0).
#' @return A C x K matrix of log probabilities.
#' @noRd
.latents_ordinal_log_probabilities <- function(intercepts, locations) {
  scores <- seq_along(c(0, intercepts)) - 1
  logits <- outer(locations, scores) + matrix(c(0, intercepts), length(locations),
                                              length(scores), byrow = TRUE)
  sweep(logits, 1L, .multilpa_log_sum_exp(logits), "-")
}

#' Row-by-profile log densities of the ordinal and count indicators
#'
#' @param extra Result of `.latents_prepare_extra()`, or `NULL`.
#' @param parameters A parameter list carrying `ordinal_intercepts` (a list of
#'   intercept vectors), `ordinal_locations` (C x J matrix) and `count_means`
#'   (C x J matrix).
#' @param n_rows,n_profiles Dimensions of the result.
#' @return An n x C matrix; zero when there is nothing to add.
#' @noRd
.latents_extra_log_density <- function(extra, parameters, n_rows, n_profiles) {
  total <- matrix(0, n_rows, n_profiles)
  if (is.null(extra)) return(total)
  ordinal <- lapply(seq_len(ncol(extra$ordinal) %||% 0L), function(j) {
    log_probability <- .latents_ordinal_log_probabilities(
      parameters$ordinal_intercepts[[j]], parameters$ordinal_locations[, j])
    code <- extra$ordinal[, j]
    contribution <- matrix(0, n_rows, n_profiles)
    observed <- !is.na(code)
    contribution[observed, ] <- t(log_probability[, code[observed], drop = FALSE])
    contribution
  })
  count <- lapply(seq_len(ncol(extra$count) %||% 0L), function(j) {
    value <- extra$count[, j]
    observed <- !is.na(value)
    contribution <- matrix(0, n_rows, n_profiles)
    means <- parameters$count_means[, j]
    contribution[observed, ] <- outer(value[observed], log(means)) -
      matrix(means, sum(observed), n_profiles, byrow = TRUE) -
      lgamma(value[observed] + 1)
    contribution
  })
  Reduce(`+`, c(ordinal, count), total)
}

#' Posterior-weighted category counts of one ordinal indicator, C x K
#' @noRd
.latents_ordinal_counts <- function(code, posteriors, n_categories) {
  observed <- !is.na(code)
  indicator <- vapply(seq_len(n_categories), function(k) {
    as.numeric(observed & code == k)
  }, numeric(length(code)))
  crossprod(posteriors, matrix(indicator, length(code), n_categories))
}

#' Maximize one ordinal indicator's expected log likelihood
#'
#' The objective `sum_ck N_ck log P(k | c)` is concave in the intercepts and
#' locations (a multinomial logit), so Newton-Raphson with the exact Hessian
#' and step halving climbs it monotonically. The free coordinates are the
#' K - 1 intercepts then the first C - 1 locations.
#'
#' @param counts C x K posterior-weighted category counts.
#' @param intercepts,locations The current point.
#' @param max_newton Most Newton steps.
#' @return A list with `intercepts` and `locations`.
#' @noRd
.latents_ordinal_maximize <- function(counts, intercepts, locations,
                                      max_newton = 100L) {
  n_profiles <- nrow(counts)
  n_categories <- ncol(counts)
  scores <- seq_len(n_categories) - 1
  # Features of category k in profile c: intercept dummies for k >= 2, then
  # the score k - 1 in the location column of c (none for the last profile).
  features <- lapply(seq_len(n_profiles), function(profile) {
    location <- matrix(0, n_categories, n_profiles - 1L)
    if (profile < n_profiles) location[, profile] <- scores
    cbind(diag(n_categories)[, -1L, drop = FALSE], location)
  })
  totals <- rowSums(counts)
  pack <- function(a, eta) c(a, eta[-n_profiles])
  unpack <- function(theta) {
    list(intercepts = theta[seq_len(n_categories - 1L)],
         locations = c(theta[n_categories - 1L + seq_len(n_profiles - 1L)], 0))
  }
  objective <- function(theta) {
    point <- unpack(theta)
    sum(counts * .latents_ordinal_log_probabilities(point$intercepts, point$locations))
  }
  derivatives <- function(theta) {
    point <- unpack(theta)
    probability <- exp(.latents_ordinal_log_probabilities(point$intercepts,
                                                          point$locations))
    pieces <- lapply(seq_len(n_profiles), function(profile) {
      feature <- features[[profile]]
      expected <- drop(crossprod(feature, probability[profile, ]))
      list(gradient = drop(crossprod(feature, counts[profile, ])) -
             totals[profile] * expected,
           information = totals[profile] *
             (crossprod(feature, feature * probability[profile, ]) -
                tcrossprod(expected)))
    })
    list(gradient = Reduce(`+`, lapply(pieces, `[[`, "gradient")),
         information = Reduce(`+`, lapply(pieces, `[[`, "information")))
  }
  theta <- pack(intercepts, locations)
  value <- objective(theta)
  iteration <- 0L
  # Newton steps refine one estimate in sequence; nothing to vectorize. A ridge
  # far below the data scale keeps an unobserved category solvable.
  while (iteration < max_newton) {
    step_parts <- derivatives(theta)
    if (max(abs(step_parts$gradient)) <= 1e-10 * (1 + sum(totals))) break
    information <- step_parts$information
    ridge <- 1e-10 * (1 + max(abs(diag(information))))
    factor <- chol(information + diag(ridge, nrow(information)))
    direction <- backsolve(factor, forwardsolve(t(factor), step_parts$gradient))
    step <- 1
    accepted <- FALSE
    while (step > 1e-10 && !accepted) {
      candidate <- theta + step * direction
      candidate_value <- objective(candidate)
      accepted <- is.finite(candidate_value) && candidate_value >= value
      if (!accepted) step <- step / 2
    }
    if (!accepted) break
    theta <- candidate
    value <- candidate_value
    iteration <- iteration + 1L
  }
  unpack(theta)
}

#' M-step of the ordinal and count indicators
#'
#' @param extra Prepared indicators.
#' @param posteriors n x C posterior weights (weighted, under sampling weights).
#' @param previous The current parameters, the Newton start for ordinal ones.
#' @return A list with `ordinal_intercepts`, `ordinal_locations` and
#'   `count_means`, each `NULL` when absent.
#' @noRd
.latents_extra_maximize <- function(extra, posteriors, previous) {
  n_profiles <- ncol(posteriors)
  ordinal <- lapply(seq_len(ncol(extra$ordinal) %||% 0L), function(j) {
    counts <- .latents_ordinal_counts(extra$ordinal[, j], posteriors,
                                      extra$ordinal_categories[[j]])
    .latents_ordinal_maximize(counts, previous$ordinal_intercepts[[j]],
                              previous$ordinal_locations[, j])
  })
  count_means <- if (is.null(extra$count)) NULL else {
    observed <- !is.na(extra$count)
    values <- extra$count
    values[!observed] <- 0
    # A Poisson mean is the posterior-weighted mean count; a profile with no
    # observed count keeps a tiny positive mean rather than an undefined one.
    crossprod(posteriors, values) / pmax(crossprod(posteriors, observed * 1),
                                         .Machine$double.xmin)
  }
  if (!is.null(count_means)) {
    count_means <- pmax(count_means, 1e-10)
    dimnames(count_means) <- list(NULL, colnames(extra$count))
  }
  list(ordinal_intercepts = if (length(ordinal) == 0L) NULL else
         stats::setNames(lapply(ordinal, `[[`, "intercepts"), colnames(extra$ordinal)),
       ordinal_locations = if (length(ordinal) == 0L) NULL else
         matrix(vapply(ordinal, `[[`, numeric(n_profiles), "locations"),
                n_profiles, length(ordinal),
                dimnames = list(NULL, colnames(extra$ordinal))),
       count_means = count_means)
}

#' Starting values from a hard assignment of the rows
#' @param extra Prepared indicators.
#' @param assignments Profile of each row.
#' @param n_profiles Number of profiles.
#' @return The parameter blocks, as `.latents_extra_maximize()` returns them.
#' @noRd
.latents_extra_initialize <- function(extra, assignments, n_profiles) {
  # Soften the partition so no profile starts with an empty category or a
  # zero mean.
  posteriors <- 0.9 * outer(assignments, seq_len(n_profiles), "==") + 0.1 / n_profiles
  start <- list(
    ordinal_intercepts = lapply(extra$ordinal_categories, function(k) numeric(k - 1L)),
    ordinal_locations = matrix(0, n_profiles, length(extra$ordinal_categories)))
  .latents_extra_maximize(extra, posteriors, start)
}

#' Standardized columns for the initial clustering
#' @return An n x J matrix (ordinal codes and log counts, mean-imputed and
#'   standardized), or `NULL`.
#' @noRd
.latents_extra_clustering <- function(extra) {
  if (is.null(extra)) return(NULL)
  columns <- cbind(extra$ordinal, if (!is.null(extra$count)) log1p(extra$count))
  if (is.null(columns) || ncol(columns) == 0L) return(NULL)
  matrix(vapply(seq_len(ncol(columns)), function(j) {
    value <- columns[, j]
    value[is.na(value)] <- mean(value, na.rm = TRUE)
    spread <- stats::sd(value)
    if (!is.finite(spread) || spread <= 0) spread <- 1
    (value - mean(value)) / spread
  }, numeric(nrow(columns))), nrow(columns))
}

#' Free parameters of the ordinal and count indicators
#' @noRd
.latents_extra_n_parameters <- function(extra, n_profiles) {
  if (is.null(extra)) return(0L)
  as.integer(sum(extra$ordinal_categories - 1L) +
               length(extra$ordinal_categories) * (n_profiles - 1L) +
               (ncol(extra$count) %||% 0L) * n_profiles)
}

#' Distinct observed rows of the ordinal and count indicators, for the
#' identification check
#' @noRd
.latents_extra_matrix <- function(extra) {
  if (is.null(extra)) return(NULL)
  cbind(extra$ordinal, extra$count)
}

#' Validate the ordinal and count indicator names and their pairings
#'
#' @param vars,categorical,ordinal,count The indicator arguments of the fit.
#' @param unsupported Named logicals: options the new indicators do not take
#'   yet; a `TRUE` one is refused.
#' @return `NULL`, invisibly.
#' @noRd
.latents_check_extra_arguments <- function(vars, categorical, ordinal, count,
                                           unsupported = list()) {
  stopifnot(
    "`ordinal` must be a character vector of indicator names" =
      is.character(ordinal) && !anyNA(ordinal),
    "`count` must be a character vector of indicator names" =
      is.character(count) && !anyNA(count))
  if (length(ordinal) == 0L && length(count) == 0L) return(invisible(NULL))
  typed <- c(categorical, ordinal, count)
  if (!all(c(ordinal, count) %in% vars) || anyDuplicated(typed)) {
    stop(errorCondition(paste(
      "`ordinal` and `count` must name indicators in `vars`, and each indicator",
      "can have one type only (categorical, ordinal or count)."),
      class = "latents_bad_argument", call = NULL))
  }
  refused <- names(unsupported)[vapply(unsupported, isTRUE, logical(1))]
  if (length(refused) > 0L) {
    stop(errorCondition(sprintf(
      "Ordinal and count indicators are not supported with %s yet.",
      paste(refused, collapse = ", ")),
      class = "latents_unsupported_indicator", call = NULL))
  }
  invisible(NULL)
}

#' Name the ordinal and count parameter blocks of a fit
#' @noRd
.latents_label_extra <- function(parameters, extra, profile_names) {
  if (is.null(extra)) return(parameters)
  if (!is.null(parameters$ordinal_intercepts)) {
    parameters$ordinal_intercepts <- stats::setNames(
      lapply(seq_along(parameters$ordinal_intercepts), function(j) {
        stats::setNames(unname(parameters$ordinal_intercepts[[j]]),
                        extra$ordinal_levels[[j]][-1L])
      }), colnames(extra$ordinal))
    dimnames(parameters$ordinal_locations) <- list(profile_names,
                                                   colnames(extra$ordinal))
  }
  if (!is.null(parameters$count_means)) {
    dimnames(parameters$count_means) <- list(profile_names, colnames(extra$count))
  }
  parameters
}

#' The measurement model a fit's indicator types make up
#' @noRd
.latents_measurement_model <- function(gaussian, categorical, ordinal, count) {
  present <- c(gaussian = gaussian, categorical = categorical,
               ordinal = length(ordinal) > 0L, count = length(count) > 0L)
  if (sum(present) == 1L) names(present)[present] else "mixed"
}

#' Inference coordinates of the ordinal and count indicators
#'
#' Per ordinal indicator: its K - 1 intercepts (shared by every profile), then
#' the first C - 1 profile locations (the last is the reference, zero); then
#' per count indicator its C means. Ordinal coordinates are already on the
#' logit scale, so the natural and unconstrained values coincide; count means
#' are estimated on the log scale.
#'
#' @param parameters A parameter list or fit carrying the blocks.
#' @param extra Prepared indicators (for names and levels).
#' @param n_profiles Number of profiles.
#' @param scale `"natural"` or `"unconstrained"`.
#' @return A named numeric vector in the package's `kind[...]` label grammar.
#' @noRd
.latents_extra_coordinates <- function(parameters, extra, n_profiles, scale) {
  if (is.null(extra)) return(stats::setNames(numeric(0), character(0)))
  names_of <- .latents_extra_names(extra)
  ordinal <- lapply(seq_along(names_of$ordinal), function(j) {
    indicator <- names_of$ordinal[[j]]
    intercepts <- unname(parameters$ordinal_intercepts[[j]])
    locations <- unname(parameters$ordinal_locations[, j])[seq_len(n_profiles - 1L)]
    c(stats::setNames(intercepts, sprintf("ordinal_intercept[shared,%s,%s]", indicator,
                                          extra$ordinal_levels[[j]][-1L])),
      stats::setNames(locations, sprintf("ordinal_location[%d,%s]",
                                         seq_len(n_profiles - 1L), indicator)))
  })
  count <- lapply(seq_along(names_of$count), function(j) {
    indicator <- names_of$count[[j]]
    means <- unname(parameters$count_means[, j])
    if (identical(scale, "natural")) {
      stats::setNames(means, sprintf("count_mean[%d,%s]", seq_len(n_profiles), indicator))
    } else {
      stats::setNames(log(means), sprintf("log_count_mean[%d,%s]",
                                          seq_len(n_profiles), indicator))
    }
  })
  unlist(c(ordinal, count))
}

#' Number of ordinal and count coordinates (the same on both scales)
#' @noRd
.latents_extra_width <- function(extra, n_profiles) {
  .latents_extra_n_parameters(extra, n_profiles)
}

#' Rebuild the ordinal and count blocks from unconstrained coordinates
#' @return A list with `ordinal_intercepts`, `ordinal_locations`, `count_means`.
#' @noRd
.latents_extra_decode <- function(theta, extra, n_profiles) {
  at <- 0L
  take <- function(n) {
    values <- theta[at + seq_len(n)]
    at <<- at + n
    values
  }
  ordinal <- lapply(extra$ordinal_categories, function(k) {
    list(intercepts = take(k - 1L), locations = c(take(n_profiles - 1L), 0))
  })
  count <- lapply(seq_len(ncol(extra$count) %||% 0L), function(j) exp(take(n_profiles)))
  stopifnot("the extra coordinates must be used exactly" = at == length(theta))
  list(ordinal_intercepts = if (length(ordinal) == 0L) NULL else
         lapply(ordinal, `[[`, "intercepts"),
       ordinal_locations = if (length(ordinal) == 0L) NULL else
         matrix(vapply(ordinal, `[[`, numeric(n_profiles), "locations"), n_profiles),
       count_means = if (length(count) == 0L) NULL else
         matrix(unlist(count), n_profiles))
}

#' Scores of the ordinal and count coordinates
#'
#' By Fisher's identity each row's score is its posterior-weighted
#' complete-data score: for an ordinal indicator the observed feature minus
#' its expectation in each profile, for a count `y - lambda` on the log mean.
#' Missing values contribute nothing.
#'
#' @param extra Prepared indicators.
#' @param posteriors n x C posteriors (weighted, under sampling weights).
#' @param parameters The blocks at which to evaluate.
#' @param group_index Group of each row, or `NULL` for the total.
#' @return Groups (or one row) by `.latents_extra_width()` columns.
#' @noRd
.latents_extra_scores <- function(extra, posteriors, parameters, group_index = NULL) {
  n_rows <- nrow(posteriors)
  n_profiles <- ncol(posteriors)
  if (is.null(extra)) {
    return(matrix(numeric(0), if (is.null(group_index)) 1L else max(group_index), 0L))
  }
  ordinal <- lapply(seq_len(ncol(extra$ordinal) %||% 0L), function(j) {
    code <- extra$ordinal[, j]
    observed <- !is.na(code)
    n_categories <- extra$ordinal_categories[[j]]
    probability <- exp(.latents_ordinal_log_probabilities(
      parameters$ordinal_intercepts[[j]], parameters$ordinal_locations[, j]))
    weight <- posteriors * observed
    intercepts <- vapply(seq_len(n_categories)[-1L], function(k) {
      rowSums(weight * (outer(code == k & observed, rep(1, n_profiles)) -
                          matrix(probability[, k], n_rows, n_profiles, byrow = TRUE)),
              na.rm = TRUE)
    }, numeric(n_rows))
    expected_score <- drop(probability %*% (seq_len(n_categories) - 1))
    locations <- vapply(seq_len(n_profiles - 1L), function(profile) {
      value <- ifelse(observed, (code - 1) - expected_score[profile], 0)
      weight[, profile] * value
    }, numeric(n_rows))
    cbind(matrix(intercepts, n_rows), matrix(locations, n_rows))
  })
  count <- lapply(seq_len(ncol(extra$count) %||% 0L), function(j) {
    value <- extra$count[, j]
    observed <- !is.na(value)
    means <- parameters$count_means[, j]
    residual <- outer(ifelse(observed, value, 0), rep(1, n_profiles)) -
      outer(observed * 1, means)
    posteriors * residual
  })
  rows <- do.call(cbind, c(ordinal, count))
  if (is.null(group_index)) matrix(colSums(rows), 1L) else
    rowsum(rows, group_index, reorder = FALSE)
}

#' Natural-by-unconstrained Jacobian of the ordinal and count block
#' @return A square diagonal matrix: one for ordinal coordinates, the mean for
#'   count ones.
#' @noRd
.latents_extra_jacobian <- function(parameters, extra, n_profiles) {
  natural <- .latents_extra_coordinates(parameters, extra, n_profiles, "natural")
  is_count <- startsWith(names(natural), "count_mean")
  diag(ifelse(is_count, natural, 1), nrow = length(natural))
}

#' The ordinal indicators as a tidy table
#'
#' @param x A fit carrying ordinal indicators.
#' @param data Optional fitting data, as for the other measurement tables.
#' @return One row per profile, ordinal indicator and category: `profile`,
#'   `indicator`, `category`, `probability`, the profile's `location` on that
#'   indicator (zero for the last, reference profile) and
#'   `location_standard_error`.
#' @noRd
.latents_ordinal_frame <- function(x, data = NULL) {
  extra <- x$extra_data
  empty <- data.frame(profile = integer(), indicator = character(),
                      category = character(), probability = numeric(),
                      location = numeric(), location_standard_error = numeric())
  if (is.null(extra$ordinal)) return(empty)
  n_profiles <- x$n_profiles
  frame <- do.call(rbind, lapply(seq_len(ncol(extra$ordinal)), function(j) {
    levels_for <- extra$ordinal_levels[[j]]
    probability <- exp(.latents_ordinal_log_probabilities(
      x$ordinal_intercepts[[j]], x$ordinal_locations[, j]))
    data.frame(profile = rep(seq_len(n_profiles), times = length(levels_for)),
               indicator = colnames(extra$ordinal)[[j]],
               category = rep(levels_for, each = n_profiles),
               probability = as.vector(probability),
               location = rep(unname(x$ordinal_locations[, j]), times = length(levels_for)),
               stringsAsFactors = FALSE)
  }))
  errors <- .multilpa_table_errors(x, data)
  frame$location_standard_error <- if (is.null(errors)) NA_real_ else
    .multilpa_match_error(errors, "ordinal_location", frame$profile, frame$indicator)
  # Indicators in the order the fit was given, then profile, then category.
  category_rank <- unlist(lapply(seq_len(ncol(extra$ordinal)), function(j) {
    rep(seq_along(extra$ordinal_levels[[j]]), each = n_profiles)
  }))
  frame <- frame[order(match(frame$indicator, colnames(extra$ordinal)),
                       frame$profile, category_rank), ]
  row.names(frame) <- NULL
  frame
}

#' The count indicators as a tidy table
#' @param x A fit carrying count indicators.
#' @param data Optional fitting data.
#' @return One row per profile and count indicator: `profile`, `indicator`,
#'   `mean` (the Poisson mean) and `standard_error`.
#' @noRd
.latents_count_frame <- function(x, data = NULL) {
  means <- x$count_means
  if (is.null(means)) {
    return(data.frame(profile = integer(), indicator = character(),
                      mean = numeric(), standard_error = numeric()))
  }
  frame <- data.frame(profile = rep(seq_len(nrow(means)), times = ncol(means)),
                      indicator = rep(colnames(means), each = nrow(means)),
                      mean = as.vector(means), stringsAsFactors = FALSE)
  errors <- .multilpa_table_errors(x, data)
  frame$standard_error <- if (is.null(errors)) NA_real_ else
    .multilpa_match_error(errors, "count_mean", frame$profile, frame$indicator)
  frame
}

#' The `get_results()` tables a fit with ordinal or count indicators adds
#' @noRd
.latents_extra_catalogue <- function(x) {
  if (is.null(x$extra_data)) return(list())
  c(if (!is.null(x$extra_data$ordinal))
      list(ordinal = .multilpa_table(function(x, data = NULL)
        .latents_ordinal_frame(x, data = data), "data")),
    if (!is.null(x$extra_data$count))
      list(count_means = .multilpa_table(function(x, data = NULL)
        .latents_count_frame(x, data = data), "data")))
}

#' Draw ordinal and count indicators from a fit, given each row's profile
#'
#' Ordinal indicators come back as ordered factors with the fit's levels, so a
#' refit re-encodes the same categories even when one is not drawn.
#'
#' @param object A fit carrying the blocks and `extra_data`.
#' @param profile Profile of each row.
#' @return A named list of columns, empty when the fit has none.
#' @noRd
.latents_draw_extra <- function(object, profile) {
  extra <- object$extra_data
  if (is.null(extra)) return(list())
  ordinal <- lapply(seq_len(ncol(extra$ordinal) %||% 0L), function(j) {
    levels_for <- extra$ordinal_levels[[j]]
    probability <- exp(.latents_ordinal_log_probabilities(
      object$ordinal_intercepts[[j]], object$ordinal_locations[, j]))
    factor(levels_for[.multilpa_draw_rows(probability[profile, , drop = FALSE])],
           levels = levels_for, ordered = TRUE)
  })
  names(ordinal) <- colnames(extra$ordinal)
  count <- lapply(seq_len(ncol(extra$count) %||% 0L), function(j) {
    as.numeric(stats::rpois(length(profile), object$count_means[profile, j]))
  })
  names(count) <- colnames(extra$count)
  c(ordinal, count)
}

#' Encode new rows' ordinal and count indicators as a fit encoded its own
#'
#' @param object A fit carrying `extra_data`.
#' @param newdata A data frame with the fit's indicator columns.
#' @return Prepared indicators on the fit's categories, or `NULL`. A category
#'   the fit never saw, or an invalid count, raises `latents_bad_data`.
#' @noRd
.latents_prepare_extra_like <- function(object, newdata) {
  extra <- object$extra_data
  if (is.null(extra)) return(NULL)
  n <- nrow(newdata)
  ordinal <- if (is.null(extra$ordinal)) NULL else
    matrix(vapply(colnames(extra$ordinal), function(name) {
      value <- newdata[[name]]
      code <- match(as.character(value), extra$ordinal_levels[[name]])
      unseen <- !is.na(value) & is.na(code)
      if (any(unseen)) {
        stop(errorCondition(sprintf(
          "`%s` has the categor%s %s, which the fit never saw.", name,
          if (length(unique(value[unseen])) > 1L) "ies" else "y",
          paste(sprintf("`%s`", unique(as.character(value[unseen]))), collapse = ", ")),
          class = "latents_bad_data", call = NULL))
      }
      as.integer(code)
    }, integer(n)), n, dimnames = list(NULL, colnames(extra$ordinal)))
  count <- if (is.null(extra$count)) NULL else
    matrix(vapply(colnames(extra$count), function(name) {
      value <- newdata[[name]]
      observed <- value[!is.na(value)]
      if (!is.numeric(value) || any(!is.finite(observed)) || any(observed < 0) ||
          any(observed != round(observed))) {
        stop(errorCondition(sprintf(
          "Count indicator `%s` must hold non-negative whole numbers.", name),
          class = "latents_bad_data", call = NULL))
      }
      as.numeric(value)
    }, numeric(n)), n, dimnames = list(NULL, colnames(extra$count)))
  list(ordinal = ordinal, ordinal_levels = extra$ordinal_levels,
       ordinal_categories = extra$ordinal_categories, count = count)
}

#' Print the ordinal and count indicators of a fit
#' @return `NULL`, invisibly; one line naming them and their tables.
#' @noRd
.latents_print_extra <- function(x) {
  pieces <- c(
    if (length(x$ordinal) > 0L) sprintf(
      "ordinal %s (adjacent-category logit; get_results(x, \"ordinal\"))",
      paste(x$ordinal, collapse = ", ")),
    if (length(x$count) > 0L) sprintf(
      "count %s (Poisson; get_results(x, \"count_means\"))",
      paste(x$count, collapse = ", ")))
  if (length(pieces) > 0L) cat(sprintf("Also %s\n", paste(pieces, collapse = "; ")))
  invisible(NULL)
}
