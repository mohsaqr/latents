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
#   count     Poisson, one mean per profile; or (count_model =
#             "negative_binomial") NB2 with mean mu and variance mu + alpha mu^2,
#             the dispersion alpha per profile or shared (Latent GOLD's
#             `poisson overdispersed`).
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
    observed <- value[!is.na(value)]
    if (!is.null(dim(value)) || !(is.factor(value) || is.numeric(value)) ||
        (is.numeric(value) && (any(!is.finite(observed)) ||
                               any(observed != round(observed))))) {
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
      if (!is.numeric(value) || !is.null(dim(value)) ||
          any(!is.finite(observed)) || any(observed < 0) ||
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
    contribution[observed, ] <- if (.latents_negative_binomial(extra)) {
      dispersion <- parameters$count_dispersion[, j]
      vapply(seq_len(n_profiles), function(profile) {
        stats::dnbinom(value[observed], size = 1 / dispersion[profile],
                       mu = means[profile], log = TRUE)
      }, numeric(sum(observed)))
    } else {
      vapply(means, function(mean) {
        stats::dpois(value[observed], lambda = mean, log = TRUE)
      }, numeric(sum(observed)))
    }
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
      # Rounding alone can lower the objective by a hair at the maximum; such
      # a step is accepted rather than halved forty times.
      accepted <- is.finite(candidate_value) &&
        candidate_value >= value - 1e-12 * (1 + abs(value))
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
  count_dispersion <- NULL
  if (.latents_negative_binomial(extra)) {
    # Negative binomial: the mean and the dispersion are solved together, from
    # the current point (or, at the start, the Poisson means and a moment
    # dispersion).
    solved <- lapply(seq_len(ncol(extra$count)), function(j) {
      start_means <- previous$count_means[, j] %||% count_means[, j]
      start_dispersion <- previous$count_dispersion[, j] %||%
        .latents_moment_dispersion(extra$count[, j], posteriors, count_means[, j])
      .latents_negative_binomial_maximize(extra$count[, j], posteriors, start_means,
                                          start_dispersion,
                                          identical(extra$count_dispersion, "equal"))
    })
    count_means <- matrix(vapply(solved, `[[`, numeric(n_profiles), "means"),
                          n_profiles, dimnames = list(NULL, colnames(extra$count)))
    count_dispersion <- matrix(vapply(solved, `[[`, numeric(n_profiles), "dispersion"),
                               n_profiles, dimnames = list(NULL, colnames(extra$count)))
  }
  list(ordinal_intercepts = if (length(ordinal) == 0L) NULL else
         stats::setNames(lapply(ordinal, `[[`, "intercepts"), colnames(extra$ordinal)),
       ordinal_locations = if (length(ordinal) == 0L) NULL else
         matrix(vapply(ordinal, `[[`, numeric(n_profiles), "locations"),
                n_profiles, length(ordinal),
                dimnames = list(NULL, colnames(extra$ordinal))),
       count_means = count_means, count_dispersion = count_dispersion)
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
               (ncol(extra$count) %||% 0L) * n_profiles +
               (ncol(extra$count) %||% 0L) * .latents_dispersion_rows(extra, n_profiles))
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
  if (!is.null(parameters$count_dispersion)) {
    dimnames(parameters$count_dispersion) <- list(profile_names, colnames(extra$count))
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
  natural <- identical(scale, "natural")
  count <- lapply(seq_along(names_of$count), function(j) {
    indicator <- names_of$count[[j]]
    means <- unname(parameters$count_means[, j])
    mean_values <- stats::setNames(
      if (natural) means else log(means),
      sprintf(if (natural) "count_mean[%d,%s]" else "log_count_mean[%d,%s]",
              seq_len(n_profiles), indicator))
    if (!.latents_negative_binomial(extra)) return(mean_values)
    rows <- .latents_dispersion_rows(extra, n_profiles)
    dispersion <- unname(parameters$count_dispersion[seq_len(rows), j])
    outcome <- if (identical(extra$count_dispersion, "equal")) "shared" else
      as.character(seq_len(rows))
    c(mean_values, stats::setNames(
      if (natural) dispersion else log(dispersion),
      sprintf(if (natural) "count_dispersion[%s,%s]" else "log_count_dispersion[%s,%s]",
              outcome, indicator)))
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
  negative_binomial <- .latents_negative_binomial(extra)
  rows <- .latents_dispersion_rows(extra, n_profiles)
  count <- lapply(seq_len(ncol(extra$count) %||% 0L), function(j) {
    list(means = exp(take(n_profiles)),
         dispersion = if (negative_binomial) rep_len(exp(take(rows)), n_profiles))
  })
  stopifnot("the extra coordinates must be used exactly" = at == length(theta))
  list(ordinal_intercepts = if (length(ordinal) == 0L) NULL else
         lapply(ordinal, `[[`, "intercepts"),
       ordinal_locations = if (length(ordinal) == 0L) NULL else
         matrix(vapply(ordinal, `[[`, numeric(n_profiles), "locations"), n_profiles),
       count_means = if (length(count) == 0L) NULL else
         matrix(vapply(count, `[[`, numeric(n_profiles), "means"), n_profiles),
       count_dispersion = if (length(count) == 0L || !negative_binomial) NULL else
         matrix(vapply(count, `[[`, numeric(n_profiles), "dispersion"), n_profiles))
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
    if (!.latents_negative_binomial(extra)) {
      residual <- outer(ifelse(observed, value, 0), rep(1, n_profiles)) -
        outer(observed * 1, means)
      return(posteriors * residual)
    }
    parts <- .latents_negative_binomial_scores(ifelse(observed, value, 0), means,
                                               parameters$count_dispersion[, j])
    mean_scores <- posteriors * observed * parts$log_mean
    dispersion_scores <- posteriors * observed * parts$log_dispersion
    if (identical(extra$count_dispersion, "equal")) {
      dispersion_scores <- matrix(rowSums(dispersion_scores), n_rows)
    }
    cbind(mean_scores, dispersion_scores)
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
  is_count <- startsWith(names(natural), "count_mean") |
    startsWith(names(natural), "count_dispersion")
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
.latents_ordinal_frame <- function(x, data = NULL, standard_errors = TRUE) {
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
  if (isTRUE(standard_errors)) {
    errors <- .multilpa_table_errors(x, data)
    frame$location_standard_error <- if (is.null(errors)) NA_real_ else
      .multilpa_match_error(errors, "ordinal_location", frame$profile, frame$indicator)
  }
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
.latents_count_frame <- function(x, data = NULL, standard_errors = TRUE) {
  means <- x$count_means
  if (is.null(means)) {
    return(data.frame(profile = integer(), indicator = character(),
                      mean = numeric(), standard_error = numeric()))
  }
  frame <- data.frame(profile = rep(seq_len(nrow(means)), times = ncol(means)),
                      indicator = rep(colnames(means), each = nrow(means)),
                      mean = as.vector(means), stringsAsFactors = FALSE)
  negative_binomial <- !is.null(x$count_dispersion)
  if (negative_binomial) frame$dispersion <- as.vector(x$count_dispersion)
  if (!isTRUE(standard_errors)) return(frame)
  errors <- .multilpa_table_errors(x, data)
  frame$standard_error <- if (is.null(errors)) NA_real_ else
    .multilpa_match_error(errors, "count_mean", frame$profile, frame$indicator)
  if (negative_binomial) {
    # A shared dispersion is reported once, under outcome `shared`.
    shared <- identical(x$extra_data$count_dispersion, "equal")
    frame$dispersion_standard_error <- if (is.null(errors)) NA_real_ else
      .multilpa_match_error(errors, "count_dispersion",
                            if (shared) rep("shared", nrow(frame)) else frame$profile,
                            frame$indicator)
  }
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
    factor(levels_for[.latents_draw_rows(probability[profile, , drop = FALSE])],
           levels = levels_for, ordered = TRUE)
  })
  names(ordinal) <- colnames(extra$ordinal)
  count <- lapply(seq_len(ncol(extra$count) %||% 0L), function(j) {
    means <- object$count_means[profile, j]
    as.numeric(if (.latents_negative_binomial(extra)) {
      stats::rnbinom(length(profile), size = 1 / object$count_dispersion[profile, j],
                     mu = means)
    } else stats::rpois(length(profile), means))
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
      if (!is.null(dim(value))) {
        stop(errorCondition(sprintf(
          "Ordinal indicator `%s` must be a vector of categories.", name),
          class = "latents_bad_data", call = NULL))
      }
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
      if (!is.numeric(value) || !is.null(dim(value)) ||
          any(!is.finite(observed)) || any(observed < 0) ||
          any(observed != round(observed))) {
        stop(errorCondition(sprintf(
          "Count indicator `%s` must hold non-negative whole numbers.", name),
          class = "latents_bad_data", call = NULL))
      }
      as.numeric(value)
    }, numeric(n)), n, dimnames = list(NULL, colnames(extra$count)))
  list(ordinal = ordinal, ordinal_levels = extra$ordinal_levels,
       ordinal_categories = extra$ordinal_categories, count = count,
       count_model = extra$count_model, count_dispersion = extra$count_dispersion)
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
      "count %s (%s; get_results(x, \"count_means\"))",
      paste(x$count, collapse = ", "),
      if (.latents_negative_binomial(x$extra_data)) sprintf(
        "negative binomial, %s dispersion", x$extra_data$count_dispersion) else "Poisson"))
  if (length(pieces) > 0L) cat(sprintf("Also %s\n", paste(pieces, collapse = "; ")))
  invisible(NULL)
}

#' The stored ordinal and count indicators as data columns
#' @param extra Prepared indicators, or `NULL`.
#' @return A named list of columns: ordered factors with the fit's levels for
#'   ordinal indicators, numbers for counts.
#' @noRd
.latents_draw_extra_columns <- function(extra) {
  if (is.null(extra)) return(list())
  ordinal <- lapply(seq_len(ncol(extra$ordinal) %||% 0L), function(j) {
    levels_for <- extra$ordinal_levels[[j]]
    factor(levels_for[extra$ordinal[, j]], levels = levels_for, ordered = TRUE)
  })
  names(ordinal) <- colnames(extra$ordinal)
  count <- lapply(seq_len(ncol(extra$count) %||% 0L), function(j) extra$count[, j])
  names(count) <- colnames(extra$count)
  c(ordinal, count)
}

#' Tidy labels of the ordinal and count coordinates
#'
#' One row per coordinate, in `.latents_extra_coordinates()` order, with the
#' level, outcome, term and parameter columns every inference table uses.
#'
#' @param object A fit carrying `extra_data`.
#' @param profiles Profile labels.
#' @param free `TRUE` names count means on the log scale they are estimated on.
#' @return A data frame, with no rows when the fit has neither type.
#' @noRd
.latents_extra_labels <- function(object, profiles, free = FALSE) {
  extra <- object$extra_data
  row <- function(outcome, term, parameter) {
    data.frame(level = rep("measurement", length(term)), outcome = outcome,
               term = term, parameter = rep(parameter, length(term)),
               stringsAsFactors = FALSE)
  }
  empty <- row(character(), character(), character())
  if (is.null(extra)) return(empty)
  ordinal <- lapply(seq_len(ncol(extra$ordinal) %||% 0L), function(j) {
    indicator <- colnames(extra$ordinal)[[j]]
    categories <- extra$ordinal_levels[[j]][-1L]
    located <- profiles[-length(profiles)]
    rbind(row(rep("shared", length(categories)),
              sprintf("%s:%s", indicator, categories), "ordinal_intercept"),
          row(located, rep(indicator, length(located)), "ordinal_location"))
  })
  dispersion_outcomes <- if (identical(extra$count_dispersion, "equal")) "shared" else
    profiles
  count <- lapply(colnames(extra$count) %||% character(), function(indicator) {
    rbind(row(profiles, rep(indicator, length(profiles)),
              if (free) "log_count_mean" else "count_mean"),
          if (.latents_negative_binomial(extra))
            row(dispersion_outcomes, rep(indicator, length(dispersion_outcomes)),
                if (free) "log_count_dispersion" else "count_dispersion"))
  })
  do.call(rbind, c(list(empty), ordinal, count))
}

#' The rows of prepared ordinal and count indicators
#' @param extra Prepared indicators, or `NULL`.
#' @param rows Row indices to keep.
#' @return The same structure on those rows only.
#' @noRd
.latents_extra_rows <- function(extra, rows) {
  if (is.null(extra)) return(NULL)
  if (!is.null(extra$ordinal)) extra$ordinal <- extra$ordinal[rows, , drop = FALSE]
  if (!is.null(extra$count)) extra$count <- extra$count[rows, , drop = FALSE]
  extra
}

#' Ordinal or count table of a transition fit, one block per measurement
#'
#' As the transition `profiles` table: estimates by profile, with an
#' `occasion` column that is `NA` when measurement is invariant.
#'
#' @param fit A `multilpa_lta` fit.
#' @param what `"ordinal"` or `"count_means"`.
#' @return A data frame.
#' @noRd
.lta_extra_table <- function(fit, what) {
  build <- if (identical(what, "ordinal")) .latents_ordinal_frame else
    .latents_count_frame
  do.call(rbind, lapply(seq_along(fit$measurement), function(b) {
    block <- c(fit$measurement[[b]],
               list(n_profiles = fit$n_profiles, extra_data = fit$extra_data))
    frame <- build(block, standard_errors = FALSE)
    cbind(occasion = rep(if (length(fit$measurement) == 1L) NA_integer_ else b,
                         nrow(frame)), frame)
  }))
}

#' Put drawn values into the type of the column they replace
#'
#' An ordinal draw is an ordered factor on the fit's levels; a column held as
#' numbers gets the numbers back, a factor column the labels.
#'
#' @param drawn Drawn values.
#' @param column The original column.
#' @return `drawn` in `column`'s type.
#' @noRd
.latents_as_column_type <- function(drawn, column) {
  if (is.factor(drawn) && is.numeric(column)) return(as.numeric(as.character(drawn)))
  if (is.factor(drawn) && is.factor(column)) {
    return(factor(as.character(drawn), levels = levels(column),
                  ordered = is.ordered(column)))
  }
  drawn
}

#' Does a fit model its counts as negative binomial?
#' @noRd
.latents_negative_binomial <- function(extra) {
  identical(extra$count_model, "negative_binomial")
}

#' How many dispersion rows each negative-binomial count indicator has
#' @return `0` for Poisson, `1` when shared, the number of profiles otherwise.
#' @noRd
.latents_dispersion_rows <- function(extra, n_profiles) {
  if (!.latents_negative_binomial(extra)) return(0L)
  if (identical(extra$count_dispersion, "equal")) 1L else as.integer(n_profiles)
}

#' Attach the count model to prepared indicators
#' @param extra Prepared indicators, or `NULL`.
#' @param count_model `"poisson"` or `"negative_binomial"`.
#' @param count_dispersion `"varying"` or `"equal"`.
#' @return `extra` with the two fields set.
#' @noRd
.latents_set_count_model <- function(extra, count_model = "poisson",
                                     count_dispersion = "varying") {
  if (is.null(extra)) return(NULL)
  extra$count_model <- count_model
  extra$count_dispersion <- count_dispersion
  extra
}

#' The negative-binomial dispersion floor; reaching it is the Poisson limit
#' @noRd
.latents_min_dispersion <- 1e-8

#' Per-row scores of a negative-binomial count on its log mean and log
#' dispersion, one column per profile
#'
#' With size r = 1 / alpha: d log f / d log mu = r (y - mu) / (r + mu), and
#' d log f / d log alpha = -r [digamma(y + r) - digamma(r) + log(r / (r + mu))
#' + (mu - y) / (r + mu)].
#'
#' @param value Counts (unobserved ones set to anything; they are masked by
#'   the caller).
#' @param means,dispersion One value per profile.
#' @return A list of two n x C matrices.
#' @noRd
.latents_negative_binomial_scores <- function(value, means, dispersion) {
  size <- 1 / dispersion
  n_profiles <- length(means)
  y <- matrix(value, length(value), n_profiles)
  mu <- matrix(means, length(value), n_profiles, byrow = TRUE)
  r <- matrix(size, length(value), n_profiles, byrow = TRUE)
  list(log_mean = r * (y - mu) / (r + mu),
       log_dispersion = .latents_negative_binomial_dispersion_derivatives(
         value, means, dispersion)$score)
}

#' Stable negative-binomial log-dispersion derivatives
#'
#' `means` is one value per profile, or a rows-by-classes matrix of means.
#' Near the Poisson limit, differences of digamma/trigamma terms lose the
#' small dispersion score to cancellation. Expand the NB log density relative
#' to Poisson in alpha through order four when alpha times the count or mean
#' is below 0.001. The coefficients follow from expanding
#' sum_j log(1 + j alpha) - (y + 1/alpha) log(1 + mu alpha).
#' Both expansion variables (alpha y and alpha mu) are then below 0.001.
#' Above that switch, large size parameters use the digamma asymptotic
#' expansion, combining its logarithm with the likelihood terms before
#' scaling. Subtracting two digamma/trigamma values there loses precision.
#' @noRd
.latents_negative_binomial_dispersion_derivatives <- function(value, means, dispersion) {
  # One mean per profile, or (for a regression) a rows-by-classes matrix.
  n_profiles <- if (is.matrix(means)) ncol(means) else length(means)
  y <- matrix(value, length(value), n_profiles)
  mu <- if (is.matrix(means)) means else
    matrix(means, length(value), n_profiles, byrow = TRUE)
  alpha <- matrix(dispersion, length(value), n_profiles, byrow = TRUE)
  r <- 1 / alpha
  d <- digamma(y + r) - digamma(r) - log1p(mu / r) + (mu - y) / (r + mu)
  d_prime <- trigamma(y + r) - trigamma(r) + mu / (r * (r + mu)) -
    (mu - y) / (r + mu)^2
  score <- -r * d
  curvature <- r * d + r^2 * d_prime
  small <- alpha * pmax(y, mu) < 1e-3
  if (any(small)) {
    count <- y[small]
    mean <- mu[small]
    a <- alpha[small]
    sums <- list(count * (count - 1) / 2,
                 count * (count - 1) * (2 * count - 1) / 6,
                 (count * (count - 1) / 2)^2,
                 count * (count - 1) * (2 * count - 1) *
                   (3 * count^2 - 3 * count - 1) / 30)
    terms <- lapply(seq_len(4L), function(order) {
      (-1)^(order + 1) * a^order *
        ((sums[[order]] - count * mean^order) / order + mean^(order + 1) / (order + 1))
    })
    score[small] <- Reduce(`+`, Map(`*`, terms, seq_len(4L)))
    curvature[small] <- Reduce(`+`, Map(`*`, terms, seq_len(4L)^2))
  }
  large <- r >= 1000 & !small
  if (any(large)) {
    size <- r[large]
    count <- y[large]
    mean <- mu[large]
    difference <- (count - mean) / (size + mean)
    log_remainder <- log1p(count / size) - log1p(mean / size) - difference
    near <- abs(difference) < 0.01
    if (any(near)) {
      # log(1 + z) - z without cancellation around zero.
      z <- difference[near]
      log_remainder[near] <- Reduce(`+`, lapply(2:12, function(order) {
        (-1)^(order + 1L) * z^order / order
      }))
    }
    score_large <- -size * log_remainder
    curvature_large <- size * (log_remainder +
      ((size / (size + count)) * difference) * difference)
    log_ratio <- log1p(count / size)
    fraction <- count / (size + count)
    powers <- c(1L, 2L, 4L, 6L, 8L)
    coefficients <- c(0.5, 1 / 12, -1 / 120, 1 / 252, -1 / 240)
    for (index in seq_along(powers)) {
      power <- powers[index]
      gap <- -expm1(-power * log_ratio)
      coefficient <- coefficients[index] * size^(1L - power)
      score_large <- score_large - coefficient * gap
      curvature_large <- curvature_large + coefficient *
        ((1L - power) * gap - power * fraction * exp(-power * log_ratio))
    }
    score[large] <- score_large
    curvature[large] <- curvature_large
  }
  list(score = score, curvature = curvature)
}

#' A moment start for the negative-binomial dispersion, per profile
#' @noRd
.latents_moment_dispersion <- function(value, posteriors, means) {
  observed <- !is.na(value)
  y <- ifelse(observed, value, 0)
  weight <- posteriors * observed
  variance <- colSums(weight * outer(y, means, "-")^2) /
    pmax(colSums(weight), .Machine$double.xmin)
  pmin(pmax((variance - means) / means^2, 0.05), 10)
}

#' Maximize one negative-binomial count indicator's expected log likelihood
#'
#' Modified Newton-Raphson on the log means and log dispersions: the analytic
#' Hessian with its eigenvalues reflected, so every direction ascends, and step
#' halving so every accepted step raises the objective (to rounding). A
#' dispersion is held at or above `.latents_min_dispersion`, the Poisson limit.
#'
#' @param value The indicator's counts (`NA` where unobserved).
#' @param posteriors n x C posterior weights.
#' @param means,dispersion The current point, one value per profile.
#' @param shared Whether one dispersion is shared by every profile.
#' @param max_newton Most Newton steps.
#' @return A list with `means` and `dispersion` (one value per profile).
#' @noRd
.latents_negative_binomial_maximize <- function(value, posteriors, means, dispersion,
                                                shared = FALSE, max_newton = 100L) {
  observed <- !is.na(value)
  y <- value[observed]
  weight <- posteriors[observed, , drop = FALSE]
  n_profiles <- ncol(posteriors)
  n_dispersion <- if (shared) 1L else n_profiles
  floor <- log(.latents_min_dispersion)
  unpack <- function(theta) {
    list(means = exp(theta[seq_len(n_profiles)]),
         dispersion = rep_len(exp(theta[n_profiles + seq_len(n_dispersion)]), n_profiles))
  }
  objective <- function(theta) {
    point <- unpack(theta)
    sum(weight * vapply(seq_len(n_profiles), function(profile) {
      stats::dnbinom(y, size = 1 / point$dispersion[profile], mu = point$means[profile],
                     log = TRUE)
    }, numeric(length(y))))
  }
  gradient <- function(theta) {
    point <- unpack(theta)
    parts <- .latents_negative_binomial_scores(y, point$means, point$dispersion)
    dispersion_gradient <- colSums(weight * parts$log_dispersion)
    c(colSums(weight * parts$log_mean),
      if (shared) sum(dispersion_gradient) else dispersion_gradient)
  }
  theta <- c(log(means), log(pmax(dispersion[seq_len(n_dispersion)],
                                  .latents_min_dispersion)))
  value_now <- objective(theta)
  iteration <- 0L
  # Newton steps refine one estimate in sequence; nothing to vectorize.
  while (iteration < max_newton) {
    slope <- gradient(theta)
    at_floor <- n_profiles + which(theta[n_profiles + seq_len(n_dispersion)] <= floor + 1e-12)
    # A dispersion on its floor whose slope points below it stays there.
    held <- at_floor[slope[at_floor] < 0]
    free <- setdiff(seq_along(theta), held)
    if (max(abs(slope[free])) <= 1e-10 * (1 + sum(weight))) break
    curvature <- -.latents_negative_binomial_hessian(unpack(theta), y, weight,
                                                     shared)[free, free, drop = FALSE]
    # Near the Poisson limit the objective is convex in a log dispersion, so
    # the curvature is indefinite there. A plain Newton step then descends in
    # that coordinate while still ascending overall through the means, and
    # stalls short of the maximum. Reflecting the curvature's eigenvalues
    # (modified Newton) ascends along every eigen-direction; the step is
    # capped because a nearly flat direction would otherwise be enormous.
    spectrum <- eigen((curvature + t(curvature)) / 2, symmetric = TRUE)
    eigenvalues <- pmax(abs(spectrum$values), 1e-10 * (1 + max(abs(spectrum$values))))
    direction <- drop(spectrum$vectors %*%
                        (crossprod(spectrum$vectors, slope[free]) / eigenvalues))
    direction <- direction * min(1, 5 / max(abs(direction)))
    step <- 1
    accepted <- FALSE
    while (step > 1e-12 && !accepted) {
      candidate <- theta
      candidate[free] <- theta[free] + step * direction
      candidate[n_profiles + seq_len(n_dispersion)] <-
        pmax(candidate[n_profiles + seq_len(n_dispersion)], floor)
      candidate_value <- objective(candidate)
      # Rounding alone can lower the objective by a hair at the maximum; such
      # a step is accepted rather than halved forty times.
      accepted <- is.finite(candidate_value) &&
        candidate_value >= value_now - 1e-12 * (1 + abs(value_now))
      if (!accepted) step <- step / 2
    }
    if (!accepted) break
    theta <- candidate
    value_now <- candidate_value
    iteration <- iteration + 1L
  }
  .latents_poisson_limit(unpack(theta), y, weight, shared)
}

#' Put a negative-binomial dispersion on its boundary when the data ask for it
#'
#' The likelihood flattens as the dispersion goes to zero, so Newton steps
#' stop short of the boundary rather than reaching it. At alpha = 0 the score
#' of the dispersion is `sum w [(y - mu)^2 - y] / 2` (the overdispersion score
#' test); when it is not positive the maximum is on the boundary (the Poisson
#' limit, Karush-Kuhn-Tucker), so the dispersion is set to the floor and the
#' mean to its Poisson maximum, the weighted mean count.
#'
#' @param point A list with `means` and `dispersion`, one value per profile.
#' @param y Observed counts; `weight` their n x C posterior weights.
#' @param shared Whether one dispersion is shared by every profile.
#' @return `point`, with boundary dispersions on the floor.
#' @noRd
.latents_poisson_limit <- function(point, y, weight, shared) {
  poisson_means <- colSums(weight * y) / pmax(colSums(weight), .Machine$double.xmin)
  score_at_zero <- colSums(weight * (outer(y, poisson_means, "-")^2 - y)) / 2
  boundary <- if (shared) rep(sum(score_at_zero) <= 0, length(poisson_means)) else
    score_at_zero <= 0
  point$dispersion[boundary] <- .latents_min_dispersion
  point$means[boundary] <- pmax(poisson_means[boundary], 1e-10)
  point
}

#' Warn when a negative-binomial dispersion is on its boundary
#' @param dispersion The fitted dispersion matrix, or `NULL`.
#' @return `NULL`, invisibly; a `latents_boundary` warning when any
#'   dispersion is at the floor.
#' @noRd
.latents_warn_poisson_limit <- function(dispersion) {
  if (length(dispersion) == 0L ||
      !any(dispersion <= .latents_min_dispersion * (1 + 1e-6))) {
    return(invisible(NULL))
  }
  warning(warningCondition(paste(
    "A negative-binomial dispersion is estimated at zero: those counts show no",
    "overdispersion within the profile, which is the Poisson limit. This is a",
    "boundary fit; consider `count_model = \"poisson\"`."),
    class = "latents_boundary", call = NULL))
  invisible(NULL)
}

#' Hessian of a negative-binomial indicator's expected log likelihood
#'
#' On the log means and log dispersions, from the closed-form second
#' derivatives of the NB2 log density (r = 1 / alpha):
#' d2/dlog mu2 = -r mu (r + y) / (r + mu)^2,
#' d2/dlog mu dlog alpha = -r mu (y - mu) / (r + mu)^2,
#' d2/dlog alpha2 = r D + r^2 dD/dr, with
#' D = digamma(y + r) - digamma(r) + log(r / (r + mu)) + (mu - y) / (r + mu) and
#' dD/dr = trigamma(y + r) - trigamma(r) + 1 / r - 1 / (r + mu) -
#' (mu - y) / (r + mu)^2.
#'
#' @param point A list with `means` and `dispersion`, one value per profile.
#' @param y Observed counts; `weight` their n x C posterior weights.
#' @param shared Whether one dispersion is shared by every profile.
#' @return A square matrix over the log means, then the log dispersion(s).
#' @noRd
.latents_negative_binomial_hessian <- function(point, y, weight, shared) {
  n_profiles <- ncol(weight)
  dispersion_derivatives <- .latents_negative_binomial_dispersion_derivatives(
    y, point$means, point$dispersion)
  terms <- lapply(seq_len(n_profiles), function(profile) {
    mu <- point$means[profile]
    r <- 1 / point$dispersion[profile]
    w <- weight[, profile]
    spread <- (r + mu)^2
    c(mean_mean = sum(w * -r * mu * (r + y) / spread),
      mean_dispersion = sum(w * -r * mu * (y - mu) / spread),
      dispersion_dispersion = sum(w * dispersion_derivatives$curvature[, profile]))
  })
  n_dispersion <- if (shared) 1L else n_profiles
  hessian <- matrix(0, n_profiles + n_dispersion, n_profiles + n_dispersion)
  invisible(lapply(seq_len(n_profiles), function(profile) {
    dispersion_at <- n_profiles + if (shared) 1L else profile
    part <- terms[[profile]]
    hessian[profile, profile] <<- part[["mean_mean"]]
    hessian[profile, dispersion_at] <<- part[["mean_dispersion"]]
    hessian[dispersion_at, profile] <<- part[["mean_dispersion"]]
    hessian[dispersion_at, dispersion_at] <<- hessian[dispersion_at, dispersion_at] +
      part[["dispersion_dispersion"]]
  }))
  hessian
}
