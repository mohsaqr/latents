#' Bivariate residuals for every pair of indicators
#'
#' The statistic per pair kind, what the test assumes and the references are
#' documented on `?get_results`, which is where a caller reaches this table from.
#'
#' @param x A fitted model of this package.
#' @param data The data the model was fitted to; `NULL` uses what it carries.
#' @param by `"profile"` assesses each profile separately, `"overall"` pools.
#' @param adjust The multiplicity correction across the indicator pairs.
#' @return A base `data.frame`, one row per indicator pair and, under
#'   `by = "profile"`, per profile.
#' @noRd
.multilpa_bivariate_residuals <- function(x, data = NULL, by = c("profile", "overall"),
                               adjust = .multilpa_p_adjust_methods) {
  stopifnot("`x` must be a fitted model of this package" = .multilpa_any_fit(x))
  .multilpa_refuse_noise(x, "Bivariate residuals")
  data <- .multilpa_resolve_data(x, data)
  stopifnot("`data` must be a data frame" = is.data.frame(data))
  .multilpa_check_row_count(x, data)
  by <- match.arg(by)
  adjust <- match.arg(adjust)
  continuous <- .multilpa_continuous_names(x)
  categorical <- x$categorical %||% character()
  extra_names <- .latents_extra_names(x$extra_data)
  stopifnot("`data` must contain the fitted indicators" =
              all(c(continuous, categorical, unlist(extra_names)) %in% names(data)))
  stopifnot("Gaussian indicators must be numeric and finite when observed" =
    all(vapply(data[continuous], function(value)
      is.numeric(value) && all(is.na(value) | is.finite(value)), logical(1))))
  # Every row is weighted by the posterior the fit holds for that position, so
  # a frame in any other order pairs each observation with someone else's
  # posterior and reports a residual for a table that was never observed.
  .multilpa_check_alignment(x, data)
  weights <- if (identical(by, "profile")) {
    stats::setNames(lapply(seq_len(x$n_profiles),
                           function(k) x$subject_posteriors[, k]),
                    paste0("profile_", seq_len(x$n_profiles)))
  } else list(overall = rep(1, x$n_observations))
  sampling <- .multilpa_level_weights(x, "individuals", x$n_observations)
  if (!is.null(sampling)) weights <- lapply(weights, `*`, sampling)

  # Categorical and ordinal indicators are both tables of categories; every
  # pair of them gets the Pearson table. Pairs with a count, or of a Gaussian
  # and a discrete indicator, get a residual correlation without a test
  # (review 2026-10-26 round 2, C2-3), so no pair of indicators is left out.
  discrete <- .multilpa_discrete_items(x, data)
  values <- .multilpa_pair_values(x, data, continuous, discrete)
  rows <- lapply(names(weights), function(label) {
    weight <- weights[[label]]
    rbind(
      .multilpa_gaussian_residuals(x, data, continuous, weight, label, by),
      .multilpa_categorical_residuals(x, discrete, weight, label),
      .multilpa_mixed_residuals(x, values, weight, label, by))
  })
  result <- do.call(rbind, rows)
  if (is.null(result) || nrow(result) == 0L) {
    return(data.frame(profile = character(), indicator_1 = character(),
                      indicator_2 = character(), kind = character(),
                      observed = numeric(), expected = numeric(),
                      residual = numeric(), effective_n = numeric(),
                      statistic = numeric(), df = integer(),
                      p_value = numeric(), p_adjusted = numeric()))
  }
  # Worst pair first. The indicator names and the profile label are explicit
  # secondary keys, so ties do not fall back on construction order.
  result <- result[order(-abs(result$residual), result$indicator_1,
                         result$indicator_2, result$profile), , drop = FALSE]
  # Every pair of indicators is a separate test, so a fit with five indicators
  # asks ten questions at `by = "overall"` and ten per profile otherwise. The
  # family is the verb's to declare, not the reader's to reconstruct, so the
  # adjusted column travels beside the raw one.
  result$p_adjusted <- stats::p.adjust(result$p_value, method = adjust)
  row.names(result) <- NULL
  result
}

#' Cramer's V of a weighted two-way table
#'
#' The association in a contingency table on the zero-to-one scale, comparable
#' with a correlation magnitude and zero for a table that is a product of its
#' own margins.
#'
#' @param counts A numeric matrix of (possibly fractional) cell counts.
#' @return A single numeric, or `NA_real_` when the table has no weight or has
#'   a dimension of one and the association is undefined.
#' @noRd
.multilpa_cramers_v <- function(counts) {
  stopifnot("`counts` must be a numeric matrix" =
              is.matrix(counts) && is.numeric(counts))
  total <- sum(counts)
  smaller <- min(nrow(counts), ncol(counts)) - 1L
  if (!is.finite(total) || total <= 0 || smaller < 1L) return(NA_real_)
  independent <- outer(rowSums(counts), colSums(counts)) / total
  # A zero cell of `independent` sits in a zero margin, so its `counts` cell is
  # zero too and contributes nothing; masking it avoids a 0/0.
  positive <- independent > 0
  chi_square <- sum((counts[positive] - independent[positive])^2 /
                      independent[positive])
  sqrt(chi_square / (total * smaller))
}

#' Residual correlation for each Gaussian pair
#' @param weight Posterior weights for the profile, or ones when pooling.
#' @param label The profile label to report.
#' @param by Whether the weights are a profile's or pooled.
#' @return One row per Gaussian pair, or `NULL` when there are fewer than two.
#' @noRd
.multilpa_gaussian_residuals <- function(object, data, continuous, weight, label,
                                         by) {
  if (length(continuous) < 2L) return(NULL)
  # On the scale the model was fitted on, which is not the caller's scale
  # when the fit centred its indicators.
  # Residual association given the class and, with direct effects, given the
  # covariates: the effects are removed before the pairs are compared.
  x <- .latents_direct_shift(.multilpa_center_like(object, as.matrix(data[continuous])),
                             object$extra_data, object)
  sampling <- .multilpa_level_weights(object, "individuals", object$n_observations)
  pairs <- utils::combn(seq_along(continuous), 2L)
  profiles <- if (identical(by, "profile"))
    as.integer(sub("^profile_", "", label)) else seq_len(object$n_profiles)
  rows <- lapply(seq_len(ncol(pairs)), function(index) {
    pair <- pairs[, index]
    observed_rows <- stats::complete.cases(x[, pair, drop = FALSE])
    effective <- sum(weight[observed_rows])
    moments <- lapply(profiles, function(profile) {
      profile_weight <- object$subject_posteriors[observed_rows, profile]
      if (!is.null(sampling)) profile_weight <- profile_weight * sampling[observed_rows]
      residual <- sweep(x[observed_rows, pair, drop = FALSE], 2L,
                        object$means[profile, pair], "-")
      covariance <- if (identical(object$covariance_model, "full"))
        matrix(object$covariances[pair, pair, profile], 2L) else
        diag(object$variances[profile, pair], 2L)
      list(observed = crossprod(residual * profile_weight, residual),
           expected = covariance * sum(profile_weight))
    })
    observed_covariance <- Reduce(`+`, lapply(moments, `[[`, "observed"))
    expected_covariance <- Reduce(`+`, lapply(moments, `[[`, "expected"))
    observed <- if (all(diag(observed_covariance) > 0))
      observed_covariance[1L, 2L] / sqrt(prod(diag(observed_covariance))) else NA_real_
    expected <- if (all(diag(expected_covariance) > 0))
      expected_covariance[1L, 2L] / sqrt(prod(diag(expected_covariance))) else NA_real_
    statistic <- if (is.null(sampling) && effective > 3 &&
                     is.finite(observed) && is.finite(expected)) {
      (atanh(pmin(pmax(observed, -0.999999), 0.999999)) -
         atanh(pmin(pmax(expected, -0.999999), 0.999999))) * sqrt(effective - 3)
    } else NA_real_
    data.frame(profile = label, indicator_1 = continuous[pair[1L]],
      indicator_2 = continuous[pair[2L]], kind = "gaussian",
      observed = observed, expected = expected, residual = observed - expected,
      effective_n = effective, statistic = statistic, df = NA_integer_,
      p_value = 2 * stats::pnorm(-abs(statistic)),
      row.names = NULL, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' The discrete indicators of a fit as category codes and class tables
#'
#' Categorical indicators with their response probabilities, then ordinal
#' indicators with the category probabilities of their adjacent-category
#' model; both are read from `data` by their category labels. A half-observed
#' pair row of a joined `local_dependence` item has no joint category and is
#' left out of that item's tables.
#' @param x A fitted model.
#' @param data The data, in fitting order.
#' @return A named list, one entry per discrete indicator: `codes` (one
#'   integer per row, `NA` when unobserved) and `probabilities` (profiles by
#'   categories).
#' @noRd
.multilpa_discrete_items <- function(x, data) {
  categorical <- x$categorical %||% character()
  nominal <- lapply(stats::setNames(categorical, categorical), function(name) {
    value <- data[[name]]
    levels_for <- x$categorical_levels[[name]]
    codes <- match(as.character(value), levels_for)
    half <- .latents_partial_levels(c(levels_for, as.character(value)))[-seq_along(levels_for)]
    stopifnot("Categorical indicators contain levels absent from the fit" =
                !any(is.na(codes) & !is.na(value) & !half))
    list(codes = codes, probabilities = x$response_probabilities[[name]])
  })
  extra <- x$extra_data
  ordinal_names <- colnames(extra$ordinal) %||% character()
  ordinal <- if (length(ordinal_names) == 0L) list() else {
    prepared <- .latents_prepare_extra_like(x, data)
    stats::setNames(lapply(seq_along(ordinal_names), function(j) {
      list(codes = prepared$ordinal[, j],
           probabilities = exp(.latents_ordinal_log_probabilities(
             x$ordinal_intercepts[[j]], x$ordinal_locations[, j])))
    }), ordinal_names)
  }
  c(nominal, ordinal)
}

#' Bivariate-residual chi-square for each pair of discrete indicators
#' @param discrete From `.multilpa_discrete_items()`.
#' @return One row per pair, or `NULL` when there are fewer than two.
#' @noRd
.multilpa_categorical_residuals <- function(object, discrete, weight, label) {
  if (length(discrete) < 2L) return(NULL)
  names_of <- names(discrete)
  posteriors <- object$subject_posteriors
  sampling <- .multilpa_level_weights(object, "individuals", object$n_observations)
  if (!is.null(sampling)) posteriors <- posteriors * sampling
  pairs <- utils::combn(seq_along(names_of), 2L)
  rows <- lapply(seq_len(ncol(pairs)), function(index) {
    first <- names_of[pairs[1L, index]]
    second <- names_of[pairs[2L, index]]
    codes_1 <- discrete[[first]]$codes
    codes_2 <- discrete[[second]]$codes
    block_1 <- discrete[[first]]$probabilities
    block_2 <- discrete[[second]]$probabilities
    complete <- !is.na(codes_1) & !is.na(codes_2)
    pair_weight <- weight[complete]
    category_1 <- factor(codes_1[complete], seq_len(ncol(block_1)))
    category_2 <- factor(codes_2[complete], seq_len(ncol(block_2)))
    observed <- as.matrix(stats::xtabs(pair_weight ~ category_1 + category_2))
    ## Conditional tables use the profile's own product distribution. The
    ## pooled table mixes those products with the complete-pair class counts.
    counts <- if (identical(label, "overall")) colSums(posteriors[complete, , drop = FALSE]) else {
      profile <- as.integer(sub("^profile_", "", label))
      replace(numeric(object$n_profiles), profile, sum(pair_weight))
    }
    expected <- Reduce(`+`, lapply(seq_len(object$n_profiles), function(k)
      counts[k] * outer(block_1[k, ], block_2[k, ])))
    positive <- expected > 0
    statistic <- if (!is.null(sampling) || sum(pair_weight) <= 0) NA_real_ else
      if (any(!positive & observed > 0)) Inf else
        sum((observed[positive] - expected[positive])^2 / expected[positive])
    degrees <- (nrow(observed) - 1L) * (ncol(observed) - 1L)
    ## Both tables are summarised by the same association measure, so that
    ## `observed`, `expected` and `residual` mean one thing across both kinds
    ## of pair. The model-implied table is a product of margins within a
    ## profile, so its V is zero there, exactly as the Gaussian `expected` is.
    observed_association <- .multilpa_cramers_v(observed)
    expected_association <- .multilpa_cramers_v(expected)
    data.frame(
      profile = label, indicator_1 = first, indicator_2 = second,
      kind = "categorical", observed = observed_association,
      expected = expected_association,
      residual = observed_association - expected_association,
      effective_n = sum(pair_weight), statistic = statistic,
      df = if (is.null(sampling)) as.integer(degrees) else NA_integer_,
      p_value = stats::pchisq(statistic, degrees, lower.tail = FALSE),
      row.names = NULL, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Numeric values of every indicator, for the residual correlations
#'
#' Gaussian indicators on the scale the model was fitted on, with direct
#' effects removed, as for the Gaussian pairs; discrete indicators as their
#' category codes (1 to K, in category order); counts as counts.
#' @param x A fitted model.
#' @param data The data, in fitting order.
#' @param continuous The Gaussian indicators.
#' @param discrete From `.multilpa_discrete_items()`.
#' @return A list with `values` (a rows-by-indicators numeric matrix, `NA`
#'   where unobserved) and `kind` (`"gaussian"`, `"discrete"` or `"count"`
#'   per column).
#' @noRd
.multilpa_pair_values <- function(x, data, continuous, discrete) {
  gaussian <- if (length(continuous) == 0L) NULL else
    .latents_direct_shift(.multilpa_center_like(x, as.matrix(data[continuous])),
                          x$extra_data, x)
  count_names <- colnames(x$extra_data$count) %||% character()
  counts <- if (length(count_names) == 0L) NULL else
    .latents_prepare_extra_like(x, data)$count
  coded <- if (length(discrete) == 0L) NULL else
    matrix(as.numeric(unlist(lapply(discrete, `[[`, "codes"), use.names = FALSE)),
           nrow(data), length(discrete), dimnames = list(NULL, names(discrete)))
  values <- cbind(gaussian, coded, counts)
  kind <- c(rep("gaussian", length(continuous)), rep("discrete", length(discrete)),
            rep("count", length(count_names)))
  names(kind) <- c(continuous, names(discrete), count_names)
  # The fit's own indicator order, so a pair reads the way the model names it.
  order_of <- order(match(names(kind), x$vars))
  list(values = values[, order_of, drop = FALSE], kind = kind[order_of])
}

#' Residual correlation for each pair that has neither a table nor a normal model
#'
#' A pair with a count indicator, or of a Gaussian and a discrete indicator,
#' is assessed by the posterior-weighted correlation of the two indicators'
#' values within the profile (pooled over the profiles under `"overall"`),
#' centred on each profile's weighted means. The model makes the indicators
#' independent within a profile, so the model-implied correlation is zero.
#' No reference distribution is attached: `statistic`, `df` and `p_value`
#' are `NA`.
#' @param values From `.multilpa_pair_values()`.
#' @return One row per such pair, or `NULL` when there is none.
#' @noRd
.multilpa_mixed_residuals <- function(object, values, weight, label, by) {
  kind <- values$kind
  if (length(kind) < 2L) return(NULL)
  pairs <- utils::combn(seq_along(kind), 2L)
  tested <- (kind[pairs[1L, ]] == "gaussian" & kind[pairs[2L, ]] == "gaussian") |
    (kind[pairs[1L, ]] == "discrete" & kind[pairs[2L, ]] == "discrete")
  pairs <- pairs[, !tested, drop = FALSE]
  if (ncol(pairs) == 0L) return(NULL)
  sampling <- .multilpa_level_weights(object, "individuals", object$n_observations)
  profiles <- if (identical(by, "profile"))
    as.integer(sub("^profile_", "", label)) else seq_len(object$n_profiles)
  rows <- lapply(seq_len(ncol(pairs)), function(index) {
    pair <- pairs[, index]
    v <- values$values[, pair, drop = FALSE]
    complete <- stats::complete.cases(v)
    v <- v[complete, , drop = FALSE]
    scatter <- Reduce(`+`, lapply(profiles, function(profile) {
      w <- object$subject_posteriors[complete, profile]
      if (!is.null(sampling)) w <- w * sampling[complete]
      if (sum(w) <= 0) return(matrix(0, 2L, 2L))
      centred <- sweep(v, 2L, colSums(v * w) / sum(w), "-")
      crossprod(centred * w, centred)
    }))
    observed <- if (all(diag(scatter) > 0))
      scatter[1L, 2L] / sqrt(prod(diag(scatter))) else NA_real_
    data.frame(profile = label, indicator_1 = colnames(values$values)[pair[1L]],
      indicator_2 = colnames(values$values)[pair[2L]],
      kind = if (any(kind[pair] == "count")) "count" else "mixed",
      observed = observed, expected = 0, residual = observed,
      effective_n = sum(weight[complete]), statistic = NA_real_,
      df = NA_integer_, p_value = NA_real_,
      row.names = NULL, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
