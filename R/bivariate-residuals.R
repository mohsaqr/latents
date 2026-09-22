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
  data <- .multilpa_resolve_data(x, data)
  stopifnot(
    "`data` must be a data frame" = is.data.frame(data),
    "`data` must have one row per observation of the fit" =
      nrow(data) == x$n_observations
  )
  by <- match.arg(by)
  adjust <- match.arg(adjust)
  continuous <- .multilpa_continuous_names(x)
  categorical <- x$categorical %||% character()
  stopifnot("`data` must contain the fitted indicators" =
              all(c(continuous, categorical) %in% names(data)))
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

  rows <- lapply(names(weights), function(label) {
    weight <- weights[[label]]
    rbind(
      .multilpa_gaussian_residuals(x, data, continuous, weight, label, by),
      .multilpa_categorical_residuals(x, data, categorical, weight, label))
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
  x <- .multilpa_center_like(object, as.matrix(data[continuous]))
  pairs <- utils::combn(seq_along(continuous), 2L)
  profiles <- if (identical(by, "profile"))
    as.integer(sub("^profile_", "", label)) else seq_len(object$n_profiles)
  rows <- lapply(seq_len(ncol(pairs)), function(index) {
    pair <- pairs[, index]
    observed_rows <- stats::complete.cases(x[, pair, drop = FALSE])
    effective <- sum(weight[observed_rows])
    moments <- lapply(profiles, function(profile) {
      profile_weight <- object$subject_posteriors[observed_rows, profile]
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
    statistic <- if (effective > 3 && is.finite(observed) && is.finite(expected)) {
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

#' Bivariate-residual chi-square for each categorical pair
#' @return One row per categorical pair, or `NULL` when there are fewer than two.
#' @noRd
.multilpa_categorical_residuals <- function(object, data, categorical, weight,
                                            label) {
  if (length(categorical) < 2L) return(NULL)
  blocks <- object$response_probabilities
  posteriors <- object$subject_posteriors
  pairs <- utils::combn(seq_along(categorical), 2L)
  rows <- lapply(seq_len(ncol(pairs)), function(index) {
    first <- categorical[pairs[1L, index]]
    second <- categorical[pairs[2L, index]]
    codes_1 <- match(as.character(data[[first]]), object$categorical_levels[[first]])
    codes_2 <- match(as.character(data[[second]]), object$categorical_levels[[second]])
    stopifnot("Categorical indicators contain levels absent from the fit" =
      !any(is.na(codes_1) & !is.na(data[[first]])) &&
      !any(is.na(codes_2) & !is.na(data[[second]])))
    complete <- !is.na(codes_1) & !is.na(codes_2)
    pair_weight <- weight[complete]
    category_1 <- factor(codes_1[complete], seq_along(object$categorical_levels[[first]]))
    category_2 <- factor(codes_2[complete], seq_along(object$categorical_levels[[second]]))
    observed <- as.matrix(stats::xtabs(pair_weight ~ category_1 + category_2))
    ## Conditional tables use the profile's own product distribution. The
    ## pooled table mixes those products with the complete-pair class counts.
    counts <- if (identical(label, "overall")) colSums(posteriors[complete, , drop = FALSE]) else {
      profile <- as.integer(sub("^profile_", "", label))
      replace(numeric(object$n_profiles), profile, sum(pair_weight))
    }
    expected <- Reduce(`+`, lapply(seq_len(object$n_profiles), function(k)
      counts[k] * outer(blocks[[first]][k, ], blocks[[second]][k, ])))
    positive <- expected > 0
    statistic <- if (sum(pair_weight) <= 0) NA_real_ else
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
      df = as.integer(degrees),
      p_value = stats::pchisq(statistic, degrees, lower.tail = FALSE),
      row.names = NULL, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
