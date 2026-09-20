#' Residual association between pairs of indicators
#'
#' Diagonal Gaussian and categorical measurement models assume indicators are
#' independent within a profile. This diagnostic compares the residual
#' association with that implied by the fitted measurement model, including
#' the estimated association under a full Gaussian covariance model.
#'
#' A large residual means the pair shares something the profiles do not capture.
#' The usual remedies are `covariance_model = "full"`, which estimates the
#' association rather than assuming it away, or another profile.
#'
#' @param x A fitted model of this package that has discrete group classes: a
#'   `multilpa`, `multilpa_covariates` or `multilpa_transitions` fit. A
#'   `multilpa_random_intercept` fit raises `multilpa_no_group_classes`.
#' @param data Optional. The data frame the model was fitted to, with one row
#'   per observation and in the order it was fitted in; when omitted it is
#'   rebuilt from the indicators, identifiers and occasions the fit stores,
#'   which round-trip exactly. Supplying it is the stronger check that the
#'   caller still holds that frame: the row count, the indicators, the group
#'   identifier and the occasion are all compared against the fit row by row,
#'   and a frame in any other order raises `multilpa_bad_inference_data` rather
#'   than pairing each row with another observation's posterior.
#' @param by `"profile"` (the default) assesses each profile separately, which
#'   is where the assumption is actually made. `"overall"` pools the profiles
#'   into one posterior-weighted table per pair.
#' @param adjust How to correct `p_value` for the number of pairs tested, one of
#'   the methods [stats::p.adjust()] accepts. Every pair of indicators is a
#'   separate test, so a five-indicator fit asks ten questions at
#'   `by = "overall"` and ten per profile otherwise. The default is `"none"`,
#'   which leaves `p_adjusted` equal to `p_value`; `"BH"` is the usual choice
#'   when the whole table is being scanned for the worst pair.
#' @return A base `data.frame`, one row per indicator pair and, for
#'   `by = "profile"`, per profile, with the columns
#'   \describe{
#'     \item{`profile`}{character: `"profile_k"`, or `"overall"` when pooled.}
#'     \item{`indicator_1`, `indicator_2`}{character: the pair, in the order the
#'       indicators were fitted.}
#'     \item{`kind`}{character: `"gaussian"` or `"categorical"`, the measurement
#'       model both indicators share.}
#'     \item{`observed`}{numeric: the observed residual association, on the
#'       correlation scale. A Gaussian pair uses the posterior-weighted residual
#'       correlation, signed; a categorical pair uses Cramer's V of the
#'       posterior-weighted contingency table, which is unsigned.}
#'     \item{`expected`}{numeric: the association the fitted measurement model
#'       implies for that pair, on the same scale as `observed`. It is zero
#'       within a profile under local independence, and nonzero when the model
#'       estimates the association (`covariance_model = "full"`) or when pooling
#'       mixes profiles whose distributions differ.}
#'     \item{`residual`}{numeric: `observed - expected`, on every row. This is
#'       the column the table is ranked by.}
#'     \item{`effective_n`}{numeric: posterior weight of the rows observed on
#'       both indicators.}
#'     \item{`statistic`}{numeric: the test statistic, which depends on `kind`.
#'       A Gaussian pair carries a Fisher z, a standard normal deviate; a
#'       categorical pair carries the bivariate-residual chi-square.}
#'     \item{`df`}{integer: degrees of freedom of `statistic`. A categorical pair
#'       carries `(categories_1 - 1)(categories_2 - 1)`. A Gaussian pair carries
#'       `NA_integer_`, because a standard normal deviate has no degrees of
#'       freedom; this is a column that does not apply to that `kind`, not a
#'       missing value.}
#'     \item{`p_value`}{numeric: the two-sided normal tail for a Gaussian pair,
#'       the upper chi-square tail for a categorical one. Descriptive only, see
#'       Details.}
#'     \item{`p_adjusted`}{numeric: `p_value` corrected for the number of pairs
#'       tested, by the method `adjust` names. Equal to `p_value` under the
#'       default `"none"`.}
#'   }
#'   Rows are ordered by decreasing `abs(residual)`, so the worst pair is the
#'   first row, with `indicator_1`, `indicator_2` and `profile` as deterministic
#'   secondary keys. `observed`, `expected` and `residual` are on one scale for
#'   both kinds, so that ordering is meaningful in a mixed fit.
#'
#'   `observed`, `expected` and therefore `residual` are `NA_real_` when the
#'   association is undefined: a Gaussian pair with no residual spread, or a
#'   pair with no effective weight. `statistic` and `p_value` are `NA_real_` for
#'   a Gaussian pair with effective size at most three or no residual spread,
#'   and for a categorical pair with no effective weight. The classic
#'   Vermunt-Magidson bivariate residual for a categorical pair is
#'   `statistic / df`.
#' @details The p-values are descriptive approximations: they do not account
#'   for estimated posterior memberships, fitted measurement parameters, or
#'   dependence within groups. They are not calibrated significance tests.
#'   Use the residuals to rank possible model misspecification; a multiplicity
#'   correction alone does not resolve these calibration limitations.
#'
#'   Pairs of different kinds are not assessed and do not appear in the table; a
#'   mixed fit therefore returns fewer rows than it has pairs.
#'   Each pair uses only rows observed on both indicators, and `effective_n`
#'   reports their posterior weight. Gaussian pooling averages within-profile
#'   second moments, excluding association explained by differences in profile
#'   means. A categorical profile is compared with its own response distribution;
#'   overall categorical tables use the posterior-weighted mixture distribution.
#'   Cramer's V of a two-way table with total `N` is
#'   `sqrt(X^2 / (N min(rows - 1, columns - 1)))`, where `X^2` is that table's
#'   own chi-square against independence; applied to the model-implied table it
#'   measures exactly the association the measurement model accounts for, and is
#'   zero within a profile, where the model is a product of margins.
#'   A model raising the error condition `multilpa_no_group_classes` has no
#'   discrete group classes and cannot be assessed this way.
#' @references Vermunt, J. K., & Magidson, J. (2004). Local dependence in
#'   latent class models. In *The Sage Encyclopedia of Social Science Research
#'   Methods*. Cramer, H. (1946). *Mathematical Methods of Statistics*.
#'   Princeton University Press.
#' @examples
#' set.seed(4)
#' school <- rep(seq_len(20), each = 10)
#' profile <- rep(rep(c(1L, 2L), length.out = 20), each = 10)
#' shared <- rnorm(200)
#' example_data <- data.frame(
#'   school = school,
#'   a = rnorm(200, ifelse(profile == 2L, 2, -2), 0.8) + shared,
#'   b = rnorm(200, ifelse(profile == 2L, 2, -2), 0.8) + shared,
#'   c = rnorm(200, ifelse(profile == 2L, 1, -1), 0.8)
#' )
#' fit <- multilpa(example_data, c("a", "b", "c"), "school", n_profiles = 2,
#'                 n_group_classes = 1, n_starts = 3, seed = 1)
#' bivariate_residuals(fit, example_data)
#' bivariate_residuals(fit, example_data, by = "overall")
#' @export
bivariate_residuals <- function(x, data = NULL, by = c("profile", "overall"),
                               adjust = .multilpa_p_adjust_methods) {
  stopifnot("`x` must be a fitted model of this package" = .multilpa_any_fit(x))
  data <- .multilpa_resolve_data(x, data)
  stopifnot(
    "`data` must be a data frame" = is.data.frame(data),
    "`data` must have one row per observation of the fit" =
      nrow(data) == x$n_observations
  )
  if (inherits(x, "multilpa_random_intercept")) {
    stop(errorCondition("Bivariate residuals require a discrete group-class model.",
                        class = "multilpa_no_group_classes", call = NULL))
  }
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
  x <- as.matrix(data[continuous])
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
