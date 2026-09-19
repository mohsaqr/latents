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
#' @param object A fitted `multilpa` or `multilpa_covariates` model.
#' @param data The data frame the model was fitted to.
#' @param by `"profile"` (the default) assesses each profile separately, which
#'   is where the assumption is actually made. `"overall"` pools the profiles
#'   into one posterior-weighted table per pair.
#' @return A base `data.frame` with one row per indicator pair and, for
#'   `by = "profile"`, per profile. Columns are `profile` (`"overall"` when
#'   pooled), `indicator_1`, `indicator_2`, `kind` (`"gaussian"` or
#'   `"categorical"`), `observed`, `expected`, `residual`, `effective_n`,
#'   `statistic`, `df` and `p_value`, ordered by decreasing absolute residual so
#'   the worst pair is the first row.
#'
#'   For a Gaussian pair, `observed` is the posterior-weighted residual
#'   correlation within the profile, `expected` is what the fitted covariance
#'   model implies for it (zero under `"diagonal"`), and the statistic is a
#'   Fisher z. For a categorical pair, `observed` and `expected` are the summed
#'   absolute cell discrepancy and the table total, and the statistic is the
#'   bivariate-residual chi-square on `(categories_1 - 1)(categories_2 - 1)`
#'   degrees of freedom.
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
#'   Gaussian p-values are unavailable when the pair has effective size at most
#'   three or zero residual spread.
#' @references Vermunt, J. K., & Magidson, J. (2004). Local dependence in
#'   latent class models. In *The Sage Encyclopedia of Social Science Research
#'   Methods*.
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
#' @export
bivariate_residuals <- function(object, data, by = c("profile", "overall")) {
  stopifnot(
    "`object` must be a fitted model of this package" = .multilpa_any_fit(object),
    "`data` must be a data frame" = is.data.frame(data),
    "`data` must have one row per observation of the fit" =
      nrow(data) == object$n_observations
  )
  if (inherits(object, "multilpa_random_intercept")) {
    stop("Bivariate residuals require a discrete group-class model.")
  }
  by <- match.arg(by)
  continuous <- .multilpa_continuous_names(object)
  categorical <- object$categorical %||% character()
  stopifnot("`data` must contain the fitted indicators" =
              all(c(continuous, categorical) %in% names(data)))
  stopifnot("Gaussian indicators must be numeric and finite when observed" =
    all(vapply(data[continuous], function(value)
      is.numeric(value) && all(is.na(value) | is.finite(value)), logical(1))))
  weights <- if (identical(by, "profile")) {
    stats::setNames(lapply(seq_len(object$n_profiles),
                           function(k) object$subject_posteriors[, k]),
                    paste0("profile_", seq_len(object$n_profiles)))
  } else list(overall = rep(1, object$n_observations))

  rows <- lapply(names(weights), function(label) {
    weight <- weights[[label]]
    rbind(
      .multilpa_gaussian_residuals(object, data, continuous, weight, label, by),
      .multilpa_categorical_residuals(object, data, categorical, weight, label))
  })
  result <- do.call(rbind, rows)
  if (is.null(result) || nrow(result) == 0L) {
    return(data.frame(profile = character(), indicator_1 = character(),
                      indicator_2 = character(), kind = character(),
                      observed = numeric(), expected = numeric(),
                      residual = numeric(), effective_n = numeric(),
                      statistic = numeric(), df = numeric(), p_value = numeric()))
  }
  result <- result[order(-abs(result$residual)), , drop = FALSE]
  row.names(result) <- NULL
  result
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
      effective_n = effective, statistic = statistic, df = NA_real_,
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
    data.frame(
      profile = label, indicator_1 = first, indicator_2 = second,
      kind = "categorical", observed = sum(abs(observed - expected)),
      expected = sum(expected), residual = statistic / max(degrees, 1L),
      effective_n = sum(pair_weight), statistic = statistic, df = degrees,
      p_value = stats::pchisq(statistic, degrees, lower.tail = FALSE),
      row.names = NULL, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
