#' Residual association between pairs of indicators
#'
#' Every model in this package assumes that indicators are independent within a
#' profile: `covariance_model = "diagonal"` says so for Gaussian indicators, and
#' the categorical measurement model says so by construction. Nothing else in
#' the package tests that assumption, so a model can fit a dataset whose
#' indicators remain strongly related inside a profile and report nothing
#' unusual. This verb looks for exactly that.
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
#' @details The tests are approximate in a way worth stating: they treat the
#'   posterior profile memberships as known rather than estimated, so the
#'   effective sample size is optimistic and the p-values are anti-conservative.
#'   Read them as a ranking of which pairs are worst, not as exact levels, and
#'   apply a multiplicity correction across the pairs before calling any one of
#'   them significant.
#'
#'   Pairs of different kinds are not assessed and do not appear in the table; a
#'   mixed fit therefore returns fewer rows than it has pairs.
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
  by <- match.arg(by)
  continuous <- .multilpa_continuous_names(object)
  categorical <- object$categorical %||% character()
  stopifnot("`data` must contain the fitted indicators" =
              all(c(continuous, categorical) %in% names(data)))
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
  ## Centre on the profile's own means when weighting by that profile, and on
  ## the posterior-weighted mixture mean when pooling, so the residual is always
  ## what the model has failed to explain rather than what it has explained.
  centre <- if (identical(by, "profile")) {
    object$means[as.integer(sub("^profile_", "", label)), ]
  } else colSums(object$subject_posteriors %*% object$means) / object$n_observations
  residual <- sweep(x, 2L, centre, "-")
  effective <- sum(weight)
  covariance <- crossprod(residual * weight, residual) / effective
  spread <- sqrt(diag(covariance))
  correlation <- covariance / outer(spread, spread)
  implied <- if (identical(object$covariance_model, "full") &&
                 identical(by, "profile")) {
    profile <- as.integer(sub("^profile_", "", label))
    block <- matrix(object$covariances[, , profile], length(continuous))
    block / outer(sqrt(diag(block)), sqrt(diag(block)))
  } else diag(length(continuous))

  pairs <- utils::combn(seq_along(continuous), 2L)
  observed <- correlation[t(pairs)]
  expected <- implied[t(pairs)]
  difference <- observed - expected
  ## Fisher's z on the difference of two correlations, with the profile's
  ## effective size. Approximate: the posteriors are treated as known.
  statistic <- (atanh(pmin(pmax(observed, -0.999999), 0.999999)) -
                  atanh(pmin(pmax(expected, -0.999999), 0.999999))) *
    sqrt(max(effective - 3, 1))
  data.frame(
    profile = label, indicator_1 = continuous[pairs[1L, ]],
    indicator_2 = continuous[pairs[2L, ]], kind = "gaussian",
    observed = observed, expected = expected, residual = difference,
    effective_n = effective, statistic = statistic, df = NA_real_,
    p_value = 2 * stats::pnorm(-abs(statistic)),
    row.names = NULL, stringsAsFactors = FALSE)
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
    observed <- as.matrix(stats::xtabs(weight ~ codes_1 + codes_2))
    ## Under local independence the expected table is the posterior-weighted
    ## product of the two response distributions.
    expected <- Reduce(`+`, lapply(seq_len(object$n_profiles), function(k) {
      sum(weight * posteriors[, k]) *
        outer(blocks[[first]][k, ], blocks[[second]][k, ])
    }))
    expected <- expected * sum(observed) / sum(expected)
    statistic <- sum((observed - expected)^2 / pmax(expected, 1e-10))
    degrees <- (nrow(observed) - 1L) * (ncol(observed) - 1L)
    data.frame(
      profile = label, indicator_1 = first, indicator_2 = second,
      kind = "categorical", observed = sum(abs(observed - expected)),
      expected = sum(expected), residual = statistic / max(degrees, 1L),
      effective_n = sum(weight), statistic = statistic, df = degrees,
      p_value = stats::pchisq(statistic, degrees, lower.tail = FALSE),
      row.names = NULL, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
