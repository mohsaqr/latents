#' Lo-Mendell-Rubin adjusted likelihood-ratio statistic
#'
#' Computes the likelihood-ratio statistic comparing two fitted models and the
#' Lo-Mendell-Rubin small-sample adjustment of it. **No p-value is returned.**
#' The reference distribution of this statistic under a class-count null is not
#' chi-square, and the package does not reproduce the Vuong-Lo-Mendell-Rubin
#' reference distribution, so reporting a chi-square tail probability here would
#' be wrong rather than approximate. Use [bootstrap_lrt()] for a
#' calibrated p-value.
#'
#' @param null_model The smaller fitted `multilpa` model.
#' @param alternative_model The larger fitted `multilpa` model, with strictly more
#'   free parameters and a log likelihood at least as large. Both fits must use
#'   the same data, measurement specification and parameter bounds. Neither
#'   class count may decrease in the larger model.
#' @param n Sample size entering the adjustment. `"individuals"` uses the
#'   individual count and matches the convention external mixture software
#'   applies; `"groups"` uses the independent group count, which is the
#'   convention this package's default BIC uses. A single positive number is
#'   also accepted.
#' @return A one-row base `data.frame` with columns `statistic` (twice the log
#'   likelihood difference), `df` (the difference in free parameters), `n`,
#'   `adjustment_factor`, `adjusted_statistic`, and `p_value`. `p_value` is
#'   always `NA_real_`, for the reason given above; the column exists so that
#'   results bind with other test tables rather than to suggest a test was
#'   performed.
#' @details The adjustment divides the statistic by `1 + 1 / (df * log(n))`,
#'   the Lo-Mendell-Rubin correction. That factor was verified against two
#'   independent genuine Mplus TECH11 runs, reproducing the reported adjusted
#'   values to the precision Mplus prints. The unadjusted statistic itself is
#'   simply `2 (L1 - L0)`.
#'
#'   Class-count comparisons violate the regularity conditions of the ordinary
#'   likelihood-ratio test, because the null places a class on the boundary of
#'   the parameter space and leaves that class's parameters unidentified. This
#'   is why the statistic is reported without a nominal p-value.
#' @references Lo, Y., Mendell, N. R., & Rubin, D. B. (2001). Testing the number
#'   of components in a normal mixture. Biometrika, 88, 767--778.
#'   Nylund, K. L., Asparouhov, T., & Muthen, B. O. (2007). Deciding on the
#'   number of classes in latent class analysis and growth mixture modeling.
#'   Structural Equation Modeling, 14, 535--569.
#' @examples
#' set.seed(7)
#' example_data <- data.frame(
#'   school = rep(seq_len(12), each = 10),
#'   score_a = rnorm(120), score_b = rnorm(120)
#' )
#' smaller <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                       n_profiles = 1, n_group_classes = 1, n_starts = 2, seed = 1)
#' larger <- multilpa(example_data, c("score_a", "score_b"), "school",
#'                      n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
#' lmr_lrt(smaller, larger)
#' @export
lmr_lrt <- function(null_model, alternative_model,
                           n = c("individuals", "groups")) {
  stopifnot(
    "`null_model` must be an `multilpa` fit" = inherits(null_model, "multilpa"),
    "`alternative_model` must be an `multilpa` fit" =
      inherits(alternative_model, "multilpa"))
  sample_size <- if (is.character(n) || missing(n)) {
    convention <- match.arg(n)
    .multilpa_check_comparable(null_model, alternative_model)
    if (convention == "individuals") null_model$n_observations else
      null_model$n_groups
  } else {
    stopifnot("`n` must be a single number greater than one" =
                is.numeric(n) && length(n) == 1L && is.finite(n) && n > 1)
    .multilpa_check_comparable(null_model, alternative_model)
    n
  }
  stopifnot("The sample size entering the adjustment must be greater than one" =
              length(sample_size) == 1L && is.finite(sample_size) && sample_size > 1)
  smaller_classes <- c(null_model$n_profiles, null_model$n_group_classes)
  larger_classes <- c(alternative_model$n_profiles, alternative_model$n_group_classes)
  if (length(smaller_classes) == 2L && length(larger_classes) == 2L &&
      (any(larger_classes < smaller_classes) || all(larger_classes == smaller_classes))) {
    stop(errorCondition("The alternative must add classes without decreasing either class count.",
                        class = "multilpa_bad_nesting", call = NULL))
  }
  df <- alternative_model$n_parameters - null_model$n_parameters
  if (df < 1L) {
    stop(errorCondition("`alternative_model` must have more free parameters than `null_model`.",
                        class = "multilpa_bad_nesting", call = NULL))
  }
  statistic <- 2 * (alternative_model$log_likelihood - null_model$log_likelihood)
  if (statistic < 0) {
    stop(errorCondition(sprintf(
      "`alternative_model` has the smaller log likelihood by %.4g; refit with more starts.",
      -statistic / 2), class = "multilpa_reversed_likelihood", call = NULL))
  }
  adjustment_factor <- 1 + 1 / (df * log(sample_size))
  data.frame(statistic = statistic, df = df, n = sample_size,
             adjustment_factor = adjustment_factor,
             adjusted_statistic = statistic / adjustment_factor,
             p_value = NA_real_)
}

#' Check that two fits describe the same observations
#' @param null_model The smaller fitted model.
#' @param alternative_model The larger fitted model.
#' @return `NULL`, invisibly; raises a classed condition when the fits differ.
#' @noRd
.multilpa_check_comparable <- function(null_model, alternative_model) {
  fields <- c("n_observations", "n_groups", "indicators", "group", "group_values",
              "group_index", "continuous", "categorical", "categorical_levels",
              "indicator_data", "categorical_data", "variance_model",
              "covariance_model", "min_variance", "min_probability")
  same <- all(vapply(fields, function(field)
    identical(null_model[[field]], alternative_model[[field]]), logical(1)))
  if (!same) {
    stop(errorCondition(paste("Both models must be fitted to the same individuals, groups,",
                              "indicators, measurement specification, and parameter bounds."),
                        class = "multilpa_incomparable_models", call = NULL))
  }
  invisible(NULL)
}
