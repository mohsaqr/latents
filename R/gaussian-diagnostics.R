#' Descriptive covariance and support diagnostics for Gaussian components
#'
#' Standardization uses the observed, sampling-weighted marginal variance
#' (divisor sum of weights), including between-component variation. Neither
#' the ratios nor membership measures classify a component as spurious or
#' establish identification. Gaussian mass uses the fit's normalized sampling
#' weights; its scale is not an independent sample count. Weight effective n
#' is Kish's concentration measure for posterior times sampling weight.
#'
#' @param parameters Gaussian means, variances and optional covariances.
#' @param data Observation-by-Gaussian-indicator matrix on the fitted scale.
#' @param posteriors Unweighted observation-by-component probabilities.
#' @param weights Nonnegative observation sampling weights, or `NULL`.
#' @param min_variance Absolute variance or covariance eigenvalue floor.
#' @param context Name of the density level or measurement block.
#' @return A data frame with one row per Gaussian component, or no rows when
#'   there are no Gaussian indicators. Missing reference variances give `NA`
#'   ratios, rather than an implied regular or irregular classification.
#' @noRd
.latents_gaussian_diagnostics <- function(parameters, data, posteriors,
                                          weights = NULL, min_variance = 1e-6,
                                          context = "profiles") {
  stopifnot(is.list(parameters), is.matrix(data), is.numeric(data),
            is.matrix(posteriors), is.numeric(posteriors),
            nrow(data) == nrow(posteriors),
            is.null(weights) || (is.numeric(weights) &&
              length(weights) == nrow(data) && all(is.finite(weights)) &&
              all(weights >= 0)),
            is.numeric(min_variance), length(min_variance) == 1L,
            is.finite(min_variance), min_variance > 0,
            is.character(context), length(context) == 1L, !is.na(context))
  empty <- data.frame(context = character(), profile = character(),
    gaussian_mass = numeric(), positive_rows = integer(),
    weight_effective_n = numeric(), min_variance_ratio = numeric(),
    min_standardized_eigenvalue = numeric(), boundary = logical())
  if (ncol(data) == 0L) return(empty)
  variances <- parameters$variances
  stopifnot(is.matrix(variances), ncol(variances) == ncol(data),
            nrow(variances) == ncol(posteriors),
            all(is.finite(posteriors)), all(posteriors >= 0))
  weights <- weights %||% rep(1, nrow(data))
  gaussian_rows <- rowSums(!is.na(data)) > 0L
  overall <- vapply(seq_len(ncol(data)), function(indicator) {
    observed <- !is.na(data[, indicator]) & weights > 0
    if (!any(observed)) return(NA_real_)
    values <- data[observed, indicator]
    weight <- weights[observed] / max(weights[observed])
    # Subtract a supported origin before means to protect large offsets and
    # exclude zero-weight outliers from both the origin and its arithmetic.
    values <- values - values[1L]
    center <- sum(values * weight) / sum(weight)
    variance <- sum(weight * (values - center)^2) / sum(weight)
    if (is.finite(variance) && variance > 0) variance else NA_real_
  }, numeric(1))
  labels <- rownames(variances) %||% paste0("profile_", seq_len(nrow(variances)))
  rows <- lapply(seq_len(nrow(variances)), function(profile) {
    contribution <- posteriors[, profile] * weights * gaussian_rows
    positive <- contribution > 0
    scaled <- if (any(positive)) contribution[positive] / max(contribution) else numeric()
    effective_n <- if (length(scaled)) sum(scaled)^2 / sum(scaled^2) else 0
    covariance <- if (is.null(parameters$covariances)) {
      diag(variances[profile, ], ncol(data))
    } else matrix(parameters$covariances[, , profile], ncol(data), ncol(data))
    ratio <- if (all(is.finite(overall))) {
      standard_deviation <- sqrt(overall)
      standardized <- sweep(sweep(covariance, 1L, standard_deviation, "/"),
                            2L, standard_deviation, "/")
      min(eigen(standardized, symmetric = TRUE, only.values = TRUE)$values)
    } else NA_real_
    data.frame(context = context, profile = labels[profile],
      gaussian_mass = sum(contribution), positive_rows = sum(positive),
      weight_effective_n = effective_n,
      min_variance_ratio = if (all(is.finite(overall)))
        min(variances[profile, ] / overall) else NA_real_,
      min_standardized_eigenvalue = ratio,
      boundary = min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values) <=
        min_variance * (1 + 1e-8), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Print the narrowest Gaussian direction without applying a cutoff
#' @param table A Gaussian diagnostic table, or `NULL` for an older fit.
#' @return `NULL`, invisibly. Called for its descriptive output.
#' @noRd
.latents_print_gaussian_diagnostics <- function(table) {
  stopifnot(is.null(table) || is.data.frame(table))
  if (is.null(table) || nrow(table) == 0L) return(invisible(NULL))
  usable <- which(is.finite(table$min_standardized_eigenvalue))
  if (length(usable) == 0L) {
    cat("Gaussian covariance ratios unavailable: insufficient observed reference variation.\n")
    return(invisible(NULL))
  }
  row <- usable[which.min(table$min_standardized_eigenvalue[usable])]
  cat(sprintf(paste0("Smallest standardized Gaussian eigenvalue: %.4g (%s, %s); ",
                     "Gaussian mass: %.4g; weight effective n: %.4g.\n"),
              table$min_standardized_eigenvalue[row], table$context[row],
              table$profile[row], table$gaussian_mass[row],
              table$weight_effective_n[row]))
  invisible(NULL)
}

#' Smallest available component diagnostic for a model comparison row
#' @param fit A fitted model, or `NULL` for a failed candidate.
#' @param column A numeric column of its Gaussian diagnostic table.
#' @return A numeric scalar; `NA` when the diagnostic is unavailable.
#' @noRd
.latents_gaussian_minimum <- function(fit, column) {
  stopifnot(is.null(fit) || is.list(fit), is.character(column), length(column) == 1L)
  values <- fit$gaussian_diagnostics[[column]]
  values <- values[is.finite(values)]
  if (length(values) == 0L) NA_real_ else min(values)
}
