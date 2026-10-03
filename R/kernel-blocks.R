# Measurement blocks: the class-conditional densities of the observed data.
#
# A block maps its parameters to a rows-by-classes matrix of log densities;
# the structure (two-level mixing, membership logits, a Markov chain over
# occasions) turns those into a likelihood and posteriors. The profile models
# measure every row with up to four blocks that are conditionally independent
# given the profile -- Gaussian indicators (diagonal or full covariance, with
# missing values integrated out), categorical, ordinal and count indicators --
# and their log densities add. This file holds that sum, once, for every
# engine that measures rows: profiles, covariate models, the Houle
# cross-level family and both transition engines. Each block's M-step,
# scores and draws stay with the block (multilpa.R, categorical.R,
# extra-indicators.R, gaussian-moments.R).

#' Diagonal Gaussian log densities of every row under every profile
#'
#' @param x Observation-by-indicator matrix, complete.
#' @param means Profiles-by-indicators matrix.
#' @param variances Profiles-by-indicators matrix.
#' @param constant How the normalizing constant is formed: `"separate"`,
#'   `log(2 pi) + log(variance)` (profile and transition models), or
#'   `"joint"`, `log(2 pi variance)` (covariate models). They are the same
#'   number up to the last bit, and each engine keeps the one its results
#'   have always been computed with, since robust standard errors near a
#'   boundary amplify that last bit.
#' @return A rows-by-profiles matrix.
#' @noRd
.latents_gaussian_log_density <- function(x, means, variances,
                                          constant = c("separate", "joint")) {
  constant <- match.arg(constant)
  matrix(vapply(seq_len(nrow(means)), function(profile) {
    residuals <- sweep(x, 2L, means[profile, ], "-")
    normalizer <- switch(constant,
      separate = log(2 * pi) + log(variances[profile, ]),
      joint = log(2 * pi * variances[profile, ]))
    -0.5 * rowSums(sweep(residuals^2, 2L, variances[profile, ], "/") +
                     matrix(normalizer, nrow(x), ncol(x), byrow = TRUE))
  }, numeric(nrow(x))), nrow(x), nrow(means))
}

#' Log densities of every row under every profile, summed over the blocks
#'
#' The Gaussian block uses the conditional moments when indicators are
#' missing or covariances are full (the M-step then reuses those moments),
#' and the diagonal closed form otherwise. With `offset`, each row's largest
#' Gaussian log density is removed before the other blocks are added: a very
#' small density (log f = -5e15) would otherwise round away the priors and the
#' small categorical contributions, and the offset is added back to the
#' likelihood by the structure.
#'
#' @param x Observation-by-indicator matrix of continuous indicators, `NA`
#'   where unobserved; it may have zero columns.
#' @param parameters A list with `means`, `variances`, optionally
#'   `covariances`, `response_probabilities` and the extra-indicator blocks.
#' @param codes Optional categorical codes.
#' @param extra Optional ordinal and count data.
#' @param noise_log_density Optional constant log density of a uniform noise
#'   component, appended as the last column before the offset.
#' @param offset Whether to remove the row maximum of the Gaussian block.
#' @param n_profiles Number of profiles (from `means` by default).
#' @param constant Passed to `.latents_gaussian_log_density()`.
#' @return A list with `log_density` (rows by profiles, offset removed),
#'   `offset` (per row; zeros when `offset = FALSE`) and `moments` (the
#'   Gaussian conditional moments, or `NULL`).
#' @noRd
.latents_measurement_log_density <- function(x, parameters, codes = NULL,
                                             extra = NULL,
                                             noise_log_density = NULL,
                                             offset = TRUE,
                                             n_profiles = nrow(parameters$means),
                                             constant = "separate") {
  gaussian <- if (ncol(x) > 0L && (anyNA(x) || !is.null(parameters$covariances))) {
    .multilpa_gaussian_moments(x, parameters)
  } else NULL
  log_density <- if (!is.null(gaussian)) gaussian$log_density else
    .latents_gaussian_log_density(x, parameters$means, parameters$variances, constant)
  if (!is.null(noise_log_density)) log_density <- cbind(log_density, noise_log_density)
  row_offset <- rep(0, nrow(x))
  if (offset) {
    row_offset <- .multilpa_row_max(log_density)
    if (any(!is.finite(row_offset))) {
      stop("All component densities vanished, or a density overflowed.")
    }
    log_density <- sweep(log_density, 1L, row_offset, "-")
  }
  # Categorical, ordinal and count indicators are conditionally independent
  # of the continuous ones given the profile.
  if (!is.null(codes)) {
    log_density <- log_density +
      .multilpa_categorical_log_density(codes, parameters$response_probabilities)
  }
  if (!is.null(extra)) {
    log_density <- log_density +
      .latents_extra_log_density(extra, parameters, nrow(x), n_profiles)
  }
  list(log_density = log_density, offset = row_offset, moments = gaussian$moments)
}
