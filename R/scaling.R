#' Column scales that make an optimization coordinate unit-free
#'
#' A multinomial logit is a parameterization of the same model under any linear
#' rescaling of a predictor: multiplying a covariate by `c` and its coefficient
#' by `1/c` leaves every fitted probability alone. A quasi-Newton search is not
#' invariant in that way, because it stops on a relative tolerance measured in
#' the coordinates it was handed. A column whose units are `1e-7` times those of
#' the rest therefore stalls the search while the score is still large, and the
#' optimizer still reports success.
#'
#' The root-mean-square of a column is exactly proportional to that column's
#' units, so dividing each column by it gives a design matrix that does not
#' change when a predictor is rescaled. Fitting in those coordinates and
#' dividing the coefficients back by the same scales returns estimates in the
#' user's units while making the search itself unit-free.
#'
#' @param design A numeric design matrix with finite entries.
#' @return A positive numeric vector with one entry per column of `design`. A
#'   column that is exactly zero, or whose root-mean-square is not finite, gets
#'   `1`, so the scaling is never itself a source of `NaN`.
#' @noRd
.multilpa_design_scale <- function(design) {
  stopifnot(
    "`design` must be a numeric matrix" =
      is.matrix(design) && is.numeric(design),
    "`design` must be finite" = all(is.finite(design)))
  if (ncol(design) == 0L) return(numeric(0))
  ## The peak comes out first so that squaring a column held in very large
  ## units cannot overflow to `Inf` before the mean is taken.
  peak <- apply(abs(design), 2L, max)
  safe <- ifelse(is.finite(peak) & peak > 0, peak, 1)
  scale <- safe * sqrt(colMeans(sweep(design, 2L, safe, "/")^2))
  unname(ifelse(is.finite(scale) & scale > 0, scale, 1))
}

#' Scales for the log-Cholesky coordinates of a residual covariance
#'
#' The unrestricted residual covariance is estimated through the lower triangle
#' of its Cholesky factor, with the diagonal logged. Those coordinates do not
#' share one unit: a logged diagonal entry is dimensionless, while an
#' off-diagonal entry `L[i, j]` carries the units of indicator `i`. Handing them
#' to a finite-difference Hessian unscaled makes the information matrix's
#' condition number grow with the square of the indicator units, and it is then
#' declared singular for no reason but the units the data arrived in.
#'
#' This is the general form of the scaling [parameter_inference()] already
#' applies to a covariate-free `multilpa` fit (R/inference.R), written so that
#' the covariate path can share it instead of rediscovering it.
#'
#' @param variances A profiles-by-indicators matrix of residual variances, the
#'   diagonal of each profile's residual covariance.
#' @param dimension Number of continuous indicators.
#' @param rows Integer indices of the rows of `variances` whose spread block is
#'   estimated separately: every profile under a varying spread, `1L` when one
#'   block is shared.
#' @return A positive numeric vector of length
#'   `length(rows) * dimension * (dimension + 1) / 2`, ordered as the
#'   coordinates themselves are: one block per entry of `rows`, and within a
#'   block the lower triangle including the diagonal, in column-major order.
#'   Diagonal entries get `1` because they are logged.
#' @noRd
.multilpa_covariance_coordinate_scale <- function(variances, dimension, rows) {
  stopifnot(
    "`variances` must be a numeric matrix" =
      is.matrix(variances) && is.numeric(variances),
    "`dimension` must be a single non-negative whole number" =
      is.numeric(dimension) && length(dimension) == 1L && is.finite(dimension) &&
      dimension >= 0 && dimension == as.integer(dimension),
    "`variances` must have at least `dimension` columns" =
      ncol(variances) >= dimension,
    "`rows` must index rows of `variances`" =
      is.numeric(rows) && length(rows) >= 1L && all(is.finite(rows)) &&
      all(rows >= 1L) && all(rows <= nrow(variances)),
    "`variances` must be positive" = all(variances > 0))
  if (dimension == 0L) return(numeric(0))
  positions <- which(lower.tri(matrix(0, dimension, dimension), diag = TRUE),
                     arr.ind = TRUE)
  unlist(lapply(rows, function(profile) {
    ifelse(positions[, 1L] == positions[, 2L], 1,
           sqrt(variances[profile, positions[, 1L]]))
  }), use.names = FALSE)
}
