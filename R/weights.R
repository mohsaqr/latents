#' Resolve sampling weights to one scaled weight per unit
#'
#' Sampling weights belong to the independent units the likelihood sums over:
#' the groups of a two-level fit, the rows of a single-level one, the persons
#' of a transition fit. A weight that varies inside a unit would be a
#' within-unit weight, which needs a different estimator, so it is refused.
#' The weights are scaled to sum to the number of units, as Mplus does, so the
#' pseudo log likelihood and the information criteria stay on the sample's
#' scale.
#'
#' @param data The data frame of the fit.
#' @param weights `NULL`, or the name of a numeric column of `data`.
#' @param group_index Integer unit of each row, `1..n_groups`.
#' @param n_groups Number of units.
#' @return `NULL` for no weights, otherwise a numeric vector of length
#'   `n_groups` that sums to `n_groups`.
#' @noRd
.latents_sampling_weights <- function(data, weights, group_index, n_groups) {
  if (is.null(weights)) return(NULL)
  if (!is.character(weights) || length(weights) != 1L || is.na(weights) ||
      !weights %in% names(data)) {
    stop(errorCondition(
      "`weights` must name one column of `data` holding a sampling weight per unit.",
      class = "latents_bad_weights", call = NULL))
  }
  values <- data[[weights]]
  if (!is.numeric(values) || anyNA(values) || any(!is.finite(values)) ||
      any(values < 0)) {
    stop(errorCondition(sprintf(
      "Sampling weights in `%s` must be finite, non-missing and non-negative.",
      weights), class = "latents_bad_weights", call = NULL))
  }
  stopifnot("one unit per row" = length(group_index) == length(values))
  unit <- as.vector(tapply(values, group_index, min))
  spread <- as.vector(tapply(values, group_index, max)) - unit
  if (any(spread > sqrt(.Machine$double.eps) * pmax(1, abs(unit)))) {
    stop(errorCondition(sprintf(paste(
      "Sampling weights in `%s` vary inside a unit. The pseudo likelihood",
      "weights the independent units (the `id` groups, or the rows of a",
      "single-level fit); a within-unit weight is not supported."), weights),
      class = "latents_bad_weights", call = NULL))
  }
  if (sum(unit) <= 0 || sum(unit > 0) < 2L) {
    stop(errorCondition(sprintf(
      "Sampling weights in `%s` must be positive for at least two units.", weights),
      class = "latents_bad_weights", call = NULL))
  }
  stopifnot("one weight per unit" = length(unit) == n_groups)
  # Normalize before summing or multiplying: finite raw weights can overflow
  # both operations even though their relative weights are well defined.
  relative <- unit / max(unit)
  relative * (n_groups / sum(relative))
}

#' Refuse an option that sampling weights do not support
#' @param what The option, as the caller should read it.
#' @return Never returns; raises `latents_unsupported_weights`.
#' @noRd
.latents_refuse_weights <- function(what) {
  stop(errorCondition(sprintf(paste(
    "Sampling weights are not supported with %s."), what),
    class = "latents_unsupported_weights", call = NULL))
}

#' Weighted fits report sandwich standard errors only
#'
#' Under pseudo maximum likelihood the inverse Hessian of the weighted log
#' likelihood is not the variance of the estimate, and neither is the outer
#' product of the weighted scores alone; the sandwich is (Skinner, 1989).
#' A weighted fit therefore defaults to `"robust"` and refuses the others.
#'
#' @param weighted Whether the fit carries sampling weights.
#' @param vcov_type The requested type, already matched.
#' @param defaulted Whether the caller left `vcov_type` at its default.
#' @return The type to use.
#' @noRd
.latents_weighted_vcov <- function(weighted, vcov_type, defaulted) {
  if (!isTRUE(weighted)) return(vcov_type)
  if (defaulted) return("robust")
  if (!identical(vcov_type, "robust")) {
    stop(errorCondition(sprintf(paste(
      "A fit with sampling weights is a pseudo likelihood fit, whose",
      "`vcov_type = \"%s\"` errors are not valid; use `vcov_type = \"robust\"`",
      "(the default for weighted fits) or `method = \"bootstrap\"`."), vcov_type),
      class = "latents_unsupported_weights", call = NULL))
  }
  vcov_type
}

#' Whether a fit carries sampling weights
#' @param fit Any fit of the package.
#' @return `TRUE` when the fit was made with `weights`.
#' @noRd
.latents_is_weighted <- function(fit) {
  is.list(fit) && !is.null(fit$sampling_weights)
}

#' Print the sampling-weight line of a weighted fit
#'
#' Kish's effective sample size, `(sum w)^2 / sum w^2`, says how much the
#' weights cost in precision; it is shown beside the column name.
#'
#' @param x A fit.
#' @return `NULL`, invisibly; prints one line for a weighted fit, nothing
#'   otherwise.
#' @noRd
.latents_print_weights <- function(x) {
  if (!.latents_is_weighted(x)) return(invisible(NULL))
  weights <- x$sampling_weights
  cat(sprintf(paste("Sampling weights `%s` (pseudo maximum likelihood): Kish",
                    "effective n %.1f of %d units\n"),
              x$weights, sum(weights)^2 / sum(weights^2), length(weights)))
  invisible(NULL)
}
