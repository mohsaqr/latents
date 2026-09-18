# Internal helpers. No Rd page: nothing here is exported.

# Polyfill `%||%` for R < 4.4. Base R 4.4 added it; this package declares
# R (>= 4.1.0), so the base version cannot be relied on. Defined
# unconditionally because the package namespace is searched first, and base's
# version on R >= 4.4 is functionally identical, so this is a harmless shadow.
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Is this any fitted model of this package?
#'
#' `multilpa_covariates` deliberately does not inherit from `multilpa`, because
#' most methods for the latter assume fixed mixing weights that a covariate
#' model does not have. The diagnostics are the exception: they read only
#' posteriors and counts, which both fits carry, so they accept either.
#'
#' @param object Any object.
#' @return `TRUE` for a fitted model of either class.
#' @noRd
.multilpa_any_fit <- function(object) {
  inherits(object, "multilpa") || inherits(object, "multilpa_covariates") ||
    inherits(object, "multilpa_random_intercept")
}
