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
#' model does not have. `multilpa_transitions` is kept separate for the same
#' reason: its mixing weights are an initial distribution and a transition
#' matrix, not one prevalence vector. The diagnostics are the exception: they
#' read only posteriors and counts, which every fit carries, so they accept any
#' of them.
#'
#' @param object Any object.
#' @return `TRUE` for a fitted model of any of these classes.
#' @noRd
.multilpa_any_fit <- function(object) {
  inherits(object, "multilpa") || inherits(object, "multilpa_covariates") ||
    inherits(object, "multilpa_transitions")
}

#' Which values are seeds that `set.seed()` takes exactly?
#'
#' A seed is a whole number within the integer range, of either sign: that is
#' what `set.seed()` accepts. A fraction is rejected because `set.seed()` would
#' truncate it without notice, so two different seeds would give one stream.
#'
#' @param seed Any vector.
#' @return A logical vector, one element per element of `seed`.
#' @noRd
.multilpa_is_seed <- function(seed) {
  if (!is.numeric(seed)) return(rep(FALSE, length(seed)))
  finite <- is.finite(seed)
  finite & abs(ifelse(finite, seed, 0)) <= .Machine$integer.max &
    ifelse(finite, seed == trunc(seed), FALSE)
}

#' Refuse a seed argument that `set.seed()` would not take exactly
#' @param seed The supplied seed, or `NULL`.
#' @return `NULL`, invisibly; raises `latents_bad_argument` otherwise.
#' @noRd
.multilpa_check_seed <- function(seed) {
  if (is.null(seed) || (length(seed) == 1L && .multilpa_is_seed(seed))) {
    return(invisible(NULL))
  }
  stop(errorCondition(sprintf(
    "`seed` must be NULL or a single whole number between -%d and %d.",
    .Machine$integer.max, .Machine$integer.max),
    class = "latents_bad_argument", call = NULL))
}
