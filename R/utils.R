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

#' Choose the start a fit reports
#'
#' Restarts of a mixture routinely reach the same optimum under different label
#' permutations, and then differ only in the last bits of the likelihood.
#' Picking by `which.max()` alone makes the reported labelling turn on
#' floating-point noise, so the maximum is taken up to a relative tolerance far
#' below any difference that could matter, a converged start is preferred among
#' those tied, and then the earliest.
#'
#' @param log_likelihood Numeric vector, one entry per start; `-Inf` for a
#'   start that failed.
#' @param converged Logical vector of the same length.
#' @param rule `"likelihood"` takes the highest likelihood whether or not that
#'   start converged. `"converged"` takes the highest likelihood among the
#'   converged starts, falling back to every start when none converged.
#' @return A single integer index.
#' @noRd
.multilpa_select_start <- function(log_likelihood, converged,
                                   rule = c("likelihood", "converged")) {
  rule <- match.arg(rule)
  stopifnot(
    "`log_likelihood` must be numeric" = is.numeric(log_likelihood),
    "`converged` must be logical and match `log_likelihood`" =
      is.logical(converged) && length(converged) == length(log_likelihood),
    "at least one start must have a finite likelihood" =
      any(is.finite(log_likelihood))
  )
  usable <- is.finite(log_likelihood)
  settled <- usable & converged
  candidates <- if (identical(rule, "converged") && any(settled)) {
    which(settled)
  } else which(usable)
  best <- max(log_likelihood[candidates])
  tied <- candidates[log_likelihood[candidates] >= best - 1e-10 * (1 + abs(best))]
  preferred <- tied[converged[tied]]
  if (length(preferred) > 0L) preferred[1L] else tied[1L]
}

#' Row maxima of a numeric matrix
#'
#' `pmax()` over the columns: the same values as `apply(x, 1L, max)`,
#' including `NA` propagation, without one R call per row. The E-steps call
#' this every iteration, where the per-row version dominated their cost.
#'
#' @param x A numeric matrix.
#' @return A numeric vector with one value per row.
#' @noRd
.multilpa_row_max <- function(x) {
  stopifnot("`x` must be a numeric matrix" = is.matrix(x) && is.numeric(x))
  if (ncol(x) == 0L) return(rep(-Inf, nrow(x)))
  do.call(pmax, lapply(seq_len(ncol(x)), function(column) x[, column]))
}

#' Raise a classed message
#'
#' Base R builds classed errors and warnings (`errorCondition()`,
#' `warningCondition()`) but has no constructor for a message, so a notice
#' that a caller may want to silence by class is built here.
#' @param text The message text.
#' @param class The condition class, first in the class vector.
#' @return `NULL`, invisibly. Called for its side effect.
#' @noRd
.multilpa_notice <- function(text, class) {
  condition <- structure(list(message = paste0(text, "\n"), call = NULL),
                         class = c(class, "message", "condition"))
  message(condition)
  invisible(NULL)
}

#' Say that a two-level verb fitted a single-level model
#'
#' `multilpa()`, `multilca()` and `enumerate_classes()` are the verbs for
#' observations nested in groups. Given `id = NULL` they fit the single-level
#' model, which is a legitimate choice, so this is a one-line message rather
#' than a warning, naming the model that was estimated. `lpa()`, `lca()`,
#' `enumerate_lpa()` and `enumerate_lca()` muffle it, since their names
#' already say it.
#'
#' @param latent_class Whether every indicator is categorical.
#' @return `NULL`, invisibly, after signalling a `latents_single_level`
#'   message.
#' @noRd
.multilpa_single_level_notice <- function(latent_class = FALSE) {
  analysis <- if (isTRUE(latent_class)) "latent class" else "latent profile"
  .multilpa_notice(sprintf("Fitted a single-level %s model.", analysis),
    class = "latents_single_level")
  invisible(NULL)
}

#' Refuse a two-level verb called without `id`
#'
#' R's own error for a missing argument names the argument but not what to do.
#' A caller who left `id` out usually has single-level data and wants the verb
#' that fits it by name.
#'
#' @param verb The two-level verb, and `single` the single-level verbs to
#'   point to.
#' @return Does not return; raises `latents_bad_argument`.
#' @noRd
.multilpa_missing_id <- function(verb, single) {
  stop(errorCondition(sprintf(
    "`%s()` needs `id` (grouping column). For single-level data use %s.",
    verb, single),
    class = "latents_bad_argument", call = NULL))
}
