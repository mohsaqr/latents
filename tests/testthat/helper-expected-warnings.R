# Muffle only the warnings the package raises on purpose, and only the ones a
# deliberately small test fixture is expected to trigger. Every other warning
# still reaches testthat and is reported.
#
# The suite previously wrapped these calls in `suppressWarnings()`, which also
# hides a warning that a change newly introduced -- precisely the signal the
# suite exists to carry. `withCallingHandlers()` observes and continues instead,
# so an unexpected warning stays visible without failing an unrelated test.

# Warnings a small fixture legitimately produces: EM that has not converged in
# the iteration budget, individual restarts that failed, a variance or
# probability resting on its boundary, and classes too small to be well
# separated.
.multilpa_expected_warnings <- c(
  "latents_unconverged",
  "latents_failed_starts",
  "latents_failed_replicates",
  "latents_boundary",
  "latents_small_classes",
  # A deliberately small or separated fixture also reaches these two.
  "latents_extreme_coefficients",
  # Not a defect: a single-level fixture is a deliberate `id = NULL`.
  "latents_single_level",
  "latents_empty_transition_row"
)
# Deliberately NOT muffled by default, because each reports a result the suite
# should not quietly accept: `multilpa_quadrature_check` (the integral is not
# accurate), `latents_unverified_alignment` (the caller's rows could not be
# checked) and anything unclassed.

#' Evaluate `expr`, muffling only expected package warnings
#'
#' @param expr Code to evaluate.
#' @param classes Character vector of condition classes to muffle.
#' @return The value of `expr`.
quietly <- function(expr, classes = .multilpa_expected_warnings) {
  stopifnot("`classes` must be a character vector" = is.character(classes))
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (inherits(w, classes)) invokeRestart("muffleWarning")
    }
  )
}
