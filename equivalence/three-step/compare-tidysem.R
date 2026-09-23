# Independent-implementation check for the BCH three-step correction.
#
# tidySEM (van Lissa) implements Bolck, Croon & Hagenaars (2004) for mixture
# models fitted with OpenMx. Its weight construction is a pure function of the
# posterior probabilities, so it can be handed this package's posteriors and the
# two implementations compared exactly rather than approximately.
#
# Run from the project root:
#   Rscript equivalence/three-step/compare-tidysem.R

suppressMessages(pkgload::load_all(".", quiet = TRUE))
# A missing external package is a skip with a reason, not a silent absence and
# not an error that reads like a defect in this package.
if (!requireNamespace("tidySEM", quietly = TRUE)) {
  cat("SKIPPED: tidySEM is not installed, so the independent BCH weight",
      "construction this comparison is made against is unavailable.\n")
  quit(status = 0L, save = "no")
}

set.seed(11)
groups <- rep(seq_len(50), each = 10L)
n <- length(groups)
truth <- 1L + as.integer(stats::runif(n) > 0.5)
data <- data.frame(
  g = groups,
  a = stats::rnorm(n, ifelse(truth == 2L, 1.1, -1.1)),
  b = stats::rnorm(n, ifelse(truth == 2L, 1.1, -1.1)),
  y = stats::rnorm(n, ifelse(truth == 2L, 10, 0))
)
fit <- multilpa(data, c("a", "b"), "g", n_profiles = 2, n_group_classes = 1,
                n_starts = 6, seed = 1)
posteriors <- fit$subject_posteriors

mine <- latents:::.multilpa_error_matrix(
  latents:::.multilpa_level_assignments(fit, "individuals"))
theirs <- as.matrix(tidySEM:::classification_probs_mostlikely(posteriors))
# `bch_weights()` returns one row per unit and class since the public returns
# were tidied; it used to return one `weight_class_k` column per class. The
# comparison is against a unit-by-class matrix, so the long frame is reshaped
# here rather than the comparison being weakened to fit the new shape.
#' The BCH weights as a unit-by-class matrix in the fit's own row order
#' @param frame The long frame returned by `bch_weights()`.
#' @return A numeric matrix with one row per unit and one column per class.
.bch_weight_matrix <- function(frame) {
  wide <- stats::reshape(frame, direction = "wide", idvar = "unit",
                         timevar = "class", v.names = "weight",
                         drop = c("level", "assigned_class"))
  ordered_rows <- wide[order(wide$unit), , drop = FALSE]
  unname(as.matrix(ordered_rows[, setdiff(names(ordered_rows), "unit"), drop = FALSE]))
}
weights_mine <- .bch_weight_matrix(get_results(fit, "bch_weights", level = "individuals"))
weights_theirs <- solve(theirs)[apply(posteriors, 1L, which.max), ]
estimate_mine <- three_step(fit, data, "y", method = "bch")$estimate
estimate_theirs <- vapply(seq_len(2), function(class) {
  sum(weights_theirs[, class] * data$y) / sum(weights_theirs[, class])
}, numeric(1))

comparison <- data.frame(
  quantity = c("classification error matrix", "BCH weights", "corrected class means"),
  max_difference = c(max(abs(unname(mine) - unname(theirs))),
                     max(abs(unname(weights_mine) - unname(weights_theirs))),
                     max(abs(estimate_mine - estimate_theirs)))
)
print(comparison, row.names = FALSE)

stopifnot(
  "the error matrices must agree to machine precision" =
    max(abs(unname(mine) - unname(theirs))) < 1e-12,
  "the weights must agree to machine precision" =
    max(abs(unname(weights_mine) - unname(weights_theirs))) < 1e-12,
  "the corrected means must agree to machine precision" =
    max(abs(estimate_mine - estimate_theirs)) < 1e-10
)
cat("\nAll three agree to machine precision.\n")
cat("tidySEM version:", format(utils::packageVersion("tidySEM")), "\n")
