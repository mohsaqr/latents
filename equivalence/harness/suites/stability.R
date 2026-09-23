# Seed stability of the fits the other suites compare.
#
# Every other suite fits its model once, under one seed, and compares the
# result against an external program. That is a claim about two
# implementations only if the fit does not depend on the seed: an agreement
# that holds under seed 1 and fails under seed 2 is measuring the seed.
#
# This suite refits three of the models the other suites rely on -- one
# single-level categorical, one Gaussian, one two-level -- under six seeds
# each, at the same start count the corresponding suite uses, and reports two
# quantities per model:
#
#   * the spread of the maximised log-likelihood across the six seeds, which is
#     zero when every seed reaches the same mode;
#   * the number of distinct optima reached, which is one when it does.
#
# A failure here is not a disagreement with another program. It says the
# corresponding equivalence row was obtained from a multi-start scheme that
# does not reliably reach the same maximum, and that the row's agreement is
# therefore conditional on its seed. That is worth knowing separately from
# whether two likelihoods match, which is why it reports as its own source.
#
# The seeds are fixed, so this is reproducible; they are not "the" seeds, they
# are six of them, which is the point.

.stability_seeds <- c(1L, 2L, 3L, 17L, 101L, 20260920L)

suite_stability <- function() {
  require_suite_packages(
    "poLCA", "multilevLCA",
    reason = "the carcinoma and dataTOY data the refitted models are fitted to")
  rbind(.stability_carcinoma(), .stability_faithful(), .stability_datatoy())
}

#' Single-level categorical: carcinoma, three classes
#'
#' The model behind the `polca` and `published` suites' carcinoma rows.
#'
#' @return A `data.frame` of compared quantities.
.stability_carcinoma <- function() {
  items <- c("A", "B", "C", "D", "E", "F", "G")
  raw <- get(utils::data(list = "carcinoma", package = "poLCA",
                         envir = environment()), envir = environment())
  frame <- raw[, items, drop = FALSE]
  frame$unit <- factor(seq_len(nrow(frame)))
  fit_one <- function(seed) {
    multilpa(frame, vars = items, id = "unit", n_profiles = 3L,
             n_group_classes = 1L, categorical = items, n_starts = 30L,
             seed = seed, tol = 1e-13, max_iter = 20000L)
  }
  seed_stability("carcinoma, 3 classes, 30 starts", fit_one, .stability_seeds,
                 tolerance = .equivalence_tolerances[["likelihood"]]) |>
    transform(source = "multilpa, refitted across seeds", dataset = "carcinoma",
              precision = "machine")
}

#' Gaussian: faithful, three components, equal full covariance
#'
#' The model behind the `mclust` suite's faithful EEE rows.
#'
#' @return A `data.frame` of compared quantities.
.stability_faithful <- function() {
  indicators <- c("eruptions", "waiting")
  frame <- as.data.frame(datasets::faithful)[, indicators, drop = FALSE]
  frame$unit <- factor(seq_len(nrow(frame)))
  fit_one <- function(seed) {
    multilpa(frame, vars = indicators, id = "unit", n_profiles = 3L,
             n_group_classes = 1L, n_starts = 60L, seed = seed, tol = 1e-13,
             max_iter = 20000L, covariance_model = "full",
             variance_model = "equal")
  }
  seed_stability("faithful, 3 components, EEE, 60 starts", fit_one,
                 .stability_seeds,
                 tolerance = .equivalence_tolerances[["likelihood"]]) |>
    transform(source = "multilpa, refitted across seeds", dataset = "faithful",
              precision = "machine")
}

#' Two-level: dataTOY, three profiles in two group classes
#'
#' The model behind the `multilevlca` suite's dataTOY rows, and the only one of
#' the three whose likelihood carries a second level.
#'
#' @return A `data.frame` of compared quantities.
.stability_datatoy <- function() {
  items <- paste0("Y_", seq_len(10))
  observed <- get(utils::data(list = "dataTOY", package = "multilevLCA",
                              envir = environment()), envir = environment())
  fit_one <- function(seed) {
    multilpa(observed, vars = items, id = "id_high", n_profiles = 3L,
             n_group_classes = 2L, categorical = items, n_starts = 30L,
             seed = seed, tol = 1e-12, max_iter = 20000L)
  }
  seed_stability("dataTOY, 3 profiles x 2 group classes, 30 starts", fit_one,
                 .stability_seeds,
                 tolerance = .equivalence_tolerances[["likelihood"]]) |>
    transform(source = "multilpa, refitted across seeds", dataset = "dataTOY",
              precision = "machine")
}
