# enumerate_lca(), enumerate_lpa() and enumerate_classes() are compositions:
# they fit the core models (lca(), lpa(), multilpa()) at each number of
# classes and add the criteria of compare_models(). The core models are
# verified against R packages, Mplus, Latent GOLD and simulab recovery on their
# own; this file checks the hand-off, so that evidence carries over:
#
# - every candidate of a grid is the direct fit: the same log likelihood and
#   the same coefficients, exactly, with the same arguments;
# - the grid's criteria are compare_models()'s on the direct fits, exactly
#   (test-enumeration-latentgold.R checks those against Latent GOLD).

criteria_columns <- c("log_likelihood", "n_parameters", "aic", "kic", "bic_groups",
                      "bic_individual", "sabic_groups", "sabic_individual", "caic_groups",
                      "caic_individual", "awe_groups", "awe_individual", "icl_groups",
                      "icl_individual", "clc_groups", "clc_individual", "profile_entropy",
                      "group_entropy", "avepp_min", "avepp_max", "avepp_mean")

#' Every candidate of `grid` against the direct fits, and the grid's criteria
#' against compare_models() on those fits
#' @noRd
check_grid <- function(grid, direct, names) {
  candidates <- as.data.frame(get_results(grid, "candidates"))
  table <- as.data.frame(do.call(compare_models, stats::setNames(direct, names)))
  invisible(lapply(seq_along(names), function(i) {
    candidate <- get_candidate(grid, names[[i]])
    expect_identical(as.numeric(stats::logLik(candidate)), as.numeric(stats::logLik(direct[[i]])),
                     label = sprintf("%s log likelihood", names[[i]]))
    expect_identical(stats::coef(candidate), stats::coef(direct[[i]]),
                     label = sprintf("%s coefficients", names[[i]]))
    row <- candidates[candidates$candidate == names[[i]], criteria_columns]
    expect_identical(unname(unlist(row)), unname(unlist(table[i, criteria_columns])),
                     label = sprintf("%s criteria", names[[i]]))
  }))
}

#' Evaluate `expr`, muffling only latents_thin_profiles: a VVI fit with one
#' more profile than generated has a profile of a few persons, which the grid
#' and the direct fit both flag; the diagnostic is expected there.
#' @noRd
without_thin_profiles <- function(expr) {
  withCallingHandlers(expr, latents_thin_profiles = function(w) invokeRestart("muffleWarning"))
}

lca_data <- function() {
  probabilities <- array(0, c(3L, 6L, 2L))
  probabilities[, , 1L] <- rbind(rep(0.85, 6L), c(0.85, 0.85, 0.85, 0.15, 0.15, 0.15),
                                 rep(0.15, 6L))
  probabilities[, , 2L] <- 1 - probabilities[, , 1L]
  data <- as.data.frame(simulab::simulate_lca(600L, probabilities = probabilities,
                                              proportions = c(0.5, 0.3, 0.2), seed = 11L))
  items <- grep("^item_", names(data), value = TRUE)
  data[items] <- lapply(data[items], factor)
  data[items]
}

test_that("every enumerate_lca() candidate is the lca() fit, with compare_models()'s criteria", {
  skip_if_not_installed("simulab")
  data <- lca_data()
  items <- names(data)
  grid <- enumerate_lca(data, items, n_classes = 1:4, n_starts = 10L, seed = 1L)
  direct <- lapply(1:4, function(k) lca(data, items, n_classes = k, n_starts = 10L, seed = 1L))
  check_grid(grid, direct, sprintf("classes%d", 1:4))
})

test_that("every enumerate_lpa() candidate is the lpa() fit, with compare_models()'s criteria", {
  skip_if_not_installed("simulab")
  data <- as.data.frame(simulab::simulate_lpa(500L, means = rbind(c(0, 0, 0), c(2, 2, 0),
                                                                  c(0, 2, 2)),
                                              proportions = c(0.4, 0.35, 0.25), seed = 12L))
  vars <- grep("^indicator_", names(data), value = TRUE)
  invisible(lapply(c("EEI", "VVI"), function(model) {
    grid <- without_thin_profiles(enumerate_lpa(data, vars, n_profiles = 1:4, model = model,
                                                n_starts = 10L, seed = 1L))
    direct <- lapply(1:4, function(k) without_thin_profiles(
      lpa(data, vars, n_profiles = k, model = model, n_starts = 10L, seed = 1L)))
    candidates <- as.data.frame(get_results(grid, "candidates"))
    check_grid(grid, direct, candidates$candidate[order(candidates$n_profiles)])
  }))
})

test_that("every enumerate_classes() candidate is the multilpa() fit, with compare_models()'s criteria", {
  skip_if_not_installed("simulab")
  data <- as.data.frame(simulab::simulate_ml_lpa(
    clusters = 40L, cluster_size = 15L, means = rbind(c(0, 0, 0), c(2, 2, 0), c(0, 2, 2)),
    profile_probabilities = rbind(c(0.7, 0.2, 0.1), c(0.1, 0.3, 0.6)),
    cluster_class_proportions = c(0.5, 0.5), seed = 13L))
  vars <- grep("^indicator_", names(data), value = TRUE)
  grid <- enumerate_classes(data, vars, id = "cluster", n_profiles = 2:3, n_group_classes = 1:2,
                            n_starts = 5L, seed = 1L)
  candidates <- as.data.frame(get_results(grid, "candidates"))
  direct <- lapply(seq_len(nrow(candidates)), function(i) {
    # The grid crosses the covariance structures of model = "basic"; each
    # candidate is the direct fit of its own structure.
    multilpa(data, vars, id = "cluster", n_profiles = candidates$n_profiles[[i]],
             n_group_classes = candidates$n_group_classes[[i]], model = candidates$model[[i]],
             n_starts = 5L, seed = 1L)
  })
  check_grid(grid, direct, candidates$candidate)
})
