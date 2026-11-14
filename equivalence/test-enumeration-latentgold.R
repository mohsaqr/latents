# The enumeration layer against Latent GOLD 6.1. enumerate_*() and
# compare_models() fit the core models (lca(), lpa(), each verified against
# Latent GOLD on its own) and add a layer: the information criteria and
# classification statistics of every candidate. This file checks that layer
# against what Latent GOLD prints for the same models: the coverage kit's two
# enumeration sequences (latentgold-coverage/, e1: six binary items at one to
# four classes, n = 800; e2: three VVI indicators at one to four profiles,
# n = 600), each latents fit started at Latent GOLD's solution
# (fit_at_reference()) so both programs describe the same maximum, then
# passed to compare_models().
#
# Latent GOLD prints, per model: LL, Npar, BIC, AIC, CAIC and SABIC "based on
# LL", the classification log-likelihood CL, the entropy EN and Entropy
# R-squared. latents' CLC is -2 CL = -2 LL + 2 EN (Biernacki, Celeux and
# Govaert 2000). Latent GOLD's Entropy R-squared is 1 - EN / (n H), with H the
# entropy of the class shares (Vermunt and Magidson's definition), while
# latents' profile_entropy is Mplus's relative entropy 1 - EN / (n log K); the
# EN both are built from is compared through CLC, and Latent GOLD's statistic
# is rebuilt from latents' EN and class shares. Tolerances, fixed before the
# comparison: the log likelihood to 6e-5 (the coverage kit's), so every
# -2LL-based criterion to 1.2e-4 plus half a unit of the fourth decimal
# (1.7e-4); the parameter count exactly; CLC, twice a sum of up to 800
# posterior terms each within 2e-5, to 2e-3; Entropy R-squared to 1e-4.
#
# test-enumeration-identity.R shows that an enumeration candidate is the direct
# fit and that enumerate_*() reports compare_models()'s criteria, so this
# comparison holds for enumerate_lca() and enumerate_lpa() as well.

coverage_dir <- file.path("latentgold-coverage")
source(file.path(coverage_dir, "compare-functions.R"), local = TRUE)

#' Latent GOLD's printed criteria of one listing
#' @noRd
lg_criteria <- function(name) {
  lines <- sub("\r$", "", iconv(readLines(file.path(coverage_dir, "returned", paste0(name, ".lst")),
                                          warn = FALSE), "latin1", "UTF-8"))
  value <- function(label) .lg_value(lines, label)
  c(log_likelihood = value("Log-likelihood \\(LL\\)"),
    n_parameters = value("Number of parameters \\(Npar\\)"),
    bic = value("BIC \\(based on LL\\)"), aic = value("AIC \\(based on LL\\)"),
    caic = value("CAIC \\(based on LL\\)"), sabic = value("SABIC \\(based on LL\\)"),
    entropy_r2 = value("Entropy R-squared"),
    classification_ll = value("Classification log-likelihood"))
}

#' compare_models() on the four fits of a sequence, against Latent GOLD
#' @noRd
check_sequence <- function(prefix) {
  cases <- readRDS(file.path(coverage_dir, "kit", "cases.rds"))
  names <- sprintf("%sk%d", prefix, 1:4)
  fits <- lapply(names, function(name) {
    fit_at_reference(cases[[name]], read_listing_any(file.path(coverage_dir, "returned",
                                                               paste0(name, ".lst"))))
  })
  table <- as.data.frame(do.call(compare_models, stats::setNames(fits, names)))
  expect_identical(nrow(table), 4L)
  invisible(lapply(seq_along(names), function(k) {
    printed <- lg_criteria(names[[k]])
    row <- table[k, ]
    label <- function(what) sprintf("%s %s", names[[k]], what)
    expect_lte(abs(row$log_likelihood - printed[["log_likelihood"]]), 6e-5,
               label = label("log likelihood"))
    expect_identical(as.integer(row$n_parameters), as.integer(printed[["n_parameters"]]),
                     label = label("parameters"))
    expect_lte(abs(row$bic_individual - printed[["bic"]]), 1.7e-4, label = label("BIC"))
    expect_lte(abs(row$aic - printed[["aic"]]), 1.7e-4, label = label("AIC"))
    expect_lte(abs(row$caic_individual - printed[["caic"]]), 1.7e-4, label = label("CAIC"))
    expect_lte(abs(row$sabic_individual - printed[["sabic"]]), 1.7e-4, label = label("SABIC"))
    if (k > 1L) {
      expect_lte(abs(row$clc_individual - (-2 * printed[["classification_ll"]])), 2e-3,
                 label = label("CLC = -2 CL"))
      entropy <- (row$clc_individual + 2 * row$log_likelihood) / 2
      shares <- get_results(fits[[k]], "profile_probabilities")$probability
      n <- nrow(cases[[names[[k]]]]$data)
      expect_lte(abs(1 - entropy / (n * -sum(shares * log(shares))) - printed[["entropy_r2"]]),
                 1e-4, label = label("entropy R-squared, Latent GOLD's definition"))
    }
  }))
}

test_that("compare_models() criteria for lca() at one to four classes match Latent GOLD", {
  check_sequence("e1")
})

test_that("compare_models() criteria for VVI lpa() at one to four profiles match Latent GOLD", {
  check_sequence("e2")
})
