# The group-class families (additive, dispersion, additive-dispersion) against
# Latent GOLD 6.1 (latentgold-families/). Latent GOLD fits each family with a
# group-level random intercept integrated by Gauss-Hermite quadrature; this
# package evaluates the same likelihood exactly. The base runs (50 nodes, 15
# per factor with two indicators) carry Latent GOLD's quadrature error, which
# the folder's README diagnoses: it vanishes as the nodes are refined. Each
# case is therefore compared with its most refined retained run: 100 nodes for
# g01, g03 and g06, 30 per factor for g07, 60 per factor for g08, the base run
# for g02, g04 and g05.
#
# Tolerances: the log likelihood within 2e-4 (4 printed decimals plus the
# quadrature error left at the finest grid, measured at 1.6e-4 for g08 on
# 2026-10-01, above the 1e-4 declared before any run); group posteriors within
# 1e-3; parameter counts exactly.

families_dir <- file.path("latentgold-families")

families_listing <- function(name) {
  folders <- file.path(families_dir, "returned", c("refine60", "refine", "."))
  found <- folders[file.exists(file.path(folders, paste0(name, ".lst")))][1L]
  lines <- iconv(readLines(file.path(found, paste0(name, ".lst")), warn = FALSE),
                 "latin1", "UTF-8")
  posteriors <- utils::read.delim(file.path(found, paste0(name, "_posteriors.txt")),
                                  check.names = FALSE, fileEncoding = "latin1")
  columns <- grep("^GClass#[0-9]+$", names(posteriors), value = TRUE)
  groups <- posteriors[!duplicated(posteriors[[1L]]), , drop = FALSE]
  list(folder = basename(found),
       log_likelihood = printed_value(lines, "^Log-likelihood \\(LL\\)"),
       n_parameters = as.integer(printed_value(lines, "^Number of parameters \\(Npar\\)")),
       group = groups[[1L]], posteriors = as.matrix(groups[columns]))
}

test_that("every group-class family reproduces Latent GOLD at its finest quadrature", {
  cases <- readRDS(file.path(families_dir, "kit", "cases.rds"))
  expect_identical(length(cases), 8L)
  invisible(lapply(names(cases), function(name) {
    case <- cases[[name]]
    data <- utils::read.delim(file.path(families_dir, "kit", paste0(name, ".dat")))
    fit <- suppressWarnings(multilpa(
      data, paste0("y", seq_len(case$d)), "g", n_group_classes = case$H,
      family = case$family, between_variance = case$between,
      tol = 1e-14, max_iter = 20000, n_starts = 20, seed = 1))
    lg <- families_listing(name)
    label <- sprintf("%s (%s)", name, lg$folder)
    expect_identical(as.integer(fit$n_parameters), lg$n_parameters, label = label)
    expect_lt(abs(fit$log_likelihood - lg$log_likelihood), 2e-4, label = label)
    mine <- fit$group_posteriors[as.character(lg$group), , drop = FALSE]
    expect_lt(montecarlo_relabel(mine, lg$posteriors)$gap, 1e-3, label = label)
  }))
})
