# Enumeration against Mplus (montecarlo/enumeration/): two Monte Carlo datasets
# fitted by genuine Mplus 9 at every number of classes. Each candidate of an
# enumeration grid is set against Mplus's run at its k; the grid's minima and
# compare_models()'s BIC weights against those Mplus's criteria imply; and a
# grid fitted on two cores against one fitted on one.
#
# Where Mplus replicated its best log likelihood, the candidate must match it.
# Where it did not (the 4- and 5-profile models, 400/80 starts), Mplus warns
# that its solution "may not be trustworthy", and nothing is asserted: the
# over-extracted likelihood has many local maxima. Recorded on 2026-10-18 for
# the 4-profile model: Mplus's best, -3105.390, has a class of 2 observations
# with a variance of 0.000 (a degenerate solution), its second, -3106.613,
# was reached twice; this package reaches -3108.709 from 30 to 300 starts and
# -3106.051 from 1000. Its starts are converged k-means partitions, which
# explore such a surface less widely than Mplus's screened random starts; see
# the roadmap (Priority E).

enumeration_dir <- file.path(montecarlo_dir, "enumeration")

enumeration_mplus <- function(name, k) {
  lines <- read_mplus_lines(file.path(enumeration_dir, sprintf("%s-k%d.out", name, k)))
  npar <- as.integer(printed_value(lines, "Number of Free Parameters"))
  values <- scan(file.path(enumeration_dir, sprintf("%s-k%d-results.dat", name, k)),
                 quiet = TRUE)
  stopifnot(isTRUE(all.equal(values[2L * npar + 1L], npar)))
  list(k = k, n_parameters = npar, log_likelihood = values[2L * npar + 2L],
       aic = values[2L * npar + 3L], bic = values[2L * npar + 4L],
       sabic = values[2L * npar + 5L],
       replicated = k == 1L || any(grepl("HAS BEEN REPLICATED", lines)))
}

enumeration_grid <- function(case, ks, cores = 1L) {
  spec <- case$spec
  if (isTRUE(spec$categorical)) {
    return(enumerate_lca(case$data, spec$vars, n_classes = ks, n_starts = 30,
                         seed = 1, tol = 1e-10, max_iter = 5000, cores = cores))
  }
  enumerate_lpa(case$data, spec$vars, n_profiles = ks, model = spec$model,
                n_starts = 30, seed = 1, tol = 1e-10, max_iter = 5000, cores = cores)
}

test_that("every candidate of a grid reproduces Mplus at its number of classes", {
  designs <- readRDS(file.path(enumeration_dir, "designs.rds"))
  cases <- montecarlo_cases()
  invisible(lapply(names(designs), function(name) {
    ks <- designs[[name]]
    grid <- suppressWarnings(enumeration_grid(cases[[name]], ks))
    table <- as.data.frame(grid)
    references <- lapply(ks, enumeration_mplus, name = name)
    invisible(lapply(seq_along(ks), function(i) {
      reference <- references[[i]]
      row <- subset(table, n_profiles == ks[i])
      label <- sprintf("%s k = %d", name, ks[i])
      expect_identical(as.integer(row$n_parameters), reference$n_parameters, label = label)
      gap <- row$log_likelihood - reference$log_likelihood
      if (reference$replicated) {
        expect_gt(gap, -5e-4, label = paste(label, "not below Mplus"))
        expect_lt(abs(gap), 5e-4, label = paste(label, "LL"))
        expect_lt(abs(row$aic - reference$aic), 1e-3, label = paste(label, "AIC"))
        expect_lt(abs(row$bic_individual - reference$bic), 1e-3, label = paste(label, "BIC"))
        expect_lt(abs(row$sabic_individual - reference$sabic), 1e-3,
                  label = paste(label, "SABIC"))
      }
    }))
    # The minimum of BIC falls where Mplus's BIC has it, among replicated runs.
    replicated <- Filter(function(reference) reference$replicated, references)
    mplus_best <- replicated[[which.min(vapply(replicated, `[[`, numeric(1), "bic"))]]$k
    best <- get_results(grid, "best")
    picked <- subset(best, criterion == "bic" & convention == "individuals")
    expect_identical(as.integer(picked$n_profiles), as.integer(mplus_best), label = name)
  }))
})

test_that("compare_models() BIC weights follow Mplus's BIC", {
  designs <- readRDS(file.path(enumeration_dir, "designs.rds"))
  cases <- montecarlo_cases()
  name <- "lca_binary_01"
  ks <- designs[[name]]
  grid <- suppressWarnings(enumeration_grid(cases[[name]], ks))
  comparison <- get_results(compare_models(grid), "comparison")
  bic <- vapply(ks, function(k) enumeration_mplus(name, k)$bic, numeric(1))
  weights <- exp(-(bic - min(bic)) / 2) / sum(exp(-(bic - min(bic)) / 2))
  expect_equal(comparison$bic_weight[order(comparison$n_profiles)], weights,
               tolerance = 1e-3)
  expect_equal(comparison$delta_bic[order(comparison$n_profiles)], bic - min(bic),
               tolerance = 1e-3)
})

test_that("a grid fitted on two cores is identical to one fitted on one", {
  skip_on_os("windows")
  case <- montecarlo_cases()$lpa_vvi_01
  strip <- function(grid) {
    grid$call <- NULL
    grid$fits <- lapply(grid$fits, function(fit) { fit$call <- NULL; fit })
    grid
  }
  one <- suppressWarnings(enumeration_grid(case, 1:4, cores = 1L))
  two <- suppressWarnings(enumeration_grid(case, 1:4, cores = 2L))
  expect_identical(strip(two), strip(one))
})
