# The comparisons kept as scripts, run as tests. Each script refits with the
# current package and stops (stopifnot() or quit(status = 1)) when it
# disagrees with its reference; the latent-transition scripts instead write
# their rows with an `agree` flag. Here each runs in its own R process from the
# project root, and a non-zero exit, or any row that does not agree, fails.
#
# The group-class families are tested in test-latentgold-families.R, against
# Latent GOLD's most refined quadrature (latentgold-families/compare.R uses the
# base runs, whose quadrature error its README documents). References here:
# Mplus TYPE = COMPLEX with WEIGHT (weights/), tidySEM's BCH weights
# (three-step/), Houle et al. (2026) arithmetic and the exhaustive likelihood
# (houle-2026/), and Mplus ex8.13 and ex8.14, Mplus occasion-specific
# measurement, depmixS4 covariate transitions and LMest occasion models
# (lta-extensions/).

run_script <- function(path) {
  root <- normalizePath("..")
  output <- system2(file.path(R.home("bin"), "Rscript"), c("-e", shQuote(sprintf(
    "setwd(%s); source(%s)", deparse(root), deparse(file.path(root, path))))),
    stdout = TRUE, stderr = TRUE)
  list(status = attr(output, "status") %||% 0L, output = output)
}

expect_script_passes <- function(path) {
  result <- run_script(path)
  expect_identical(as.integer(result$status), 0L,
                   label = path,
                   info = paste(utils::tail(result$output, 25L), collapse = "\n"))
}

test_that("weights agree with Mplus TYPE = COMPLEX and the weighted closed form", {
  skip_if_not(dir.exists(file.path("..", "..", "JStats", "validation", "r-reference")))
  expect_script_passes(file.path("equivalence", "weights", "compare-mplus-complex.R"))
})

test_that("BCH weights and corrected means agree with tidySEM", {
  skip_if_not_installed("tidySEM")
  expect_script_passes(file.path("equivalence", "three-step", "compare-tidysem.R"))
})

test_that("Houle et al. (2026) likelihood and fit-index arithmetic agree", {
  expect_script_passes(file.path("equivalence", "houle-2026", "check-likelihood.R"))
  expect_script_passes(file.path("equivalence", "houle-2026", "check-fit-indices.R"))
})

test_that("latent transition extensions agree with Mplus, depmixS4 and LMest", {
  skip_if_not(dir.exists(file.path("..", "..", "JStats", "validation", "r-reference")))
  scripts <- c("compare-ex814.R", "compare-ex813.R", "compare-mplus-measurement.R",
               "compare-depmix-covariate.R", "compare-lmest-occasion.R")
  results <- do.call(rbind, lapply(scripts, function(script) {
    stem <- sub("^compare-(.*)[.]R$", "result-\\1.rds", script)
    path <- file.path("lta-extensions", stem)
    unlink(path)
    run <- run_script(file.path("equivalence", "lta-extensions", script))
    expect_identical(as.integer(run$status), 0L, label = script,
                     info = paste(utils::tail(run$output, 25L), collapse = "\n"))
    # A script that ran but wrote no result must fail, not be skipped.
    expect_true(file.exists(path), label = stem)
    if (file.exists(path)) readRDS(path) else NULL
  }))
  disagreeing <- subset(results, !agree)
  expect_identical(nrow(disagreeing), 0L,
                   info = paste(capture.output(print(disagreeing)), collapse = "\n"))
})

test_that("ordinal and count indicators agree with Latent GOLD", {
  # latentgold-ordinal/compare.R writes its table and has no failure exit of
  # its own, so the table is judged here after a fresh run. Latent GOLD prints
  # the log likelihood to 4 decimals and the profile table to about as many.
  path <- file.path("latentgold-ordinal", "comparison.csv")
  unlink(path)
  expect_script_passes(file.path("equivalence", "latentgold-ordinal", "compare.R"))
  expect_true(file.exists(path))
  comparison <- utils::read.csv(path)
  expect_identical(nrow(comparison), 7L)
  expect_true(all(abs(comparison$ll_difference) < 5e-5 + 1e-4),
              info = paste(capture.output(print(comparison)), collapse = "\n"))
  expect_true(all(comparison$max_parameter_difference < 2e-4),
              info = paste(capture.output(print(comparison)), collapse = "\n"))
})
