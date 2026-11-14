# three_step() on a categorical distal outcome against Latent GOLD 6.1's
# step-three model of a nominal dependent variable, for the four methods
# (`step3 modal bch`, `modal none`, `proportional none`, `modal ml`). Latent
# GOLD read latents' step-one posteriors from the kit; its listings are kept in
# returned/, so this runs offline. The class probabilities of every category
# are compared for every method; the robust standard errors and the Wald test
# of the class effect for BCH, modal and ML, with `vcov_type = "independent"`.
# See latentgold-categorical-distal/README.md for what is not comparable.

base <- "latentgold-categorical-distal"
source(file.path("latentgold", "compare-functions.R"), local = TRUE)
source(file.path(base, "common.R"), local = TRUE)
source("package-kit-runner.R", local = TRUE)

test_that("the kit is the one Latent GOLD was given", {
  recorded <- readLines(file.path(base, "returned", "INPUTS.md5"))
  references <- lapply(distal_cases, function(case) distal_case(case))
  # The data files latents writes now reproduce the ones Latent GOLD read.
  invisible(Map(function(case, made) {
    written <- utils::read.delim(file.path(base, "kit", sprintf("%s.dat", case)))
    # The step-one fit is re-estimated here, so its posteriors match the file
    # to the optimizer's reproducibility across platforms, not to the last
    # bit: 1e-12 held under macOS Accelerate and failed under Linux's
    # reference BLAS, where posteriors differ by 1e-10.
    expect_equal(written, distal_data_file(made), tolerance = 1e-8, ignore_attr = TRUE,
                 label = case)
  }, distal_cases, references))
  files <- sort(list.files(file.path(base, "kit"), pattern = "[.](dat|lgs)$", full.names = TRUE),
                method = "radix")
  expect_identical(sprintf("%s  %s", unname(tools::md5sum(files)), basename(files)), recorded)
})

test_that("categorical distal outcomes match Latent GOLD's step-three model", {
  snapshot <- matrix_capture_source("..")
  references <- lapply(stats::setNames(distal_cases, distal_cases), function(case) {
    distal_reference(distal_case(case))
  })
  comparison <- distal_comparison(references, file.path(base, "returned"))
  kit_save_comparison(base, comparison, snapshot)
  differing <- subset(comparison, status != "agree" & status != "not comparable")
  expect_identical(nrow(differing), 0L, label = paste(differing$quantity, collapse = "; "))
  # Every probability of every method, and the robust standard errors and
  # Wald tests of BCH, modal and ML, are compared.
  expect_identical(sum(comparison$status == "agree"), 144L)
  probabilities <- subset(comparison, startsWith(quantity, "probability:"))
  expect_true(all(probabilities$status == "agree"))
})
