# three_step(method = "ml", step_one = ) against Latent GOLD 6.1. Latent GOLD
# fitted the step-three model of a numeric and a categorical outcome with the
# classification-error and class-size logits read as fixed parameters, once as
# plain values and once with their Jacobian and the step-one covariance from
# an independent likelihood (the first-order correction of Bakk, Oberski and
# Vermunt 2014). Its own step-one run wrote its classification-error file. The
# listings are kept in returned/, so this runs offline. See
# latentgold-step-one/README.md.

base <- "latentgold-step-one"
source(file.path("latentgold", "compare-functions.R"), local = TRUE)
source(file.path(base, "common.R"), local = TRUE)
source("package-kit-runner.R", local = TRUE)

test_that("the kit is the one Latent GOLD was given", {
  recorded <- readLines(file.path(base, "returned", "INPUTS.md5"))
  invisible(lapply(step_cases, function(case) {
    written <- utils::read.delim(file.path(base, "kit", sprintf("%s.dat", case)))
    expect_equal(written, step_data_file(step_case(case)), tolerance = 1e-12,
                 ignore_attr = TRUE, label = case)
  }))
  files <- sort(list.files(file.path(base, "kit"), pattern = "[.](dat|lgs|txt)$",
                           full.names = TRUE), method = "radix")
  expect_identical(sprintf("%s  %s", unname(tools::md5sum(files)), basename(files)), recorded)
})

test_that("step-one uncertainty matches Latent GOLD's correction", {
  snapshot <- matrix_capture_source("..")
  comparison <- step_comparison(file.path(base, "returned"),
                                readRDS(file.path(base, "kit", "independent.rds")))
  kit_save_comparison(base, comparison, snapshot)
  differing <- subset(comparison, status != "agree" & status != "not comparable")
  expect_identical(nrow(differing), 0L, label = paste(differing$quantity, collapse = "; "))
  # 38 rows per step_one setting and 4 of Latent GOLD's own step one (s02's
  # own step one is at another maximum and is not comparable).
  expect_identical(sum(comparison$status == "agree"), 156L)
  corrected <- subset(comparison, step_one == "estimated" &
                        startsWith(quantity, "robust standard error"))
  plain <- subset(comparison, step_one == "fixed" &
                    startsWith(quantity, "robust standard error"))
  # The correction is not vacuous: it raises the standard errors.
  expect_true(all(corrected$latents >= plain$latents))
  expect_gt(stats::median(corrected$latents / plain$latents), 1.1)
})
