# The Mplus 8.8 User's Guide examples latents fits, against Mplus's own output
# files (mplus-ug/listings/, read by mplus-ug/compare-functions.R, where the
# tolerances are fixed): chapter 7 -- mixture regression (ex7.1), latent class
# models (ex7.3, 7.4, 7.6, 7.7), latent profile models (ex7.9, 7.10, 7.22) --
# and chapter 8 -- latent class growth analysis (ex8.9), hidden Markov
# (ex8.12), latent transition with a covariate (ex8.14) and mover-stayer
# (ex8.15). Disagreements are allowed only where
# mplus-ug/known-differences.csv explains them.

mplus_ug_dir <- file.path("mplus-ug")
source(file.path(mplus_ug_dir, "compare-functions.R"), local = TRUE)

source("package-kit-runner.R", local = TRUE)

test_that("the Mplus User's Guide listings are the files statmodel.com posted", {
  manifest <- utils::read.table(file.path(mplus_ug_dir, "LISTINGS.md5"),
                                col.names = c("md5", "file"))
  expect_identical(unname(tools::md5sum(file.path(mplus_ug_dir, "listings", manifest$file))),
                   manifest$md5)
  expect_setequal(manifest$file, list.files(file.path(mplus_ug_dir, "listings")))
})

kit_reference_tests(
  "mplus-ug", function() compare_all_ug(mplus_ug_dir), n_rows = 516L, groups = "Mplus")
