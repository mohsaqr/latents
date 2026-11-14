# Growth (ex8.9, a latent class growth analysis of a binary outcome) and
# latent transition models (ex8.12 hidden Markov, ex8.14 with a covariate,
# ex8.15 mover-stayer) against the Mplus 8.8 User's Guide, chapter 8: Mplus's
# own output files in mplus-ug8/listings/, read by
# mplus-ug8/compare-functions.R, where the tolerances are fixed.
# Disagreements are allowed only where mplus-ug8/known-differences.csv
# explains them.

mplus_ug8_dir <- file.path("mplus-ug8")
source(file.path(mplus_ug8_dir, "compare-functions.R"), local = TRUE)

source("package-kit-runner.R", local = TRUE)

test_that("the Mplus User's Guide listings are the files statmodel.com posted", {
  manifest <- utils::read.table(file.path(mplus_ug8_dir, "LISTINGS.md5"),
                                col.names = c("md5", "file"))
  expect_identical(unname(tools::md5sum(file.path(mplus_ug8_dir, "listings", manifest$file))),
                   manifest$md5)
  expect_setequal(manifest$file, list.files(file.path(mplus_ug8_dir, "listings")))
})

kit_reference_tests(
  "mplus-ug8", function() compare_all_ug8(mplus_ug8_dir), n_rows = 114L, groups = "Mplus")
