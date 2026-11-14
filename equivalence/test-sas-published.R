# lca(), lpa(), enumerate_lpa() and mixture_regression() against SAS output:
# PROC LCA as printed by UCLA OARC (the latent class data analysis example and
# McCutcheon's textbook examples) and in SAS and R, Example 8.21, and PROC FMM
# in the SAS/STAT 14.2 User's Guide (sas-published/). Tolerances are fixed in
# sas-published/compare-functions.R. Disagreements are allowed only where
# sas-published/known-differences.csv explains them.

sas_published_dir <- file.path("sas-published")
source(file.path(sas_published_dir, "compare-functions.R"), local = TRUE)

source("package-kit-runner.R", local = TRUE)

test_that("the SAS published kit's data are the files posted and transcribed", {
  manifest <- utils::read.table(file.path(sas_published_dir, "DATA.md5"),
                                col.names = c("md5", "file"))
  expect_identical(unname(tools::md5sum(file.path(sas_published_dir, "data", manifest$file))),
                   manifest$md5)
  expect_setequal(manifest$file, list.files(file.path(sas_published_dir, "data")))
})

kit_reference_tests(
  "sas-published", function() compare_all_sas(sas_published_dir), n_rows = 322L, groups = "SAS")
