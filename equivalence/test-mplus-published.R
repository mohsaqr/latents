# lca(), lpa(), lta() and the Vuong-Lo-Mendell-Rubin statistics against
# genuine Mplus output printed on the UCLA OARC Mplus pages (mplus-published/):
# McCutcheon (2002), Vermunt and Magidson (2002) and Langeheine and van de Pol
# (2002) in Applied Latent Class Analysis, Dayton (1998), and the Mplus data
# analysis example of latent class analysis with its TECH11. Tolerances are
# fixed in mplus-published/compare-functions.R. Disagreements are allowed only
# where mplus-published/known-differences.csv explains them.

mplus_published_dir <- file.path("mplus-published")
source(file.path(mplus_published_dir, "compare-functions.R"), local = TRUE)

source("package-kit-runner.R", local = TRUE)

test_that("the Mplus published kit's data are the files UCLA posted", {
  manifest <- utils::read.table(file.path(mplus_published_dir, "DATA.md5"),
                                col.names = c("md5", "file"))
  expect_identical(unname(tools::md5sum(file.path(mplus_published_dir, "data", manifest$file))),
                   manifest$md5)
  expect_setequal(manifest$file, list.files(file.path(mplus_published_dir, "data")))
})

kit_reference_tests(
  "mplus-published", function() compare_all_mplus(mplus_published_dir), n_rows = 286L,
  groups = "Mplus")
