# The paquid growth mixtures of Proust-Lima, Philipps and Liquet (2017), the
# current lcmm on paquid and on its data_hlme (growth-published/). Tolerances
# are fixed in growth-published/compare-functions.R.
# Disagreements are allowed only where growth-published/known-differences.csv explains
# them, and every entry there must still match one.

invisible(lapply(c("lcmm", "withr"), function(p) skip_if_not_installed(p)))

growth_published_dir <- file.path("growth-published")
source(file.path(growth_published_dir, "compare-functions.R"), local = TRUE)

source("package-kit-runner.R", local = TRUE)

test_that("growth-published agrees with its references", {
  known <- utils::read.csv(file.path(growth_published_dir, "known-differences.csv"),
                           stringsAsFactors = FALSE)
  comparison <- compare_all_growth_published(growth_published_dir)
  registered <- apply_register(comparison, known)
  expect_identical(nrow(comparison), 37L)
  unexplained <- subset(registered$comparison, status == "disagrees")
  expect_identical(nrow(unexplained), 0L,
                   info = paste(unexplained$case, unexplained$reference, unexplained$quantity,
                                collapse = "; "))
  expect_identical(nrow(registered$stale), 0L)
})
