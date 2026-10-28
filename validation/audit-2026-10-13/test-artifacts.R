library(testthat)
out <- "validation/audit-2026-10-13"
report <- readLines(file.path(out, "REVIEW.md"))
test_that("report links point to existing artifacts and sources", {
  text <- paste(report, collapse = "\n")
  targets <- regmatches(text, gregexpr("\\]\\([^)]*\\)", text))[[1L]]
  targets <- sub("^\\]\\(", "", sub("\\)$", "", targets))
  local <- targets[!grepl("^https?://", targets)]
  expect_true(all(file.exists(file.path(out, local))))
})
test_that("review destinations exist for every function", {
  x <- read.csv(file.path(out, "function-inventory.csv"))
  expect_false(any(grepl("(expressions|coverage_percent|unit_suite)[.][0-9]+$", names(x))))
  paths <- unique(x$specialist_report)
  paths[paths == "REVIEW.md shared services and coverage"] <- "REVIEW.md"
  expect_true(all(file.exists(file.path(out, paths))))
})
test_that("report counts match saved source suite results", {
  x <- read.csv(file.path(out, "unit-tests.csv"))
  y <- read.csv(file.path(out, "equivalence-tests.csv"))
  expect_equal(sum(x$passed) + sum(y$passed), 6349)
  expect_equal(nrow(x) + nrow(y), 958L)
  expect_equal(sum(x$failed) + sum(y$failed), 0)
  expect_false(any(x$error) || any(y$error))
  expect_equal(sum(x$warning) + sum(y$warning), 0)
  expect_equal(sum(x$skipped) + sum(y$skipped), 0)
})
test_that("numbered findings comprise 27 distinct audit entries", {
  ids <- unique(unlist(regmatches(report, gregexpr("(SV|MC|RA|PD)[0-9]{2}", report))))
  expect_equal(length(ids), 27L)
})
