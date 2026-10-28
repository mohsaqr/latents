library(testthat)
out <- "validation/audit-2026-10-13"
inventory <- read.csv(file.path(out, "function-inventory.csv"), stringsAsFactors = FALSE)
source <- read.csv(file.path(out, "source-files.csv"), stringsAsFactors = FALSE)
docs <- read.csv(file.path(out, "documentation-inventory.csv"), stringsAsFactors = FALSE)
test_that("inventory contains every source file and valid function spans", {
  expect_setequal(source$file, list.files("R", pattern = "[.]R$", full.names = TRUE))
  expect_equal(nrow(source), 80L)
  expect_equal(sum(inventory$definition == "function"), 1023L)
  expect_equal(anyDuplicated(paste(inventory$file, inventory$name)), 0L)
  expect_true(all(inventory$line > 0 & inventory$end_line >= inventory$line))
  expect_true(all(inventory$end_line <= source$lines[match(inventory$file, source$file)]))
  expect_identical(source$md5, unname(tools::md5sum(source$file)))
})
test_that("semantic recency ignores source formatting and recognizes new functions", {
  expect_false(any(inventory$body_changed_since_head[inventory$file == "R/centering.R"]))
  expect_true(all(inventory$body_changed_since_head[inventory$file == "R/selection.R"]))
  expect_match(inventory$arguments[inventory$name == "pool_imputations"], "imputed")
  expect_match(inventory$validation_calls[inventory$name == "pool_imputations"], "stopifnot")
  expect_match(inventory$called_functions[inventory$name == "pool_imputations"], "parameter_inference")
})
test_that("public documentation inventory includes every authored Rd and Rmd source", {
  expected <- c(list.files("man", pattern = "[.]Rd$", full.names = TRUE),
                list.files("vignettes", pattern = "[.]Rmd$", full.names = TRUE, recursive = TRUE))
  expect_true(all(expected %in% docs$file))
  expect_identical(docs$md5, unname(tools::md5sum(docs$file)))
})
