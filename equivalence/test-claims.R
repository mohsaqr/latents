# The claims registry (CLAIMS.csv): every sentence in the documentation that
# says this package agrees with another program or a published result, with
# the test that backs it. A claim whose wording changes, or whose test is
# renamed or removed, fails here, so the registry cannot drift from either.
# A test file named with a path is read from the project root (tests outside
# this folder); a bare name is a file in this folder.

claims_root <- ".."
claims <- utils::read.csv("CLAIMS.csv")

claims_source_text <- function(source) {
  lines <- readLines(file.path(claims_root, source), warn = FALSE)
  # Roxygen lines are compared without their comment marker.
  if (startsWith(source, "R/")) lines <- sub("^#' ?", "", lines)
  lines
}

claims_test_path <- function(test_file) {
  if (grepl("/", test_file, fixed = TRUE)) file.path(claims_root, test_file) else test_file
}

claims_test_present <- function(test_file, test_name) {
  text <- paste(readLines(claims_test_path(test_file), warn = FALSE), collapse = "\n")
  # Tests generated in a loop are registered under their full name; the file
  # holds the stem before the last ": ".
  stem <- sub(": [^:]*$", "", test_name)
  grepl(test_name, text, fixed = TRUE) || grepl(stem, text, fixed = TRUE)
}

test_that("the registry has its columns and no duplicate claims", {
  expect_identical(names(claims),
                   c("source", "phrase", "reference", "test_file", "test_name", "note"))
  expect_false(anyNA(claims))
  expect_false(anyDuplicated(claims[c("source", "phrase", "reference")]) > 0L)
})

test_that("every claim is still worded as registered in its source", {
  invisible(Map(function(source, phrase) {
    expect_true(file.exists(file.path(claims_root, source)), label = source)
    found <- any(grepl(phrase, claims_source_text(source), fixed = TRUE))
    expect_true(found, label = sprintf("%s: %s", source, phrase))
  }, claims$source, claims$phrase))
})

test_that("every backed claim names a test that exists", {
  backed <- subset(claims, test_file != "NONE")
  invisible(Map(function(test_file, test_name) {
    expect_true(file.exists(claims_test_path(test_file)), label = test_file)
    expect_true(claims_test_present(test_file, test_name),
                label = sprintf("%s: %s", test_file, test_name))
  }, backed$test_file, backed$test_name))
})

test_that("unbacked claims are explained and do not grow unnoticed", {
  unbacked <- subset(claims, test_file == "NONE")
  expect_true(all(nzchar(unbacked$note)))
  # One: the `boundary = "fix"` convention described as the one Mplus
  # follows, which states a convention rather than a numerical agreement.
  expect_identical(nrow(unbacked), 1L)
})
