# Attribution rules shared by the verification matrix (build.R) and the CVS
# reporter (../cvs/cvs-report.R): which reference, of which kind, a harness
# quantity or a testthat block compares with. One definition, so the matrix and
# CVS can never disagree about what a check is evidence for.
#
# The rule tables live beside this file: harness-references.csv (a regex on a
# harness `source`), harness-functions.csv (suite plus a regex on the quantity)
# and test-attribution.csv (directory, file and a regex on the test name; `;`
# lists several functions, references and kinds). The first matching rule wins.

#' Read one rule table
#' @param name File name of the table, e.g. "test-attribution.csv".
#' @param dir Directory holding the tables.
#' @return A data.frame, one row per rule.
#' @noRd
matrix_read_rules <- function(name, dir) {
  path <- file.path(dir, name)
  stopifnot("the rule table exists" = file.exists(path))
  utils::read.csv(path, stringsAsFactors = FALSE)
}

#' Index of the first pattern each text matches (NA when none does)
#' @noRd
.matrix_first_match <- function(text, patterns) {
  vapply(text, function(value) {
    hits <- which(vapply(patterns, function(pattern) grepl(pattern, value, perl = TRUE),
                         logical(1)))
    if (length(hits)) hits[[1L]] else NA_integer_
  }, integer(1), USE.NAMES = FALSE)
}

#' The references and kinds one rule names
#'
#' A rule lists several references separated by `;`, with one kind each or one
#' kind for all of them.
#' @param reference,kind One rule's `reference` and `kind` cells.
#' @return A data.frame with columns `reference` and `kind`, one row per reference.
#' @noRd
matrix_rule_references <- function(reference, kind) {
  stopifnot("one rule at a time" = length(reference) == 1L && length(kind) == 1L)
  refs <- strsplit(reference, ";", fixed = TRUE)[[1L]]
  kinds <- strsplit(kind, ";", fixed = TRUE)[[1L]]
  stopifnot("one kind per reference" = length(kinds) %in% c(1L, length(refs)))
  data.frame(reference = refs, kind = rep_len(kinds, length(refs)), stringsAsFactors = FALSE)
}

#' Which rule of harness-references.csv each harness source matches
#' @param source Harness `source` strings (e.g. "poLCA 1.6.0.2").
#' @param rules The harness-references.csv table.
#' @return Integer row indices into `rules`, NA where no rule matches.
#' @noRd
matrix_harness_reference <- function(source, rules) {
  distinct <- unique(source)
  .matrix_first_match(distinct, rules$pattern)[match(source, distinct)]
}

#' Which rule of harness-functions.csv each harness quantity matches
#' @param suite,quantity Harness suite names and quantities, of equal length.
#' @param rules The harness-functions.csv table.
#' @return Integer row indices into `rules`, NA where no rule of the suite matches.
#' @noRd
matrix_harness_function <- function(suite, quantity, rules) {
  stopifnot("`suite` and `quantity` have equal length" = length(suite) == length(quantity))
  vapply(seq_along(suite), function(i) {
    candidates <- which(rules$suite == suite[[i]])
    hit <- .matrix_first_match(quantity[[i]], rules$pattern[candidates])
    if (is.na(hit)) NA_integer_ else candidates[[hit]]
  }, integer(1))
}

#' Which rule of test-attribution.csv each test matches
#' @param dir,file Directory ("equivalence", "tests/equivalence") and file name
#'   of each test; recycled against `test`.
#' @param test Test names (the `test_that()` description).
#' @param rules The test-attribution.csv table.
#' @return Integer row indices into `rules`, NA where no rule of the file matches.
#' @noRd
matrix_test_rule <- function(dir, file, test, rules) {
  n <- length(test)
  dir <- rep_len(dir, n)
  file <- rep_len(file, n)
  vapply(seq_len(n), function(i) {
    own <- which(rules$dir == dir[[i]] & rules$file == file[[i]])
    hit <- .matrix_first_match(test[[i]], rules$pattern[own])
    if (is.na(hit)) NA_integer_ else own[[hit]]
  }, integer(1))
}
