# Attribution rules shared by the verification matrix (build.R) and the CVS
# reporter (../cvs/cvs-report.R): which reference, of which kind, a harness
# quantity, a testthat block or an argument-sweep check compares with. One
# definition, so the matrix and CVS can never disagree about what a check is
# evidence for.
#
# The rule tables live beside this file: harness-references.csv (a regex on a
# harness `source`), harness-functions.csv (suite plus a regex on the quantity),
# test-attribution.csv (directory, file and a regex on the test name; `;`
# lists several functions, references and kinds), and sweep-functions.csv and
# sweep-references.csv (module plus a regex on "<case>: <quantity>"; the sweep
# writes the result into its results.csv). The first matching rule wins.

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

#' Which rule of a sweep table each sweep check matches
#'
#' The argument sweep's tables (sweep-functions.csv, sweep-references.csv)
#' match a regex on the check's title, "<case>: <quantity>", within its module
#' (`module` "*" applies to every module). The first matching rule wins.
#' @param module,title Module names and check titles, of equal length.
#' @param rules One sweep rule table.
#' @return Integer row indices into `rules`, NA where no rule matches.
#' @noRd
matrix_sweep_rule <- function(module, title, rules) {
  stopifnot("`module` and `title` have equal length" = length(module) == length(title),
            "the rules have module and pattern columns" = all(c("module", "pattern") %in% names(rules)))
  key <- paste(module, title, sep = "\r")
  distinct <- !duplicated(key)
  hits <- vapply(which(distinct), function(i) {
    own <- which(rules$module %in% c(module[[i]], "*"))
    hit <- .matrix_first_match(title[[i]], rules$pattern[own])
    if (is.na(hit)) NA_integer_ else own[[hit]]
  }, integer(1))
  hits[match(key, key[distinct])]
}

#' The function, reference and kind of every argument-sweep check
#'
#' @param module,case,quantity The sweep's columns of the same names.
#' @param dir Directory holding sweep-functions.csv and sweep-references.csv.
#' @return A data.frame, one row per check: `fn` ("" for a check of no matrix
#'   function, such as a plot), `reference`, `kind`.
#' @section Errors: a check no rule of either table matches raises
#'   `latents_sweep_unattributed`, so a new check cannot silently drop out of
#'   the matrix or reach CVS untagged.
#' @noRd
matrix_sweep_attribution <- function(module, case, quantity, dir) {
  title <- paste0(case, ": ", quantity)
  functions <- matrix_read_rules("sweep-functions.csv", dir)
  references <- matrix_read_rules("sweep-references.csv", dir)
  functions$fn[is.na(functions$fn)] <- ""
  fn_hit <- matrix_sweep_rule(module, title, functions)
  ref_hit <- matrix_sweep_rule(module, title, references)
  missing <- is.na(fn_hit) | is.na(ref_hit)
  if (any(missing)) {
    stop(errorCondition(sprintf(
      "%d argument-sweep check(s) match no rule of sweep-functions.csv / sweep-references.csv, e.g. %s: %s",
      sum(missing), module[missing][[1L]], title[missing][[1L]]),
      class = "latents_sweep_unattributed", call = NULL))
  }
  data.frame(fn = functions$fn[fn_hit], reference = references$reference[ref_hit],
             kind = references$kind[ref_hit], stringsAsFactors = FALSE)
}
