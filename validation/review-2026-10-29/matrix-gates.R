# Read-only verification-matrix outcome-import reproduction.
# Run from the latents repository root.
matrix_environment <- new.env(parent = globalenv())
build_source <- readLines('equivalence/verification-matrix/build.R')
stopifnot(grepl('^print\\(build_verification_matrix', tail(build_source, 1)))
eval(parse(text = head(build_source, -1L)), envir = matrix_environment)

test_path <- tempfile(fileext = '.R')
writeLines(c(
  'testthat::test_that("reference crashed after structural check", { testthat::expect_true(TRUE); stop("reference unavailable") })',
  'testthat::test_that("reference skipped", { testthat::skip("missing dependency") })'
), test_path)
actual_results <- as.data.frame(testthat::test_file(test_path, reporter = 'silent'))
cat('Actual testthat outcomes:\n')
print(actual_results[c('test', 'nb', 'failed', 'error', 'skipped', 'passed')], row.names = FALSE)

results_path <- tempfile(fileext = '.csv')
utils::write.csv(actual_results[c('file', 'test', 'nb', 'failed', 'error', 'skipped', 'passed')], results_path, row.names = FALSE)
manifest <- data.frame(role = 'tests', dir = 'tests/equivalence', path = results_path,
                       run = '2026-10-29', code = '1.0.0', current = TRUE)
imported <- matrix_environment$.matrix_test_results(manifest)
cat('\nImported outcomes (error and skipped are absent):\n')
print(imported, row.names = FALSE)

evidence <- data.frame(
  evidence_id = paste0('outcome-probe-', seq_len(nrow(imported))), origin = 'testthat',
  location = 'tests/equivalence/test-weighted-logits.R', item = imported$test,
  fn = 'multilpa', reference = c('crashed reference', 'skipped reference'),
  kind = 'R package', checks = imported$nb, values = NA_real_, failed = imported$failed,
  run = imported$run, code = imported$code, current = imported$current,
  dataset = NA_character_, shared = FALSE, known = FALSE, one_class = FALSE)
cells <- matrix_environment$.matrix_cells(
  evidence, matrix_environment$.matrix_rules('candidates.csv'))
probe_cells <- cells[cells$reference %in% evidence$reference,
                     c('reference', 'current_checks', 'current_failed', 'status')]
cat('\nMatrix treatment of the actual error and skip:\n')
print(probe_cells, row.names = FALSE)
stopifnot(all(probe_cells$current_checks == 1L),
          all(probe_cells$current_failed == 0L), all(probe_cells$status == 'current'))
cat('\nCONFIRMED: both errored and skipped tests receive current source credit with zero failed checks.\n')
unlink(c(test_path, results_path))
