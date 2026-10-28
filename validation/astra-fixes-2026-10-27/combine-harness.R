# Consolidate the full run with complete replacements for affected suites.
# Preserve every measured row; do not alter tolerances, values or status.
suppressMessages(pkgload::load_all('.', quiet = TRUE))
source('equivalence/harness/registry.R')
base <- 'validation/astra-fixes-2026-10-27/harness'
replacements <- c('depmix-repair', 'invariance-repair', 'refusals-repair')
read_evidence <- function(directory, filename) {
  out <- utils::read.csv(file.path(base, directory, filename), stringsAsFactors = FALSE)
  out$validation_run <- rep(directory, nrow(out))
  out
}
combine <- function(filename) {
  original <- read_evidence('full', filename)
  newer <- do.call(rbind, lapply(replacements, read_evidence, filename = filename))
  rbind(original[!original$suite %in% unique(newer$suite), , drop = FALSE], newer)
}
run <- list(comparisons = combine('report.csv'), suites = combine('suites.csv'),
            coverage = combine('coverage.csv'))
marked <- .mark_known_open(run$comparisons, .known_open(), run$suites$suite)
run$comparisons <- marked$comparisons
run$known_open <- marked$register
available <- sub('[.]R$', '', basename(list.files('equivalence/harness/suites',
                                                  pattern = '[.]R$', full.names = TRUE)))
stopifnot(!anyDuplicated(run$suites$suite), setequal(run$suites$suite, available),
          .equivalence_passes(run, require_all = TRUE))
directory <- file.path(base, 'final')
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
invisible(Map(function(value, name) {
  utils::write.csv(value, file.path(directory, name), row.names = FALSE)
}, run, c('report.csv', 'suites.csv', 'coverage.csv', 'known-open-status.csv')))
saveRDS(run, file.path(directory, 'results.rds'))
summary <- c(suites = nrow(run$suites), quantities = nrow(run$comparisons),
              values = sum(run$comparisons$n_values),
              disagreements = sum(!run$comparisons$agrees),
              incomplete_outputs = sum(run$coverage$status %in% c('partial', 'not covered')),
              warnings = sum(run$suites$warnings))
print(summary)
writeLines(c('# Final combined full-output validation', '',
  'This is a combined validation matrix, not a second uninterrupted full run.',
  'The original full run executed all 15 suites. Four affected suites were rerun:',
  '`depmixs4` and `hard-depmixs4` after completing mixed residual references,',
  '`invariance` after completing excluded-observation references, and `refusals`',
  'after restoring the classed empty-data error. The other 11 suites retain their',
  'original measured comparisons. Each CSV row identifies its validation_run.', '',
  'The strict acceptance gate passes; all 15 registered suites are present.', '',
  paste(sprintf('- %s: %s', names(summary), format(summary, scientific = FALSE)), collapse = '\n'), '',
  'Original and rerun REPORT.md/SESSION.txt files remain in the sibling folders.',
  'No reference value, tolerance, status or exemption was changed during combination.'),
  file.path(directory, 'REPORT.md'))
