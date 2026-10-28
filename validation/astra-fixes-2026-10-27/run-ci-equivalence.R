suppressMessages(pkgload::load_all('.', quiet = TRUE))
print(sessionInfo())
r <- testthat::test_dir('tests/equivalence', load_package = 'none',
                         reporter = 'summary', stop_on_failure = FALSE)
d <- as.data.frame(r)
saveRDS(r, 'validation/astra-fixes-2026-10-27/ci-equivalence-results.rds')
write.csv(d[setdiff(names(d), 'result')],
          'validation/astra-fixes-2026-10-27/ci-equivalence-results.csv', row.names = FALSE)
print(c(files = length(unique(d$file)), tests = nrow(d), passed = sum(d$passed),
         failed = sum(d$failed), errors = sum(d$error), warnings = sum(d$warning),
         skipped = sum(d$skipped)))
quit(status = as.integer(any(d$failed > 0 | d$error)), save = 'no')
