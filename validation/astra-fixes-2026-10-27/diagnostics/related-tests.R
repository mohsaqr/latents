results <- testthat::test_local('.',
  filter = '^(review-1027-diagnostics|diagnostics|tidy-diagnostics|grid-contract|report-surface|plot-named)$',
  reporter = 'summary', stop_on_failure = FALSE)
frame <- as.data.frame(results)
print(data.frame(tests=nrow(frame), passed=sum(frame$passed), failed=sum(frame$failed),
                 errors=sum(frame$error), skipped=sum(frame$skipped)))
saveRDS(results, 'validation/astra-fixes-2026-10-27/diagnostics/related-results.rds')
