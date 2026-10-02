# Retain and validate the native package-check result after R CMD check exits.
evidence <- file.path('validation', 'review-2026-10-01')
checked <- file.path('tmp', 'review-check', 'latents.Rcheck')
log_path <- file.path(checked, '00check.log')
lines <- readLines(log_path)
status <- lines[startsWith(lines, 'Status:')]
stopifnot(length(status) == 1L, identical(status, 'Status: OK'))
stopifnot(!any(grepl('(ERROR|WARNING|NOTE)$', lines)))
stopifnot(file.copy(log_path, file.path(evidence, '00check.log'), overwrite = TRUE))
stopifnot(file.copy(file.path(checked, 'tests', 'testthat.Rout'),
                   file.path(evidence, 'package-tests.log'), overwrite = TRUE))
stopifnot(file.copy(file.path(checked, 'latents-Ex.Rout'),
                   file.path(evidence, 'package-examples.log'), overwrite = TRUE))
write.csv(data.frame(status = 'OK', errors = 0L, warnings = 0L, notes = 0L),
          file.path(evidence, 'package-check-summary.csv'), row.names = FALSE)
cat('Native package check: OK; 0 errors, 0 warnings, 0 notes.\n')
