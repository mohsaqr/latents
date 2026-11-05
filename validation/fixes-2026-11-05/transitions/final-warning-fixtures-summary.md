# Final transition warning fixture cleanup

Only tests/testthat/test-review-1015-transitions.R changed. The time-ordering test generates sequences that all start in the low profile and finish in the high profile. Its numeric, padded-text and date fits correctly emit latents_boundary for a probability at its numerical bound. Each valid-encoding fit now explicitly requires that warning; the likelihood-equivalence and invalid-ordering assertions remain intact. No production code or broad warning suppression changed.

Whole-file validation uses the final installed latents1.1.17 from /private/tmp/latents-fixes-1105-vtcym3fq/final-check/latents.Rcheck, NOT_CRAN=true, _R_CHECK_LIMIT_CORES_=true and single-thread BLAS settings. No package was reinstalled.

- Before:38 passed,0 failed,0 errors,3 warnings,0 skipped.
- After:41 passed,0 failed,0 errors,0 warnings,0 skipped.

The three added expectations require the known boundary qualification rather than hiding it. final-warning-fixtures.rds/csv/log retain the complete final file results; final-warning-fixtures-before.rds/log retain the warning sources. Fixture parsing passed.
