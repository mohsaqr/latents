# Review evidence

See [the review](../REVISION_REVIEW_2026-09-28.md) for findings and scope.

- `tests.csv`: full source suite, 3,919 assertions, before the final R3STEP
  empty-label fix.
- `controls-followup.csv`: affected controls file after that fix, 56 passing
  assertions, six of them new. Do not add these counts to the full-suite count:
  50 assertions overlap.
- `harness-before.csv`: initial runnable comparison after repairing the stale
  namespace guard. Five comparisons fail: two mclust probabilities, two poLCA
  probabilities and the published inferior-mode search.
- `harness-after.csv` / `harness-suites.csv`: latest successful execution of
  each suite, 1,018 agreements. Seven suites come from the complete rerun; the
  poLCA suite comes from its subsequent rerun with an explicitly refined
  external reference. No agreement thresholds were relaxed.
- `latentgold.csv`: comparison against original retained external output with
  all 12 package targets freshly estimated. The Step-3 target explicitly uses
  observed-information covariance to match the listing, after a follow-up
  regeneration of that target. Seven quantities are non-comparable, not passes.
- `R-CMD-check.log`: full build/check, including vignette rebuilding. Network
  and time-verification notes are environmental limitations, not successful
  online checks.
- `R-CMD-check-followup.log`: final package after the R3STEP fix, with already
  rendered vignettes retained and `--no-vignettes`. Zero errors/warnings, same
  two notes; installed tests pass (3,000 assertions, 60 CRAN skips).
- `houle-arithmetic.txt`: checks of published rounded arithmetic only.
- `SOURCE.md5`: final package R-source file checksums.
- `SESSION.txt`: local R environment used to collect the evidence.

The older committed external reports and Monte Carlo CSVs were preserved.
The large Monte Carlo studies were not all rerun. No proprietary executable
was run; Latent GOLD and Mplus comparisons use retained output.
