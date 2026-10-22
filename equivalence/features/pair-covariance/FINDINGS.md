# Findings: covariances for chosen pairs (`local_dependence`, continuous)

## 1. Covariance matrices were not exactly symmetric (fixed)

- Test: `equivalence/test-features-pair-covariance.R`, "the covariance
  pattern holds exactly".
- Repro: `Rscript equivalence/features/pair-covariance/findings/asymmetric-covariance.R`
- The M-step's weighted `crossprod(residuals, residuals * responsibility)`
  can differ between its two triangles in the last bit. In 65 of the 120
  reference-study fits, and in VVV fits too, the reported covariance matrix
  was symmetric only to 1 ulp, so `"y1:y2"` and `"y2:y1"` could report
  covariances 1.1e-16 apart.
- Fix (2026-10-21): `.multilpa_bound_covariance()` (`R/gaussian-moments.R`)
  symmetrizes its input, `(S + t(S)) / 2`, before the eigenvalue check and
  returns that matrix whether or not the bound acts. Re-run: all 120 study
  fits exactly symmetric; the suite passes (1242 expectations).
- Shared code path: every full-covariance structure goes through this
  function. On 48 VVV/EEE/VEV/EVE fits before and after the change, the
  same start won and the same number of iterations was run in all 48; the
  largest change was 5e-10 in the log likelihood and 4.4e-8 in a mean (VVV),
  about 1e-12 elsewhere. The existing equivalence suites
  (`Rscript equivalence/run.R`, 39 files, 5419 expectations) pass unchanged.
