# Random-starts study (2026-10-18)

Do 10 or 20 random starts reach the best maximum, at the true number of
classes and above it, and does the number of starts change the number of
classes BIC selects? Stopped at 29 of 216 planned datasets (machine load):
16 LPA, 13 LCA, one replicate per condition. Rates are rough.

Files: `simulate.R` (design, data, 500-start fits), `check-emulation.R`
(blocks of a long run equal shorter runs), `analyse.R` (emulated 10/20-start
fits, BIC selection), `make-mplus.R` + `mplus/` (Mplus STARTS = 400 80 at K
and K+1, replicate 1), `compare-mplus.R`, `codex/REVIEW.md` (independent
review that recomputed every number).

## Findings

- **True K:** 10 starts reached the 500-start best in 93.8% of runs, 20 starts
  in 95.0%. All misses were low separation with n = 300; in
  `lpa_k4_low_n300_r1` 1 of 500 single starts found the best, and even that
  is 4.76 below Mplus. Earlier claim "10-20 starts are reliable at the true
  K" holds only for moderate or better separation.
- **LCA** agrees with Mplus within 8e-4 at K and K+1 in all 13 datasets.
- **LPA at K+1:** latents' 500-start best ranges from 12.56 below to 12.21
  above Mplus. Every case above Mplus is a one-observation class with all
  variances at `min_variance` (1e-6): a likelihood spike. Verified on
  `lpa_k4_high_n300_r1`, k = 5: class counts 89, 36, 114, 1, 60; the fit
  raises `latents_boundary` and `latents_small_classes`. Mplus at K+1 is not
  clean either (9 of 12 print a variance of 0.000, 11 of 12 a
  non-positive-definite warning).
- **Selection:** 10/20-start runs never chose more classes than the 500-start
  reference. Where they chose fewer (LPA, true K = 4), the reference had
  chosen K+1 because of a singleton spike and the short runs chose the true
  K. Dropping boundary starts makes all 2,349 comparisons agree. The claim
  "optimisation error can only bias BIC toward fewer classes" is false in
  general (BIC shifts by 2 x the likelihood shortfall at each k, which can
  favour either side).
- **Cause** (k-means starts vs Mplus's screened starts) is untested.

## Implemented: `boundary = "exclude"` (opt-in, 0.95.1)

Takes the best start with no variance at `min_variance`. On
`lpa_k4_high_n300_r1`, `enumerate_lpa(n_profiles = 3:5, n_starts = 500)`:
default k = 5 log likelihood -2048.678 (boundary), BIC picks 5; exclude
k = 5 -2067.124 (no boundary), BIC picks 4 (4353.073), as Mplus does. Tests:
`tests/testthat/test-boundary-exclude.R`.
