# Follow-up re-inspection of the 2026-10-01/02 review

An independent pass over the fixes recorded in
[../review-2026-10-01/README.md](../review-2026-10-01/README.md): rerun its
gates, re-derive its riskiest fixes, and check that archived simulations still
reproduce. Package version unchanged (0.9.7).

## Rerun of the earlier review's gates (pre-edit source)

| Gate | Result |
| --- | --- |
| Source identity (`../review-2026-10-01/verify.R --verify-saved`) | 273 files match; saved results 5,295 assertions / 813 cases, 0 failures |
| Lint | No lints |
| Native `R CMD build` + `R CMD check --no-manual` | Status: OK (0 errors, warnings, notes), recorded by `record-package-check.R` |
| Archive identity (`verify-archive.R`) | 272 reviewed files match the built tarball |
| Full source rerun of `verify.R` | Stopped by a 1-hour background limit before finishing (machine shared with four review agents and the check); superseded by the final run below |

## Independent re-derivation of the fixes

| Area | Verdict | Evidence |
| --- | --- | --- |
| Ordinal reference rebasing | Correct | Algebra re-derived; all 6/24 permutations of 3-/4-profile mixed, ordinal-only, NB, weighted and covariate fits preserve probabilities (<= 2.3e-15) and LL (<= 9.1e-13); naive location permutation changes LL by -579.76 |
| NB dispersion derivatives, labels, bootstrap | Correct | Series vs 600-bit Rmpfr finite sums: 2.1e-12 (score), 1.05e-11 (curvature); SEs equal `solve(-numDeriv::hessian)`; labels right for one/many profiles and equal/varying dispersion |
| LTA second-order coordinates | Correct | n_parameters equals hand counts in 11 configurations (HEAD over-counted: 27 vs 19, 59 vs 43, 61 vs 45); zero-effect second-order LL equals first-order exactly |
| LTA covariance charts | Correct | All 14 structures: coordinate count, round trip <= 8.9e-16, decoded covariance <= 2.9e-15, LL <= 9e-13 |
| Sampling weights | Correct for the claimed fixes | Integer weights equal row replication across lpa/lca/multilpa/covariates/FIML/additive/cross-level; scale invariance from 1e-300 to near xmax |

## Defects found and fixed here (all present before the earlier review)

1. NB M-step stalled near the Poisson limit and reported convergence
   (dispersion 4.19e-7 vs profile maximum 1.435e-5; Wald singular). Modified
   Newton in `.latents_negative_binomial_maximize`.
2. `bootstrap_lrt()` on a full-covariance transition null swallowed the classed
   refusal per replicate. Now refused at entry.
3. Transition `parameter_inference()` silently ignored `method = "bootstrap"`,
   `boundary = "fix"` and (general `lta()`) `adjust`. Now refused / applied.
4. `bootstrap_lrt(data = )` failed for single-level fits.
5. Weighted classification `estimated_n`/`estimated_proportion`/OCC were unweighted.
6. Weighted restricted cross-level composition and class counts were unweighted.

Regression tests: `tests/testthat/test-review-followup.R` (8 cases, 37
assertions). Run against an untouched `git archive HEAD` export, every case
fails or errors; with the fixes all pass.

## Archived simulations

Sampled replicates (`simulation-repro-*.R`, run from the repository root):

| Study | Rows compared | Max abs difference in estimates |
| --- | --- | --- |
| Additive (`additive-simulation-full`) | 500, all 11 cells | 5.1e-15 (no change in inference availability) |
| Families | 184, all 4 cells | 4.9e-15 |
| LTA extensions | 182, all 4 cells | 4.3e-8 — identical on HEAD, so it predates this and the earlier review; ~3e-7 of an SE |

The archived simulations do not exercise sampling weights, ordinal/NB
indicators or occasion-varying second-order transitions, so their conclusions
neither cover nor are affected by those fixes.

## Final gates (edited source)

- `verify.R` here: **5,332 assertions in 821 cases, 0 failures, errors,
  warnings or skips**; source identity verified over 274 files
  (`tests-final.log`, `SOURCE.csv`, `SESSION.txt`).
- Lint: no lints.
- Native `R CMD build` + `R CMD check --no-manual` of the edited source:
  **Status: OK, 0 errors, 0 warnings, 0 notes** (`00check.log`); installed
  tests `[ FAIL 0 | WARN 0 | SKIP 4 | PASS 5305 ]` (`package-tests.log`; the
  four skips need source files and pass in the source run); examples in
  `package-examples.log`. A first check flagged an undocumented `...` in
  `parameter_inference.multilpa_transitions.Rd` (introduced by this pass);
  it was restored and both gates were rerun on the final files.

## Open

- NB closed-form dispersion derivatives lose ~1e-6 relative accuracy just
  above the series switch (alpha * max(y, mu) = 1e-3) for counts ~500,
  growing with the count; a stable large-r form is needed.
- Not reviewed: weighted `additive_dispersion`, robust/OPG LTA covariance,
  second-order fits with occasion gaps.
- Cosmetic: single-level fits print "1 group classes"; the NB equal-vs-varying
  `bootstrap_lrt()` refusal does not name the dispersion configuration.

## Release

Released as 0.9.8. The version bump (`DESCRIPTION` `Version`, `NEWS.md`
heading) was made after the gates above; those are the only two files that
differ from `SOURCE.csv` (272 of 274 unchanged), and the release-item tests
pass on the bumped files (29/29). The `tests-final.rds` result object is kept
locally and not committed; `tests-final.csv/log` carry the same results.

## Post-0.9.8: transition full-covariance simulation and bootstrap

Implemented after release 0.9.8 (version "development" until the next bump):

- `.lta_simulate()` draws full-covariance profiles through the Cholesky factor
  of each profile's covariance (`.lta_draw_continuous()`); diagonal fits keep
  their exact random stream. `bootstrap_lrt()` therefore compares
  full-covariance transition models.
- `parameter_inference(method = "bootstrap", data = )` for
  `multilpa_transitions` and `multilpa_lta` (`R/lta-bootstrap.R`): persons are
  resampled, refitted with the same specification (weights included), and
  matched to the original on measurement (profiles) and on matched
  initial/transition parameters (group classes; stayer only to stayer).
  Initial logits are rebased when the reference profile moves; stay-referenced
  transition logits are renamed only.
- Regression fixed: 0.9.8 refused `boundary = "fix"` on transition fits, which
  `get_results(, "responses")`/`summary()` request, so categorical transition
  response tables lost their standard errors (0.020/0.023 in 0.9.7, `NA` in
  0.9.8).

Evidence:

- `tests/testthat/test-lta-bootstrap.R` (100 assertions): diagonal stream
  identical; full-covariance draws recover each covariance (2e5 draws, error
  < 0.05); relabelling preserves every implied initial, transition and
  second-order probability (K = 3, two classes, covariate); permuted
  homogeneous fit aligns back to 1e-12; person-shuffled refits that switched
  labels align to 1e-3; bootstrap/Wald SE ratio within 0.6-1.7 and identical
  layout; reproducibility and seed restoration; weighted fits; classed errors.
  Sabotage (no rebasing; diagonal factor) fails the targeted tests.
- `transition-bootstrap-validation.R/.log/.txt` and CSVs: 200 datasets from a
  VVV occasion-transition truth (240 persons x 4 occasions).
  - Percentile coverage (iter = 99): pooled 0.924; per parameter 0.890-0.940
    (MCSE about 0.02). Label-switched datasets (94) covered 0.933 and the rest
    0.915, so matching is not the cause. ML logit bias up to 0.23 SD.
  - `bootstrap_lrt()` size, first vs second order with a VVV null: 0.065 at
    0.05 (MCSE 0.017), 0.115 at 0.10; KS uniformity p = 0.281.
  - Sensitivity (`TB_ITER=199`, same 200 datasets):
    `transition-bootstrap-validation-iter199.log`. Pooled coverage 0.928
    (per parameter 0.890-0.945); label-switched 0.930 vs 0.926. The number of
    replicates is not the cause; the mild undercoverage is a finite-sample
    property of percentile intervals for ML logits at 240 persons, not of the
    label matching. Read intervals for weakly informed moves with that in mind.
- Final gates on the post-0.9.8 source: source suite **5,436 assertions in
  832 cases, 0 failures/errors/warnings/skips**, identity over 276 files
  (`tests-final.log`); lint clean; native `R CMD check --no-manual`
  **Status: OK (0/0/0)**, installed tests `[ FAIL 0 | WARN 0 | SKIP 4 | PASS
  5409 ]` (`00check.log`, `package-tests.log`, `package-examples.log`).

Released as 0.9.9: the version bump (`DESCRIPTION`, `NEWS.md` heading) was
made after the gates above and is the only difference from `SOURCE.csv`
(274 of 276 files unchanged); release-item tests pass on the bumped files (29/29).

## Post-0.9.10: false non-convergence in transition fits

CI for 0.9.8-0.9.10 failed on macOS and R 4.1 at one negative-binomial LTA
test. Cause: `.lta_quasi_newton()` counted only L-BFGS-B code 0 as converged;
code 52 (ABNORMAL_TERMINATION_IN_LNSRCH) occurs at the maximum on some
platforms. Locally 5 of 120 test-design datasets were flagged unconverged with
scores of about 1e-6 and a restart gain below 1e-12; 0 after the fix (restart
once, accept when the gain is within `tol`). The 0.9.10 `max_iter` test change
was a misdiagnosis and is reverted; a regression test with the five datasets
fails 5 times on 0.9.10 and passes now. Gates: source suite 5,446 assertions
in 833 cases, 0 failures/errors/warnings/skips, identity over 276 files; lint
clean; native `R CMD check --no-manual` Status: OK (0/0/0), installed tests
`[ FAIL 0 | WARN 0 | SKIP 4 | PASS 5419 ]`.

Released as 0.9.11: the version bump (`DESCRIPTION`, `NEWS.md` heading) was
made after the gates above and is the only difference from `SOURCE.csv`
(274 of 276 files unchanged); release-item tests pass (29/29).
