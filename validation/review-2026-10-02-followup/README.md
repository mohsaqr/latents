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
