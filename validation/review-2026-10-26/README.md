# Review of 2026-10-26: what kind of defects were found

This folder holds the review of commit 8ac3576 and the nine-area mathematical
audit that followed. `REPORT.md` is the consolidated report, with every defect,
its fix and the verification runs. `area-reports/` holds one report per area,
with reproductions and computed checks.

## Summary

The core estimators agree with independent implementations. These include
brute-force enumeration, mclust, depmixS4, poLCA, nlme, lme4, glm, polr,
sandwich and numDeriv (see "Verified correct" in `REPORT.md`). The defects were
confined to edge conditions:

- degenerate or constant data;
- extreme measurement scales;
- unusual labels;
- empty or near-empty classes.

None arose on a standard fit to typical, well-scaled data. Several did,
however, give a wrong answer silently instead of raising an error. Those are
fixed, and each fix has a regression test that fails without it.

"Unhappy paths" describes only part of this. Sorted by what went wrong, the
findings fall into three groups.

### 1. Invalid input accepted instead of refused

| Id | Defect |
|---|---|
| R1 | Exactly constant indicators and covariates not detected, because their computed SD is about 1e-17 rather than 0 (three places) |
| CAT-F1 | Direct-effect covariate constant where its indicator is observed; a coefficient of −3.6e13 was reported as converged |
| ADD-A1 | Zero within-group scatter (1.9e-31) not flagged |
| MR-F1 | Empty or thin class reaching the QR fast path in mixture regression |
| LTA-F1 | Rank-deficient `lta()` covariate designs accepted, and their unidentified logits counted in AIC/BIC (**open**) |
| ENUM-1 | `enumerate_lpa()`, `enumerate_lca()` and `enumerate_classes()` returned a grid of failed cells, with no condition, for a column that does not exist or a constant, missing or text indicator, all of which `multilpa()` refuses (found by the refusal-parity suite) |
| MR-F7 | A constant Gaussian outcome in `mixture_regression()` was refused as "every start degenerated", with advice (fewer classes, more starts) that cannot help |

### 2. Valid input refused, or output lost

| Id | Defect |
|---|---|
| R2 | Category labels containing `:`, `,`, `[` or `]` lost their standard errors |
| R3 | Re-reading the original data from CSV refused for person-centred fits |
| SEL-F2 | Homogeneous and occasion `lta()` fits of the same data refused for comparison (row names) |
| CAT-F2, CAT-F3 | Unobserved joint-category labels missing; 17-digit ordinal labels merged on redraw |
| ADD-A2 | Identical data with rows reordered within groups was refused by additive inference and `bootstrap_lrt()`: the sufficient statistics were compared exactly, and summation order moves them by a few ulps. Now compared to 64 ulps; a changed rating is still refused |

### 3. Valid input, silently wrong answer

These are not unhappy paths, and they matter most.

| Id | Defect |
|---|---|
| MR-F4 | Negative-binomial dispersion forced to its floor, lowering the likelihood |
| K1 | EVI/VEI variances biased by 12% for a near-empty profile (common when many classes are enumerated) |
| G1 | Growth models with outcomes in small units broke EM; log likelihood 380 below nlme |
| MR-F3 | Ordinal regressions at large predictor scales stopped short (thresholds 4.4e-4 off polr) while reporting converged |
| MR-F5 | Standard errors depended on predictor units; a calendar-year predictor gave all-`NA` standard errors even with one class. With MR-F6, the differencing step was not scale-aware either (0.2% error under a shift of origin) |
| ML-1 | Zero-weight units changed the log likelihood (by 19/15) and BIC (622.99 → 784.74), so they could change the class count BIC selects; in the additive, regression and growth families too |
| G2 | At extreme outcome scales a growth fit stops short while reporting converged. Measured in standard errors (y x s against the s = 1 fit carried through exactly): negligible for s from 1e-4 to 1e4 (at most 0.009 SE); at s = 1e-5 a random-slope SD is 0.06 SE off and one standard error is 6.2 times too large; at s = 1e5 and 1e6 the random correlation was 0.085 SE off. Fixed by fitting on a standardized outcome and mapping back exactly, with standard-error-scaled differencing steps for the information and the delta method: every scale from 1e-5 to 1e6 now agrees to 8e-7 standard errors, with exact standard errors |
| TS-1 | Rows with no observed indicator were counted in the classification error matrix, the BCH weights and every three-step estimator, though not in entropy or the criteria. Their posterior is the class prior, so P(assigned 1 \| true 2) was 0.025 instead of 5e-5 on a fixture with 3 empty rows of 120 |
| INF-7 | In a binary item, the complement of a probability on the floor (1 - 1e-10) had a standard error of exactly 0 and the interval [1, 1] when it was the last category, but `NA` when it was the first |
| INF-6 | When an item's last category is on the floor, the whole row was held, interior probabilities included, so the other standard errors became conditional and too small: class-share errors 0.0396/0.0307/0.0348 against 0.0521/0.0462/0.0348 on gss82, and the errors depended on the order of the categories |
| GD-1 | `gaussian_diagnostics$positive_rows` counted posteriors that did not underflow (87 of 300 for a component with SD 0.01), not rows with a value; now the exact count |

MR-F5 (calendar-year predictors), ML-1 (zero weights) and INF-6 (an item
category one class never uses) were the most likely to arise in real
analyses; all three are now fixed, as is G2. LTA-F1 and ADD-A2 remain open.

## Round two: every output compared

Every finding added after the first round (TS-1, INF-7, INF-6, GD-1, ENUM-1,
MR-F7) came from a rebuilt equivalence harness. None was visible before,
because each suite compared a few chosen numbers: likelihood, criteria and
point estimates. None compared posteriors, the classification error matrix,
the BCH weights or most standard errors. The harness now compares every
output column a fit reports, as listed by `get_results(fit, "all")`. A column
it cannot compare needs a stated exemption, or the run fails. See
`equivalence/harness/README.md`.

| | Before | After |
|---|---|---|
| Quantities compared | 1,018 | 4,172 |
| Individual values compared | about 1,018 | 2,704,355 |
| Output columns neither compared nor exempted | not measured | 0 |
| Suites | 8 | 10 (adds `invariance` and `refusals`) |

- **`invariance`** refits under transformations with known effects: rows
  permuted; groups, categories and indicators renamed with `: , [ ]`;
  y -> a + c y at c up to 1e6; data duplicated; zero-weight units; a
  predictor shifted by 2000. It confirmed ML-1 and MR-F5, both now fixed,
  and found nothing new.
- **`refusals`** checks 50 pairs of degenerate input and entry point; it found
  ENUM-1.

Three further suites (`hard-mclust`, `hard-polca`, `hard-depmixs4`; 373,673
values) run the same references on data built to be difficult. The package
matched every reference at the reference's own mode, and its search never
finished below a reference's maximum. They found one problem: GD-1 is worse
than first measured (87 of 300 for a component with SD 0.01). They also
measured that `lta()` stops 5.7e-5 from the maximum on a flat likelihood. That
is at most 9.4e-5 of the parameter's own standard error (about 0.01%), and
every reference program stops the
same way, on a change in the log likelihood, so it is recorded as an
observation, not a defect. The search claim compares estimates in standard
errors, the yardstick for "close enough to the maximum"; machine precision
remains the yardstick for "the same computation".

Open defects are listed in `equivalence/harness/known-open.csv`. They are
reported by id instead of failing the run, and an entry that no longer
matches anything fails the run. TS-1, INF-7, INF-6, ENUM-1, MR-F7, ML-1, MR-F5 and G2 are fixed, each
with a regression test that fails without the fix
(`tests/testthat/test-review-1026-empty-rows.R`,
`test-review-1026-zero-weight-likelihood.R`, `test-review-1026-mixreg.R`,
`test-review-1026-refusals.R`, `test-review-1026-categorical.R`).

## Verification

All results are from a fresh install of the final tree.

| Run | Result |
|---|---|
| Full suite | 1,629 tests, 11,327 passes, 0 failures, errors or warnings, 15 fork skips (Accelerate BLAS) |
| New review tests (`tests/testthat/test-review-1026-*.R`) | 111 tests, 729 passes |
| External oracles (`equivalence/test-review-1026-oracles.R`) | 9 tests, 197 passes |

The 15 fork tests were not run under the reference BLAS.
