# Review 2026-10-26: growth mixture models (slug `growth`)

## Scope

Files: `R/growth-mixture.R`, `R/growth-mixture-engine.R`,
`R/growth-mixture-correlated.R`, `R/growth-mixture-inference.R`,
`R/growth-mixture-methods.R`, `R/growth-mixture-parameter-inference.R`,
`R/growth-mixture-predict.R`, `R/growth-mixture-simulate.R`,
`R/growth-mixture-plot.R`, `R/trajectory-views.R`. I reviewed Codex's diff
(8ac3576) for these files, then the algorithms: the marginal likelihood
through Woodbury and the determinant lemma, the empirical Bayes posteriors,
ECM M-steps (beta, sigma^2, G varying/equal/proportional/diagonal), the AR
and Brownian residual processes, analytic scores, the Wald covariance,
BIC/ICL, simulate, predict, the bootstrap class matching and relabelling. I
checked each by computation. All runs used `OMP_NUM_THREADS=1
VECLIB_MAXIMUM_THREADS=1` and `devtools::load_all()`. Scratch scripts are in
`/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/caf01cef-53f3-419a-95fe-ebd08f111cb5/scratchpad/growth/`
(`v1.R` to `v15.R`).

## Verdict on each Codex change

| Change | Verdict | Evidence |
|---|---|---|
| `.growth_patterns`: Gram-matrix keys `signif(., 12)` replaced by `sprintf("%a")` hex keys | Correct. It can fragment patterns but never merges them, and fragmentation costs speed only | Hex strings are exact, so two different Gram matrices never share a key. Same-time persons whose rows are shuffled and whose times are fractional (`time/3 + 0.1`) give 27 patterns for 18 distinct time sets (1-ulp summation-order differences). With integer times: 18 = 18. The log likelihood depends only on each person's own Z'Z, so it is unaffected (shuffled vs sorted LL difference 3.0e-12). |
| `.growth_dense_patterns`: `format(digits = 17)` replaced by `%a` | Correct | Keys are built from raw times and design rows, so both forms are exact. The only difference is that `-0` and `0` now split into separate patterns (speed only). The correlated log density equals a dense `dmvnorm` to 3.6e-15. |
| `.growth_random_design()`: newdata evaluated with the stored `terms` (predvars), `xlevels` and contrasts | Correct | `random = ~ poly(time, 2)`: `predict(newdata = data)` equals `predict()` to 2.7e-15. The design of a 3-person subset with time <= 2 equals the fitted rows to 2.3e-16 (before the change, `poly()` was refitted on the subset). `spec$random_terms` and `spec$random_xlevels` are set in `R/mixture-regression.R:696-699`. Older objects fall back through `%||%`. |
| Brownian newdata times must be >= 0 | Correct | min(t, s) is not PSD for negative times. `predict()` on a row with time -0.5 raises `latents_bad_data` (computed). This matches the fit-time check in `.growth_residual_time()`. |
| Plot: person facet labels keyed on native ids, with `make.unique` | Correct (read only, not rendered) | `factor(match(person, people), levels = ordering, labels = labels[ordering])` maps each person to their own label. Alphabetical ordering is kept as before. No visual output was rendered (project rule: plot issues are reported, not changed). |
| Plot: spread band uses `.growth_random_design()` | Correct | Same function as predict. The `random = ~ 1` path (`data.frame(row = ...)`) is exercised by the existing "every growth plot draws" test (results below). |
| `.trajectory_check_data`: every `model_data` column must match, and the sampling weights must match | Correct, one style note | The fitted data passes. A changed `y` raises `latents_bad_inference_data`. Changed weights raise `latents_bad_inference_data`. The right weights pass (all computed). Style note: `tryCatch(..., error = function(condition) FALSE)` is a blanket handler (CLAUDE.md Rule 5), but it is turned into a classed error at once, so nothing fails silently. |

I searched the partition for the lead's `sd <= 0` constancy pattern. Codex
made no such change in these files. The only spread test is in
`.growth_trajectory_seed` (`R/growth-mixture.R:219-220`, `spread > 0`).
There, a column is constant only when its values are bitwise copies of the
pooled coefficient, for example a time-invariant covariate. R's two-pass
`mean()` makes `sd()` of identical copies exactly 0: 2000 random values
repeated 90 times gave 0 nonzero sds, and `sd(rep(0.1, 3))` is 0. The
`time + age` example gave sd exactly 0 for `age`. Not affected.

## Findings

| id | severity | file:line | problem | reproduction (computed) | status |
|---|---|---|---|---|---|
| G1 | medium | `R/growth-mixture-engine.R:216-219` (`.growth_factor`) | The Cholesky jitter was `1e-12 * max(1, max(diag(G)))`: an absolute floor of 1e-12. With an outcome in small units (variances near 1e-10 or smaller), it swamps G. The likelihood is biased and EM's monotonicity breaks. | `v3b.R`: the one-class data of `v1.R`, `y * s`. At s = 1e-5: **error** "EM likelihood decreased by 0.00983 ... This is a defect". At s = 1e-7: LL is 379.9 below `nlme::lme` ML and G11 is inflated 90.3x (expected 1x). After the fix (jitter relative to the covariance's scale, falling back to 1e-12 for a zero matrix): s = 1e-5 fits, LL within 3.95e-3 of nlme. s = 1e-7: within 6.3e-3, G11 ratio 0.988. The rest is G2. Fits with max diag(G) >= 1 are bit-for-bit unchanged. With max diag(G) < 1 the change is below 1e-12 relative. | **fixed** (regression tests fail 3/28 before, pass 28/28 after) |
| G2 | medium | `R/growth-mixture-inference.R:278-302` (`.growth_quasi_newton`), with the relative EM stop in `R/kernel-em.R` (outside the partition) | The fit is not scale-free in the outcome's units. EM stops on a log-likelihood change relative to \|LL\|, which shifts by -N log s with the units. The L-BFGS-B coordinates (coefficients, off-diagonal Cholesky) are in raw units, so the problem is badly conditioned for very large or very small y. The fit stops short of the maximum but reports `converged = TRUE`. | `v12.R`, same data, `y * s`, fitted LL minus lme4 ML LL (the latents density at lme4's MLE matches lme4 to 1e-12): s = 1e-4: -6.7e-05; s = 1e-5: -3.95e-03; s = 1, 1e3, 1e4: +1.96e-09; s = 1e5 and 1e6: -6.57e-03. All `converged = TRUE`. At s = 1e-5, G11 is off by 0.9%, and the Wald information is singular (step rule `plus1`, h = 1e-5 absolute, close to the coefficient itself), so the SEs are NA. Wald SEs are scale-equivariant for s in [1e-3, 1e3] (ratio 1 +/- 1.1e-4). | reported, not fixed (fixing it means internal standardization of y, or relative LL/step rules in the kernel: a design change that touches ordinary fits and files outside the partition). Workaround: rescale y to O(1). |
| G3 | low | `R/growth-mixture-methods.R:282-300` (fit table) with weight normalization in the kernel | Persons with sampling weight 0 still count in `n_persons`, `nobs()` and log(n) of BIC. Weights are normalized to sum to all units, zeros included. | `v10.R`: 5 of 90 persons at weight 0. Coefficients equal the fit without them (max diff 8.2e-8). LL = (90/85) x LL of the dropped fit (-675.53 vs -638.00). `nobs` = 90 vs 85, BIC 1409.56 vs 1333.76. Consistent with the documented "scaled to sum to the number of units", but debatable. | reported (convention; the normalization lives in the kernel) |
| G4 | low | `R/growth-mixture-engine.R:70` | Hex keys fragment patterns that are equal up to an ulp. | 27 vs 18 patterns (see the verdict table). Speed only. | reported (no action needed) |
| G5 | info | `tests/testthat/test-review-1026-growth.R` | The new calibration tests use `nlme` under `skip_if_not_installed("nlme")`. `nlme` is not in DESCRIPTION Suggests (nor is `mvtnorm`, which other reviewers' tests use). `R CMD check` will flag undeclared test dependencies. | | reported: add `nlme` to Suggests (DESCRIPTION is outside the partition) |

## Algorithms verified correct

| Algorithm | Check | Max abs difference |
|---|---|---|
| One-class marginal ML (LMM) | vs `nlme::lme(method = "ML")`, unbalanced (15% of rows dropped), random intercept + slope, unstructured G | LL 2.4e-8; beta 1e-6; G 5e-5; sigma^2 3e-6; vs `lme4::lmer(REML = FALSE)` LL 2.0e-9 |
| Woodbury/determinant-lemma log density, 2 classes, common covariate, unbalanced | vs `mvtnorm::dmvnorm` on V = Z G Z' + sigma^2 I | 2.4e-10 (log density); total LL 4.5e-13 |
| Empirical Bayes means G Z'V^-1 r and covariances (G^-1 + Z'Z/sigma^2)^-1 | dense solve | means 4.0e-11; covariances 7.3e-12 |
| One-class BLUPs | vs `nlme::ranef` | 3.6e-5 (nlme tolerance) |
| Penalized-residual fallback for the unstable Woodbury difference | read: r'V^-1 r = min_a \|\|r - Z L a\|\|^2/sigma^2 + \|\|a\|\|^2 (identity) | n/a |
| ECM M-steps (beta given E[b], sigma^2 with tr(Z'Z S), G = E[bb'] / n_k; proportional closed forms) | read (Laird and Ware 1982) and confirmed by matching optima: equal-G 2-class vs `lcmm::hlme(nwg = FALSE)` | LL 5.2e-6 |
| Proportional G (`nwg = TRUE`) | latents LL -590.1379 vs lcmm's own fit -600.5419 (lcmm stuck in a local maximum). lcmm started at the latents estimate converges to the same point | 5.2e-6 |
| AR process W = s^2 exp(-rho\|dt\|), Brownian s^2 min(t, s) | density vs dense `dmvnorm`; full fit vs `lcmm::hlme(cor = AR(time))` / `BM(time)` | density 3.6e-15; LL 4.5e-6 (AR), 4.5e-6 (BM) |
| Analytic scores (varying, equal, proportional, diagonal, equal residual variance) | vs `numDeriv::grad` / central differences, at a perturbed point | relative 2.5e-8, 1.8e-8, 3.4e-9, 1.2e-8, 7.7e-9 |
| Scores in the high variance-ratio branch (kernel > 1e6, 2.3e6 here) | central differences | ~1e-3 absolute on gradients of order 1e-2 to 1e3 (within finite-difference noise) |
| Observed information / Wald SEs | vs inverse of `numDeriv::hessian` | max relative SE difference 6.1e-9 |
| BIC with persons | `get_results(fit, "fit")$bic` = `BIC(fit)` = -2LL + log(n_persons) k | exact (1114.491) |
| Class relabelling | permuted parameters | LL difference 0; posterior swap 0 |
| Time origin shift (intercept + slope, unstructured G) | parameters mapped by T = [[1, 2.5], [0, 1]] | 3.7e-11; refits 4.1e-5 apart (class 2 on a boundary in that draw) |
| One observation per person (20 singletons added) | dense density | 1.4e-9; converged; SEs finite |
| simulate() (AR, 1 class) | 20000 draws: sample covariance vs V; mean vs X beta | 1.2% of max V (sampling SE ~1%); mean 0.6% of sd |
| Wald SE scale equivariance | y x 1e-3, 1e-2, 1e2, 1e3 | SE ratios 1 +/- 1.1e-4 |

## Tests added

`tests/testthat/test-review-1026-growth.R` (local seeds via `withr`; about
8 s with NOT_CRAN=true):

- "the per-person log density and posteriors equal a dense computation"
- "one class reproduces the linear mixed model of nlme::lme (ML)" (`skip_if_not_installed("nlme")`; also BIC with persons)
- "analytic scores equal central differences of the log likelihood"
- "AR and Brownian residual covariances follow their formulas and are positive definite"
- "the correlated-residual log density equals a dense computation"
- "relabelling the classes leaves the likelihood unchanged and swaps posteriors"
- "shifting the time origin leaves an intercept-and-slope model's likelihood unchanged"
- "a random-effect covariance in small units is factorized without an absolute jitter" (regression for G1)
- "an outcome in small units fits like the same outcome in ordinary units" (regression for G1)

Result with the fix: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 28 ]`. With the fix
reverted by hand: `[ FAIL 3 | WARN 0 | SKIP 0 | PASS 25 ]`. The failures
were the factor and density assertions at lines 212 and 228, plus an error
at line 240 ("EM likelihood decreased by 0.0104"). Under CRAN settings (no
NOT_CRAN): `[ FAIL 0 | WARN 0 | SKIP 7 | PASS 9 ]`.

## Related existing tests

These were run with the fix in place (NOT_CRAN=true, `testthat::test_file`,
silent reporter; script `runall.R`):

| file | failed | errors | skipped | warnings | expectations passed |
|---|---|---|---|---|---|
| test-growth-mixture.R | 0 | 0 | 0 | 0 | 122 |
| test-growth-correlated.R | 0 | 0 | 0 | 0 | 25 |
| test-growth-predict.R | 0 | 0 | 0 | 0 | 21 |
| test-growth-parameter-inference.R | 0 | 0 | 0 | 0 | 39 |
| test-growth-r3step.R | 0 | 0 | 0 | 0 | 16 |
| test-multilevel-growth.R | 0 | 0 | 0 | 0 | 33 |
| test-audit-extension-repairs.R | 0 | 0 | 0 | 0 | 67 |
| test-engine-fixes.R | 0 | 0 | 0 | 0 | 26 |
| test-review-1015-regression.R | 0 | 0 | 0 | 1 | 30 |
| test-review-1020-zero-weights.R | 0 | 0 | 0 | 0 | 18 |
| test-api-compat-0921.R | 0 | 0 | 0 | 0 | 44 |
| test-review-1015-accessors.R | 0 | 0 | 0 | 0 | 34 |

The one warning, at test-review-1015-regression.R:122, is a
mixture-regression separation warning ("A fitted logit or log mean is
extreme") from a model without random effects. That model does not call
`.growth_factor`, so the warning is unrelated to this change.

## Could not check

- Plots were not rendered: no HTML was made under `./tmp/`. The label and band changes were checked by reading the code only.
- Multilevel growth (cluster) scores and the two-level M-step were not checked independently. The existing `test-multilevel-growth.R` was run (above).
- The parametric bootstrap LRT driver is in `R/mixture-regression-enumerate.R` (another reviewer's partition). Only its growth draw (`.growth_draw`) was checked here, through simulate moments.
- The person bootstrap of `parameter_inference(method = "bootstrap")` was read, not run. Class matching is by posterior membership, and `.growth_permute_params` leaves the proportional scales un-rebased, which is harmless because the natural values read `.growth_covariances()`.
