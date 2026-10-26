# Review: finite mixtures of regressions (slug `mixreg`), 2026-10-26

## Scope

Files reviewed: `R/mixture-regression.R`, `-engine.R`, `-ordinal.R`, `-predict.R`,
`-starts.R`, `-inference.R`, `-methods.R`, `-enumerate.R` (and, read-only,
`R/kernel-structures.R`, `R/inference-service.R`). I reviewed Codex's commit
8ac3576 for these files and the algorithms beneath them: E-step, Newton/IRLS
M-steps, NB2 dispersion, ordinal thresholds, parameter counts, observed /
robust covariance, prediction, plot helper and identifiability refusals.
Every number below comes from a run in this session (R 4.5.2, `OMP_NUM_THREADS=1`,
`load_all`). Scripts are in
`/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/caf01cef-53f3-419a-95fe-ebd08f111cb5/scratchpad/mixreg/`.
For before/after comparisons, the pre-Codex tree (HEAD~1) was exported with
`git archive` to `.../scratchpad/prev`.

## Verdict on each Codex change (my partition)

| Change | Verdict | Evidence |
|---|---|---|
| `.mixture_linear_predictors`: centre columns with \|mean\| > 1e6 before forming eta | **Correct** | This is an exact reparametrisation, since the intercept absorbs c'beta. Shifting x by 1e7 leaves fitted eta equal to 2.7e-9 (gaussian) and 1.0e-9 (poisson). |
| M-step: centre large-origin columns, convert the intercept back with b0_raw = b0_c - c'beta - c_z'gamma | **Correct algebra, with a constancy regression** | The forward and back maps are inverses (eta = b0_c + (x-c)'b). With a 1e7 shift, slopes agree to 1.3e-10 and LL to 2.7e-8. Before the fix, HEAD~1 gave dLL = -276 at a 1e4 shift. **The `spread == 0` constancy test was a regression (F2).** |
| Gaussian QR shortcut when q == 0 and a column is centred | **Regression (F1)** | `sqrt(tau/max(tau))` is NaN for an empty class, and `qr.solve` stops on a class with rank < p. Both raised unclassed errors where the ridge path stayed solvable. Sampling weights *are* included, because `tau` is already weighted (`.mixture_weigh`), and the weighted least-squares solution does not depend on the scale of the weights. |
| Newton solve in information-standardised coordinates when the ridge exceeds 1e-7 x some diagonal entry | **Correct (Jacobi preconditioning, same solution)**, but **incomplete**: the ordinal M-step kept the raw ridge (F3) | An empty class (zero diagonal) falls back to the 1e-10 ridge with direction 0, which is checked. |
| Separation flag: \|eta\| > 30 on supported rows (GLM), \|t_c - eta\| > 30 (ordinal) | **Sensible** | This is unit-invariant, unlike \|beta\| > 30. For logit it means a fitted probability within ~1e-13 of 0 or 1, close to `glm.fit`'s "numerically 0 or 1" rule (\|eta\| > ~33.7). For log links, eta < -30 means mu < 1e-13, which matches glm's "fitted rates numerically 0" (eta < ~-36). eta > 30 means mu > 1e13, which is never plausible. The ordinal `outer(thresholds, eta, "-")` is (C-1) x n per class and is right: P(Y<=c) = F(t_c - eta). Caveat (low): the check is not weighted by tau, so a class can be flagged for rows it does not own. |
| Bernoulli and two-category ordinal refusal rewording | **Correct and more accurate** | The old text ("a mixture of Bernoullis is a Bernoulli, so it is not identified") holds only at one fixed covariate value. With a continuous covariate, sum_k pi_k F(a_k + b_k x) is a finite mixture of logistic CDFs in x, which is identifiable (Follmann & Lambert 1991; Grün & Leisch 2008, J. Classification 25:225-247, give sufficient conditions on trials and covariate design). The refusal is still a valid support restriction. It keeps class `latents_not_identified`, so the released API is unchanged. |
| `.mixture_prediction_regression` split; plot uses it (no membership groups) | **Correct** | Prediction tests pass: hand-computed priors x means, and hand-computed posteriors, both within 1e-10. |
| Typical-row helper keeps `ordered` and the `contrasts` attribute | **Correct** | — |
| `random_terms` / `random_xlevels` in the spec | Correct (no math) | — |
| Seed documentation, warning text | No math | — |

## Findings

| id | sev | file:line | problem | reproducing input; actual vs expected | status |
|---|---|---|---|---|---|
| F1 | high | R/mixture-regression-engine.R:237 (QR shortcut) | The large-origin Gaussian least-squares shortcut errors on an empty class (`max(tau)=0` makes the weights NaN) or on a class with rank < p (`qr.solve` stops). The earlier ridge path kept both solvable. | (a) `.mixture_update_regression(spec, null_params, cbind(1, 0))`, spec from `y ~ I(x+1e7)`, n=300: error "NA/NaN/Inf in foreign function call"; with an unshifted x it returns finite betas (class 2 = 0). (b) `set.seed(17); x <- round(rnorm(20),1)+1e7; y <- round(1+0.5*(x-1e7)+rnorm(20)); mixture_regression(y ~ x, d, 4, n_starts=10, seed=17)`: unclassed error "singular matrix 'a' in solve". Across seeds 1:40, 2 errored; HEAD~1 and the unshifted data gave none (22 no_valid_start / 18 ok). | **fixed**: per-class `qr()` with a rank check; an empty or rank-deficient class falls back to the ridge-stabilised Newton step. Well-posed fits are unchanged because `qr.coef` on the same decomposition is what `qr.solve` does. |
| F2 | high | R/mixture-regression-engine.R:216 (centre) | The same constancy pattern the lead found: `spread == 0` misses an exactly constant column whose colMeans rounds. `colMeans(rep(0.1, 99))` != 0.1 and the spread is 1.9e-16. The column is then centred into rounding noise next to the intercept. | `y ~ x + k`, `k = 0.1`, n=99, K=1. Gaussian: error "singular matrix 'a' in solve". Poisson: "EM likelihood decreased by 0.0661 ... This is a defect". HEAD~1 fitted both (slope 0.97494 / 0.234936, as without `k`). | **fixed**: constancy is now read on the raw values (`max == min`). Both fits reproduce HEAD~1, and the LL equals the fit without `k`. |
| F3 | medium | R/mixture-regression-ordinal.R:350 | The ordinal Newton step kept the raw ridge 1e-10(1+max diag). With large predictor units it swamps the threshold curvature, so the inner Newton stalls and EM reports convergence at a non-maximum. Codex's GLM fix had not been carried over to this path. | `o ~ I(x*1e6)`, K=1, n=300. Against polr: dLL = -7.8e-6, thresholds off by 4.4e-4, `converged = TRUE`. One M-step from null ended 0.039 below the max. K=2: LL -320.0676572 vs -320.0675199 at unit scale. | **fixed**: the same standardised solve as the GLM path, gated identically. Now, across scale 1, 40, 1e3, 1e6, LL = -330.704007777818 (K=1) and -320.067519935509 (K=2) at every scale. Scale 1 is bit-identical to before. |
| F4 | medium | R/mixture-regression-engine.R:431 (NB dispersion) | The KKT override "score at alpha=0 <= 0, so set alpha to the floor" is only a *local* condition. The weighted alpha profile can be bimodal, so the step can move from a higher interior point to the floor and *lower* the M-step objective, which breaks GEM monotonicity. | y=(0,0,16,0,0,6,0); mu=(0.261227,0.908894,15.731,1.18125,2.08141,3.39144,0.352174); tau=(0.0263183,0.570606,0.205777,0.938473,0.189234,0.0251619,0.88033); start alpha=20. Returned alpha=1e-8 with objective -2.880383, against -1.712351 at the start and an interior max of -1.706305 at alpha=26.81. | **fixed**: a class moves to the floor only when that does not lower its own objective. Now returns alpha 26.81364 with objective -1.706305. glm.nb calibration is unchanged (dLL 4e-13). |
| F5 | high (pre-existing) | R/mixture-regression-inference.R:215 (`.mixture_invert`) | The singularity rule (min eig <= 1e-10 max \|eig\| of the raw information) depends on units and origin. An ordinary, well-identified fit with a calendar-year predictor gets **all SEs NA** plus `latents_singular_information`, even with one class, where glm() is fine. | `year` in 2000:2020, n=300, K=1 or 2: every SE is NA (eigenvalue ratio 2.1e-12, min eigenvalue 5.5e-3 > 0). With `t = year - 2010` the SEs are finite (slope 0.0063). Shifts of 1e3 and above, and scaling a predictor by 1e6, also give NA. | **reported, not fixed**: `.mixture_invert` is shared with growth-mixture and lta-three-step, and the rule is a judgement call (it also guards near-empty classes). Suggested fix: test and invert the Jacobi-scaled information D^-1/2 I D^-1/2, which is exactly scale-invariant. For origin invariance, compute the information in the centred chart the M-step already builds and map it back with V = J V_c J'. |
| F6 | medium (pre-existing) | R/mixture-regression-inference.R:148 (`h_rule = "plus1"`) | The finite-difference step h = 1e-5(1+\|theta\|) ignores each parameter's scale. When a predictor is in large units its slope (and SE) is tiny, so h is many SEs wide and the observed SEs carry truncation error. | x*1e4 (estimates are invariant to 1e-14): the slope SE ratio against the unit-scale fit is off by 3.6e-4 (gaussian), 7.6e-3 (poisson) and 5.4e-3 (ordinal). At unit scale it is off by < 2e-10. | reported, not fixed (shared inference service; changes SE numbers). |
| F7 | low | R/mixture-regression-engine.R:309-325 | The separation flag is not weighted by tau, so a class can be flagged at rows it does not own. It is only a warning. | Reasoning, not reproduced as a false positive | reported |
| F8 | low | EM default `tol = 1e-8` (relative LL) | At the default tolerance the score at the selected optimum is ~1e-2. The estimates match flexmix to ~1e-4, against SEs of ~0.05. With `tol = 1e-13` they match to 1e-7. | Gaussian K=2, n=300: max\|grad\| = 6.0e-3. Poisson vs flexmix at default tol: dcoef 8.7e-5; at tol 1e-13: 1.4e-7, dLL 1.9e-11. | reported (note for users who need tight estimates) |
| F9 | low | ordinal, large origin | Ordinal x has no intercept to absorb a large origin. Thresholds near 1e7 cancel against eta, so the slope is accurate to about 1e-5 at a 1e7 shift (1e-11 at 1e4). | `o ~ I(x+1e7)`, K=2: dslope 1.9e-5, dLL 4e-8 | reported |

## Algorithms verified correct (checks run)

- **One-class limit equals the GLM.**
  - Gaussian vs glm: coefficients 1.8e-16, LL 1.1e-13; SE vs glm x sqrt((n-p)/n) to 4.4e-13.
  - Binomial counts: coefficients 3.2e-12, SE 8.6e-8.
  - Poisson with offset: coefficients 1.4e-11, SE 1.8e-7.
  - NB vs MASS::glm.nb: coefficients 1.3e-8, alpha 8.8e-10, LL 4e-13. The beta SE differs by 5.7e-3 because glm.nb uses the expected information conditional on theta; against numDeriv the full observed SEs agree to 8.2e-10.
  - Ordinal vs MASS::polr: slope 4.2e-10, thresholds 7.1e-10, LL 1.1e-13; slope and first-threshold SEs equal polr's to 7 digits.
- **Observed SEs vs the numDeriv Hessian of the LL** (max relative difference): Gaussian K=2 with membership 1.8e-7; equal variance 2.6e-10; Poisson K=2 8.6e-10; NB K=1 7.4e-9; NB K=2 3.1e-9; ordinal K=2 1.2e-7; Gaussian with `common` 7.9e-7; grouped binomial (`class_level = "group"`) 3.3e-9. In every case `n_parameters` equals the length of the packed theta.
- **Robust covariance**: one-class Poisson equals the GLM sandwich bread x meat x bread x n/(n-1) (test, tolerance 1e-5).
- **flexmix with K=2**, started from our posterior, tol 1e-12, our tol 1e-13: Poisson dLL 1.9e-11, coefficients 1.4e-7; grouped binomial (`| id`) dLL 7e-12, coefficients 6e-7. For the Gaussian, our LL is higher by 3.06e-3 because flexmix's sigma is df-corrected, sqrt(sum w r^2/mean(w)/(n-p)), not the MLE (read from `flexmix:::FLXMRglm`).
- **NB2 score at alpha=0** equals ((y-mu)^2 - y)/2; derived by hand, it matches the code.
- **Ordinal derivatives** are confirmed indirectly by the numDeriv agreement of the SEs and by polr agreement.
- **Invariants**:
  - Class relabelling (`.mixture_permute`) leaves LL equal to 1e-12, and tau and means permute.
  - A 1e7 translation leaves slopes, eta and LL unchanged (gaussian and poisson).
  - Scaling a predictor by 1e6 multiplies its slopes by 1e-6 and leaves LL unchanged (gaussian, poisson, and ordinal after F3).
- **Prediction**: `type = "response"` = sum_k prior_k(w) mean_k(x), with priors from the membership softmax (prior, not posterior, as it should be for new data). `"class_response"` priors match. `"posterior"` on new rows equals the hand-computed Bayes posterior to 1e-10, and on the training rows reproduces `fit$expectation$tau`.

## Tests added

`tests/testthat/test-review-1026-mixreg.R` (~6 s):

- one-class gaussian, binomial and poisson fits are glm()
- one-class negative-binomial and ordinal fits are glm.nb() and polr()
- observed-information errors equal a numerical Hessian of the likelihood
- the robust covariance of a one-class Poisson fit is the GLM sandwich
- relabelling the classes leaves the likelihood and the fitted means unchanged
- translating a predictor by 1e7 changes only the intercepts
- rescaling a predictor rescales its slopes and nothing else (regression test for F3)
- the large-origin least-squares step keeps an empty or tiny class solvable (F1)
- an overfitted large-origin gaussian mixture fails only with a classed condition (F1)
- an aliased constant column with an inexact mean adds nothing (F2)
- predictions use prior class probabilities and posteriors reproduce the fit
- the negative-binomial dispersion step never falls back to a lower Poisson limit (F4)

Result: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 54 ]`.

Each fix was reverted by hand and the matching tests then failed:

- F1 and F3 reverted together: 5 failures/errors (ordinal rescaling x2, empty-class NaN, seed-17 singular, aliased gaussian).
- F2 alone: 2 failures (gaussian LL; poisson "EM likelihood decreased").
- F4 alone: 2 failures.

After each check the files were restored, and `diff` against the backup was empty.

## Related existing tests

These were run after all four fixes, each file through `testthat::test_file`.

| file | passed | failed | errors | skipped | warnings | seconds |
|---|---|---|---|---|---|---|
| test-review-1026-mixreg.R | 54 | 0 | 0 | 0 | 0 | 5.0 |
| test-mixture-regression.R | 78 | 0 | 0 | 0 | 0 | 13.1 |
| test-mixture-regression-labels.R | 7 | 0 | 0 | 0 | 0 | 0.6 |
| test-mixture-regression-prediction-design.R | 12 | 0 | 0 | 0 | 0 | 0.2 |
| test-negative-binomial-regression.R | 17 | 0 | 0 | 0 | 0 | 1.0 |
| test-ordinal-regression.R | 38 | 0 | 0 | 0 | 0 | 11.1 |
| test-review-1015-regression.R | 30 | 0 | 0 | 0 | **1** | 4.0 |
| test-review-regression.R | 50 | 0 | 0 | 0 | 0 | 0.2 |
| test-growth-mixture.R | 122 | 0 | 0 | 0 | 0 | 24.3 |
| test-growth-predict.R | 21 | 0 | 0 | 0 | 0 | 3.9 |
| test-multilevel-growth.R | 33 | 0 | 0 | 0 | 0 | 10.8 |

**The new warning comes from F3's fix, and it is a true warning.** It is raised at `test-review-1015-regression.R:122`, in `grade ~ hours`, ordinal, K=2, on `study_hours`.

- Class 2 has 0 rows in the top category (the class-by-category table shows 97 rows in class 1 and 0 in class 2). Its top threshold MLE is therefore +Inf: genuine separation.
- The old raw ridge damped that low-curvature coordinate, stopping it at 24.73 with no warning.
- The standardised step lets it run to 30.28, so `latents_separation` now fires. The LL moves from -747.567194612276 to -747.567194604271, so the fit is the same.
- The test still passes. The test file is outside my partition, so I left it alone. The lead may want to wrap that call in `expect_warning(..., class = "latents_separation")`.

## Could not check

- No external NB-mixture or ordinal-mixture reference is installed, so K >= 2 for those families is checked only through numDeriv SEs, invariants and the one-class limits.
- The Bernoulli identifiability statement is checked against the literature, not by computation: identifiability is a property of the parameter map and cannot be settled by a finite run.
- F5 and F6 are not fixed (shared inference code; they change SE numbers).
- The full test suite and R CMD check were not run, only the related files listed above.
