# Review: additive group-class, cross-level and structure-inference math (slug: additive)

Date: 2026-10-26. Reviewer partition: R/additive.R, R/additive-inference.R, R/additive-compare.R,
R/additive-methods.R, R/cross-level.R, R/structure-inference.R.
Every check ran with OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 under `devtools::load_all()`. Scripts are in
`/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/caf01cef-53f3-419a-95fe-ebd08f111cb5/scratchpad/check*.R`.

## Scope
- Codex diff in 8ac3576 for my files: R/additive-inference.R (centred, unit-scaled information chart; Jacobi-conditioned
  inverse, OPG, rank check and Newton step), R/additive-compare.R (min_variance nesting rule; `tolerance = 0` data
  identity; `search`/`n_screen` allowed in enumeration), R/additive-methods.R (weights passed to `.additive_prepare`;
  `tolerance = 0` data identity; sampling-weight comparison). Codex did not change additive.R, cross-level.R or
  structure-inference.R.
- I also reviewed every likelihood, E/M step, score, information, packing, df count, KKT boundary pass, Newton polish,
  cross-level likelihood/BIC, and the structure-chart map, Jacobian and gradients.

## Verdict on each Codex change

| Change | Verdict | Evidence |
|---|---|---|
| `.additive_information`: centre the averages and means, scale the mean coordinates by sqrt(W+T), absolute step | **Correct, and it fixes a real defect** | With the HEAD~1 code, moving a converged 2-class additive_dispersion fit by +1e6 made inference refuse ("observed information is not positive definite (eigenvalue ratio -1)"), for observed, robust and opg alike. With the new code the SEs are invariant: max relative SE change 2.1e-11 (observed), 5.6e-11 (robust), 6.2e-11 (opg) at +1e6, and 8.3e-9 / 1.8e-8 / 4.8e-9 at +1e8. Scale x10: means SE x10 and variances SE x100 to 9.7e-12 (observed). SEs agree with the inverse of a numDeriv Hessian of the total log likelihood to 4.4e-8 relative. The information agrees with -numDeriv::hessian to <= 1.7e-8 relative in all 5 structures. |
| Jacobi conditioning of the inverse, the refusal eigen-ratio and the Newton solve | **Correct** | S C^-1 S = I^-1 exactly. Diagonal scaling cannot hide genuine singularity. The reported `condition_ratio` is now scale-free (0.1167464 at the original data, at +1e6 and at x10). Behaviour change: `condition_ratio` and the 1e-10 refusal now apply to the conditioned matrix, not the raw one. That is a deliberate improvement, but the reported attribute values differ from 0.9.21. |
| OPG `chol2inv(chol(meat*SS'))*SS'`, rank check on column-scaled scores | **Correct** | Same mathematical quantity. Rank is unchanged by column scaling, and the scaling makes the qr tolerance meaningful. |
| `.additive_nested`: alternative `min_variance` > null `min_variance` means not nested | **Correct direction** | The alternative's space must contain the null's. A higher floor in the alternative excludes the null's low within variances. Tested: `.additive_nested()` is FALSE and `bootstrap_lrt()` raises `latents_bad_argument`. The rule is conservative (it refuses even when neither floor is active), which is acceptable. |
| `same_statistics(..., tolerance = 0)` in bootstrap_lrt and parameter_inference(data=) | **Regression (low)** | See A2. The same data with rows reordered within groups is refused. |
| parameter_inference(data=) now prepares with `x$weights` and compares sampling weights | **Correct** | Weights were previously ignored by the check. The 1e-12 tolerance for weights next to tolerance 0 for averages is inconsistent but harmless, because preparation is deterministic. |
| enumeration accepts `search`, `n_screen` | **Correct** | Both are in `.additive_arguments()` and are passed through to `multilpa()`. |

## Findings

| id | severity | file:line | problem | reproducing input; actual vs expected | status |
|---|---|---|---|---|---|
| A1 | low | R/additive.R:808 (now :817; flag computed at :76-82) | `zero_within_scatter = stats$pooled_scatter <= 0` is the same exactness pattern the lead found. If an indicator is exactly constant within every group but the group means are not representable, the scatter is about 1e-31, not 0, so the flag is FALSE. | 20 groups of 3. y2 = 0.1/0.7/1.3/2.9 repeated within each group, y1 ~ N(0,1). `pooled_scatter` = c(29.94, 1.877782e-31); flag = c(FALSE, FALSE); expected c(FALSE, TRUE). The fit does hit within_floor, so only this diagnostic field was wrong. | **fixed**: `.additive_prepare()` now records an exact per-indicator test on the raw ratings (`nrow(unique(cbind(index, x[, r]))) == n_groups`, i.e. one distinct value per group), and the fit reports that field. No numerical result changes (scatter is untouched). Regression test fails with the old line (actual FALSE FALSE) and passes with the fix. |
| A2 | low | R/additive-methods.R:309-310; R/additive-compare.R:74-75 | `all.equal(..., tolerance = 0)` requires bit-identical group means. The same data with rows permuted within groups gives means that differ by 1 ulp, so the data check refuses. The fit's own likelihood is row-order invariant. | Draw 60 groups, then permute rows within each group (group order kept). `ids` identical; max abs average difference 4.44e-16. `parameter_inference(fit, data = permuted)` and `bootstrap_lrt(f1, f2, data = permuted)` both raise `latents_bad_inference_data`. With HEAD~1 (default all.equal tolerance) these were accepted. | **reported, not fixed**: choosing the tolerance is a judgement call. Suggestion: an elementwise tolerance of a few ulps of each column's magnitude, e.g. `abs(a - b) <= 64 * .Machine$double.eps * (abs(a) + 1) * n_j` for averages (scatter similarly), shared by one helper `.additive_same_statistics()` used in both places. |
| A3 | info | R/additive.R (min_variance) | `min_variance` is an absolute floor, so boundary status is not scale-equivariant. | Data x 1e-4: every within variance equals 1e-6 and `within_floor` is all TRUE, with `latents_boundary`. | documented argument; not a defect. Users with small-scale data need a smaller `min_variance`. |

No other defects found.

## Algorithms verified correct (with the check used)
- **Sufficient-statistic likelihood = row-level Gaussian likelihood.** Compared with `mvtnorm::dmvnorm` of each group's ratings under W I + T 11' (unbalanced sizes 2-7, H = 2, class-specific W and T): max |log density diff| 1.42e-14. Mixture LL -913.9202 both ways. (Searle, Casella & McCulloch 1992, the balanced one-way random effects determinant/quadratic decomposition.)
- **Scores (Fisher identity) = numDeriv gradient** in all 5 structures (additive varying/equal, dispersion, additive_dispersion varying/equal): max abs diff 6.7e-8 to 4.1e-7 (numDeriv noise; relative <= 1e-6 in the test).
- **pack/unpack** round trip 1.4e-16. **n_parameters = length(theta)** for all 5 structures (11, 9, 9, 13, 11).
- **The EM/Newton maximum is a maximum.** `optim(BFGS, reltol 1e-14)` from the fitted theta gains at most 3.4e-13 LL in all 5 families. Plain EM is monotone (smallest step >= 6e-12) for dispersion (ECM), additive_dispersion and additive.
- **KKT boundary pass.** Data with zero between variance on y2: the fit holds both y2 between variances at 0 with kkt = TRUE. Interior optim from T = 1e-3 drives them to 0 and gains -8.6e-7 (no improvement).
- **Sampling weights.** Integer weights reproduce duplicated groups: means 2.5e-13, within 6.9e-14, between 2.2e-14, class weights 1.4e-14. LL ratio 1.86 = mean weight 1.86.
- **Invariance.** Refit on data + 1e6: LL identical (-902.8096), non-mean estimates 6e-11, SE 5.9e-11 relative. Data x 10: estimates 9.1e-15, SE 1.1e-11 relative. Row permutation: LL identical. Class relabelling: natural SEs permute to 7.9e-13.
- **Cross-level.** The full and restricted LLs equal hand-written brute-force likelihoods, sum_j log sum_h omega_h f(m_j|h) prod_i sum_k pi_{k|h} f(x_ij|k) and the factorized product: diffs 0, 6.8e-13, 2.3e-13, 1.1e-13 (K = 2; H = 2 and H = 1). With weights the weighted pseudo-LL equals sum_j w_j l_j (w normalized to J): diffs 2.0e-12 (restricted) and 1.1e-13 (full). Full H = 1 equals restricted H = 1 (diff 7.1e-10). Restricted equals the single-level profile fit of the ratings plus the single-level fit of the group means (diff 8.8e-10). BIC uses log(J) and bic_individual uses log(N) (exact). df: full 19, restricted 18 (hand count).
- **Structure charts.** `.multilpa_chart_map` Jacobian equals the numDeriv Jacobian for all 10 constrained structures at d = 3, k = 3 (max diff 5.8e-12 to 2.3e-9). Chart widths equal `mclust::nVarParams()` (1, 3, 5, 7, 8, 10, 12, 12, 14, 16). The mapped Sigma are symmetric PD and satisfy their constraint (shared volume/shape spread <= 1.3e-15). The gradient (P S P - w P)/2 and the mean gradient equal numDeriv derivatives of the weighted Gaussian log density (tolerance 1e-8 passes).

## Tests added
File `tests/testthat/test-review-1026-additive.R` (10 tests, about 3.5 s):
1. sufficient-statistic density equals the dense Gaussian density of the ratings
2. analytic scores and information match numerical derivatives in every family
3. standard errors are translation invariant and scale equivariant (fails on the HEAD~1 inference code with "The observed information is not positive definite (eigenvalue ratio -0.0281)")
4. relabelling the group classes permutes the natural standard errors
5. fitted maxima cannot be improved and integer weights equal duplicated groups
6. an exactly constant within-group indicator is flagged as zero scatter (A1 regression; fails before the fix: actual FALSE FALSE)
7. a higher alternative variance floor breaks nesting, by class (error path `latents_bad_argument`)
8. cross-level likelihoods equal their brute-force definitions
9. covariance-structure charts have analytic Jacobians and mclust's counts
10. the covariance gradient is the derivative of the weighted Gaussian log density

Result: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 81 ]`

## Related existing test files run (NOT_CRAN=true, after the fix)
tests/testthat: additive-compare 29, additive-engine 30, additive-families 24, additive-fit 43, additive-methods 82,
cross-level-thin-profiles 51, cross-level 14, structure-inference 42, review-structures 84, covariance-structure 167,
review-1015-accessors 34. Every file had FAIL 0, ERR 0, WARN 0, SKIP 0 (the counts are passes).
tests/equivalence (helpers sourced): additive-engine 7, additive-families 20, additive-fit 13, additive-methods 4,
structure-inference 56. Every file had FAIL 0, ERR 0.

## Files changed
- R/additive.R: `.additive_prepare()` adds the exact `zero_within_scatter` field (and its @return doc). `.additive_fit()` reports it. No numerical result changes.
- tests/testthat/test-review-1026-additive.R: new.

## Not checked, and why
- Wald coverage calibration (`.additive_weak_class_threshold = 50`) depends on Codex's simulation study. I did not rerun it (too long for this pass).
- The SQUAREM accelerator's validity and restore functions: checked only indirectly, since accelerated fits reach the optim-verified maximum.
- Cross-level families have no inference (refused by design), so there are no SE checks.
- The structure-chart encode/decode round trip was checked only by the existing test at d = 2, not at d = 3.
- Latent GOLD / external oracle comparison: not run.
