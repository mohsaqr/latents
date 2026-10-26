# Review 2026-10-26: multilevel (two-level LPA/LCA engine)

## Scope

Files reviewed: `R/multilpa.R`, `R/multilpa-starts.R`, `R/multilpa-pool.R`,
`R/covariates.R`, `R/covariate-inference.R`, `R/lpa.R`, `R/start-pool.R`. I read the
helpers they call in other files (`R/kernel-structures.R` nested E-step
`.multilpa_nested_structure()` / `.multilpa_logit_structure()` / `.multilpa_weigh()`,
`R/weights.R`, `R/noise.R` `.multilpa_group_sample_size()`, `R/scaling.R`) but did not
edit them.

Codex commit 8ac3576 touched three of my files: `R/covariates.R` (+35/-5),
`R/covariate-inference.R` (+16/-3) and `R/multilpa.R` (roxygen only).

All checks were run in this session with
`OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 Rscript`, under `devtools::load_all()`, with
`cores = 1`. The scripts are in
`/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/caf01cef-53f3-419a-95fe-ebd08f111cb5/scratchpad/ml/` (`c1.R` to `c9.R`).

I did not edit any file under `R/`. The only repository file I added is
`tests/testthat/test-review-1026-multilevel.R`.

## Verdict on each Codex change

| # | Change | Verdict | Evidence |
|---|---|---|---|
| C1 | `.multilpa_weighted_logits()` drops rows whose expected counts are all zero before building the QR coordinate basis | **Correct** | Zero-count rows add nothing to the objective, the score or the information, so removing them leaves the problem unchanged; only the basis moves. I appended 3 zero-count rows with covariates of about 1e8 to 1e9 (c4.R). The new solver returns coefficients identical to the unpadded ones (max abs diff 0) and `converged = TRUE`. With the filter removed by hand (the old code), the solver stops with coefficients off by 2.49 and `converged = FALSE`, scaled score 4.85e-6. So the old code failed to converge rather than falsely converging; the fix is still right. `!any(active)` returns `initial` and is harmless. |
| C2 | `.multilpa_cov_designs()` rank check uses only positive-weight rows (stacked profile design and group design) | **Correct** | These are exactly the rows that enter the likelihood, so a predictor that is constant among contributing rows is now refused (`latents_bad_covariate`) instead of reaching a singular Newton step. Under `weights = NULL` it is the old check. Inference rebuilds the designs without `row_weights` and checks all rows. If the active rows have full rank then all rows do too, so no regression. |
| C3 | `.multilpa_cov_extreme()` gains `profile_weights` / `group_weights` and drops zero-weight rows from the separation basis; used in `.multilpa_cov_assemble()` and in `.multilpa_cov_covariance()` | **Correct** | Reproduced: one zero-weight row with covariate 1e7 makes the old verdict `TRUE` (spurious separation, which refuses Wald inference). With the weights the verdict is `FALSE`. This is now a regression test. With `sampling_weights = NULL`, `NULL[group_index]` is `NULL`, so the unweighted path is unchanged. |
| C4 | Wald `transform` uses positive-weight rows for the membership coordinate blocks | **Correct (conditioning only)** | The final covariance is `moves %*% inverse(I_moves) %*% t(moves)`, which does not depend on which invertible basis is used. The change only affects conditioning. Row bookkeeping `rep(active[group_index], G)` matches the `rbind` order of `profile_design`. Robust vcov equals a numDeriv sandwich to 1.1e-9 (see below). |
| C5 | `logLik()` / `nobs()` for `multilpa_covariates` use `.multilpa_group_sample_size()` instead of `n_groups` | **Correct (it fixes an inconsistency)** | `$bic` already used the informative-group count. Before the change, `BIC(fit)` disagreed with `fit$bic` whenever a group observed nothing. With an all-missing group under `missing = "fiml"` (c5.R): `nobs = 20`, `BIC() = $bic = 785.0454`. With `drop_empty = FALSE`: `nobs = 21`, `BIC() = $bic = 785.6796`. This matches `nobs.multilpa` and `nobs.multilpa_lta`. |
| C6 | `multilpa()` roxygen for `search`: shared covariance does not make every degenerate fit impossible | **Correct** | Documentation only. The new wording is accurate: an equal-covariance model can still reach the floor on atomic data. |

None of Codex's changes in my files repeats the lead's `sd <= 0` constancy regression.
The constancy tests in my partition are exact:

- `.multilpa_prepare_indicators()` uses `length(unique(.)) < 2`.
- The group-covariate constancy check uses an exact `!=` against each group's first row.

## Findings

| id | severity | file:line | problem | reproducing input; actual vs expected | status |
|---|---|---|---|---|---|
| ML-1 | medium | `R/weights.R:56` (normalization, outside my partition), `R/multilpa-starts.R:166` and `R/noise.R:181` (group sample size), and the same in `R/transitions.R:981` and `R/lta-engine.R:895` | Zero-weight groups are counted as units. The weights are normalized to sum to `n_groups` including zero-weight groups, and `n_informative_groups` counts them too. Padding a data set with zero-weight groups therefore leaves every estimate unchanged but inflates the pseudo log likelihood by `J_all / J_positive` and raises `log(n)` in BIC. Since −2·LL and the penalty scale differently, this can change which number of classes BIC prefers. Codex's own tests assert that zero-weight groups do not change coefficients or robust vcov, but nothing checks LL, nobs or BIC. | c3.R: 15 groups with weights in {1,2,3}, plus 4 zero-weight copies. Means identical (max diff 2.2e-15). LL −289.7777 vs −228.7719 (ratio 1.266667 = 19/15). `nobs` 19 vs 15. BIC 611.9442 vs 487.3323. Expected: LL, nobs and BIC equal to the unpadded fit. | **reported, not fixed.** It needs a coordinated change to `R/weights.R` (normalize over positive-weight units) and to the group sample size in every engine. That is outside my partition, and it changes the numbers of existing weighted fits that have zero weights. |
| ML-2 | low | `R/multilpa.R:1963` (`.multilpa_prepare_indicators`) | The constant-indicator check runs on all rows, including zero-weight rows. An indicator that is constant among the contributing rows but varies in zero-weight rows is not refused. The fit goes to the variance floor and only warns `latents_boundary`. Dropping the zero-weight rows instead gives the `latents_bad_data` error. This breaks the "zero-weight rows influence nothing" contract that the audit established. | c8.R: `y2 = 0.1` in the 10 weight-1 groups, `rnorm` in the 10 weight-0 groups. Result: `latents_boundary` warning, `y2` variances = 1e-06 in both profiles. Expected: `latents_bad_data` "Constant indicators cannot identify Gaussian profiles". | **reported, not fixed.** The check runs before the weights are resolved (`.multilpa_setup` line 1703 vs line 1723), so a fix means reordering setup, and the call would change from returning a fit to erroring. That is a judgement call. |
| ML-3 | low | `R/multilpa.R:452,458,505` (`.multilpa_validate_start`) | Symmetry and equal-across-profiles checks of a supplied start use an absolute 1e-12 tolerance, which depends on scale. | Not reproduced with realistic input. `starting_values(fit)` round-trips cleanly at indicator scales 1, 1e4 and 1e6, for both equal and varying full covariance (c7.R; the asymmetry of produced covariances is exactly 0, refit LL diff ≤ 5.9e-9). It could only affect hand-built starts on very large or very small scales. | reported (observation), not fixed |
| ML-4 | info | `R/covariates.R:174` | The `value > baseline + 1e-8` guard is absolute, while the accept rule allows `1e-12·(1+|value|)` per step. Only a concern for objectives of about 1e7 or more. | Not reproduced. | observation only |

No critical or high defect was found in my partition.

## Algorithms verified correct

| Algorithm | Check | Result (max abs difference) |
|---|---|---|
| Two-level likelihood L = ∏_j Σ_h ω_h ∏_i Σ_k π_{k\|h} f_k(y_ij) (Vermunt 2003), computed in log space | Brute-force enumeration of all (h, z_1..z_4) configurations for 12 groups of 4 (c1.R) | LL diff 0. Group posteriors 5.6e-16. Subject posteriors 5.6e-16. |
| M-step (closed form): ω_h = mean_j P(h\|y_j); π_{k\|h} = Σ P(h, z=k) / Σ n_j P(h); posterior-weighted means and variances | Fitted parameters re-derived from the brute-force posteriors | ω 5.2e-11, π 2.1e-10, μ 9.5e-10, σ² 8.4e-10 |
| Parameter count (H−1) + H(K−1) + 2KJ; BIC on the group count; `bic_individual` on the informative rows | Hand formula, plus `BIC()` vs `$bic` | 11 = 11; BIC 305.3125 = 305.3125; `nobs = 12` |
| Single-level limit (one group class = ordinary mixture) | G = 1 two-level vs `lpa()`; vs mclust `meVVI` at tol 1e-12; vs mclust df | LL identical (−230.8535 both); vs meVVI −3.1e-9; df 9 = 9. Default `Mclust` VVI is 3.4e-4 lower (looser tolerance); EEE, 3 profiles: latents −227.0919 vs mclust −227.0934 |
| Label permutation of profiles and group classes | `.multilpa_expectation()` with permuted parameters | LL 0, posteriors 0 |
| Sampling weights (pseudo ML) | Integer weights vs physically replicated groups; weights × 1e6 | means 5.5e-8, variances 1.4e-7; LL ratio = mean raw weight exactly (2.066667); ×1e6: LL and means diff 0 |
| Membership Newton solver | 3-category weighted multinomial vs `optim` BFGS at reltol 1e-16; binary vs `glm(quasibinomial)`; vs `nnet::multinom` (exploratory, not in Suggests) | Objective lower than BFGS by 3.3e-8 (BFGS stops earlier). Matches glm to the 1e-8 test tolerance. nnet 3.5e-5 (nnet's own tolerance). |
| Covariate two-level likelihood (multinomial-logit ω_h(w_j), π_{k\|h}(x_ij)) | Independent R implementation from the reported coefficients | LL diff 0 |
| Covariate analytic scores (means, log variances, log-Cholesky, β, γ; weighted) | `numDeriv::grad` at a perturbed θ, for diag-varying, diag-equal, full, weighted and group-class slopes | 3.4e-7, 3.4e-8, 7.3e-8, 2.1e-8, 2.1e-7 (max \|grad\| about 7 to 11) |
| Observed-information vcov | `solve(-numDeriv::hessian(LL))` | 6.4e-10 (diag varying), 8.6e-8 (equal), 8.1e-9 (full), 7.7e-9 (group-class slopes); relative SE ≤ 8.7e-8 |
| Robust (CR1 sandwich) vcov, weighted | numDeriv H⁻¹ · (G/(G−1))·ΣS_jS_jᵀ · H⁻¹ | 1.1e-9 |
| Empty group handling (`drop_empty`) | Adding an all-missing group under FIML (covariate-free and covariate models) | LL equal; ω equal to 3e-7; nobs unchanged (20); covariate coefficients within 4e-5 (flat EM optimum at tol 1e-10) |
| Start reproducibility | Same seed twice (kmeans and emEM; covariate model); caller's `.Random.seed` before and after | `identical()` TRUE for starts tables, means and coefficients; global RNG untouched |

## Tests added

The new file is `tests/testthat/test-review-1026-multilevel.R`, with 13 `test_that` blocks:

1. the nested likelihood and both posteriors equal brute-force enumeration
2. a converged fit is a fixed point of the closed-form M-step
3. one group class is the single-level mixture, at least as high as mclust (skip_if_not_installed mclust)
4. the E-step is invariant to relabelling profiles and group classes
5. integer group weights equal replicated groups and weights are scale free
6. the parameter count and BIC use the group sample size
7. an all-missing group changes neither the fit nor the group sample size
8. the membership Newton step matches glm and an independent optimizer
9. rows without expected mass do not move the membership step (regression test for C1 and C3)
10. the covariate likelihood equals an independent two-level computation
11. covariate scores and observed information match numerical derivatives (skip_if_not_installed numDeriv)
12. the covariate structure is invariant to relabelling the group classes
13. seeded two-level fits are reproducible and leave the caller's stream alone

Result of `testthat::test_file("tests/testthat/test-review-1026-multilevel.R")`:
`[ FAIL 0 | WARN 0 | SKIP 0 | PASS 55 ]`, in about 5 s.

**Sanity check that the tests can fail.** I temporarily commented out Codex's two
zero-mass filter lines in `R/covariates.R` by hand. The file then gave
`[ FAIL 1 | WARN 0 | SKIP 0 | PASS 54 ]`, failing at line 268 (test 9). I then restored
`R/covariates.R` from a byte copy, and `git diff --stat R/covariates.R` is empty again.

## Related existing test files

I ran these with `NOT_CRAN=true`, one `testthat::test_file()` call per file under
`load_all()`. The runner script is `ml/run_related.R`; its output is in `ml/related.txt`.

| file | test blocks | failed | errors | skipped | warnings | expectations passed |
|---|---|---|---|---|---|---|
| test-multilpa.R | 17 | 0 | 0 | 0 | 0 | 163 |
| test-covariates.R | 5 | 0 | 0 | 0 | 0 | 23 |
| test-covariate-inference.R | 8 | 0 | 0 | 0 | 0 | 40 |
| test-covariate-separation.R | 4 | 0 | 0 | 0 | 0 | 12 |
| test-covariate-slopes.R | 7 | 0 | 0 | 0 | 0 | 20 |
| test-covariate-missing.R | 6 | 0 | 0 | 0 | 0 | 24 |
| test-weights.R | 14 | 0 | 0 | 0 | 0 | 76 |
| test-review-1015-covariates.R | 11 | 0 | 0 | 0 | 0 | 33 |
| test-review-1020-zero-weights.R | 2 | 0 | 0 | 0 | 0 | 18 |
| test-audit-core-repairs.R | 13 | 0 | 0 | 0 | 0 | 64 |
| test-zero-weight-initialization.R | 6 | 0 | 0 | 0 | 0 | 40 |
| test-review-1026-multilevel.R (new) | 13 | 0 | 0 | 0 | 0 | 55 |

No regressions. The full `devtools::test()` suite was not run.

## What I could not check, and why

- **Results with `cores > 1`.** Under `devtools::load_all()` the PSOCK workers do not
  have the development namespace (`could not find function ".latents_collect_conditions"`).
  Parallel equivalence needs an installed build.
- **Mplus / Latent GOLD two-level comparisons.** These were not run here; they are covered
  by the repository's equivalence harness, which I did not re-run.
- **Ward / k-means start quality.** This was reviewed for reproducibility only, as scoped.
- **Categorical (LCA) and ordinal/count blocks inside the two-level model.** These live
  in kernel files outside my partition and were not separately calibrated here.
