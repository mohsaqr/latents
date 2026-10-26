# Review 2026-10-26: shared EM kernel and Gaussian measurement model (slug: kernel)

## Scope
Files: R/kernel-em.R, R/kernel-blocks.R, R/kernel-structures.R, R/gaussian-moments.R,
R/covariance-structure.R, R/acceleration.R, R/centering.R, R/scaling.R, R/fixed-measurement.R,
R/gaussian-diagnostics.R, R/prior.R, R/noise.R, R/robust.R.

Functions reviewed and exercised: `.latents_em`, `.latents_em_step`, `.latents_squarem_cycle`,
`.latents_quasi_newton`, the start runners; `.latents_gaussian_log_density`,
`.latents_measurement_log_density`; `.mixture_row_lse`, `.mixture_log_softmax`,
`.mixture_normalize_log_weights`, `.mixture_structure_expectation` (observation / group / two-level),
`.mixture_multinomial`, `.mixture_update_mixing`, `.mixture_structure_scores`,
`.latents_forward_backward`, `.latents_scaled_forward_backward`, `.multilpa_nested_structure`,
`.multilpa_logit_structure`, `.multilpa_weigh`; `.multilpa_gaussian_moments`,
`.multilpa_maximize_moments`, `.multilpa_bound_covariance`; all 14 structure M-steps,
`.multilpa_structure_parameters`, `.multilpa_bounded_variances`, `.multilpa_bounded_covariances`,
the EVE/VVE MM orientation step; the SQUAREM validity map; the prior M-steps; the noise hypervolume;
`.multilpa_indicator_scales`, `.multilpa_check_centered`; `.multilpa_group_scores` and the full-covariance
group scores.

Of these files, Codex's commit 8ac3576 changed only kernel-blocks.R, kernel-structures.R and centering.R.
The other ten files have no diff in that commit.

## Verdict on each Codex change

| Change | Verdict | Evidence |
|---|---|---|
| kernel-blocks.R: offsets for the ordinal/count blocks, plus a second row-max centring when extra blocks are present | **Correct** | Gaussian (with NA) + 2 categorical + ordinal + 2 count indicators, Poisson and NB. Max relative \|(log_density + offset) - raw\| = 2.9e-16. The nested log-likelihood matches path-free brute force (prod over rows, sum over group classes) to all 12 printed digits: -440.997769377347 (Poisson) and -440.831331975574 (NB) by all three routes. Posteriors with and without offset differ by at most 2.6e-15. On extreme data (log f = -5e15 for both profiles, categorical evidence only) the posteriors equal prior x P(code) to 1.1e-16. |
| kernel-structures.R: row-max centring in `.mixture_log_softmax` | **Correct** (algebraically identical; differs only in rounding) | Covered by the brute-force structure checks below. |
| kernel-structures.R: `.mixture_normalize_log_weights` and relative densities for tau / group_tau / rho / tau_by_group_class | **Correct, and fixes a real defect** | Ordinary data, against brute-force probability arithmetic: observation tau 2.2e-16, LL 1.4e-14; group tau 8.0e-16; two-level rho 2.2e-16, tau 2.2e-16, LL 0. Per-row offsets of -5e15 x U(0.5, 1.5): observation tau 1.1e-16, group 1.1e-16, two-level rho 2.2e-16, tau 3.3e-16. The old arithmetic `exp(joint - unit)` was off by **0.437** on that input. The LL is still computed from the un-offset joint, as before. A row of all -Inf gives NaN tau and LL -Inf, the same as before the commit. |
| centering.R: relative-magnitude constancy test replaced by `sd <= 0` | **Regression** (found by the lead) | Three values of 0.1 give sd about 1e-17 and passed the test. The lead's uncommitted range-based fix is correct: `ranges` is length(present) x ncol, empty cells give NA and `%in% TRUE` turns them into FALSE, and test-review-1026-constancy.R passes (12). I did not modify it. |

## Findings

| id | severity | file:line | what is wrong | reproduction (numbers computed here) | status |
|---|---|---|---|---|---|
| K1 | medium | R/covariance-structure.R:201 (EVI), :225 (VEI) | The closed-form EVI and VEI M-steps floored the **scatter** diagonals at `min_variance` (`pmax(diagonals, min_variance)`). They are reached only when every `diagonals / weights > min_variance`. A profile with effective membership below 1 can therefore have a scatter diagonal below the floor while its variance is above it, and the floor then moved the maximizer. | weights = c(50, 0.5), scatter = rbind(c(40,60), c(8e-7,1)), floor 1e-6, EVI. **Before:** variances [2,] = (9.70e-4, 970.11), objective 97.9353613322. **Exact closed form / KKT solver:** (8.68e-4, 1084.62), objective 97.9353498430. The floor does not bind at the maximum (min 8.68e-4). VEI on the same input: [2,2] 1.000002 vs 1.000001. | **Fixed.** The floor is now `.Machine$double.xmin`; the guard already guarantees diagonals > 0. A no-op whenever every diagonal is at least `min_variance`, which covers every profile with effective membership of at least 1, so ordinary fits are unchanged. The regression test failed before the fix (3 expectations) and passes after. All three mclust equivalence suites still pass (41 / 14 / 51). |
| K2 | low | R/centering.R:110-128 (`.multilpa_indicator_scales`) | The sd is formed from raw squares, so data that does vary is refused as "zero or undefined" when its spread underflows or overflows. | `scaling = "grand"` refuses (`latents_bad_data`) both c(1e-170, 2e-170, 3e-170) (sd underflows to 0; base `sd()` also returns 0) and c(1e160, 2e160, 3e160) (squares overflow to Inf). Before 8ac3576 the overflow case was **not** refused: it divided by Inf and silently produced zeros. So the change is an improvement, but the refusal message is wrong. | reported, not fixed. Only extreme magnitudes; lead's uncommitted code; judgement call (scale before squaring, as `.multilpa_design_scale` does). |
| K3 | low | R/kernel-blocks.R:85,112; R/multilpa.R:11 | "All component densities vanished, or a density overflowed." is a plain `stop()`. A user cannot catch it by class. | style only | reported |
| K4 | info | R/covariance-structure.R `.multilpa_bounded_covariances` (VVE) | With a binding floor, the result is 1.4e-7 above the constrained optimum (115.37780037 vs 115.37780023, relative 1.2e-9). | Inside its documented stopping rule tol*(1+\|obj\|) = 1.7e-6. A generalized-EM step, so not a defect. | verified-acceptable |
| K5 | info | (mclust) | `mclust::mstep(prior = priorControl())` differs from `mclust::me` for EEI (0.066), EVI (0.089) and EEV (0.024, mean 3.7e-4). latents matches `me` (the engine `Mclust` uses) to <= 1.4e-14. | equivalence/test-mclust-prior.R with NOT_CRAN, 51 pass | verified-correct (an inconsistency inside mclust, not in latents) |

## Algorithms verified correct

- **14 covariance M-steps vs `mclust::mstep`** (n = 300, d = 3, G = 3, soft z). Max \|Sigma diff\|:
  - closed forms: EII 4.4e-16, VII 4.4e-16, EEI 4.9e-15, EVI 4.0e-15, VVI 8.9e-16, EEE 5.3e-15, EEV 8.9e-15, EVV 2.7e-15, VVV 7.1e-15;
  - fixed points: VEI 1.3e-10, VEV 3.9e-10, VEE 1.1e-6.
  - EVE and VVE reach a **lower** M-step objective than mclust (by 9.5e-6 and 4.5e-2), so mclust's inner iteration stops earlier.
  - All 14 parameter counts equal `mclust::nVarParams`.
  - The existing equivalence/test-mclust.R (41), test-mclust-noise.R (14) and test-mclust-prior.R (51) all pass, before and after K1.
- **Bounded (floored) M-steps vs `constrOptim`** in the volume/shape/angle chart:
  - EVI (floor binding): latents objective 154.6316378946 vs 154.6316435264, so latents is lower; equal volume holds to 3.3e-16.
  - VEI: 4.1e-8 apart.
  - Ellipsoidal d = 2: EVE, VEV and EEV equal to 8 decimals; VVE as K4.
- **FIML:**
  - Observed-pattern log densities vs `mvtnorm::dmvnorm`: 1.8e-15.
  - Conditional means and covariances vs the hand formulas mu_m + S_mo S_oo^-1 (x_o - mu_o) and S_mm - S_mo S_oo^-1 S_om: 4.4e-16.
- **Scores:**
  - Analytic `-.multilpa_score` vs `numDeriv::grad` of the observed LL, with 2 group classes and missing data: full/varying 5.3e-8, full/equal 3.3e-8, diagonal 6.9e-9.
  - The numDeriv self-difference between two step settings is of the same order: 2.0e-7, 9.1e-8 and 1.4e-8.
  - Robust `.multilpa_group_scores` column sums give the same differences.
- **Forward-backward** (scaled and log-space, unequal sequence lengths) vs enumeration of every path: LL 4.4e-16, occasion posteriors 3.3e-16.
- **SQUAREM** (Varadhan and Roland 2008, S3 with backtracking; an extrapolation is accepted only if it does not lower the LL). Final LL vs plain EM: VVV -775.65908186 vs -775.65908355 (72 vs 363 iterations); VEI, EVE and VVI agree to <= 6.4e-7. Every likelihood history has min diff >= 0.
- **Invariants:**
  - Scaling by 1000 and by 0.01 gives dLL = -n d log c to <= 9.3e-11 (VVV, VEI, EVE).
  - Translating by 1e7 changes the LL by <= 8.7e-9.
  - Relabelling the profiles leaves the LL unchanged (diff 0) and permutes the posteriors exactly.
- **Prior M-steps vs `mclust::me`:** <= 1.4e-14 per structure; the whole MAP EM to 1.4e-10.
- **Noise hypervolume vs `mclust::hypvol`:** the existing suite passes. A constant indicator is refused upstream (exact `unique() < 2` test in multilpa.R), so log V = -Inf cannot reach the E-step.
- **Read and judged correct:**
  - the multinomial Newton step (gradient X'(c - Tp), information X' diag(T p_a (delta_ab - p_b)) X);
  - local-dependence masking (the pairs are refused if they overlap, so the pattern is block-diagonal with 2x2 blocks, and the masked scatter with blockwise eigenvalue clipping is the exact MLE);
  - VVV/EEE eigenvalue clipping (exact under an eigenvalue floor);
  - the Celeux and Govaert (1995) EEV / EVV / EVI / VEI forms;
  - the Browne and McNicholas (2014) MM step for EVE/VVE.

## Tests added
`tests/testthat/test-review-1026-kernel.R`, 11 tests:
1. measurement offsets leave the nested likelihood and posteriors unchanged
2. an extreme Gaussian density does not round away priors and categorical evidence
3. structure posteriors match brute force, and survive a -5e15 offset per row
4. the forward-backward recursions equal path enumeration
5. FIML densities and conditional moments are the multivariate normal ones
6. the analytic FIML score is the numerical gradient of the likelihood
7. each covariance M-step reproduces mclust's, or attains a lower objective
8. a floored axis-parallel M-step stays the constrained maximizer (the K1 regression test)
9. SQUAREM reaches plain EM's maximum with a monotone history
10. the profile likelihood is equivariant under rescaling and translation
11. the noise hypervolume is mclust's, in logs

Result: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 278 ]`, about 3 s. Before the K1 fix the same file failed 3 expectations in test 8.

Existing tests run after the fix, all with 0 failures:
- default environment: test-kernel 48, test-acceleration 26 (1 skip), test-prior 25 (1 skip), test-noise 22, test-centering 34, test-scaling 31, test-scaling-invariance 37, test-gaussian-diagnostics 41, test-gaussian-diagnostic-tables 76, test-robust 35, test-missing-covariance 50, test-review-1026-constancy 12, test-lta-acceleration 59, test-lta-covariance-boundary 52;
- NOT_CRAN=true: test-covariance-structure 167, test-review-structures 84, test-fixed-measurement 71, test-structure-inference 42, test-acceleration 29, test-prior 25;
- equivalence: test-mclust 41, test-mclust-noise 14, test-mclust-prior 51.

## Not checked
- I did not check at the level of a whole EM run that K1 causes a likelihood decrease. That needs a fit in which a profile collapses to an effective membership below 1. The defect is shown at the M-step level only.
- I did not trace every caller of the measurement offsets in the transition engines (lta-engine.R and transitions.R are outside this partition). I checked only the shared block and the structures.
- I did not verify the Gaussian diagnostics table numerically beyond its existing tests (41 + 76 pass). I read it and judged it correct, including the origin shift that protects large offsets.
- The full package test suite was not run.
