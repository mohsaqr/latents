# Round 2 review: kernel (EM driver, Gaussian measurement, covariance structures, noise, prior, fixed blocks, centring/scaling)

Date: 2026-10-26. Reviewer: round-2 agent, area **kernel**. No file under `R/`, `tests/`,
`equivalence/` or `man/` was edited. Every number below was produced in this session with
`NOT_CRAN=true OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 Rscript`, under
`pkgload::load_all()` (R 4.5.2, mclust 6.1.2, mvtnorm 1.3.6, numDeriv 2016.8.1.1). Scripts and
their outputs are in
`/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/caf01cef-53f3-419a-95fe-ebd08f111cb5/scratchpad/round2-kernel/`
(written `<S>/` below). Each script `NN-name.R` has its output in `NN-name.out`.

## Scope (files)

R/kernel-em.R, R/kernel-blocks.R, R/kernel-structures.R, R/acceleration.R,
R/covariance-structure.R, R/gaussian-moments.R, R/gaussian-diagnostics.R,
R/fixed-measurement.R, R/noise.R, R/prior.R, R/robust.R, R/scaling.R, R/centering.R. I also read
the callers in R/multilpa.R (`.multilpa_expectation`, `.multilpa_maximization`, `.multilpa_em`,
the noise, prior and centring set-up in `multilpa()`) and R/multilpa-starts.R (parameter count).

Round 1 (`area-reports/kernel.md`: K1 fixed; K2, K3 reported) is not re-reported. I re-ran its
suites (see Verified correct). The structural results of round 1 still hold. The findings below
concern **units**: absolute constants compared with quantities whose units depend on the data,
and transformations whose outcome depends on rows that should not matter.

## Findings

| id | sev. | file:line | what the code does | what it should do (reference) | reproduction | measured effect | proposed fix |
|---|---|---|---|---|---|---|---|
| K2-1 | **medium** (wrong refusal and wrong M-step on plausible data) | R/covariance-structure.R:239 (also :231), `.multilpa_vei_variances` | The VEI fixed point floors the **unnormalised shape sum** `colSums(W_k / lambda_k)` at `min_variance`. That sum is about `n * B_j`, a count times a determinant-one shape, so it has no units. `min_variance` is in squared data units. When the indicators' scales differ widely, the sum for the smallest-scale indicator falls below 1e-6 while every variance is far above the floor. The floor then moves the fixed point, so the M-step is no longer a maximizer. | Celeux & Govaert (1995) `[lambda_k B]`: `B ∝ sum_k W_k / lambda_k`, with no floor. The caller's guard (`diagonals / weights > min_variance`) already guarantees positivity. This is the same defect K1 fixed in the EVI/VEI scatter, but in a different quantity. | `<S>/07-vei-shape-floor.R`, `<S>/10-tiny-units-lowered-floor.R` | (a) M-step on diagonals with variances (1e-4, 1e14) and (3e-4, 2e14), weights (60, 40): the package's first-indicator variances are 12.98x and 10.74x the free fixed point's (relative differences 11.98 and 9.74), and the second indicator's are 44% and 54% low. Its objective is worse by 183.27. Against `mclust::mstep("VEI")` on simulated data of these scales, the maximum relative difference is 13.45. (b) End to end: a rate (SD 0.010 to 0.015) and a population (SD 1e7 to 1.5e7), SD ratio 1.16e-9, n = 250, VEI. `multilpa()` fails with "All 10 starts failed: EM likelihood decreased beyond numerical roundoff." (also with mclust's start). `mclust::me` reaches LL −3803.703 and `mclust::Mclust` −3803.711. At SD ratio 1.16e-8 both packages agree (difference −1.8e-12). (c) Iris ×1e9 with `min_variance` scaled to the units (1e-6·c²): the VEI log likelihood is 44.85 below the unit-mapped value. Every other structure is within 1.7e-4. | Replace `min_variance` by `.Machine$double.xmin` in both `pmax()` calls (lines 231 and 239), as was done for K1. The guard at :169 keeps the diagonals positive. Re-run `<S>/07` and `<S>/10`: VEI should match `mclust::me`, and the gap should be 0. |
| K2-2 | **medium** (wrong answer on plausible data; only a generic warning) | R/multilpa.R:1381 (`min_variance = 1e-6`), applied at R/multilpa.R:124, R/gaussian-moments.R:65-68, R/covariance-structure.R:184 | The variance and eigenvalue floor is an absolute 1e-6 in squared input units, and its docs say so. Data recorded in small units hit it during ordinary fits. The only signal is `latents_boundary`, worded "A covariance eigenvalue reached min_variance; this is a bound-active constrained fit". The documentation of `boundary` attributes that to a profile shrinking around one observation, not to the units of the data. | A maximum-likelihood fit should not depend on the units the data are recorded in. mclust has no absolute floor. G1 (round 1) fixed the same class of defect for growth models. | `<S>/01-small-units.R` | Iris, 3 profiles, LL mapped back to centimetres (should be 0): **VVV in metres (×0.01): −1.00**, ×0.001: −541.0, ×1e-4: −1756.1. VVI ×0.001: −414.4, ×1e-4: −1629.5. EEE ×0.001: −464.9. Profile means move by up to 0.62 cm at ×0.001 and 1.25 cm at ×1e-4. Each of these fits warned only `latents_boundary`. | Make the default floor scale-equivariant without breaking the argument: for example `min_variance = NULL`, meaning `1e-6 × min_j var(x_j)` computed on the fitted indicators. A floor proportional to a common variance keeps all 14 structures equivariant under a common rescaling, and a numeric value keeps its present meaning. At minimum, have the boundary warning report the floor relative to the smallest indicator variance and suggest rescaling. |
| K2-3 | **medium** (degenerate data from a documented transform, silent in mixed designs) | R/centering.R:100-167 (`.multilpa_indicator_scales`, refusal at :153) and :54-79 | Under `centering = "person", scaling = "person"`, a group with exactly two observations of an indicator becomes exactly ±1/√2 on it. Two-wave units therefore turn into 2^d atoms. The refusal covers only groups with one observation or a zero SD. | A Gaussian mixture likelihood is unbounded on atoms. The transform should refuse these groups, as it refuses singletons, or warn with a class. | `<S>/05-person-centring-atoms.R` | 150 groups of 2, d = 2: there are 4 distinct transformed rows among 300. The fit reaches LL **+1270.88** with every variance at the floor (1e-6). With 100 groups of 5 plus 30 groups of 2 and K = 3, there is **no warning**. Profile 2 puts 28.9% of its mass on the size-2 rows (10.7% of all rows). Its mean (−0.848, 0.675) is pulled to the atom (−0.707, 0.707), and its variances are (0.077, 0.028), against (0.128, 0.142) when those groups are removed. Person centring alone with singletons did not collapse (LL −971.5, min variance 0.595). | Under `scaling = "person"` with `centering = "person"`, refuse with `latents_bad_data`, naming the count, any group with fewer than three observations of an indicator. Alternatively, raise a classed warning and document that such groups contribute only sign patterns. |
| K2-4 | low (needs a user-lowered floor and units of about 1e-9) | R/covariance-structure.R:593 (`shape_and_volume` in `.multilpa_eve_covariances`) | The EVE/VVE conditional shape floors the rotated scatter diagonals at the absolute `.Machine$double.eps`, which is in squared units × counts. | The M-step must satisfy M(c²W) = c²M(W). | `<S>/11-eve-absolute-floor.R`, `<S>/10-tiny-units-lowered-floor.R` | M(c²W)/c² against M(W): relative error ≤ 2.9e-15 for c² down to 1e-16, but **0.531 (EVE) and 0.491 (VVE) at c² = 1e-18**. End to end, on iris ×1e-9 with `min_variance = 1e-18`: EVE fails ("All 3 starts failed: EM likelihood decreased"), and VVE's Jacobian-mapped LL is 170.3 below the cm fit. | Use a relative floor, `pmax(values, max(values) * .Machine$double.eps)`, or `.Machine$double.xmin`. |
| K2-5 | low (documentation contradicts the math) | R/multilpa.R:1037-1041 (`@param scaling`), man/multilpa.Rd | The text says grand scaling changes only units for EEI, VVI, EEE, VVV, EVI and VEI, and changes the solution for EII and VII. It omits VEE and EVV, which are also equivariant, and implies that EEV, VEV, EVE and VVE are equivariant, which they are not. Under Σ → DΣD a shared spectrum or a shared eigenbasis is not preserved. | Equivariant: EEI, VEI, EVI, VVI, EEE, VEE, EVV, VVV. Solution changes: EII, VII, EEV, VEV, EVE, VVE. | `<S>/02-grand-scaling-equivariance.R` | Scaled LL minus the unit-mapped raw LL, columns rescaled by (1, 5, 0.2, 2): EII 282.1, VII 292.2, **EVE −3.94, VVE −13.12, EEV 7.41, VEV 11.46**. The eight equivariant structures are within 5.7e-13. | Correct both lists in the roxygen and regenerate the Rd. |
| K2-6 | low | R/kernel-em.R:176 (`settled <- abs(gain) <= tol * (1 + abs(LL))`) | The relative stopping rule depends on units, because rescaling y → cy shifts LL by −n·d·log c and so loosens or tightens the tolerance. This is the rule G2 changed for growth models only. | A unit-free rule, e.g. per observation: `|gain| <= tol * max(1, n_effective)`. | `<S>/03-stopping-rule.R` | Overlapping 2-profile VVI, n = 3000, distance at stop from a 1e-15 reference, in Wald SE: default SQUAREM 0.0063 (×1), 0.022 (×0.01), 0.017 (×1e3, ×1e6). Plain EM (`acceleration = "none"`): 0.062, 0.082, 0.137 and **0.189** at ×1, ×0.01, ×1e3 and ×1e6. Both report `converged = TRUE`. | Either adopt a unit-free tolerance or document that `tol` is relative to |LL| and so depends on units. The default (SQUAREM) path is at most 0.022 SE off. |
| K2-7 | low (cosmetic, contradicts docs) | R/kernel-em.R:158-161 | Every accelerated cycle adds 3 to `iterations`. A cycle that falls back to plain EM (no valid extrapolation, `norm_v = 0`) took two M-steps. `?multilpa` says "`iterations` and `max_iter` count EM steps under either setting". | Count the M-steps actually taken. | `<S>/14-squarem-iteration-count.R` | A toy accelerator that is never valid: 20 M-steps taken, `iterations = 30`, history length − 1 = 30 (the history duplicates the second LL). | Return the number of steps from `.latents_squarem_cycle()` and add that. |

Already reported elsewhere this round, confirmed independently here, and **not counted again**: multilevel
ZW-4. `centering`/`scaling = "grand"` compute the offset and SD over every row, including rows of
weight zero (R/centering.R:145-147, :174). In `<S>/13-zero-weight-grand-scaling.R`, adding 4
zero-weight groups (12 outlying rows) to a weighted VVI fit under `scaling = "grand"` changed the
reported LL from −666.997 to −0.994 and BIC by −1332.0. Under EII the same addition changed BIC
by −1326.6 and the means by 0.68. With `centering = "grand"` alone the LL is unchanged, but the
reported means move by 1.29. Person centring and scaling are unaffected (difference 0).

Info, not a defect: under `fixed = "means"` the held means return with a drift of 2.8e-17
(one ulp, from the internal translation round trip, `<S>/12-fixed-means.out`).
`?fit_staged` says "exactly the first stage's".

## Verified correct

| algorithm | check | max discrepancy |
|---|---|---|
| 14 covariance M-steps, complete data, vs `mclust::me` from a common z (`equivalence/test-mclust.R`) | re-run this session | 4 tests, 41 passes, 0 fail |
| Noise component: EM, hypervolume, **parameter count** (= `mclust::nMclustParams(noise = TRUE)`), BIC, for all 14 structures (`equivalence/test-mclust-noise.R`) | re-run | 14 passes. LL ≤ 2.05e-12 for the 12 closed and fixed-point structures, with `parameters` difference 0 for all 14. EVE −3.30e-11; VVE +2.17e-3 (higher than mclust). BIC ≤ 4.1e-12 for the exact 12 |
| Noise + MAP prior together (10 prior structures) | same file | LL ≤ 9.09e-13; means, covariances and posteriors ≤ 5.2e-14 |
| MAP prior M-steps and EM vs `mclust::me(prior = priorControl())` (`equivalence/test-mclust-prior.R`) | re-run | 51 passes. One M-step: covariances ≤ 1.42e-14. Whole EM: LL ≤ 1.41e-10 (EEV), posteriors ≤ 8.4e-10 |
| Prior bounded fallback (EVI/VEI pseudo-weights) | derivation: the pseudo-weights reproduce the prior denominators `dof + n + 1 + bump` (EVI) and `dof + n_k d + 2 + bump` (VEI) in the volume stationarity equations | exact (algebra) |
| Two-level E-step with FIML + categorical + sampling weights + VEV, H = 2 | direct formula `L_g = Σ_h ω_h Π_i Σ_k π_{k|h} f_k(y_i,obs) P_k(c_i)` with `mvtnorm` marginals (`<S>/15`) | LL 1.1e-13; group posteriors 6.7e-16; row posteriors 3.3e-16 |
| FIML M-step under constrained structures EVI, VEI, VEV | direct `optim` (BFGS, 4 starts) of the observed-data LL in each structure's own chart (`<S>/08`; the first run's EVI row used a wrong chart and is superseded by `08-fiml-constrained-EVI.out`) | EVI −3.1e-12, VEI 2.0e-12, VEV 4.9e-12 |
| Integer sampling weights = duplicated groups (FIML, two-level, VEV/EVE/EEV, VVV) | weighted solution evaluated and re-iterated on the duplicated data (`<S>/06`) | LL ratio = Σw/G = 1.916667 exactly. Duplicated-data EM from the weighted solution: LL gain ≤ 2.3e-10; means ≤ 4.8e-8; covariances ≤ 2.7e-8 |
| SQUAREM monotonicity and its maximum, on 11 paths the harness does not reach (noise, FIML VVV/EVI/VEV, EVE, VVE, fixed means, two-level mixed categorical, two-level weighted, person centring+scaling, local dependence) | history and tight plain EM from the same start (`<S>/04`) | min step gain ≥ 0 on every path. Accelerated minus plain-tight LL ≥ −2.24e-4 (VEV-FIML, VVE), all converged |
| SQUAREM S3 step (Varadhan & Roland 2008) | read: α = −‖r‖/‖v‖, α = −1 gives θ₂ exactly, backtracking (α−1)/2, acceptance only if the landed LL ≥ the second step's | — |
| `fixed = "means"`: constrained M-step around held means, parameter count | numDeriv gradient of an independent LL in the free coordinates at the fit (`<S>/12`) | |∇| ≤ 3.4e-6 (VVI), 1.7e-7 (EVI), 5.9e-6 (VVV); count reduced by K·d = 4 |
| Grand scaling: equivariant structures | `<S>/02` | EEI, VEI, EVI, VVI, EEE, VEE, EVV, VVV within 5.7e-13 |
| Equivariance of all M-steps at units ×1e-9 and ×1e9 with `min_variance` scaled (excluding K2-1 and K2-4) | `<S>/10` | ≤ 1.65e-4 (VVE at ×1e9; EEV 1.44e-4 at ×1e-9 and 4.1e-5 at ×1e9, all GEM or stopping-point differences). All others ≤ 4.4e-6 |
| EVE/VVE M-step equivariance down to c² = 1e-16 | `<S>/11` | 2.9e-15 |
| Round-1 kernel regression file `tests/testthat/test-review-1026-kernel.R` | re-run (NOT_CRAN) | 10 tests, 102 passes, 0 fail |
| Read and judged correct | `.latents_em` stopping and decrease guard (matches `@param tol`); EM checkpoint/resume; FIML conditional moments and the expected-scatter M-step (`E[(x−μ)(x−μ)'] = (e−μ)(e−μ)' + C`); eigenvalue clipping (exact under an eigenvalue floor); local-dependence block masking; robust diagonal log-variance scores `0.5(r²/v − 1)`; noise split and the noise-only Gaussian denominators (`sum(weights)` in the equal path); `.multilpa_nested_structure` offsets | — |

## Not checked

- Inference for the structure charts (R/structure-inference.R) and the full-covariance robust
  scores are outside these files. Round 1 checked the scores against numDeriv; I did not re-run that.
- The prior's bounded VEV fallback and the EVI/VEI prior fallback were checked by algebra
  only, with no numerical case where the floor binds under a prior.
- The numbers in `.latents_gaussian_diagnostics` were checked by reading only (GD-1 was fixed in round 2).
- The emEM screening and ranking logic (`.latents_run_screened`) and `.latents_quasi_newton`
  were read, not exercised. Only growth and LTA use the latter, and those belong to other areas.
- I did not measure K2-6 for engines other than `multilpa()`. The rule is shared, so the
  same unit dependence applies wherever |LL| moves with units.
- I did not re-run the full package suite.

## Fix status (same day, after the review)

All seven findings were fixed. Regression tests are in
`tests/testthat/test-review-1026-kernel-r2.R` (12 tests, 54 passes). To check that each test
fails on the old code, a copy of the package was made outside git, the old lines of one fix were
restored in it (`<S>/sabotage.py`), and the file was run against the copy:

| Fix restored to the old code | Tests that failed |
|---|---|
| K2-1 | 3 |
| K2-4 | 2 |
| K2-6 kernel rule | 2 |
| K2-7 | 1 |
| K2-2 | 4 |
| K2-3 | 1 |

| id | fix | measured after |
|---|---|---|
| K2-1 | shape sums floored at `.Machine$double.xmin` | VEI M-step = free fixed point (objective difference 0); rate/population VEI fits, LL = `mclust::me` to 4.5e-13; VEI = `mclust::mstep` to 4.3e-10 at ×1e-9, ×1, ×1e9 |
| K2-4 | EVE/VVE floor `.Machine$double.xmin` | M(c²W)/c² vs M(W) ≤ 3.4e-15 at c² = 1e-18; objective ≤ mclust's at ×1 and ×1e9 (mclust non-finite at ×1e-9) |
| K2-6 | `|ΔLL| ≤ tol · max(1, observed values in rows of positive weight)`, passed by every engine (`stop_scale`) | plain EM: 126 steps and 0.048 SE from the maximum at ×1, ×1e-6, ×1e6 alike (before 0.062 to 0.189) |
| K2-7 | `iterations` = M-steps kept | toy: 20 = 20 |
| K2-2 | `min_variance = NULL` → 1e-6 × smallest indicator variance (`multilpa`, `lta`, `fit_staged`, additive and cross-level families, covariate path); mixture regression and growth already relative | iris VVV/VVI/EEE at ×0.01, ×1e-3, ×1e-4 within 4.3e-8 of the cm fit (before up to −1756) |
| K2-3 | `latents_bad_data` when every group has two observations; `latents_atomic_groups` warning otherwise | 150 groups of 2 refused; 30 of 130 groups (10.7% of rows) warned |
| K2-5 | `@param scaling` lists the equivariant and non-equivariant structures | — |

Observation, not fixed: the SQUAREM steplength `α = −‖r‖/‖v‖` is a norm over coordinates in
different units (means, variances and probabilities). The accelerated path, and therefore its
stopping distance, still depends on units: 0.0063 SE (×1), 0.018 SE (×1e-6) and 0.00027 SE (×1e6)
from the maximum. A unit-free norm, with each block in the scale of its indicators, would remove
this.
