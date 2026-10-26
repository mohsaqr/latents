# Review of the Codex audit commit and a mathematical audit — 2026-10-26

Scope: commit 8ac3576 ("Fix numerical, inference, and validation defects from the
audit"), then a nine-way review of every algorithm in R/ by partition. Area reports
with reproductions and computed checks: `area-reports/*.md`.

## Verification
| Run | Result |
|---|---|
| Codex commit as delivered, fresh install, full suite | 1,527 tests, 10,794 passes, 0 fail/error/warning, 15 fork skips (Accelerate BLAS) |
| Final tree, fresh install, full suite | 1,629 tests, 11,327 passes, 0 fail/error/warning, 15 fork skips |
| New review tests (10 files under tests/testthat) | 111 tests, 729 passes, 0 fail/warning |
| `equivalence/test-review-1026-oracles.R` (depmixS4, nlme, poLCA, Rmpfr, sandwich, mvtnorm) | 9 tests, 197 passes |

Fork tests were not run under the reference BLAS in this session.

## Codex regressions (fixed)
| Id | Where | Problem | Fix |
|---|---|---|---|
| R1 | centering.R, lta-three-step.R, mixture-regression-engine.R | `sd <= 0` / `spread == 0` constancy tests miss exact constants whose mean rounds (three 0.1 -> sd 1.7e-17; LTA basis entries 7.2e16) | exact range test (`max == min`) |
| R2 | covariate-inference.R:967 | category labels unescaped while accessor keys escaped: SEs NA for categories containing `: , [ ]` | escape with `.latents_parameter_component()` |
| R3 | accessors.R | CSV round trip of person-centred data refused (reconstruction rounds once more) | accept within 0.5 unit of the 15th digit + 4 eps (|r|+|o|); changed or swapped values still refused |

## Further defects (fixed, each with a regression test that fails before the fix)
| Id | Sev. | Where | Problem |
|---|---|---|---|
| CAT-F1 | high | direct-effects.R | covariate constant in the observing rows accepted; effect -3.6e13 reported converged |
| MR-F1 | high | mixture-regression-engine.R | Codex QR shortcut: NaN / singular error for empty or thin class |
| K1 | medium | covariance-structure.R | EVI/VEI floored scatter diagonals at min_variance: 12% variance bias for near-empty profile |
| G1 | medium | growth-mixture-engine.R | absolute Cholesky jitter: small-unit outcomes broke EM; LL 380 below nlme |
| MR-F3 | medium | mixture-regression-ordinal.R | raw ridge at large predictor scales; thresholds 4.4e-4 off polr while "converged" |
| MR-F4 | medium | mixture-regression-engine.R | NB dispersion forced to floor on a local-only KKT condition, lowering the objective |
| SEL-F1 | medium | three-step-ml.R | blanket `tryCatch(error = NULL)`; narrowed to chol() rejection, results identical on 61 cases |
| SEL-F2 | medium | compare-models.R | homogeneous vs occasion lta fits of the same data refused (row names) |
| ADD-A1 | low | additive.R | zero within-group scatter (1.9e-31) not flagged |
| CAT-F2 | low | local-dependence.R | NA labels for unobserved joint categories |
| CAT-F3 | low | extra-indicators.R | 17-digit ordinal labels merged categories on redraw |
| SEL-F4 | low | selection.R | non-Okabe-Ito fill `#BDBDBD` -> `#56B4E9` (render: tmp/trust-view-colour.html) |

## Reported, not fixed (decisions needed)
1. ML-1 (medium): zero-weight units are counted in weight normalisation and nobs. Reproduced: +4 zero-weight groups to 15 leaves estimates unchanged (3.9e-16) but scales logLik by 19/15 and BIC 622.99 -> 784.74. Affects weights.R and group counts in every family.
2. MR-F5 (high, pre-existing): `.mixture_invert` singularity rule is unit/origin dependent; a calendar-year predictor gives all-NA SEs with one class. Shared with growth and LTA. MR-F6: finite-difference step not scale-aware (0.76% SE error at x1e4).
3. LTA-F1 (medium): `lta()` accepts rank-deficient covariate designs and counts unidentified logits in AIC/BIC (17 vs 15). LTA-F2: logits of all-missing occasions counted as parameters. LTA-F3: tna.R uses prior-modal states for all-missing occasions.
4. G2 (medium): EM stop relative to |LL| and raw-unit quasi-Newton; at y x1e-5 or x1e6 fits stop ~5e-3 short while converged = TRUE.
5. ADD-A2 (low): exact identity refuses the same data with rows reordered within groups (additive inference, bootstrap_lrt).
6. Low: ML-2 constant-indicator check counts zero-weight rows; INF-4 escaping incomplete for ordinal levels and continuous names; INF-5 escaped text visible in coef(); K3 unclassed "densities vanished" stop(); BVR p-values ignore estimation (documented as approximate); ICL uses soft entropy; `.multilpa_r3step_frame` turns a negative variance into SE 0; growth `.trajectory_check_data` tryCatch converts to a classed error.

## Verified correct (independent checks; see area reports for maxima)
Two-level likelihood and posteriors vs brute-force enumeration (5.6e-16); 14 covariance M-steps vs mclust::mstep; FIML moments vs mvtnorm; forward-backward vs path enumeration and depmixS4; LCA vs poLCA; NB derivatives vs Rmpfr; observed/robust/OPG covariances vs numDeriv; CR1 vs sandwich::vcovCL; Rubin's rules vs mice; IC/entropy/AvePP/OCC vs hand formulas; BCH/ML three-step vs hand computation; VLMR vs Vuong W; ICC vs anova; one-class mixture regression vs glm/glm.nb/polr; growth vs nlme/lme4/lcmm.
