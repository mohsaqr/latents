# Transition fixes for the November review

The transition slice fixes the three major and six minor findings in the 4 November transition report, plus the explicit empty-row BIC and covariance-nesting inconsistencies identified during the help review. Existing estimation defaults are preserved. Source is ready for the parent agent to regenerate documentation and run integrated validation.

## Behavior and numerical changes

| Finding | Change | Effect on results |
|---|---|---|
| Near-empty origins | Transition tables expose origin_expected_count and supported. Fewer than one expected outgoing move warns with latents_weak_transition_support and withholds Wald inference. | Nonzero rows retain their formally free coordinates and IC penalty. This threshold is a weak-support diagnostic, not a declaration of structural non-identification. Exact-zero origins retain the uniform/count-reduction path. |
| Short-panel identification | Known no-covariate two-/three-occasion ridges refuse at setup with latents_unidentified. A necessary K^T-1 latent-path dimension bound also refuses excess chain coordinates. | Previously arbitrary class/stayer parameters are no longer returned as fitted estimates. Covariate-dependent identification is outside this fixed-rule guard. |
| Ragged bootstrap | Coefficient names are checked before label alignment; extraction errors/nonfinite or differently sized estimates become counted failed replicates. | Previously crashing bootstrap calls complete or raise a classed insufficient-valid-replicates error. Warnings explicitly state intervals condition on retained resamples; absent coefficient reasons remain attached. |
| Second-order first move | Its probability table labels occasion 2; first-order origin support reads only that first move. | Table label corrected. Exactly empty first-move origins are no longer incorrectly counted because of support in later second-order moves. |
| Second-order averaging | Public lta help states that the second-order table averages across groups and applicable occasions, weighted by move frequency. | Documentation only. |
| Empty interior occasions | Second-order warnings state that remaining coordinates stay counted and the parameter count is an upper bound. | Documentation/warning only. |
| Homogeneous probability boundaries | Homogeneous fits mark initial/transition/group probabilities at the effective numerical bound and warn, matching the general engine convention. Structural one-class probability one is excluded. | Boundary metadata changes; bootstrap_lrt newly refuses bound-active homogeneous inputs that the general engine already refused. Categorical measurement-only flags were not widened. |
| Extreme covariate row | Boundary checks use a category near zero across all contributing rows or extreme coefficients in the original fitting chart. | Finite interior fits with one extreme covariate row regain valid Wald output. Checking the fitting chart preserves invariance under shifts such as 1e8. |
| Trailing dot terms | Coefficient parsing/relabeling preserves complete terms, including terminal dots. | A factor level such as Zed. keeps its name, estimated status and standard errors; fitted likelihood unchanged. |
| Evaluate-only fitting | max_iter=0 general fits do not emit a false non-convergence warning. | Warning only. |
| Empty-row BIC | Both engines count all positive-weight rows for individual BIC under explicit drop_empty=FALSE. n_informative remains the number of observed rows. | Corrected homogeneous BIC changes by k log(n_all/n_observed); general BIC unchanged. Raw n_observations remains distinct. |
| Covariance nesting | Removed redundant variance_model rank check; covariance-structure nesting governs EII within EEI/EEE. Enumeration boundary diagnostics use the shared reporting predicate. | Valid nested comparisons are accepted; categorical bound status is reported consistently without widening stored flags. |

## Fresh numerical checks

- The final-only-state counterexample has 2.922e-7 outgoing expected moves in the homogeneous engine and approximately 1e-17 to 1e-24 in the occasion-specific engine. Both tables now show supported=FALSE; Wald inference is withheld. Log likelihood differences from the retained pre-fix fits are exactly zero, and the formal counts remain 38 and 50 respectively.
- The ragged 50-person case completes 15 bootstrap draws with eight usable replicates and seven recorded absent-occasion/design-term failures. The call returns a bootstrap table and a conditional-interval warning, replacing the previous unclassed signature crash. This is not a bootstrap-coverage claim.
- A finite transition slope with one destination at -6.5 SD now has boundary=FALSE and finite slope SEs 0.2551345 and 0.2748485. A factor term titleZed. has finite SEs 0.3818928 and 0.4460480, and every counted coordinate is identified.

## Files and tests

Production source: R/lta-validation.R (new helpers), R/transitions.R, R/lta-engine.R, R/lta-inference.R, R/transition-inference.R, R/lta-bootstrap.R, R/lta-methods.R, R/lta-compare.R, R/simulation-service.R.

Regression tests: tests/testthat/test-fixes-1105-transitions.R; the existing explicit transition-table schema and selective expected-warning fixtures were updated in test-transitions.R, test-review-lta.R and test-lta-compare.R.

- New regression file: 58 passed, zero failures/errors/warnings/skips; final-new-tests.log and new-tests.rds retain the results.
- Affected-file rerun: 215 passed, zero failures/errors/warnings/skips (earlier new-file count 52 plus existing transitions and review-1026-lta-open).
- Explicit-warning fixture rerun: 336 passed, zero failures/errors/warnings/skips (lta-compare and review-lta).
- Initial broader transition slice: 757 passed, one stale schema assertion, one chart-check error, 26 warnings. The two defects and all warning fixtures were repaired and verified in the targeted reruns above. Those runs overlap and are not additive.

Reproduction scripts/logs: near-empty-after.R/log/rds, near-empty-comparison.log, ragged-after.R/log/rds, covariate-after.R/log/rds, before.log and smoke.log in this directory. The parent performs final current-source integration and package checks. No man, NEWS, claims registry or root learning/handoff files were edited by this slice.
