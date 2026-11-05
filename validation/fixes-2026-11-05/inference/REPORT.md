# Inference and equivalence repairs

The four inference findings reproduced in the 5 November review are repaired. Unsupported combinations now give explicit classed refusals, and the OPG grand-transform correction uses the derivative of the estimating score. Related documentation and regression tests are included. Existing Benjamini-Hochberg residual defaults and the earlier held-categorical-information fix are preserved.

## Changes in inference

| Area | Current behavior | Source |
|---|---|---|
| Staged or partly held measurement with estimated step-one uncertainty | `three_step(step_one = "estimated")` refuses with `latents_unsupported_three_step`. `step_one = "fixed"` remains available. The correction cannot reconstruct measurement uncertainty from a covariance that omits the held parameters. | `R/three-step.R` |
| Separated step-three logits | Pooled, group-class and transition regressions refuse Wald inference with `latents_boundary_fit`. Checks cover coefficients above 20 in orthonormal design coordinates and a category's average conditional probability at or below 1e-6. A single extreme covariate row is not itself a boundary. | `R/three-step.R`, `R/lta-three-step.R` |
| External-measurement staged bootstrap | Both stages must contain the same multiset of indicator observations. External or unverifiable measurement is refused before resampling with `latents_unsupported_inference`. Different row order, identifier names, factor-level order, grouping and sequence metadata do not alone disqualify the same observations. | `R/bootstrap-inference.R` |
| OPG with estimated grand transforms | The requested OPG covariance is retained, but transform sensitivity uses inverse observed information. The original probability-SE counterexample now has ratios of one with and without grand centring. The same repair covers covariate fits. | `R/inference.R`, `R/covariate-inference.R` |
| Weighted grand transforms | Mean and population-SD nuisance equations use the fitted row weights and their observed weight totals, matching the centring/scaling changes integrated by the root agent. Weighted robust probability SEs remain invariant to grand transformation. | `R/inference-service.R` |

The staged bootstrap comparison uses exact observed values, including missing-value positions within each indicator row and duplicate row multiplicity. Categorical columns are compared by their category labels, so reversing factor-level declarations does not change sample identity. The first stage of `fit_staged()` has one group class and unconditional profile measurement, so identifier and sequence metadata do not change those measurement observations.

The original rare-transition case required an average-probability check: its tiny destination probability exceeded 1e-6 on a few covariate rows, so a check requiring every row to be tiny missed it. The regression now refuses this case under both observed and robust covariance, rather than returning the earlier misleading significant sandwich test.

## Additional verified repairs

- Binary categorical outcome contrasts count each complementary pair of category rows once for multiplicity adjustment. All six checked corrections now agree with the equivalent zero/one numeric outcome, including Holm, Bonferroni and BY. Outcomes with more than two categories still adjust across every reported category difference, as documented.
- A modal categorical overall test warns when an absent class category makes its logit undefined; its unavailable test remains `NA`.
- The three-step chain's homogeneous, first-order structure and pooled classification-error matrix are documented. Occasion-specific measurement is refused because one pooled error matrix cannot represent changing occasion-specific errors.
- Step-one documentation explicitly states that its contribution uses model-based observed information even when the step-three sandwich clusters the observations. The numeric-ML `variance` attribute description is limited to the result that actually carries it.
- Covariate inference warns when the remaining score corresponds to more than one percent of a standard error, matching the existing profile-path warning. Tight fits in the regression control do not warn.
- Zero estimated Wald variance now produces unavailable test statistics and p-values instead of an infinite statistic arising from division by zero.
- R3step QR factorization centres and standardizes covariates before checking rank. Pooled and group-class slopes, SEs and likelihoods remain stable after shifting a covariate by 1e8. This also prevents the boundary guard from mistaking a large raw intercept for separation.

## Equivalence registry and standardization

All 44 rows in `equivalence/CLAIMS.csv` are retained. Twenty-three stale documentation phrases were replaced with surviving statements of the same model behavior. Every reference, backing test, test name and explanatory note is unchanged; the single explained `NONE` convention remains. No removed citation text was restored.

`equivalence/test-claims.R` and the equivalence README now explain that a documentation anchor can describe behavior without naming its reference. The registry continues to check exact source text and backing-test presence; it does not establish completeness or prove that a named test adequately tests the behavior.

`equivalence/test-mplus-standardize.R` now checks the package against retained Mplus STANDARDIZE output directly using population SD. The old n-versus-n-minus-one likelihood offset is removed. Both the hand calculation and built-in grand transform agree with the Mplus likelihood, and the explicit scale-factor check matches the population SD to 1e-12. The equivalence README describes the same convention.

## Verification

- `tests/testthat/test-fixes-1105-inference.R`: **75 assertions passed**, zero failures, warnings or skips. These cover independent invariance checks, live refusals before inference/resampling, same-sample positive bootstrap controls, categorical label/missingness provenance, affine covariate recoding, and multiplicity equivalence.
- The related inference/three-step/staged-bootstrap/centring/scaling suite, including grand-transform calibration: **1,610 assertions across 161 cases passed**, zero failures, warnings or skips. That run preceded the final private QR refinement; the new regression file was rerun afterward.
- Filtered equivalence (`claims` and `mplus-standardize`): **186 assertions passed**, zero failures, warnings or skips.
- Final focused QR integration (`fixes-1105-inference`, `latents-controls`, `lta-three-step`, `three-step`): **223 assertions across 47 cases passed**, zero failures, warnings or skips. This caught and repaired the one-group-class case, where there is no free group-logit block to check.
- Fresh before/after probes are retained in the review and repair directories. The after-probe observes the intended refusals and repaired OPG invariance; a classed refusal caught by its harness is an expected result, not a test failure.

Evidence: `tests.log`, `test-results.rds`, `related-tests.log`, `related-test-results.rds`, `claims-standardize-tests.log`, `claims-reconciliation.log`, `probes-after.log`, and `reconcile-claims.R`. The final focused QR integration results are in `final-qr-tests.log` and `final-qr-tests.rds`.

All fitting used one core per fit with `VECLIB_MAXIMUM_THREADS=1`, `OMP_NUM_THREADS=1` and `NOT_CRAN=true`. No git mutations were performed. Root integration owns generated Rd files, NEWS, the condition catalog, package-wide verification and the overall report.
