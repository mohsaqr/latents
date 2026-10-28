# Inference and model release audit — 2026-10-28

Reviewed working tree at e907bff (0.98.0), including recent changes since d5d4ae8. Production code was not changed. This is a targeted numerical/static review, not an exhaustive proof of the statistical implementation.

## Findings

### INF-A1 — P1: general LTA inference depends on outcome measurement units

Location: `R/lta-inference.R:523-525` (raw mean finite-difference steps and raw information eigenvalue test); underlying rules in `R/inference-service.R:70-80,117-129`.

The general LTA inference chart standardizes transition/initial covariates but leaves continuous measurement means in their outcome units. Its fourth-order score differences use `h = 1e-4 * max(1, abs(phi))`. Changing an indicator from units near 1 to units near 1e-4 makes that step approximately an entire within-profile SD. The subsequent positive-definiteness check compares eigenvalues of differently scaled coordinates without standardizing the information.

Reproduction: `units-probe.R` simulates 80 persons x 5 occasions, two continuous indicators, two states, and a transition covariate. Fits at scale 1, 1e-3, 1e-5, 1e3, and 1e5 all converge to equivalent log likelihoods after the exact Jacobian correction (maximum difference 3.64e-12). `units-inference.R` then exactly rescales the same fitted object to remove optimization as a confound.

Observed:

| Outcome multiplier | Maximum relative SE error | Public inference |
| --- | ---: | --- |
| 1e-2 | 1.13e-9 | Returns |
| 1e-3 | 9.16e-6 | Returns |
| 3e-4 | 0.000996 | Returns |
| 1e-4 | 0.027695 | Returns incorrect covariance |
| 3e-5, 1e-5 | — | `vcov()` and `parameter_inference()` raise `latents_singular_information` |
| 1e6 | — | Same false refusal |

The transition/initial logit SEs also change at scale 1e-4 (maximum relative change 0.000971); those coefficients do not change units at all. At scale 1e-5, using a safe derivative step and correlation-scaled inversion yields SEs agreeing with the original-unit reference within 1.18e-8 relative. The raw information eigenvalues range from 6.409 to 2.988e12, while the correlation-scaled eigenvalues are 0.493 to 1.570: the model is well identified.

Expected: covariance transforms by a diagonal multiplier on mean coordinates, with all log-variance and transition/initial covariance entries invariant; any well-identified fit remains well identified under unit conversion.

Recommendation: apply outcome reference scaling or score-derived coordinate scales to the numerical differences, and condition/invert information in a standardized coordinate system. The same coordinate scaling must cover robust and OPG paths. Existing unit tests compare fitting paths and replication, not this inference behavior.

Status: newly confirmed in this review, not listed in the handoff's open issues. The problematic derivative/inversion rules predate the latest commit, so this is an existing release defect rather than a regression from the new `coef()` fix.

### INF-A2 — P2: zero-weight rows can remove mixture-regression inference

Location: `R/mixture-regression-inference.R:183-187` and calls at 202-205.

`.mixture_chart()` computes predictor and membership design centers/spreads over all rows/groups, including units with zero sampling weight. The fit now excludes these units from start search, but restores them before calling inference. An irrelevant zero-weight extreme predictor can therefore make the informative design almost constant in inference coordinates, causing a false singularity.

Reproduction: `zero-weight-inference.R` fits a one-class Gaussian regression to 100 synthetic observations with weights 1. Baseline coefficient SEs are 0.1021741 and 0.1019461. Appending a row `x = 1e8, y = 1e8, w = 0` leaves the log likelihood and coefficients exactly unchanged but returns NA for both coefficient SEs, p-values, and confidence limits, with `latents_singular_information`. At zero-weight predictor values 1e4, 1e5, and 1e6, maximum relative SE changes are 3.75e-10, 1.08e-7, and 8.22e-6 respectively; these are numerically negligible. Values 1e7 and 1e8 both cause complete loss of inference. The explicit row-deletion reference is the baseline fit; coefficients and log likelihood are exactly equal at every tested value. Recomputing inference on `.mixture_carrying_spec(fit$spec)` restores the entire baseline covariance exactly (maximum absolute difference 0).

Expected: zero-weight units have no effect on estimation or inference, as documented in the current release.

Recommendation: derive chart centers/spreads from positive-weight rows or positive-weight groups according to each design's level, or run inference on the carrying specification consistently.

Status: newly confirmed here, not in the handoff's open issues. The chart itself predates the latest commit; the new fitting exclusion does not close its inference counterpart.

## Checks that passed

- General LTA `coef()` means match the reported profile means exactly after the latest centering fix, and coefficient names match covariance dimnames. The constant translation between internal and public mean coordinates leaves covariance unchanged. Internal call sites use `.lta_theta()` or `.lta_parameters()` where appropriate.
- Independent balanced Gaussian random-intercept reference (`growth-reference.R`): one-class growth-mixture fixed effects agree with closed-form maximum likelihood within 8e-15; residual and random-intercept variances within 6.5e-10; fixed-effect covariance within 7.7e-12.
- Independent modal three-step reference (`three-step-reference.R`): class outcome means and CR1 sandwich SEs exactly equal direct calculations from modal assignments; the pairwise contrast SE agrees within 1e-12.
- Read recent changes to LTA empty-origin handling, bootstrap covariance masking and parameter extraction, growth covariance-bound reporting, mixture-regression carrying specifications, three-step thin-class handling, and inference argument validation. No additional confirmed issue in these inspected changes.
- Coordinated the equal-dispersion coefficient-name regression with the equivalence reviewer; their report owns that finding. The same naming pattern affects Gaussian `log_sigma.all` and negative-binomial `log_dispersion.all`; growth packs explicit names and does not use that fragile concatenation.

## Artifacts and limits

Scripts and logs: `units-probe.*`, `units-inference.*`, `zero-weight-inference.*`, `growth-reference.*`, and `three-step-reference.*` in this directory. Scripts load the current repository through `pkgload::load_all()` and run with one BLAS thread. No external service or reference program was needed: checks use exact reparameterizations and closed-form likelihoods.

The parent owns the full test suite and cross-platform/release checks. This review did not independently rerun every bootstrap distribution, every growth residual structure, or every three-step equivalence case. Existing handoff issues remain separate; passing the checks here does not close them.
