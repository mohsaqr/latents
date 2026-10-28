# ML three-step / VLMR review, 2026-10-20

Reviewed HEAD e29a310: R/three-step-ml.R, ML dispatch and variance/contrast integration in R/three-step.R, R/lmr.R, inference attributes consumed by VLMR, public get_results integration, and feature tests. No source changes. Read LEARNINGS.md and HANDOFF.md. Targeted feature tests: 85 expectations pass (`tests-output.txt`). Numerical/public probes add cases beyond the passing suite.

## Confirmed public defects

### P2: VLMR discards the inference engine's scale protection

Location: R/lmr.R:48-49 (`solve(alternative$information)` and `solve(null$information)`); propagation at line 261.

Generate 300 observations from two separated Gaussian classes (seed 21), two indicators, fit EEI 1:2 profiles with n_starts=10, seed=1, tol=1e-12. `get_results(grid, "lmr")` succeeds: LR=336.1164, reference mean=4.418811, SD=3.286323. Refit exactly the same indicators multiplied by 1e8: both candidates converge, their robust `parameter_inference()` succeeds, but `get_results(grid,"lmr")` aborts with `system is computationally singular: reciprocal condition number = 7.87581e-17`.

This is an unscaled coordinate inversion, not nonidentification. Whitening each information matrix with D=diag(1/sqrt(diag(I))) and applying the same coordinate transform to its score columns produces the original reference mean/SD and LR from the rescaled public fits. A separate positive-definite synthetic regular nested model reproduces the raw solve failure after changing one coordinate's units by 1e8.

Impact: a legitimate change of indicator or direct-effect units makes the entire comparison table fail, including other valid comparisons, despite successful model fitting and inference. The error is also not caught for the documented NA/note row behavior.

Fix direction: compute scores/information and inverses in well-scaled coordinates; use the inference service's coordinate scales or consistently whiten the matrices and scores before forming W. Handle per-comparison computation errors where appropriate. Pin public unit invariance, including both large and small indicator units and covariance/direct-effect coordinates.

Evidence: probes.R, probes-output.txt, additional-probes.R, additional-output.txt, models.rds.

### P2: ML distal-outcome SEs depend catastrophically on the outcome's origin

Location: R/three-step-ml.R:227-233: `up[j] <- up[j] + 1e-5*scale[j]` and analogous subtraction on the absolute class-mean coordinates. The scale at lines 207-210 tracks spread, but not representability relative to the coordinates' location.

On the same 300-observation public fit, `three_step(fit,data,"distal",method="ml",vcov_type="independent")` returns SEs approximately 0.090 and 0.084. Adding 1e12 to distal preserves its spread and estimates (up to the representable precision of the shifted values), but returns SEs 2177.97 and 4165.42. At offset 1e10 SE ratios already differ by approximately 10%; at 1e11 the maximum ratio is approximately 0.687. The code neither rejects nor warns about the false uncertainty.

A known-answer hard-classification example, where the correct class-mean covariance is simply the sample normal-mean covariance, yields sums of squared influence 0.02101978/0.39638209 at origin zero versus approximately 1.006e13/1.094e19 after the 1e12 shift.

Cause: the finite-difference perturbation becomes rounded, distorted or exactly zero when the mean is large relative to its standard deviation; symmetrizing and inverting this corrupted Hessian can still produce a nominal result.

Impact: confidence intervals and pairwise p-values can be grossly wrong solely from changing the origin of a continuous outcome. Outcome magnitude can correspond to epoch-style timestamps while its variation is small. This extends the unit-invariance issue already addressed by 8d6e782.

Fix direction: fit and infer the outcome in centered coordinates, preferably scaled too; translate means and intervals back for reporting. Check the conditioned information for positive definiteness rather than relying only on solve() failing. Pin both an exact hard-assignment reference and public location invariance over large origin/spread ratios.

Evidence: probes.R/probes-output.txt; additional-probes.R/additional-output.txt.

## Confirmed numerical edge cases; public trigger not reproduced

### Negative LR reaches an unsupported negative-q integration regime

Location: R/lmr.R:125 and 133; LR is passed through unchanged at lines 188 and 193.

`.latents_imhof_upper(-1,1)`, `.latents_imhof_upper(-1,-1)` and `.latents_imhof_upper(-1,c(1,-1))` all throw `latents_inaccurate_probability` while integrating [10000,100000]. Their reference distributions are well defined; for lambda=1 the exact upper tail at -1 is 1. The q<=0 branch assumes there is no oscillation, which is true only at q=0. Negative q needs symmetry/complement handling, such as P(Q>q)=1-P(-Q>-q) for continuous nondegenerate Q, or an explicit invalid-comparison result if the API disallows negative LR.

A nested model's global maximum cannot have a negative LR, but the table admits any two EM-converged neighboring fits; a local optimum can be worse than its neighbor and convergence alone does not establish otherwise. Fifty ordinary synthetic 2:3-profile one-start grids did not yield such a public example, so report the helper defect as confirmed and the public exposure as a risk, not a reproduced user failure. Per-pair `latents_inaccurate_probability` is not caught, so one such case would abort the table.

Evidence: first three entries of probes-output.txt.

### P3: the relative weight cutoff can destroy more than the advertised probability accuracy

Location: R/lmr.R:88 drops every weight <=1e-10*max(abs(lambda)).

Known answer: for independent X,Y~chi-square(1), P(X-bY>0)=P(F(1,1)>b). At b=1e10, `.latents_imhof_upper(0,c(1,-1e10))` returns 2.24875673638e-13; exact `pf(1e10,1,1,lower.tail=FALSE)` is 6.36619772346e-6. Just below the cutoff (b=1e9), the helper agrees to 1e-16. A tiny weight's probability contribution near zero is of order sqrt(abs(weight)), so small eigenvalue magnitude does not imply sub-1e-9 tail error. No public grid generating this specific weight imbalance and near-zero LR was reproduced.

Fix direction: distinguish numerical zero eigenvalues from legitimate small weights using a distribution-error bound or a stable rank treatment; include exact F cases crossing the filtering threshold. If deliberately excluding these cases, signal that limitation instead of silently reporting an inaccurate tail.

Evidence: extreme-weights.R and extreme-weights-output.txt.

## Positive checks and limits

The normal-distal score gradient and fixed-error likelihood formulas are coherent on inspection and their existing known-answer/Latent GOLD tests pass. ML pair contrasts use the joint influence matrix, preserve covariance between means, and follow the existing CR1 clustering path. Existing ML tests cover moderate shifts, multiplicative unit changes, both variance models and a Latent GOLD comparison; the missing stress dimension is origin divided by spread. Existing Imhof tests cover q>=0 and mixed weights with substantial magnitude; both identified integration edges lie outside that grid. Baseline and coordinate-whitened public VLMR results agree, giving independent evidence for the scale-inversion diagnosis. Root reviewer owns the full package suite and other feature families.
