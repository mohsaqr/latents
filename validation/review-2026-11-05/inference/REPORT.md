# Inference and three step review

Reviewed 5 November 2026 on the working tree reporting version 1.1.16, with committed HEAD `a5551081`. Four important inference issues remain reproducible. The held categorical measurement fix in the current uncommitted changes is correct on an independent synthetic case.

## Findings

### High priority Three step regression silently reports significance at separation

**Location:** `R/three-step.R:1125` through `1141`. Related transition path: `R/lta-three-step.R:264` through `280`.

An almost separating predictor is accepted by `r3step()` without a boundary warning. Convergence code zero and an invertible numerical Hessian are insufficient to establish regular inference when a logit tends to infinity. The sandwich can shrink while the inverse information grows because the individual scores vanish at separation.

Fresh synthetic case, seed 15055, 400 persons and two well separated Gaussian profiles: the slope is **302.381**. The default robust SE is **10.029**, giving **p = 1.03e-199**. Requesting observed information on the identical fit gives SE **1,998,670.74**, **p = 0.999879**. Neither call emits a warning. These are unreliable Wald results at a boundary; the observed result is a diagnostic comparison, not an asserted valid reference distribution.

Require a boundary or separation check before reporting Wald inference. Withhold the affected tests and issue a classed condition when finite interior estimates cannot be established. The same guard is absent in the transition path by source inspection; that path was not independently rerun here.

**History:** the pooled path and absence of this guard predate both 1.1.9 and the 0.9.21 submission. It is a remaining defect, not a newly demonstrated regression in the November changes.

### High priority Estimated first stage uncertainty omits held measurement

**Location:** `R/three-step.R:593`, `605` and `606`.

`three_step(step_one = "estimated")` accepts a `fit_staged()` result but obtains its covariance and derivative coordinates only for the parameters free in the second stage. The measurement parameters whose uncertainty the option should propagate remain held. The result is still labelled `step_one = "estimated"`, with no warning.

Fresh synthetic case, seed 15053, 100 groups of 10 persons, two profiles and two group classes:

| Fit | Unconstrained covariance size | Class 1 SE fixed | Class 1 SE estimated | Class 2 SE fixed | Class 2 SE estimated |
|---|---:|---:|---:|---:|---:|
| Joint | 15 by 15 | 0.055828 | 0.057419 | 0.058068 | 0.060349 |
| Staged | 3 by 3 | 0.056501 | 0.056599 | 0.056826 | 0.056964 |

The two estimators differ, so the joint SE is not presented as an exact reference for the staged estimator. The decisive evidence is that all held measurement coordinates are excluded from both the derivative and its covariance. The small increase in the staged rows propagates structural uncertainty only.

Refuse this combination until the first-stage covariance and required propagation are available, or explicitly offer a correction whose name and documentation describe the narrower uncertainty it includes.

**History:** introduced with the option in commit `4c8f5aa7`, 1 November, version 1.1.11; this is after the 1.1.9 baseline.

### High priority Staged bootstrap changes the source population of external measurement

**Location:** `R/bootstrap-inference.R:137` and `154`; the held-fit exemption is at `199`.

`fit_staged(measurement = external_fit)` permits measurement fitted on other data. Bootstrap inference retains the externally held coefficients as the point estimates, but refits the measurement on each resample of the current data. The resulting replicates describe a different estimator and can be far from the supplied measurement.

Fresh synthetic case, seed 15054: measurement was fitted to a separate sample with indicator `a` shifted by 3. The returned profile mean is **4.4851**, while its percentile interval is **[1.2876, 1.7242]**. The other profile's mean is **1.5413**, with interval **[-1.4389, -1.2649]**. All **16 of 16** replicates are accepted; 10 of 18 point estimates lie outside their intervals. Sixteen replicates do not establish interval coverage, but the population shift and the source code establish the estimator mismatch.

Reject external-measurement staging in this bootstrap, or define separate resampling of the external sample. A conditional bootstrap should retain the held measurement throughout. Any provenance comparison needs to tolerate equivalent row order and explicitly handle categorical data.

**History:** staged refitting and the exemption already exist in 0.9.21. The supplied measurement replay changed on 22 October (`f0a8f9ca`), before 1.1.9. This is a remaining defect in the reviewed surface, not a new November regression.

### High priority OPG uncertainty changes under an irrelevant indicator shift

**Location:** `R/inference.R:1165`; the same choice is made in `R/covariate-inference.R:595`.

The estimated grand-transform correction uses the inverse outer product of scores as its sensitivity matrix when `vcov_type = "opg"`. The derivative of the estimator with respect to the transform requires the inverse observed information. Mixing probabilities are unchanged by grand centring of the Gaussian indicators, yet their OPG SEs change.

Fresh synthetic case, seed 15052, 70 groups of eight observations:

| Probability | OPG SE without centring | OPG SE with grand centring | Ratio |
|---|---:|---:|---:|
| Profile given group class 1 | 0.042708 | 0.047657 | 1.115876 |
| Profile given group class 2 | 0.045854 | 0.050904 | 1.110130 |
| Group class | 0.098905 | 0.091273 | 0.922835 |

Observed and robust SE ratios equal one to printed precision in the same experiment. Retain the selected covariance estimator for its own variance, but use inverse observed information for the transform sensitivity. The covariate path has the same source pattern; its quantitative effect was not rerun here.

**History:** `R/inference.R:1165` dates to 26 October (`d5d4ae82`), before 1.1.9 and after 0.9.21.

## Confirmed fixes and controls

- **Held categorical information is fixed in the current tree.** On seed 15051, five binary indicators, two profiles, and a held measurement, independently reconstructed two-logit likelihood and `numDeriv::hessian()` give SEs 0.1275603 and 0.1401586. Package SEs agree; unconstrained `vcov()` is correctly 2 by 2, and `confint()` returns two rows. Relevant fix: `R/covariate-inference.R:790`, adding `response_logit` to the held labels. This bug was introduced with the expanded fixed-measurement covariate path in 1.1.10 and is repaired by the uncommitted 1.1.16 changes.
- **Bootstrap replicate identity is protected.** The fresh bootstrap converts to 288 replicate rows (16 times 18 parameters). Reversing the parameter table then requesting replicates is refused with `latents_bad_argument`, as intended. No plot was generated or altered.
- **Observed and robust mixing-probability SEs respect grand-centering invariance** on the OPG counterexample, providing controls for the transformation and fit matching.

## Scope and reproducibility

The review traced current inference, held-parameter, staged bootstrap and step-one propagation paths, compared their history against `a10a919b` (0.9.21) and `8518dffe` (1.1.9), and ran new synthetic probes with seeds 15051 through 15055. Earlier reports supplied leads; every numerical result above comes from the fresh script in this directory.

Run from the repository root:

```sh
VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 NOT_CRAN=true Rscript validation/review-2026-11-05/inference/probes.R
```

`probes.R` contains the simulations, an independent likelihood/Hessian check, working controls and observations of the remaining defects. `probes.log` retains the output, including conditions and actual output dimensions. Assertions validate the repaired conditional inference and replicate conversion/refusal. No production source, generated documentation, existing tests or git state were modified by this review slice.

This slice did not rerun coverage simulations, external software comparisons, all covariance structures, Rubin pooling, transition regression separation, or every minor finding in the 4 November reports. The parent review runs package tests separately. The numerical counterexamples establish reproducible failures in these combinations; they do not establish the incidence of those failures across ordinary applications.
