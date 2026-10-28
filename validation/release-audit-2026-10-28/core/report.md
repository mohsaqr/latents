# Core correctness and numerical accuracy review — 2026-10-28

Reviewed current 0.98.0 production source (`e907bff`) and recent core changes since `d5d4ae8`. No production files changed. The existing uncommitted release-engine test edit was left intact.

## Confirmed finding

**CORE-AUDIT-1 — P2: zero-weight groups can make an additive model fail.**

Location: `R/additive.R:274-275`, with the consequential assertion at `R/additive.R:300-302`.

The additive E-step computes posterior probabilities as `exp(log_joint - group_log_likelihood)`. For excluded groups far from the estimated classes, the log densities have sufficiently large negative magnitude that subtracting the log-sum-exp loses the small normalizing term. The posterior row then no longer sums to one within the asserted tolerance. The assertion covers excluded groups, so a group of weight zero can abort an otherwise valid fit. This violates the documented exclusion rule even though every supplied indicator is finite.

Reproduction: `reproduce-additive-zero-weight.R`, output in the adjacent `.log`. Seed 271 generates 30 groups of four observations, ordinary continuous values, all weights one. `multilpa(..., family="additive", n_group_classes=2, n_starts=2, seed=3, weights="w")` converges at log likelihood **-178.4893059150279**. Add two groups with weight zero and observations around a finite sentinel value:

| Excluded values centred at | Result |
| --- | --- |
| 1,000 | Same log likelihood |
| 10,000 | Same log likelihood |
| 100,000 | `simpleError`: `posterior rows must sum to one` |
| 999,999 | Same error |
| 1,000,000 | Same error |

These are fixture-specific tested points, not a universal failure threshold. The practical trigger is a large finite missing-value/sentinel code retained on intentionally excluded groups. Users can work around it by removing zero-weight groups before fitting. The same probe through the profile and full cross-level engines preserves the baseline likelihood through excluded values of 1e9.

Recommended repair: normalize posterior probabilities using row-max-shifted exponentials and their row sums, while retaining the separate stable log-sum-exp likelihood calculation. Carry the exclusion convention through any remaining zero-times-nonfinite calculations. Add a regression test comparing the fitted positive-weight sample before/after excluded sentinel groups. This audit did not apply or validate a repair.

## Independent validation passed

1. **All fourteen Gaussian covariance families, 42 fits and inference calls.** `probe-inference-units.R` generates a synthetic two-class, two-indicator sample (450 observations, seed 842), then fits each of EII, VII, EEI, VEI, EVI, VVI, EEE, VEE, EVE, VVE, EEV, VEV, EVV and VVV at scales 1, 1e-6 and 1e6. All converged and all supplied inference. `summarize-inference.R` transforms estimates and standard errors back and compares each estimate using its own baseline SE: maximum discrepancy **7.78e-7 SE**; maximum relative standard-error discrepancy **1.18e-7**. Assertions require both below 0.001. Full per-family results: `inference-units-summary.csv`.
2. **Weighted FIML against an independent analytical reference.** `probe-analytic-reference.R` uses 500 synthetic observations, unequal positive weights, ten zero weights, and missingness in both Gaussian indicators. The one-class diagonal fit agrees with direct weighted Gaussian ML formulas: means and variances within 1e-10, summed observed-data log likelihood within 1e-9. Assertions passed.
3. **Six multilevel families under small units.** `probe-multilevel-units.R` fits profile, additive, dispersion, additive-dispersion, full cross-level and restricted cross-level models to the same 60-group synthetic sample at scales 1 and 1e-6. All twelve fits converged. `summarize-multilevel.R` checks the density Jacobian, counting both ratings and derived group means for cross-level working likelihoods. Profile and interior dispersion inference succeeded at both scales. Additive and additive-dispersion estimates were at variance boundaries and consistently refused Wald inference; cross-level inference was consistently refused as documented.

## Reproduction

The scripts use the fresh audit install at `/private/tmp/latents-release-audit-20261028/lib` (or `pkgload::load_all()` in the exploratory zero-weight probe). Run from the repository root, with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1`:

```sh
Rscript validation/release-audit-2026-10-28/core/reproduce-additive-zero-weight.R
Rscript validation/release-audit-2026-10-28/core/probe-inference-units.R
Rscript validation/release-audit-2026-10-28/core/summarize-inference.R
Rscript validation/release-audit-2026-10-28/core/probe-analytic-reference.R
Rscript validation/release-audit-2026-10-28/core/probe-multilevel-units.R
Rscript validation/release-audit-2026-10-28/core/summarize-multilevel.R
```

## Coverage and limits

Read the new acceleration charts, EM and start-selection rules, zero-weight handling, core setup/finalisation, additive and cross-level engines, and relevant observed/robust inference code. The numerical panel specifically exercises common changes of units; it does not certify arbitrary unequal indicator rescaling, all missingness patterns, every boundary, or external-package equivalence. Parent reviewer owns the complete test suite and platform/release checks; another reviewer owns external equivalence and general LTA inference. No separate parallel-process suite was run here to avoid duplicating the parent's complete suite. No additional substantive core defect was confirmed.
