# Latest-model review, 1–2 October 2026

Base revision: `38450580fefbf16c0eceb80958707e861f4678c5`, version 0.9.7.
The review covers the eight releases 0.9.0–0.9.7, their changes since
`d7890cc`, and integration with the older engines they modified. Findings
were reproduced with synthetic data before fixes. Changes are uncommitted.

## Release inventory

| Release | Changes inspected and exercised |
| --- | --- |
| 0.9.0 | Additive, dispersion, additive-dispersion and restricted/full cross-level families; inference/refusals, comparison, enumeration and weak-class contracts. |
| 0.9.1 | General LTA covariates, occasion transitions/measurement, second order and covariance structures; path likelihoods, scores and retained references. |
| 0.9.2 | LTA enumeration/BLRT, mover-stayer, FIML and probability boundaries; nesting, original-data identity and simulated refits. |
| 0.9.3 | Sampling weights throughout every engine and inference/reporting; replication equivalence, sandwich checks, safe normalization and weighted effective counts. |
| 0.9.4 | Newton membership-logit M-step and convergence/scaling suites. CI configuration was inspected; this review runs locally on R 4.5.2 and does not claim a fresh R 4.1 or hosted CI run. |
| 0.9.5 | Ordinal and Poisson measurement, wrappers, mixed indicators, inference, prediction, matching and seven retained measurement references (including later NB additions). |
| 0.9.6 | Explicit unsupported-covariance refusal for membership-covariate fits; retained guard tests and supported full-covariance checks. |
| 0.9.7 | Ordinal/count covariate and LTA plumbing; NB equal/varying dispersion, derivatives, boundaries, simulation/refits and singleton reporting. |

## Model coverage

| Model or extension | Authoritative review and new evidence |
| --- | --- |
| Single-/two-level Gaussian profiles: EII, VII, EEI, VEI, EVI, VVI, EEE, VEE, EVE, VVE, EEV, VEV, EVV, VVV | `test-review-structures.R`: direct multivariate Gaussian densities for all 14 weighted multilevel likelihoods and integer-replication M-steps. Existing structure/formula, inference-chart, boundary, missing-covariance, prior and noise tests are retained in the full run. |
| Categorical/mixed, ordinal, Poisson, NB with shared/varying dispersion | [Indicator review](../review-2026-10-01-extra.md); `test-review-extra.R` and `test-review-integration.R`: stable derivatives, categorical identity, wrapper parity, prediction/relabel invariance, NB bootstrap. Seven bundled Latent GOLD fixtures remain in `test-extra-indicators.R`. |
| Membership covariates at individual/group levels, shared/group-specific slopes | Indicator review; weighted full-covariance/FIML analytic-score and Jacobian checks, singleton inference; all existing covariate fixture tests. |
| Additive, dispersion, additive-dispersion, equal/varying between variance | [Family review](../review-2026-10-01-families.md); `test-review-families.R`: all five distinct restrictions, independent dense covariance likelihood/score/Hessian, EM-only stationarity, sandwich covariance, data/nesting validation. |
| Restricted/full cross-level, equal/varying within and between variance | Family review: all eight restrictions independently checked against raw-rating working likelihoods and full numerical gradients. Ordinary Wald inference is deliberately refused. |
| Homogeneous/occasion/covariate transitions, initial covariates, invariant/occasion measurement, second order, mover-stayer, grid gaps/FIML, structures, extra indicators, weights | [Transition review](../review-2026-10-01-lta.md); `test-review-lta.R`: exact path sums, per-group/total scores, used occasion coordinates, reporting, boundary and data/nesting validation. Existing external LTA fixtures are enabled in the full run. |
| Gaussian/Poisson/binomial regression mixtures at row/group/two-level nestings | `test-review-regression.R`: nine independent weighted direct likelihood/score checks away from maxima; weighted class/composition tables and invalid-control regressions. Existing regression likelihood, lm/glm, prediction, label and inference tests remain. |
| Shared weights, comparison, bootstrap, prediction, sensitivity and imputation label matching | `test-review-integration.R`: two-/three-class ordinal reference rebasing and NB block permutation; original bootstrap/pooling/sensitivity suites; stable extreme-weight scaling. |

## Integration fixes

- `lca()`, `multilca()` and `enumerate_lca()` partition ordinal/count columns
  from the remaining categorical indicators, implementing the documented support.
- Profile alignment includes ordinal probabilities and count means in its matching
  signature. It moves NB means/dispersions with profiles and rebases ordinal
  locations/intercepts to the new last-profile reference. Both ordinary bootstrap
  and covariate imputation alignment use the correction.
- Row alignment verifies ordinal/count input values. Cross-sectional bootstrap
  compares indicator specifications and reconstructs NB data with its fitted
  configuration, instead of evaluating an NB fit with a Poisson density.
- One-profile varying NB dispersion coordinates retain profile labels so the
  count table can match and report their standard errors.
- Extended transition coefficient vectors use complete measurement covariance
  charts for all fourteen structures and reconstruct the fitted likelihood.
- Single-level covariate dispatch retains its `single_level` identity.
- Sampling-weight normalization handles finite weights as large as `1e308`
  without overflowing their sum or multiplication.
- Regression mixture class counts/shares, group-class composition/shares and
  smallest-share diagnostics use weighted posterior totals. Unit posterior
  tables remain normalized; modal counts and classification quality remain
  descriptive per-unit quantities.
- Regression counts/iteration/bootstrap controls and seeds now validate before
  integer conversion, instead of silently truncating fractions or overflowing.

## Evidence files

`tests-current.log/.rds/.csv` is the initial full-suite baseline (4625
assertions, 766 cases, zero failures/warnings/skips); it was started before
all fixes and is **not** the final integrated-source gate. The first attempted
combined review run (`review-tests.log`) recorded a test-reference vapply shape
error while the transition review was still in progress; this was repaired.

Final results are recorded in `tests-final.log/.rds/.csv`, `SOURCE.csv`,
`SESSION.txt`, `lint-final.log`, `package-build.log`, `package-check.log`,
the copied native `00check.log`, `package-check-summary.csv`,
`package-tests.log` and `package-examples.log`.
The six new review test files add **670 assertions across 47 cases**.
The combined source suite passed **5,295 assertions in 813 cases**, with zero
failures, errors, warnings or skips. Its original verifier then failed because
it had included test-generated `tests/testthat/Rplots.pdf` in the source
manifest. A direct comparison confirmed that every one of the 273 actual
source files was unchanged. The corrected verifier excludes that artifact and
also explicitly checks test errors; `source-verification.log` records its
successful check of the original results and source manifest. No production
source change occurred between these checks. The original console failure is
retained in `tests-final.log` so its cause remains visible.

Use the final results and source manifest when evaluating completion. Per-agent
records distinguish earlier checks from final integrated checks.

Reproduce the source suite:

```sh
NOT_CRAN=true Rscript validation/review-2026-10-01/verify.R
Rscript -e 'x <- lintr::lint_package(); print(x); stopifnot(length(x) == 0L)'
```

`verify-archive.R`, `ARCHIVE-SOURCE.csv` and `archive-identity.log` confirm that
272 built-package files (R, Rd, tests, NAMESPACE and NEWS) match the reviewed
source exactly. DESCRIPTION is excluded from byte comparison because R adds
build metadata.

The final native package check reports **Status: OK — zero errors, warnings
and notes**. It covers installation, examples, tests and vignette rebuilding.
`record-package-check.R` validates that status and retains its structured
summary and logs. Package build/check commands are recorded in `HANDOFF.md`.
Generated package/check working directories are under ignored `tmp/`.
The processx-based wrapper failed to load under the sandbox (`Operation not
permitted`); `package-check-wrapper.log` records that environment issue. The
native R commands are the package gate. The installed-package test run reports
5,268 passes, zero failures/warnings and four expected skips for unavailable
R/vignette source files. Those four source-dependent cases are covered by the
full source run above, where nothing was skipped.

`pure-extra-probes.R/.log/.rds` additionally check ten core/covariate combinations
(pure ordinal, Poisson, NB, categorical and mixed discrete) for typed all-table
results, complete unconstrained coefficients, reproducible refits and simulated
rows; core prediction densities also pass. Natural core coefficients include
reference probabilities and thus have more entries than `n_parameters`.
Covariate fits have their own class and refitting surface; they do not currently
export a `predict()` method. Some deliberately simple pure-item fits produce
expected singular-information or NB Poisson-limit notices/refusals, retained in
the probe log; those are distinct from full-suite warning failures.

## Limits of the evidence

No new external-program run or large recovery study is claimed. The review
reruns retained executable reference fixtures and adds independent formula,
likelihood, score, curvature, stationarity and integration checks. Historical
simulation limitations remain those documented in `ADDITIVE_SIMULATION.md`
and `LTA_SIMULATION.md`: weak/small classes, variance ML bias, and unstudied
coverage combinations. The original Latent GOLD family comparison had quadrature
mismatches; its refinement records explain them. Cross-level likelihoods reuse
ratings through manifest means and are comparable only within that family on
the same data. Unsupported features/refusals are separate from defects fixed
here. A fresh agent should inspect the fixes and rerun the listed commands,
then assess these limitations before making broader statistical claims.
