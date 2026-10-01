# Validation evidence that is not a comparison

Everything that compares multilpa with other software, published output or an
independent implementation lives in [`equivalence/`](../equivalence/README.md).
What stays here is evidence of a different kind. Like `equivalence/`, this
folder is excluded from the build (`.Rbuildignore`) and kept under version
control.

| Path | What it is | Run from the project root |
|---|---|---|
| `MATH_AUDIT.md` | The mathematical audit: numerical fixes, independent regression checks and what they established. | — |
| `ADDITIVE_MODEL_DESIGN.md`, `additive-design-check.R` | Proposed raw-rating additive model and validation gates; checks exact likelihood, intercept moments, interior EM, and balanced one-class boundary/interior MLEs. The estimator now ships as `multilpa(family = "additive")` (experimental). | `Rscript validation/additive-design-check.R` |
| `ADDITIVE_SIMULATION.md`, `additive-simulation.R`, `additive-simulation-*.csv/.txt` | Predeclared 11-cell simulation of `multilpa(family = "additive")`, 1000 datasets per cell: recovery, SE calibration and coverage (observed and robust), refusals; registry written before running. | `Rscript validation/additive-simulation.R full` (~27 min on 9 cores) |
| `simulab-recovery/` | A Monte Carlo recovery study. It measures bias and coverage against the generating values, not agreement with a fixed number. | `Rscript validation/simulab-recovery/recovery-study.R` |
| `mixture-regression-recovery.R` | Monte Carlo recovery and 95% interval coverage for the two-level mixture regression, the one `mixture_regression()` model with no external implementation to compare against. 200 replications, labels aligned to the truth; writes `mixture-regression-recovery.csv`. | `Rscript validation/mixture-regression-recovery.R` |
| `synthetic-demo.R` | The synthetic demonstration data, also the input to the generated Mplus runs in `equivalence/mplus/generated/`. | `Rscript validation/synthetic-demo.R` |

The harness's API check (`equivalence/harness/registry.R`,
`check_validation_api()`) scans the scripts here as well as those in
`equivalence/`, so a renamed argument is caught in both.
