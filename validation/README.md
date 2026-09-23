# Validation evidence that is not a comparison

Everything that compares multilpa with other software, published output or an
independent implementation lives in [`equivalence/`](../equivalence/README.md).
What stays here is evidence of a different kind. Like `equivalence/`, this
folder is excluded from the build (`.Rbuildignore`) and kept under version
control.

| Path | What it is | Run from the project root |
|---|---|---|
| `MATH_AUDIT.md` | The mathematical audit: numerical fixes, independent regression checks and what they established. | — |
| `simulab-recovery/` | A Monte Carlo recovery study. It measures bias and coverage against the generating values, not agreement with a fixed number. | `Rscript validation/simulab-recovery/recovery-study.R` |
| `synthetic-demo.R` | The synthetic demonstration data, also the input to the generated Mplus runs in `equivalence/mplus/generated/`. | `Rscript validation/synthetic-demo.R` |

The harness's API check (`equivalence/harness/registry.R`,
`check_validation_api()`) scans the scripts here as well as those in
`equivalence/`, so a renamed argument is caught in both.
