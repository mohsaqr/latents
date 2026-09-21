# Equivalence tests

These tests compare multilpa against other software or published output. They
are not shipped with the package: the folder is listed in `.Rbuildignore`.
The package's own `tests/testthat/` keeps only unit tests and invariant
checks, which run with nothing more than the package and testthat installed.

Run from the project root:

```
Rscript equivalence/run.R                  # every file
Rscript equivalence/run.R jstats mclust    # only files whose names match
```

`run.R` loads the package from source and then runs this folder with
testthat. It stops with a non-zero exit status if any test fails.

## What is here

| File | Reference |
|---|---|
| `test-mplus-public.R`, `test-mplus-public-parser.R` | Mplus User's Guide examples 7.9 and 7.10 |
| `test-mplus-twolevel-*.R`, `test-mplus-missing-full.R` | Two-level Mplus runs (public and synthetic data, full covariance with FIML) |
| `test-mplus-covariates.R`, `test-mplus-inference.R`, `test-mplus-robust.R` | Mplus covariate estimates and their standard errors: ML Hessian and MLR sandwich |
| `test-mplus-categorical.R`, `test-mplus-random-intercept.R`, `test-mplus-lmr.R` | Mplus categorical two-level LCA, the one-profile random intercept model, and TECH11 |
| `test-mclust.R` | mclust: the M-steps, parameter counts and likelihood for all 14 covariance structures, and the VVV/EEE limits of the full-covariance model |
| `test-tidysem-bch.R` | tidySEM's classification error matrix and BCH estimate |
| `test-numerical-equivalence.R` | The in-house exhaustive-enumeration likelihood (`tests/testthat/helper-independent-likelihood.R`) |
| `test-jstats-lta.R` | Latent transition references pinned by JStats: depmixS4 (Gaussian and binary), two Mplus LTA runs, and Table 5 of Muthen & Asparouhov (2022) |

## Fixtures

Each fixture is stored in exactly one place, and `equivalence_fixture()`
stops if it finds a name in neither place or in both.

- `equivalence/fixtures/` holds fixtures that only these tests read. This
  includes `jstats/`; see `jstats/PROVENANCE.md` for its source and checksums.
- `tests/fixtures/mplus/` keeps the three Mplus fixtures whose data the
  shipped unit tests also use as input.

## Relation to `validation/`

`validation/equivalence/` is a separate harness. It writes a tidy report
(`report.csv`, `REPORT.md`) of every compared quantity, and stays as it is.
This folder holds testthat assertions: each one passes or fails.
