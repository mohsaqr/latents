# Equivalence testing

Everything that compares multilpa with other software, published output or an
independent implementation lives in this folder. None of it ships: the folder
is listed in `.Rbuildignore`. The package's own `tests/testthat/` keeps only
unit tests and invariant checks, which run with nothing more than the package
and testthat installed. Evidence that is not a comparison, such as the
mathematical audit and the Monte Carlo recovery study, stays in
[`validation/`](../validation/README.md).

Run everything from the project root.

| Path | What it is | Run |
|---|---|---|
| `test-*.R`, `helper-*.R`, `fixtures/` | testthat assertions: each one passes or fails. Listed below. | `Rscript equivalence/run.R` |
| `harness/` | Runs every external comparison and writes one tidy report of every compared quantity (`report.csv`, `REPORT.md`). See its [README](harness/README.md). | `Rscript equivalence/harness/run.R` |
| `mplus/` | Retained Mplus runs, each beside the `compare.R` that checks the package against it. | `Rscript equivalence/mplus/<run>/compare.R` |
| `latentgold/` | The Latent GOLD 6.1 kit, its retained output and the comparison. See its [README](latentgold/README.md). | `Rscript equivalence/latentgold/compare.R` |
| `houle-2026/` | Checks against the article and supplement of Houle et al. (2026). | `Rscript equivalence/houle-2026/check-fit-indices.R` |
| `three-step/` | The BCH correction against tidySEM. | `Rscript equivalence/three-step/compare-tidysem.R` |
| `deferred.R`, `fixtures.R`, `independent-likelihood.R` | Helpers the scripts above source. | — |

Scripts other than the testthat files rewrite committed reports. To check
only that they still run, run them in a scratch copy of the repository.

## The testthat assertions

```
Rscript equivalence/run.R                  # every file
Rscript equivalence/run.R jstats mclust    # only files whose names match
```

`run.R` loads the package from source and then runs this folder's test files
with testthat. It does not descend into the subfolders. It stops with a
non-zero exit status if any test fails.

### Files

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
