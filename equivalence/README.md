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
| `test-mixture-regression-flexmix.R` | flexmix: `mixture_regression()`'s likelihood at flexmix's estimates, its optimum, binomial and Poisson estimates, shared coefficients and parameter counts, for row-level and group-level classes with and without a concomitant model |
| `test-tidysem-bch.R` | tidySEM's classification error matrix and BCH estimate |
| `test-mplus-classification.R` (`helper-mplus-classification.R`) | The package's own entropy and AvePP tables (`get_results(fit, "entropy")`, `"classification"`) against the values Mplus prints, to its three decimals: User's Guide examples 7.9 and 7.10, and both levels (`CW`, `CB`) of the two-level FIML runs. Classes are matched by Mplus's posterior-based counts |
| `test-missing-patterns.R` (`missing-patterns/`) | Five missingness patterns (MCAR cells, empty rows, MAR, one indicator 70% missing, a two-level set with empty groups) fitted by genuine Mplus 9 Demo and Latent GOLD 6.1: likelihood, sample size, BIC, entropy, AvePP and row-by-row posteriors. The programs' differing conventions for rows and groups with nothing observed are recorded in its README |
| `test-montecarlo.R`, `test-montecarlo-extended.R`, `test-enumeration-mplus.R` (`montecarlo/`, `helper-montecarlo.R`) | Monte Carlo equivalence (Equivalence 2.0): 160 simulated datasets with known parameters (simulab), each fitted by genuine Mplus 9 Demo and Latent GOLD 6.1: likelihood, parameters, N, criteria, posteriors, entropy, AvePP matrix, OCC, counts, implied means, group level, `predict()`, enumeration grids, `compare_models()` and `cores`; plus a retained global-maximum audit. See `montecarlo/README.md` |
| `test-jstats-mplus.R` (`jstats-mplus/`) | 42 genuine Mplus mixture runs kept in JStats: LPA EEI/VVI, LCA binary, polytomous, nominal and count, `c ON x`, R3STEP, a mixture regression and an occasion-specific LTA |
| `test-mplus-standardize.R`, `test-mplus-person-scaling.R` (`mplus/standardize/`, `mplus/person-scaling/`) | `scaling = "grand"` against Mplus DEFINE: STANDARDIZE (Mplus divides by the SD with denominator n, this package by n - 1: equal up to that constant), and within-person z-scores against Mplus on pre-standardized data |
| `test-mplus-mixture-weights.R` (`mplus/mixture-weights/`) | a weighted three-profile mixture against Mplus WEIGHT (MLR): pseudo likelihood, posteriors, estimates and robust standard errors (equal up to the documented CR1 factor) |
| `test-mplus-tech14.R` (`mplus/tech14-null/`) | `bootstrap_lrt()` against TECH14: statistic, parameter difference and null likelihood on three runs, and the p-value where the null is true |
| `test-latentgold.R` | the Latent GOLD kit (`latentgold/`) against the current package: targets are rebuilt from the package, not read from the kit |
| `test-latentgold-families.R` | the group-class families against Latent GOLD at its most refined quadrature |
| `test-external-scripts.R` | the comparisons kept as scripts (weights, BCH against tidySEM, Houle et al., latent-transition extensions, ordinal and count against Latent GOLD), each run in its own process; a non-zero exit or a disagreeing row fails |
| `test-claims.R` (`CLAIMS.csv`) | the registry of every documented claim of agreement with another program or a published result, each with the test that backs it: a claim reworded, or a test renamed or removed, fails; the one unbacked claim is a convention (`boundary = "fix"`) |
| `test-tidylpa-classification.R` | Relative entropy and AvePP against tidyLPA's formulas on the same posteriors (`srl`, VVI, 3 to 5 profiles). An independent tidyLPA fit is not compared, because mclust can stop at a different maximum |
| `test-numerical-equivalence.R` | The in-house exhaustive-enumeration likelihood (`tests/testthat/helper-independent-likelihood.R`) |
| `test-jstats-lta.R` | Latent transition references pinned by JStats: depmixS4 (Gaussian and binary), two Mplus LTA runs, and Table 5 of Muthen & Asparouhov (2022) |

## Fixtures

Each fixture is stored in exactly one place, and `equivalence_fixture()`
stops if it finds a name in neither place or in both.

- `equivalence/fixtures/` holds fixtures that only these tests read. This
  includes `jstats/`; see `jstats/PROVENANCE.md` for its source and checksums.
- `tests/fixtures/mplus/` keeps the three Mplus fixtures whose data the
  shipped unit tests also use as input.

## Sabotage battery

`sabotage/run.R` injects small realistic defects into the loaded package (in
memory; no source file changes): a wrong entropy denominator, empty rows or
groups counted, a transposed AvePP matrix, one parameter too many, a scaled
Gaussian density, inflated M-step variances, a wrong mixing update, dropped
weights, a wrong scaling denominator and a wrong bootstrap p-value. For each,
the tests named for it are run, and the defect counts as detected when one
expectation fails. A control run (the target replaced by itself) must come
out undetected, or the battery stops. Results in `sabotage/results.csv`.

Run of 2026-10-18: the control passed (0 of 112 failed) and all 11 defects
were detected. Expectations failing per defect: entropy denominator 151,
empty rows counted 109, empty groups counted 1 (the m5 group level in
`test-missing-patterns.R`, the only reference with an empty group), AvePP
transposed 119, one parameter too many 1,274, Gaussian density 828, M-step
variances 56 (and 6 errors), mixing update 10, weights dropped 6, scaling
denominator 2, bootstrap p-value without +1 1 (a unit test: against TECH14
the p-values agree only within Monte Carlo error, which this defect stays
inside). Each filter names tests that run in this process; the scripts run by
`test-external-scripts.R` start their own R and never see a defect.
