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
| `latentgold-coverage/` | Growth mixtures, `lca()`, `lpa()`, `lta()`, negative-binomial and ordinal mixture regression, including individual `lca()`/`lpa()` fits at one to four classes, against Latent GOLD 6.1 (retained listings, offline). | `Rscript equivalence/latentgold-coverage/compare.R` |
| `lca-packages/`, `lpa-packages/`, `regression-packages/`, `lta-packages/`, `pool-packages/` | Each function against independent R implementations (glca, randomLCA, BayesLCA, multilevLCA, depmixS4; Rmixmod, EMCluster, mixtools; mixtools, MoEClust, flexmix; LMest, seqHMM; mitml, mitools), plus published values where they exist. Shared runner: `package-kit-runner.R`. | `Rscript equivalence/<kit>/compare.R` |
| `enumeration-packages/` | `enumerate_lpa()` against `mclustBIC()`, `enumerate_lca()` against poLCA, `bootstrap_lrt()` against `mclustBootstrapLRT()` and mclust's vignette (mclust run to 1e-12; its default 1e-5 stops short). | `Rscript equivalence/enumeration-packages/compare.R` |
| `growth-published/`, `stata-published/`, `mplus-growth/` | Published analyses reproduced (Proust-Lima et al. 2017 on paquid; Stata [SEM] Examples 52g and 50g and Masyn 2013 Table 25.12) and genuine Mplus 9 growth runs from JStats. | `Rscript equivalence/<kit>/compare.R` |
| `argument-sweep/` | Every argument of the fitting and post-fit functions at each of its options, on several datasets: the likelihood rebuilt from the result tables (mixture, hidden Markov, mixture of regressions, random effects), mclust, poLCA and flexmix continued from the package's solution, one-class limits (`lm`, `glm`, `glm.nb`, `polr`, `lmer`), `optim` and `numDeriv` for the step-three and Wald results, invariances and classed refusals. See its [README](argument-sweep/README.md). Writes `tmp/argument-sweep/`. | `Rscript equivalence/argument-sweep/run.R` |
| `cvs/` | Reports every run of the testthat equivalence files, the package's own tests, the harness and the argument sweep to CVS, the validation service in `../validation` (one JSON file per run; see [Reporting to CVS](#reporting-to-cvs)). | used by the runners; `Rscript equivalence/cvs/test-local.R` for the package's tests |
| `REFERENCE-CANDIDATES.md` | Published and software references not yet wired in, by function, with data access and the numbers to match. | — |
| `verification/` | CVS manifests: one `modules/<fn>.yaml` per verified function and `references.bib`, generated from the matrix for CVS's Verification tab (see [Reporting to CVS](#reporting-to-cvs)). | `Rscript equivalence/verification/build-manifests.R` |
| `verification-matrix/` | **The living verification record**: every major function against every reference that checks it (strength as one + per source and one P per paper, datasets, checks), built from saved results. Read [`VERIFICATION-MATRIX.md`](verification-matrix/VERIFICATION-MATRIX.md); rebuild after any equivalence run ([how](verification-matrix/README.md#keeping-it-current)). | `Rscript equivalence/verification-matrix/build.R` |

Scripts other than the testthat files rewrite committed reports. To check
only that they still run, run them in a scratch copy of the repository.

The ten new comparison kits, and the suites `lrt`,
`latentgold-categorical-distal`, `latentgold-step-one` and
`moeclust-covariate-structures`, also run as separate CI jobs. Run one locally
with `Rscript equivalence/run-package-kit-tests.R pool-packages` (substitute
the kit or suite name). It writes fresh test outcomes and their source
identity under `tmp/equivalence-kits/` and fails on skipped tests, errors or
disagreements.
The fast validation regressions run with
`Rscript equivalence/verification-matrix/tests/run.R`; they inject missing
results, incorrect counts, skipped tests and changed source files without
refitting models.

## Reporting to CVS

`equivalence/run.R`, `harness/run.R`, `argument-sweep/run.R` and
`cvs/test-local.R` (the ci job `test`) each end by writing one CVS report,
failing runs included: a JSON file in the vitest layout CVS ingests, named
`latents-<run>-<time>.json` (the prefix assigns it to the CVS target
`latents`). One entry per source file, one assertion per item:

| Run | Entry (file) | Assertion | Status |
|---|---|---|---|
| testthat (`equivalence/`, `tests/testthat/`) | the test file | one `test_that()` block, `checkKind` "testthat" | failed when an expectation failed or errored, skipped when it skipped, else passed |
| harness | `harness/suites/<suite>.R` | one compared quantity, with difference, tolerance, reference and output | passed when it agrees, else failed (also when a known open defect explains it); a suite that did not run is one "(suite ran)" assertion |
| argument sweep | `argument-sweep/modules/<module>.R` | one check, with difference, tolerance and the reference the run recorded for it | agrees: passed; disagrees or error: failed; known difference, no reference: skipped with the reason |

Each file declares its CVS tier (`_cvs.tier`), since CVS otherwise
classifies by file name and an R file name carries none: `equiv` for the
equivalence files, harness suites and sweep modules, `unit` for
`tests/testthat`. CVS honours the declaration from validation PR #23 on.

Each assertion also names the reference it is evidence for (`_cvs` `axis`,
`authority`, `reference`), by the verification matrix's own rules
(`verification-matrix/attribution.R`, shared with `build.R`): R package →
`r-equivalence`, Latent GOLD / Mplus / Stata → `latentgold` / `mplus` /
`stata`, published → `paper`, independent code → `independent`, and internal,
component, provenance and infrastructure checks → `internal` (shown, never
counted). Sweep checks carry the reference their run recorded in
`results.csv` (`fn`, `reference`, `kind`, from `verification-matrix/sweep-functions.csv`
and `sweep-references.csv`; the sweep stops with `latents_sweep_unattributed`
on a check no rule names). Every test block in this folder checks one
reference axis: a comparison with several programs is one block per program,
and a package kit is one block per reference group (`kit_reference_tests()` in
`package-kit-runner.R`); `test-cvs-report.R` fails on a rule that names
references on two axes. A block whose references still fall on several axes
keeps `reference` but gets no axis, and the package's own tests name no
reference and stay untagged. A commercial
program the map does not know stops the reporter
(`latents_cvs_unknown_program`); the reporter's tests check every program in
the rule tables first. `cvs_retag_report()` tags a report saved before
tagging, keeping its name, run, start and version.

The manifests CVS's per-target Verification tab reads
(`../validation/docs/TARGET-VERIFICATION.md`) are generated from the matrix:
`Rscript equivalence/verification/build-manifests.R` writes
`verification/modules/<fn>.yaml` and `verification/references.bib`; rebuild
them after the matrix (`test-verification-manifests.R` fails when they are
stale). A manifest never claims a block that has no single axis.

The report goes to the CVS inbox when reporting is on (`CVS_INBOX_DIR`, else
`../validation/data/inbox`; `CVS_DISABLED=1` turns it off) and a copy always
goes to `tmp/cvs/`. On ci, mini-ci turns reporting on for the named jobs of a
`.mini-ci` that says `report cvs`, and fetches `tmp/cvs` back otherwise. A
filtered run (`equivalence/run.R mclust`, a shard) is labelled by its filter.
The reporter is tested by `test-cvs-report.R`.

### Evidence for the matrix

Each run also saves its results with a provenance sidecar
(`verification-matrix/provenance.R`: the runtime source captured before the
run, bound to the result file after it), which the verification matrix needs
to count a result as current:

| Run | Result |
|---|---|
| `equivalence/run.R` | `tmp/equivalence-tests/<label>-results.csv` (one row per test) |
| `equivalence/run.R --dir=tests/equivalence` | `tmp/equivalence-tests/tests-equivalence-results.csv` |
| a package kit's or a Latent GOLD kit's testthat file | `tmp/kits/<kit>/comparison.csv` (one row per compared quantity, from the current package; `kit_save_comparison()` in `package-kit-runner.R`) |
| `harness/run.R` | `report.csv` in its output folder |
| `argument-sweep/run.R` | `tmp/argument-sweep/results.csv` |

`.mini-ci` fetches all of them. To refresh the matrix and the CVS manifests,
run the jobs `harness`, `sweep`, `equivalence-a`, `-b`, `-c` and `-t`, copy the
fetched results into a dated `validation/matrix-refresh-<date>/` folder (and
each kit's `comparison.csv` with its sidecar into its kit folder, setting its
`run` and `code` in `verification-matrix/kit-evidence.csv`), point
`.matrix_evidence_files()` in `verification-matrix/build.R` at them, run
`build.R` and `verification/build-manifests.R`, then the job `equivalence-v`
(the verification files check the rebuilt matrix and manifests) and copy the
manifests to the CVS server (`../validation/docs/TARGET-VERIFICATION.md`).

## The testthat assertions

```
Rscript equivalence/run.R                  # every file
Rscript equivalence/run.R jstats mclust    # only files whose names match
Rscript equivalence/run.R --dir=tests/equivalence   # the tests/equivalence/ files
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
| `test-mplus-standardize.R`, `test-mplus-person-scaling.R` (`mplus/standardize/`, `mplus/person-scaling/`) | `scaling = "grand"` against Mplus DEFINE: STANDARDIZE, both using the population SD with denominator n and matching likelihoods directly; within-person z-scores against Mplus on pre-standardized data |
| `test-mplus-mixture-weights.R` (`mplus/mixture-weights/`) | a weighted three-profile mixture against Mplus WEIGHT (MLR): pseudo likelihood, posteriors, estimates and robust standard errors (equal up to the documented CR1 factor) |
| `test-mplus-tech14.R` (`mplus/tech14-null/`) | `bootstrap_lrt()` against TECH14: statistic, parameter difference and null likelihood on three runs, and the p-value where the null is true |
| `test-latentgold.R` | the Latent GOLD kit (`latentgold/`) against the current package: targets are rebuilt from the package, not read from the kit |
| `test-latentgold-families.R` | the group-class families against Latent GOLD at its most refined quadrature |
| `test-external-scripts.R` | the comparisons kept as scripts (weights, BCH against tidySEM, Houle et al., latent-transition extensions, ordinal and count against Latent GOLD), each run in its own process; a non-zero exit or a disagreeing row fails |
| `test-claims.R` (`CLAIMS.csv`) | documented behavior and conventions checked against software or published references. Each row preserves a documentation anchor, its reference and the backing test, including when help prose omits the reference name. Changed anchors or renamed tests fail; the one unbacked row explains the `boundary = "fix"` convention. This checks registered text and test presence; it does not discover every claim or prove that a named test is sufficient. |
| `test-boundary-polca.R` | The categorical boundary flag (`get_candidates()`'s `boundary`, the trust view) against poLCA on five binary items, with and without items that copy the class: which candidates reach a probability bound (poLCA at the package's 1e-6 floor), likelihood to 1e-6, parameter counts, and the interior two-class estimates to 1e-6 |
| `test-cvs-report.R` | The CVS reporter (`cvs/cvs-report.R`): statuses of real testthat results, harness and sweep rows, inbox and local copy, totals equal to the assertions in the file; the reference tags (axis per kind of reference, untagged when several axes or none) |
| `test-verification-manifests.R` | The CVS manifests (`verification/`): current with the matrix, valid YAML with CVS's required fields, claimed files are report entry names, needles claim exactly as CVS matches |
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
