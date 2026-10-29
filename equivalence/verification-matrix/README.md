# Verification matrix

**[VERIFICATION-MATRIX.md](VERIFICATION-MATRIX.md) is the project's living
record of how each function is verified**: read it on GitHub, or open
`VERIFICATION-MATRIX.html` / `.pdf` for the full detail. All three are
generated; rebuild them (below) whenever the evidence changes, and never edit
them by hand.

Every major function against every reference it is checked by, built from
saved evidence only: nothing is fitted. References are separated into model
comparisons, published analyses, component checks, independent reference code
and internal checks. Source strength counts complete, agreeing model comparisons.

```
Rscript equivalence/verification-matrix/build.R      # from the project root
```

## Keeping it current

Rebuild after anything that changes the evidence, and commit the regenerated
files with the evidence that changed them:

```
Rscript equivalence/verification-matrix/test-datasets.R   # only if a testthat kit's cases changed
Rscript equivalence/verification-matrix/build.R
```

| When | Also update |
|---|---|
| A kit's `compare.R` was rerun | preserve its `comparison.csv.provenance.csv` sidecar with the result; shared kit runners write it automatically |
| A new kit | `kit-evidence.csv` (path, format, kind, run date, wrapper), `kit-datasets.csv`, and `kit-functions.csv` if its rows have no `fn` |
| The harness was rerun, or equivalence tests were rerun into a new results file | that file in `.matrix_evidence_files()` in `build.R`, with its run date, package version and source/result sidecar |
| A new `test-*.R` file | `test-attribution.csv` (one rule per test when its tests compare different references) |
| A testthat kit gained or lost cases | `test-datasets.R`, then rerun it |
| A new reference found, or one wired in | add it to `candidates.csv`, or remove it there |

The build stops if any harness quantity, harness source, kit case, test file or
kit dataset is left unattributed, so a new piece of evidence cannot silently
drop out of the counts.

## Inputs

| File | What it holds |
|---|---|
| `harness-functions.csv` | Which function each harness quantity belongs to: suite plus a regex on the quantity; first match wins |
| `harness-references.csv` | Which reference each harness `source` is, and its kind |
| `test-attribution.csv` | For every `test-*.R` in `equivalence/` and `tests/equivalence/`: function, reference and kind, by a regex on the test name. `;` lists several; such tests count once per function |
| `candidates.csv` | References found but not yet wired in (see `../REFERENCE-CANDIDATES.md`) |
| `kit-evidence.csv` | Per-quantity comparison files of the kits (`latentgold*/`, `*-packages/`, `*-published/`): path, format, reference, kind, run date, and the testthat wrapper whose expectations the kit replaces (so nothing is counted twice) |
| `kit-functions.csv` | For kit files without a function column: a regex on the case name, and the function it checks |
| `kit-datasets.csv` | For every kit: a regex on the case name and the dataset it names (`sub()` replacement); first match wins |
| `test-datasets.csv` | One row per simulated dataset that a testthat file's kit sent to another program (file, function, dataset); written by `test-datasets.R` from the kits' `cases.rds` |

The saved results it reads are listed in `.matrix_evidence_files()` in
`build.R` and `kit-evidence.csv`. A result counts as current only when its
`<result>.provenance.csv` sidecar matches both the result's contents and every
runtime file in `R/` plus `NAMESPACE`. Missing sidecars, changed results,
changed source and added/deleted runtime files make evidence historical.
A documentation-only version bump does not change the runtime identity.

For new result-producing scripts, capture source identity before running and
record it only after writing the result:

```r
source("equivalence/verification-matrix/provenance.R")
snapshot <- matrix_capture_source()
# Run comparisons and write result_path here.
matrix_record_evidence(result_path, snapshot)
```

The recorder refuses source changes during a run. Shared package-kit runners
and `run-package-kit-tests.R` use this protocol. MD5 detects accidental file
changes; these manifests are not security signatures. Preserve sidecars with
their results instead of recording fresh source hashes onto old results.

The initial 2026-10-28/29 evidence was migrated with explicit retrospective
source-equivalence attestations, not represented as a fresh rerun. Their
sidecars document the comparison against retained source manifests and
commits 344dc7c/7d0e0df. The older 0.96.1 results remain historical. See
`attest-historical.py` for the migration checks and scope.

The latest current result for each test file is preferred; historical results
are shown separately and earn no source strength. Test errors invalidate the
reference, skipped tests are incomplete, and neither earns source credit.
The build stops if any harness quantity, source or test file is unattributed.

## Outputs

`evidence.csv` (one row per check, function and reference), `matrix.csv` (one
row per function and reference), `summary.csv` (one row per function), and
`VERIFICATION-MATRIX.html` / `.pdf`. The PDF is printed by headless Chrome; set
`LATENTS_CHROME` when it is not at the macOS default.

## Equivalence strength

A source is an independent program or package (another R package, Mplus,
Latent GOLD, Stata, Python) with complete, agreeing model comparisons on the
current runtime code for two or more classes. A failed or incomplete current
reference earns no source strength. Component formulas (for example tidyLPA
classification statistics on latents posteriors) and project-written EM
oracles using MASS/VGAM remain visible under their own evidence types and
do not count as independent model implementations. Agreement on one-class reductions only (one
class equals `lm()`, `glm.nb()`, `polr()`, lme4, a linear mixed model) is shown
as "(1 class)" and does not make a source. The number of testthat checks is
not used as a threshold: a single expectation can compare a whole table (the
wider Latent GOLD kit compares 172 quantities in 3 expectations).

Strength is one `+` per source and one `P` per published paper whose numbers
are reproduced: `++++ PP` is four sources and two papers, `0` none. The badge
colour marks none (grey), one source (orange), two (blue), three or more
(green).

**Datasets** counts the distinct datasets on which a function was compared
with a source or a paper, on current code, for models with two or more
classes: the harness's `dataset` column, the kits' datasets
(`kit-datasets.csv`) and the simulated datasets of the testthat kits that the
other program ran (`test-datasets.csv`). Names match case-insensitively and
without a trailing description in parentheses, so a dataset used by several
kits counts once. Tests that build small data inline are not counted: the
number is a lower bound. Component checks, independent code
(numDeriv, dense likelihoods) and internal checks (invariance, seed stability,
refusal parity) are listed but do not count. Results last run under 0.96.1
are shown but do not count either; the "to strengthen" column names what a
rerun would add.

## Reading the result

Read the generated [summary and per-function references](VERIFICATION-MATRIX.md)
for current counts. The HTML/PDF additionally show the reference-by-function
grid. Source counts, comparison rows, individual values and datasets describe
different things; none is a count of independent replications.

Known differences are shown separately and do not count as agreements.
Zero unexplained failures therefore does not mean that all compared quantities
agree. A comparison started from a reference solution verifies fitting and
reporting at that solution; consult the kit's search record before drawing
conclusions about independent default searches.

The external-check headline and per-function check totals exclude internal
checks, provenance and test infrastructure. Component checks and independent reference code are included
in that headline but kept separate from independent-model source strength.

Fast regression checks (source/result tampering, errors/skips, incomplete
comparisons and overbroad exceptions) run without model fitting:

```
Rscript equivalence/verification-matrix/tests/run.R
Rscript equivalence/run.R '^verification-attribution$'
python3 -B -m unittest discover -s equivalence/verification-matrix/tests -p 'test_*.py'
```
