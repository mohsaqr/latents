# Verification matrix

**[VERIFICATION-MATRIX.md](VERIFICATION-MATRIX.md) is the project's living
record of how each function is verified**: read it on GitHub, or open
`VERIFICATION-MATRIX.html` / `.pdf` for the full detail. All three are
generated; rebuild them (below) whenever the evidence changes, and never edit
them by hand.

Every major function against every reference it is checked by, built from
saved evidence only: nothing is fitted. Modelled on the JStats verification
system: one status per function, references by kind, and a 2-of-3 trust policy.

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
| A kit's `compare.R` was rerun | nothing (its `comparison.csv` is read directly) |
| A new kit | `kit-evidence.csv` (path, format, kind, run date, wrapper), `kit-datasets.csv`, and `kit-functions.csv` if its rows have no `fn` |
| The harness was rerun, or equivalence tests were rerun into a new results file | that file in `.matrix_evidence_files()` in `build.R`, with its run date and package version |
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
`build.R` (the harness report, the 2026-10-28 testthat results, the
2026-10-29 rerun of every equivalence file that had no current result, and
the 0.96.1 results as a last resort) and in `kit-evidence.csv`. For each test file the current (2026-10-28 or 2026-10-29) result is used where it exists;
otherwise the 0.96.1 result (2026-10-20), shown in grey and never counted
toward status. The build stops if any harness quantity, harness source or test
file on disk is left unattributed.

## Outputs

`evidence.csv` (one row per check, function and reference), `matrix.csv` (one
row per function and reference), `summary.csv` (one row per function), and
`VERIFICATION-MATRIX.html` / `.pdf`. The PDF is printed by headless Chrome; set
`LATENTS_CHROME` when it is not at the macOS default.

## Equivalence strength

A source is an independent program or package (another R package, Mplus,
Latent GOLD, Stata, Python) whose results agree on the current code for a
model with two or more classes. Agreement on one-class reductions only (one
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
number is a lower bound. Independent code
(numDeriv, dense likelihoods) and internal checks (invariance, seed stability,
refusal parity) are listed but do not count. Results last run under 0.96.1
are shown but do not count either; the "to strengthen" column names what a
rerun would add.

## Result (2026-10-29)

Every function has three or more sources and no failed current check:

| function | strength | sources | papers | datasets | checks |
|---|---|---|---|---|---|
| `lpa()` | ++++++++ PPP | 8 | 3 | 163 | 8,049 |
| `multilpa()` | +++++++ | 7 | 0 | 91 | 7,367 |
| `lca()` | +++++++++ PPPPPPP | 9 | 7 | 161 | 13,072 |
| `multilca()` | ++++ P | 4 | 1 | 25 | 5,376 |
| `lta()` | +++++ PPP | 5 | 3 | 31 | 2,456 |
| `mixture_regression()` | +++++ PPP | 5 | 3 | 21 | 719 |
| growth (`random =`) | +++ P | 3 | 1 | 14 | 827 |
| `three_step()` | +++++ | 5 | 0 | 77 | 2,350 |
| `bootstrap_lrt()` / VLMR | +++ | 3 | 0 | 63 | 5,485 |
| enumeration | ++++ | 4 | 0 | 10 | 245 |
| `pool_imputations()` | +++ | 3 | 0 | 2 | 136 |

Attribution fixes in this build: `test-review-1026-oracles.R` and
`test-features-ml-three-step.R` are attributed per test (a file-wide rule had
credited nlme, poLCA and sandwich to tests that compare with depmixS4 or Latent
GOLD, and "sandwich" in the ML three-step file is the estimator, not the
package), and the first test of `tests/equivalence/test-cluster-robust-cr1.R`
(the sandwich arithmetic on `lm()`) is independent code, not a mixture source.
