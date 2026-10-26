# The equivalence harness

This harness runs every external comparison and writes one tidy table. It
sits in `equivalence/`, which is excluded from the build (`.Rbuildignore`)
but kept under version control, because the artifacts are the evidence: a
claim that a number matches Mplus is worth nothing without the Mplus output
beside it. Run it from the project root:

```
Rscript equivalence/harness/run.R              # every suite
Rscript equivalence/harness/run.R glca mclust  # named suites only
```

It writes six generated files, so no number in the documentation is
transcribed by hand:

| File | Contents |
|---|---|
| `report.csv` | One row per compared quantity: reference value, the value multilpa produced, the difference, the tolerance, and the multilpa version the row was produced under. |
| `suites.csv` | One row per suite: whether it ran, was skipped or failed, with the reason, how many warnings it raised, and how long it took. |
| `coverage.csv` | One row per output column of every fit a suite compared: how many values the fit reported, how many were compared, and, for anything not compared, the reason. |
| `known-open-status.csv` | The register of open defects (`known-open.csv`), with how many disagreements each explains in this run. |
| `REPORT.md` | The summary tables, output coverage and known open defects. |
| `SESSION.txt` | The full `sessionInfo()` the run was produced under. |

The run exits non-zero when a suite fails, a quantity disagrees without a
registered open defect to explain it, an output is neither compared nor
exempted, or a register entry no longer matches anything.

## Every output, not one or two measures

Until 2026-10-26 each suite chose by hand what to compare, and the choices were
narrow: the depmixS4 suite compared the first indicator's means and the
transition matrix, mclust the means and mixing proportions, poLCA ran without
standard errors. No suite compared posteriors, classifications, the
classification error matrix or the BCH weights. That is how a defect in the
error matrix (rows with nothing observed were counted; fixed 2026-10-26) went
unnoticed. `outputs.R` makes the comparison complete and its completeness visible.

- **The inventory comes from the package.** `fit_outputs()` reads every table
  of `get_results(fit, "all")` and lists each value column with its number of
  cells. A column added to the package later appears in the inventory without
  anyone listing it.
- **Derived tables are rebuilt from the reference.** `harness_reference()`
  holds another program's estimates and posteriors in one shape.
  `reference_tables()` rebuilds every table from them, writing each derived
  table from its definition in harness code: counts, classification, average
  posteriors, odds of correct classification, the error matrix, BCH weights,
  entropy, all sixteen information criteria, bivariate residuals and Gaussian
  diagnostics. Every comparison of a derived table therefore tests the
  package's posteriors and the package's formula at once. ICL, AWE and CLC
  charge the complete-data entropy, `H(group classes) + sum_i E_C[H(X_i | C,
  y_i)]`, under both conventions (a transition model keeps each level's own);
  a two-level reference supplies the rows' within-class entropies as
  `harness_reference(row_conditional_entropy = )`, computed from its own
  estimates (`two_level_row_conditional_entropy()` in
  `two-level-lca-reference.R`). A reference with several group classes that
  lacks them has the three criteria exempted by name, never filled in.
- **One row per column, every cell checked.** `compare_tables()` joins on key
  columns and `compare_block()` compares every cell of a column. Missing values
  must be missing in both, and infinite values infinite with the same sign. The
  row shows the worst cell and `n_values`, so ten thousand posteriors take one
  row of `report.csv` and are still all checked.
- **Nothing is left out without a reason.** `output_ledger()` gives every
  output column one status:
  - `compared`: all its values were compared;
  - `partial`: some were;
  - `exempt`: a reason is recorded, either a default one for things that are
    not estimates (`data`, the random-start trace, iteration counts) or the
    suite's own;
  - `not covered`: neither, which fails the run.

  An exemption for an output the fit does not report, or a comparison naming
  one, is refused, so the ledger cannot go stale.
- **Standard errors have an independent reference.** Where the reference
  program reports none, or reports a different estimator, the suite takes a
  numerical Hessian (numDeriv) of a likelihood written from scratch in base R,
  at the reference's estimates, in the natural free parameters. The mclust
  suite is the worked example.
- **References are refined, not tolerances loosened.** mclust and depmixS4 stop
  on a relative change in the log-likelihood. Where the surface is flat, that
  leaves their estimates short of the fixed point. Each suite continues the
  program's own EM to the fixed point before comparing. A tolerance that is
  wider than the default says why next to it: OCC, for example, divides by
  1 - AvePP.

Two suites ask a different question from "does another program agree":

| Suite | What it checks |
|---|---|
| `invariance` | A transformation whose effect on every output is known exactly is applied, the model refitted, and every output table compared with the original's carried through by hand: rows permuted; groups, categories and indicators renamed to labels containing `: , [ ]`; y -> a + c y at moderate and extreme c (the variance floor carried into the new units); every unit duplicated (estimates fixed, SEs / sqrt 2, likelihood and counts doubled); zero-weight units added; a predictor's origin moved by 2000. |
| `refusals` | One degenerate input through every entry point (`multilpa`, `lpa`, `lca`, `lta`, the enumerations, `mixture_regression`) must give one classed answer. An unclassed error, another class or a silent fit fails. |

### Known open defects

`known-open.csv` lists defects the harness finds that await a decision, each
with a suite, a pattern on the quantity, and the decision needed. A matching
disagreement is reported under its id in `REPORT.md` instead of failing the
run, so an open item does not turn every run red and hide the next new
failure. An entry that no longer explains any disagreement fails the run: the
defect was fixed, or the check that found it is gone, and the entry has to be
removed or corrected.

### What the harness refuses to let you miss

The external review of 2026-09-20 found that the saved report of 1,012
agreements could not be reproduced: multilpa 0.8.0 had renamed `indicators` to
`vars` and `group` to `id`, no validation script had been updated, and a fresh
run stopped at suite load with a bare `unused arguments` naming neither file
nor line. The saved report kept being read as evidence for 0.10.0, a version it
had never been run against. Five things now stand between that and a reader.

1. **A stale call cannot reach a suite.** `check_validation_api()` parses every
   `.R` file under `equivalence/` and `validation/` and checks each call to an exported multilpa
   verb against that verb's current formals, before anything is fitted. An
   argument that no longer exists is reported with its file, its line and its
   replacement from the 0.8.0 migration table. `run.R` runs this first and
   stops if it finds anything.
2. **A missing external package is a skip with a reason.** A suite calls
   `require_suite_packages()`, which raises `multilpa_suite_skipped`; the
   harness records the suite as skipped, says which package was absent and why
   the suite needed it, and carries on with the rest. `require_suite_files()`
   does the same for retained program output. Nothing disappears silently and
   one absent dependency does not abort a run.
3. **A report carries its own provenance.** The multilpa version, the R
   version, the platform and every external package version are recorded by the
   harness, not typed in. The version and an md5 fingerprint of `R/` are taken
   both before the first suite and after the last: if the package was edited
   while the run was in progress, the report says so in bold rather than
   quietly claiming the closing version.
4. **Every tolerance is explicit.** `compare_values()` has no default
   tolerance, so a comparison cannot inherit one by accident, and the pinned
   values live in one place, `.equivalence_tolerances`. The summary table
   reports the loosest tolerance claimed by each source alongside the worst
   difference seen.
5. **A single-seed agreement is not reported as a finding.** The `stability`
   suite refits three of the models the other suites rely on -- one
   single-level categorical, one Gaussian, one two-level -- under six seeds
   each, and reports the spread of the maximised log-likelihood and the number
   of distinct optima reached. If a comparison's agreement depends on its seed,
   that shows up here rather than being invisible.

### Read the two kinds of difference separately

Every row carries a `precision` column, and the two kinds must not be pooled
into one headline number.

`machine` means the reference was computed here, in this session, at full double
precision: another package's fitted object. Nothing bounds the agreement except
the two implementations, so these differences are the real equivalence
evidence and they run at 1e-6 or below.

`printed` means the reference was read from published digits -- a journal table,
an Mplus `.out`, a package vignette. Rounding limits the available precision,
but a discrepancy can also reflect optimization or model differences. Rows
built by `compare_printed()` express their gap in half-unit rounding bounds;
other rows, including the Mplus comparisons, use explicit absolute or relative
tolerances. Reference origin and difference scale are separate concepts.
The generated summary compares the gap divided by its stated tolerance, so
one is the acceptance boundary regardless of the row's scale.

Pooling the two is how a table ends up reporting 2e-02 for mclust. That figure
was the mclust vignette's BIC of -2314.316, which mclust reproduces at its own
default tolerance to 3e-04 and which sits 0.021 below the converged optimum
because its EM stops early. Against a converged mclust, multilpa agrees to
3e-07. One number was hiding the other.

**The sibling directories** under `equivalence/` hold retained program output
— genuine Mplus and Latent GOLD runs and the fixtures built from them —
together with a `compare.R` that checks the package against each.

## What is being compared, and why each source is here

| Suite | Reference | What only this source establishes |
|---|---|---|
| `published` | Linzer & Lewis (2011), Agresti & Lang (1993), Zhou & Lange (2010), Goodman (2002), Dayton (1998), McCutcheon (1987) | Agreement with numbers the field already accepted, rather than with another implementation. Includes a published *pair* of likelihood modes, which tests the optimizer rather than the likelihood. |
| `mplus-twolevel-lca` | Mplus 8.8, User's Guide example 10.7 | The only public, fully specified two-level model with latent classes at **both** levels and printed output for all 99 parameters. Mplus stops on a log-likelihood change of 0.001 (TECH8 shows it still gaining 0.00096 per step), so printed values are compared at Mplus's own estimates, with half-unit plus estimate-rounding bounds, and the converged fit against those estimates carried to the maximum by Newton steps on a from-scratch likelihood. Mplus's MLR errors check that likelihood's sandwich; the package's robust errors add the CR1 factor sqrt(G/(G-1)). |
| `glca` | glca 1.4.2 | An independent implementation of the same Vermunt (2003) nonparametric model, on real data **with missing responses**. glca's own (outer-product) errors are not a reference: its scores are wrong in groups with a missing response (summing to 2.69 at its maximum), and on gss08 its generalised inverse returns zeros. Errors are checked against a from-scratch numDeriv information instead. |
| `multilevlca` | multilevLCA 2.1.6 and Lyrvall et al. (2025) | A third independent two-level implementation, plus peer-reviewed published values on CRAN-bundled public data. Fixes the correspondence between its `BIClow`/`BIChigh` and multilpa's two sample-size conventions. Its errors are an individual-level outer product with no package counterpart, so they are not compared. Its reference is run 2,000 EM steps past its own stopping rule, which halts on a likelihood ridge while the group probabilities still move. |
| `polca` | poLCA 1.6.0.1 | The single-level categorical limit across five datasets, and a likelihood written here from the definition as a third opinion. poLCA's errors are outer-product (BHHH, read from `poLCA.se`) and are compared with `vcov_type = "opg"`; the observed-information column is compared with a numDeriv Hessian of the textbook likelihood in baseline-category logits. |
| `mclust` | mclust 6.1.2 and Scrucca et al. (2023) | The Gaussian limit across the four basic covariance structures in this harness. The parent folder's `test-mclust.R` and `test-structure-inference.R` cover all 14 fitted structures. |
| `depmixs4` | depmixS4 1.5.1 | Latent transitions: forward-backward posteriors, two-slice transition counts and per-sequence likelihoods. Standard errors use depmixS4's constrained inverse with numDeriv in place of its `fdHess`, which is 1-15% off on these fits (0.287 against 0.072 for a mean near zero). |
| `invariance`, `refusals` | The package against itself under known transformations, and its entry points against each other | Equivariance and refusal parity, which no external program can establish. |
| `hard-mclust`, `hard-polca`, `hard-depmixs4` | The same three references on data built to be difficult | Each case makes two claims. (A) Started at the reference's mode, the package reproduces every output there. (B) The package's own search reaches at least the reference's maximum. The fixtures cover near-empty and rare classes, heavy overlap, r = 0.99, variance ratios of 100, a 1e6 offset, units of 1e-4 and 1e5, 40 rows for 62 parameters, ties, sparse and floored categories, 6 to 8 categories, near-saturated models, reversed categories, near-absorbing and structurally zero transitions, two-occasion and unequal sequences, and an extra state. On hard data the references needed more care than the package: mclust's and depmixS4's EM stop short of their fixed points; the numerical Hessian's 10% first step leaves the parameter space; mclust refuses a covariance with condition number 1.7e18; depmixS4 hits `maxit` without a warning. In five cases the package found a higher maximum than the reference, which the reference's own EM then confirmed. |

## Three claims that a converged-to-converged comparison cannot separate

Comparing two fitted models conflates agreement of the likelihood, agreement of
the optimizer, and agreement of the parameterization. The suites separate them:

1. **The likelihood function.** `multilpa(..., max_iter = 0)` performs no update
   and returns the model evaluated at `start`, so `logLik()` scores a parameter
   set produced elsewhere. The Mplus and poLCA suites both use this to evaluate
   the other program's own published estimates.
2. **The optimizer.** Zhou & Lange (2010) publish the dominant mode of the
   carcinoma four-class likelihood *and* an inferior local mode. The `published`
   suite runs a deterministic and a random start under each of 120 seeds,
   with plain EM, retains each start's likelihood and checks that both appear.
3. **The parameterization.** Free-parameter counts are compared everywhere. The
   Mplus example 10.7 count of 99 decomposes as 4 x 10 x 2 thresholds, 4
   group-class logits and 5 x 3 prevalence logits, which pins the model
   structure independently of any fitted value.

## Deliberate non-equivalences

A validation suite that only ever agrees is not measuring anything. These are
recorded because they are real and understood.

- **poLCA under missing data is not an equivalence target.** Its `na.rm = FALSE`
  path does not report the observed-data log-likelihood. On `election` it
  returns -10848.839695 where both multilpa and a likelihood written here from
  the definition give -10834.784505 at poLCA's own estimates, and its reported
  class shares sum to 0.99216 rather than one. The `polca` suite therefore
  compares complete cases only; missing-data equivalence is established against
  glca instead, which agrees to 6e-08.
- **Mplus example 10.6 is not used, deliberately.** It is the parametric
  counterpart of 10.7, fitting a *continuous* between-level factor
  (`f BY c#1 c#2; f ON w`) rather than between-level latent classes, and its
  published log-likelihood of -3888.076 is a value multilpa's nonparametric
  model should not reproduce. It is not included as a negative control because
  it would not be a clean one: the model also puts a within-level covariate on
  class membership (`c ON x`), which multilpa cannot combine with a between-level
  factor, so any gap would confound the parametric-versus-nonparametric
  difference with the omitted covariates. `fit_random_intercept()`, the
  parametric alternative in this package, takes continuous indicators only and
  has no covariate arguments, so it cannot fit this example either. **The
  continuous-indicator two-level model consequently has no published external
  reference**; its evidence is the retained Mplus runs under `mplus/` and the
  Monte Carlo recovery study, not a citation.
- **Near-saturated models are compared on fit, not on parameters.** `cheating`
  with three classes has 15 degrees of freedom against 14 free parameters. poLCA
  and multilpa agree on its log-likelihood to 3e-08 while reporting class shares
  of 0.035/0.111/0.854 and 0.045/0.078/0.877, because the likelihood has a flat
  ridge there. The suite compares parameters only where at least three residual
  degrees of freedom remain.
- **mclust's `diabetes` data are excluded.** mclust's own help page flags them as
  flawed, and three mutually inconsistent published analyses exist.

## Provenance of the retained artifacts

| Path | Origin |
|---|---|
| `equivalence/mplus/twolevel-lca/` | `https://www.statmodel.com/usersguide/chap10/chapter10.zip`, sha256 `ae29930be639669311d8652af790c531a861d8aa5414783766cd175f7c32b1d0`. Contains example 10.7's input, data and Mplus 8.8 output, gzipped. |
| `equivalence/mplus/*/` (other) | Runs of Mplus retained with their `.out`, checked for normal termination and a replicated best log-likelihood before any value is read. |
| `equivalence/latentgold/` | Latent GOLD 6.1 comparison, **run 2026-09-22: 372 quantities compared, 372 agree, 0 disagree** across twelve cases (two-level continuous, full covariance, categorical, mixed, FIML, covariates, bivariate residuals, three-step, latent transitions, and the bundled `course_engagement` data the README reports). Likelihoods, parameter counts, posteriors at both levels, means and variances recomputed from Latent GOLD's posteriors, and the information criteria in both sample-size conventions. Latent GOLD's own listings and posterior files are retained in `equivalence/latentgold/returned/` with `PROVENANCE.md` and an input fingerprint, so `compare.R` runs offline and Latent GOLD is not needed again unless the kit changes. Standard errors, the listing's own estimates, the ranking of bivariate residuals and the Step-3 results are read from the listing's tables and agree too; seven further rows are shown as `not comparable`, each with its reason. |
| `equivalence/deferred.R`, `equivalence/fixtures.R` | Shared by the Mplus scripts: `deferred_verb()` loads a feature moved to `future/` so the comparison that would restore it still runs; `mplus_fixture()` resolves a fixture to `tests/fixtures/` or `equivalence/fixtures/` by who reads it. |
| `equivalence/three-step/` | The BCH correction against tidySEM's independent construction of the same weights. Its README also records the Monte Carlo checks of the three-step methods. |
| `validation/simulab-recovery/` | Not equivalence evidence and not moved: a Monte Carlo recovery study, which measures bias and coverage rather than agreement with a fixed number. |
