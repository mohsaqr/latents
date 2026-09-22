# Validation evidence for multilpa

Everything in this directory is evidence that the package computes what it
claims to compute. It is excluded from the build (`.Rbuildignore`) but kept
under version control, because the artifacts are the evidence: a claim that a
number matches Mplus is worth nothing without the Mplus output beside it.

Two kinds of material live here.

**`equivalence/`** is a harness that runs every external comparison and writes
one tidy table. Run it from the project root:

```
Rscript validation/equivalence/run.R              # every suite
Rscript validation/equivalence/run.R glca mclust  # named suites only
```

It writes four generated files, so no number in the documentation is
transcribed by hand:

| File | Contents |
|---|---|
| `equivalence/report.csv` | One row per compared quantity: reference value, the value multilpa produced, the difference, the tolerance, and the multilpa version the row was produced under. |
| `equivalence/suites.csv` | One row per suite: whether it ran, was skipped or failed, with the reason, how many warnings it raised, and how long it took. |
| `equivalence/REPORT.md` | The summary table. |
| `equivalence/SESSION.txt` | The full `sessionInfo()` the run was produced under. |

The run exits non-zero when any quantity disagrees or any suite fails.

### What the harness refuses to let you miss

The external review of 2026-09-20 found that the saved report of 1,012
agreements could not be reproduced: multilpa 0.8.0 had renamed `indicators` to
`vars` and `group` to `id`, no validation script had been updated, and a fresh
run stopped at suite load with a bare `unused arguments` naming neither file
nor line. The saved report kept being read as evidence for 0.10.0, a version it
had never been run against. Five things now stand between that and a reader.

1. **A stale call cannot reach a suite.** `check_validation_api()` parses every
   `.R` file under `validation/` and checks each call to an exported multilpa
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
an Mplus `.out`, a package vignette. The printed precision bounds how closely
*any* implementation can agree, so these differences measure the reference, not
the package. The largest of them is Agresti & Lang's G-squared, printed to one
decimal, where a gap of up to 0.05 is arithmetically unavoidable.

Pooling the two is how a table ends up reporting 2e-02 for mclust. That figure
was the mclust vignette's BIC of -2314.316, which mclust reproduces at its own
default tolerance to 3e-04 and which sits 0.021 below the converged optimum
because its EM stops early. Against a converged mclust, multilpa agrees to
3e-07. One number was hiding the other.

**The remaining directories** hold retained program output — genuine Mplus
runs and the fixtures built from them — together with a `compare.R` that checks
the package against each.

## What is being compared, and why each source is here

| Suite | Reference | What only this source establishes |
|---|---|---|
| `published` | Linzer & Lewis (2011), Agresti & Lang (1993), Zhou & Lange (2010), Goodman (2002), Dayton (1998), McCutcheon (1987) | Agreement with numbers the field already accepted, rather than with another implementation. Includes a published *pair* of likelihood modes, which tests the optimizer rather than the likelihood. |
| `mplus-twolevel-lca` | Mplus 8.8, User's Guide example 10.7 | The only public, fully specified two-level model with latent classes at **both** levels and printed output for all 99 parameters. |
| `glca` | glca 1.4.2 | An independent implementation of the same Vermunt (2003) nonparametric model, on real data **with missing responses**. |
| `multilevlca` | multilevLCA 2.1.6 and Lyrvall et al. (2025) | A third independent two-level implementation, plus peer-reviewed published values on CRAN-bundled public data. Fixes the correspondence between its `BIClow`/`BIChigh` and multilpa's two sample-size conventions. |
| `polca` | poLCA 1.6.0.1 | The single-level categorical limit across five datasets, and a likelihood written here from the definition as a third opinion. |
| `mclust` | mclust 6.1.2 and Scrucca et al. (2023) | The Gaussian limit, across all four covariance parameterizations multilpa can express. |

## Three claims that a converged-to-converged comparison cannot separate

Comparing two fitted models conflates agreement of the likelihood, agreement of
the optimizer, and agreement of the parameterization. The suites separate them:

1. **The likelihood function.** `multilpa(..., max_iter = 0)` performs no update
   and returns the model evaluated at `start`, so `logLik()` scores a parameter
   set produced elsewhere. The Mplus and poLCA suites both use this to evaluate
   the other program's own published estimates.
2. **The optimizer.** Zhou & Lange (2010) publish the dominant mode of the
   carcinoma four-class likelihood *and* an inferior local mode. The `published`
   suite runs 120 single starts and checks that both appear.
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
| `mplus/twolevel-lca/` | `https://www.statmodel.com/usersguide/chap10/chapter10.zip`, sha256 `ae29930be639669311d8652af790c531a861d8aa5414783766cd175f7c32b1d0`. Contains example 10.7's input, data and Mplus 8.8 output, gzipped. |
| `mplus/*/` (other) | Runs of Mplus retained with their `.out`, checked for normal termination and a replicated best log-likelihood before any value is read. |
| `latentgold/` | Latent GOLD 6.1 comparison, **run 2026-09-22: 260 quantities compared, 260 agree, 0 disagree** across twelve cases (two-level continuous, full covariance, categorical, mixed, FIML, covariates, bivariate residuals, three-step, latent transitions, and the bundled `course_engagement` data the README reports). Likelihoods, parameter counts, posteriors at both levels, means and variances recomputed from Latent GOLD's posteriors, and the information criteria in both sample-size conventions. Latent GOLD's own listings and posterior files are retained in `latentgold/returned/` with `PROVENANCE.md` and an input fingerprint, so `compare.R` runs offline and Latent GOLD is not needed again unless the kit changes. Standard errors, bivariate residuals and the Step-3 tables are still `awaiting parser`. |
| `deferred.R`, `fixtures.R` | Shared by the Mplus scripts: `deferred_verb()` loads a feature moved to `future/` so the comparison that would restore it still runs; `mplus_fixture()` resolves a fixture to `tests/fixtures/` or `equivalence/fixtures/` by who reads it. |
| `simulab-recovery/`, `three-step/` | Monte Carlo recovery studies, which measure bias and coverage rather than agreement with a fixed number. |
