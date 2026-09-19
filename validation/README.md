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

It writes `equivalence/report.csv`, one row per compared quantity with its
reference value, the value multilpa produced, the difference and the tolerance,
and `equivalence/REPORT.md`, the summary table. Both are generated, so no
number in the documentation is transcribed by hand.

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
| `latentgold/` | Fixtures and targets prepared for a Latent GOLD cross-check. Latent GOLD 6.1 is free for academic use but Windows-only, so no output exists yet; `compare.R` reports that rather than implying otherwise. |
| `simulab-recovery/`, `three-step/` | Monte Carlo recovery studies, which measure bias and coverage rather than agreement with a fixed number. |
