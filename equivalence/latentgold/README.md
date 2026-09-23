# Latent GOLD comparison kit

Latent GOLD is the reference implementation of the multilevel latent class
model this package fits (Vermunt 2003) and of the three-step corrections
(Vermunt 2010; Bakk, Tekle & Vermunt 2013). It is a Windows program, but runs
on this machine under Wine; see "Pilot" below. This directory builds the data
and syntax, has Latent GOLD estimate them, and compares the output here.

**Status: run, 2026-09-22. 372 quantities compared, 372 agree, 0 disagree;
7 more shown but not comparable, each with its reason.**
Latent GOLD's own output is retained in `returned/`, so `compare.R` runs
offline: **Latent GOLD is not needed again** unless the kit's data or syntax
change, which `compare.R` detects (`multilpa_stale_latentgold_output`).

| What was compared | Rows | Result |
|---|---|---|
| Log likelihood and parameter count | 24 | likelihoods within Latent GOLD's printed precision (4.8e-5 at worst); counts exact |
| Posteriors, over every unit and class | 22 | largest difference 5.6e-6 |
| Means, variances, group-class proportions, recomputed from Latent GOLD's posteriors | 118 | largest difference 9.7e-7 |
| Information criteria, in both sample-size conventions | 75 | largest difference 4.8e-5, which is twice the likelihood's rounding |
| Sample sizes (cases, groups) | 21 | exact |
| Estimates read from the listing's tables: means, variances, class probabilities, logits, Step-3 means and R3STEP logits | 58 | largest difference 4.9e-5, inside the printed rounding of 5e-5 |
| Standard errors read from the listing, including R3STEP's | 51 | largest difference 5.0e-5 (4.95e-5), inside the printed rounding |
| Bivariate residuals, as the ranking of indicator pairs | 3 | identical ranking |

Seven rows are `not comparable`, and show both values with the reason:

- the six Step-3 distal-mean standard errors. This package's are
  cluster-robust, which on this single-level case is the Huber-White
  sandwich. Latent GOLD's are its `standard` errors, the option its listing
  records, and for the modal and proportional models they assume one variance
  shared by the classes. Its BCH errors are also well below this package's.
  Settling that would need Latent GOLD rerun with robust standard errors.
- one group-class logit in `c07`, whose Latent GOLD coding is a sum of two
  estimated terms. Its standard error needs their covariance, which the
  listing does not print.

Two findings came out of the first full run, both fixed:

- **The kit wrote invalid syntax** for every case without covariates:
  `GClass <- 1 + ;`, because `paste0(" + ", character(0), collapse = "")` is
  `" + "`, not `""`. Only the one case with covariates escaped it. Every
  generated file is now checked in `test-compare.R`.
- **This package was stopping short of the optimum.** Eight posterior
  comparisons sat at 1.1e-4 to 1.7e-4 against a 1e-4 tolerance. Tightening
  Latent GOLD changed nothing; tightening this package from `tol = 1e-10` to
  `1e-14` dropped them to about 2e-6. The targets are now fitted at `1e-14`,
  rather than the tolerance being widened to fit the result.

## What it checks that nothing else does

The package is already checked against Mplus (continuous two-level
likelihoods), depmixS4 (single-level latent transitions) and tidySEM (BCH
arithmetic). Latent GOLD adds:

| Case | Model | Anchors |
|---|---|---|
| `c01_continuous_varying` | 2 profiles x 2 group classes, 3 continuous, class-specific variances | likelihood, posteriors at both levels, means, variances, group proportions, standard errors |
| `c02_continuous_equal` | as c01, one variance per indicator | as c01 |
| `c03_full_covariance` | 2 correlated continuous, class-specific full covariance | likelihood, posteriors |
| `c04_categorical` | 4 nominal items, 3 categories, two-level | likelihood, posteriors |
| `c05_mixed` | 2 continuous + 2 nominal, two-level | likelihood, posteriors, means, variances |
| `c06_missing_fiml` | c01 with 10% of values missing at random, FIML | likelihood, posteriors |
| `c07_covariates` | c01 with a profile predictor and a group-class predictor | likelihood, posteriors, standard errors |
| `c08_bivariate_residuals` | y1 and y2 locally dependent | likelihood, posteriors, residual ranking |
| `c09_three_step` | single-level measurement, then four Step-3 analyses | likelihood, posteriors, BCH / modal / proportional means, R3STEP slope |
| `c10_lta` | latent transitions, 2 states, 4 occasions | likelihood, state posteriors |
| `c11_lta_mixture` | c10 with two sequence classes | likelihood, posteriors at both levels |
| `c12_course_engagement` | **the package's own bundled data**: 1,422 enrolments in 106 students, five indicators, 2 x 2 classes — the model the README reports | likelihood, posteriors at both levels, means, variances, information criteria |

Eleven cases are simulated from a fixed seed; c12 is the bundled
`course_engagement` data, so the package's published example is checked too.
Every case is fitted here with 50 starts, and
refused as a target unless its best likelihood is reached by at least two
starts (every case reaches it with all 50).

## Running it

Only steps 1 and 4 are needed to reproduce the comparison from the retained
output. Steps 2 and 3 are needed only after the kit's data or syntax change.

1. **Build** (development machine, ~7 minutes):
   `Rscript equivalence/latentgold/make-kit.R`
   Writes `kit/` (data, syntax, `run-all.bat`, `MANIFEST.txt`) and `targets/`.
2. **Estimate** with a licensed Latent GOLD 6.1. On Windows, set the `LG` path
   at the top of `run-all.bat` and run it. Under Wine on this machine, run each
   model as `wine "C:\\Program Files\\LatentGOLD6.1\\lg61.exe" <case>.lgs /b /o
   <case>.lst` with `WINEPREFIX=~/.wine-latentgold`, from a directory inside
   the Wine drive. Either way the models run in the order `MANIFEST.txt` lists:
   `c09_step1` must precede the four `c09_*` Step-3 models, which read its
   posterior file.
3. **Collect** every `.lst` listing and every `*_posteriors.txt` file into
   `equivalence/latentgold/returned/`.
4. **Compare**: `Rscript equivalence/latentgold/compare.R`
   Writes `comparison.csv` and `REPORT.md`, and exits with status 1 if any
   quantity disagrees. After a fresh Latent GOLD run, record it with
   `Rscript equivalence/latentgold/record-run.R`, which fingerprints the inputs
   and writes `returned/PROVENANCE.md`.

`Rscript equivalence/latentgold/test-compare.R` tests the comparison itself on
imitated Latent GOLD output: rounded, label-swapped posteriors must agree, and
a moved posterior, a wrong likelihood, a wrong parameter count, a renamed
column, an absent file and a blank listing must each be caught.

## Pilot: done, except the licence

**Latent GOLD 6.1 runs on this machine**, under Wine (Wine Devel 11.17 in
`~/.wine-latentgold`, installed 2026-09-22), so the kit no longer needs a
Windows machine. The pilot of 2026-09-22 confirmed, against Latent GOLD itself:

| Assumption | Result |
|---|---|
| Batch form `lg61.exe file.lgs /b /o file.lst` | works |
| `//LG6.1//` header, `options`/`variables`/`equations` blocks | accepted |
| `bayes categorical=0 variances=0 latent=0 poisson=0` | accepted; the listing reports `Log-prior 0.0000`, so it is maximum likelihood |
| `groupid`, `latent GClass group nominal 2`, `Cluster <- 1 + GClass` | accepted (same forms as the bundled `peetmis_smmr.lgs`) |
| `y1 \| Cluster` for class-specific variances | accepted: 15 parameters against 12 without |
| **omitting** the variance line for equal variances (c02) | identical to writing the bare `y1;` line the bundled examples use: 12 parameters, same likelihood |
| `y1 <-> y2 \| Cluster` | accepted; adds one covariance per class |
| nominal indicators (c04, c05) | accepted |
| `independent x, w` with no scale type | accepted; one slope, as this package fits it |
| `missing includeall` (c06) | accepted |
| `caseid`, `State nominal dynamic`, `State[=0]`, `State <- 1 \| State[-1]` | accepted, and **exactly equivalent** to the bundled `DFG.lgs` form `(~tra) 1 \| State[-1]`: same parameter count, same likelihood to four decimals |
| `State <- 1 \| State[-1] GClass` for class-specific transitions (c11) | likewise equivalent to the bundled multilevel Markov form |
| `step3 modal bch` with `latent Cluster nominal posterior = ( ... )` | accepted |

Two findings were fixed in `compare-functions.R`, both confirmed against a real
listing rather than assumed:

- The parameter count is labelled **`Number of parameters (Npar)`**. `Npar`
  alone appears only as a column heading of the summary table, on a line with
  no value; the parser looked for that and would have failed on every run.
- The listing is **Latin-1**, not UTF-8 (it prints `Entropy R²` as one byte).
  Read as UTF-8 it yields invalid strings that match no pattern.

Both are covered by `test-compare.R`, whose imitated listing now copies the
real labels, separators and encoding.

### What is still blocked: the licence

The installed Latent GOLD runs in **Demonstration Mode**, which opens only its
own bundled sample data and writes no output files. With the kit's data it
reports `Can't open Input data`, and `outfile` silently writes nothing. That is
a licence restriction, not a syntax error: byte-identical copies of a bundled
file fail under any other name.

So the two things still unverified are exactly the two the licence gates:

1. the models on this package's data, and
2. the **posterior file** `outfile ... classification keep id g` writes, whose
   column names (`Cluster#1`, `GClass#1`, ...) `compare.R` expects.

Latent GOLD is free for academic use. Register from its Registration dialog:

```
WINEPREFIX=~/.wine-latentgold ~/Applications/"Wine Devel.app"/Contents/Resources/wine/bin/wine "C:\Program Files\LatentGOLD6.1\lg61.exe"
```

Once it is licensed, run the kit with `kit/run-all.bat`'s commands through
Wine, put the output in `returned/`, and run `compare.R`.

## How agreement is judged

- **Log likelihood**: within half a unit of the last digit Latent GOLD printed,
  plus 1e-5 for two independent searches stopping at their own criteria.
- **Parameter count**: exactly.
- **Posteriors**, at each level: the largest absolute difference over every
  unit and class, after matching Latent GOLD's arbitrary class labels to this
  package's by the ordering that brings them closest. Tolerance: the file's
  printed rounding plus 1e-4.
- **Means and variances**: recomputed from Latent GOLD's own posteriors. At
  convergence EM is at a fixed point, so its estimates are exactly the
  posterior-weighted means and variances; this recovers them without scraping
  the listing's parameter tables. The tolerance is the bound the posteriors'
  printed rounding puts on a weighted moment, derived in
  `.lg_moment_rows()`. Only for complete data with a diagonal covariance, where
  that fixed point is the plain weighted moment.
- **Group-class proportions**: the mean of Latent GOLD's group posteriors.
- **Information criteria**: each is matched to the column computed on the same
  sample size, read from the listing, not by its label. Latent GOLD's "case" is
  an observation in a two-level model but a whole sequence in a `caseid` model,
  which is this package's *group*; matching by label alone compares a sequence
  model's BIC against `bic_individual`. A criterion is `-2 * LL` plus a penalty
  in the parameter count and the sample size, both compared separately and
  exactly, so the tolerance is twice the likelihood's plus the printed
  rounding. Latent GOLD's `AIC3`, `CLC`, `AWE`, entropy R-squared and
  classification errors are **not** compared: those definitions differ from
  this package's, so a comparison would measure the definition.

These are all invariant to how either program parameterises the model, which
is why they are compared first: Latent GOLD's logit coding differs from this
package's, so its raw coefficients are compared only after re-coding (see
"The listing's tables").

## The listing's tables

Standard errors, bivariate residuals and the Step-3 results are printed only
in the listing's tables. `compare-functions.R` reads them. The parsers were
written against the real listings in `returned/`, and `test-compare.R` tests
them on those listings, with deliberate breaks that must disagree:

- **Parameters** come from the `Parameters` section. A name spans three
  fields in one listing and five in another, so the coefficient is found by
  content, not position. Natural-scale values and their standard errors come
  from `EstimatedValues-Model`.
- **Logits** are re-coded before comparison. Latent GOLD codes against its
  first class and this package against its last. A log odds between two
  classes is the same in any coding, so each of this package's coefficients
  is a difference of Latent GOLD's, taken through the class alignment found
  from the posteriors.
- **Bivariate residuals** are a different statistic in the two programs.
  Latent GOLD prints one per indicator pair, pooled over classes. This
  package tests each pair within each profile. What is compared is the
  ranking of the pairs.
- **Step-3** maps as follows.

  | this package | Latent GOLD |
  |---|---|
  | `three_step(method = "bch")` | `step3 modal bch` |
  | `three_step(method = "modal")` | `step3 modal none` |
  | `three_step(method = "proportional")` | `step3 proportional none` |
  | `r3step()` | `step3 modal ml` |
