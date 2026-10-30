# latents 1.1.3

## Zero-weight units keep their classification (2026-10-30)

* 1.1.2 (a GitHub release) reported the class, uncertainty and posterior
  probabilities of units with sampling weight zero as `NA` in assignment,
  posterior and prediction tables. These units are again reported as in 1.1.1.
  A unit of weight zero stays out of every count, criterion and starting
  value, but keeps its row, posterior and assignment in output. Only rows and
  groups with no observed indicator are left without a class.

# latents 1.1.2

## Rows without data are no longer assigned a class (2026-10-29)

Estimates, log likelihoods and posterior probabilities used in estimation
are unchanged.

* A row whose every indicator is missing contributes nothing to the fit, so
  its posterior equals the class proportions. With the default
  `drop_empty = TRUE`, such rows already counted for nothing in the sample
  sizes and classification summaries. The assignment table then gave them
  the largest class, however. Assignment and posterior tables now keep these
  rows in their input order, but report their class, uncertainty and
  probabilities as `NA`. The same applies to `predict()` on fitted or new
  data, `lta()` assignments, `descriptives()` by class, `sensitivity()`
  agreement and sequence summaries. A group whose rows are all empty has no
  group class. An empty row in a group with data keeps the group's class.
  `drop_empty = FALSE` keeps the earlier output.
* A categorical `local_dependence` pair joins its two items into one. A row
  could answer only one item of the pair, with no other indicator observed.
  Such a row was wrongly dropped from the sample sizes, information criteria
  and classification summaries, and is now counted as observed. Likelihoods
  are unchanged.
* `predict()` for a model with a noise component now refuses a row with only
  some indicators observed, with condition class `latents_unsupported_noise`.
  The fitted noise volume does not define a density for a subset of
  indicators. A row with every indicator missing is accepted and, like any
  empty row, reported without a class.

## Documentation

* Rewrote the introductions to the help pages, the bundled dataset
  descriptions, the package overview (`?latents`), the DESCRIPTION summary and
  the opening of the LPA vignette.

# latents 1.1.1

## Help pages (2026-10-29)

* S3 methods (`print()`, `summary()`, `as.data.frame()`, `plot()`,
  `get_results()`, `parameter_inference()` and others for each class) no
  longer have help entries of their own: the help index has 57 entries instead
  of 222, and the page for each function says which methods its objects have.
  The arguments only `plot()` takes are listed in a section of
  `?"latents-plots"`, and `?parameter_inference` covers every model family.
  The topics `latents-print`, `latents-summary` and `latents-as-data-frame`
  are gone; their content is on `?get_results`. No function changed.
* Clearer descriptions of `diagnostics()`, `get_results()`, `report()`,
  `pool_imputations()`, `mixture_regression()` and `bootstrap_lrt()`, and of
  two conditions in `?"latents-conditions"`.
* Tests: 112 tests that compare results with other packages, numerical
  routines, Latent GOLD, published values or independent re-implementations
  moved from the package tests to `tests/equivalence/` in the source
  repository, where continuous integration runs them. No test was removed.

# latents 1.1.0

## CRAN release (2026-10-29)

This release brings the changes since CRAN 0.9.21 (the GitHub releases
0.9.22 to 1.0.0) to CRAN. The estimation code is the same as in 1.0.0 and
0.98.1, so every result agrees with those versions. Some results differ
from 0.9.21: "Upgrading from CRAN 0.9.21" under 0.98.1 lists when to refit.

* Tests: the boundary-fit test that uses numDeriv for its standard errors is
  skipped when numDeriv, a suggested package, is not installed.

# latents 1.0.0

## Verification against independent programs and published analyses (2026-10-29)

No estimation code has changed since 0.98.1; every result is the same. This
release records how each function is verified.

* Comparisons cover independent model fits, published analyses, component
  formulas and independent reference implementations. The verification
  matrix distinguishes these evidence types and counts only complete,
  agreeing model comparisons toward its independent-source strength. The
  record, rebuilt from saved results, is
  `equivalence/verification-matrix/VERIFICATION-MATRIX.md` in the source
  repository.
* New comparisons: Latent GOLD 6.1 (growth mixtures, latent class, profile
  and transition models, negative-binomial and ordinal mixture regression,
  and individual fits at one to four classes), glca, randomLCA, BayesLCA,
  multilevLCA, depmixS4, Rmixmod,
  EMCluster, mixtools, MoEClust, flexmix, LMest, seqHMM, lcmm, mitml,
  mitools, and mclust and poLCA for enumeration and the bootstrap LRT; and the
  published values of SAS PROC FMM, Stata [SEM] Examples 50g and 52g, Masyn
  (2013), Proust-Lima, Philipps and Liquet (2017), Bartolucci, Pandolfi and
  Pennoni (2017), Pennoni, Pandolfi and Bartolucci (2024), Benaglia et al.
  (2009), Grün and Leisch (2008) and the SAS and R blog's latent class example.
* Where a reference reaches a different maximum the difference is recorded
  with its cause. mclust's default EM tolerance stops up to 0.04 short of the
  maximum, and the statistics printed in its vignette are that stopping point.
  The default kmeans start search can return a lower maximum for VVV models
  with two profiles (banknote, crabs), where `search = "emEM"` reaches the
  higher one.
* Tests: the release-engine comparison of grand-scaled EVV and VVV fits now
  converges both fits and compares posteriors up to label order (it failed on
  Linux and Windows), and the poLCA residual statistic is compared to 1e-5.
* Verification safeguards reject missing comparisons and incomplete result
  tables. Known-difference exceptions cannot waive parameter-count errors;
  the LCA wrapper includes the added published HELP reference. The matrix
  gives no source credit to skipped or errored tests, verifies saved results
  against runtime-source and result hashes, and excludes internal checks
  from its external-check headline. All ten new comparison kits have CI
  jobs that fail on unexpected skips as well as errors and disagreements.

# latents 0.98.1

## Further release audit repairs (2026-10-28)

* Equal-variance Gaussian and equal-dispersion negative-binomial mixture
  regressions again report their residual standard deviations or dispersions
  and standard errors. Class names on the internal variances had changed the
  shared parameter's name, leaving the result-table lookup empty, even with
  one class.
* General-engine LTA covariance calculations now differentiate and invert in
  standardized parameter coordinates for observed, robust and OPG inference.
  Small outcome units could distort standard errors or cause a false
  singular-information refusal; large units could also cause that refusal.
* Zero-weight regression units are excluded from inference transformations
  and scores, as they already were from fitting. Extreme values on excluded
  units could otherwise erase standard errors. Additive models normalize
  their posteriors more stably, so excluded groups with large finite values
  no longer abort an otherwise valid fit.
* `compare_models()` retains every model and its ordering through wrappers
  forwarding `...`, including models named `n_profiles`, `n_group_classes`
  or `cutoff`. These labels could silently drop a model and change BIC weights.
* Diagnostics count the noise component as a classification alternative.
  One Gaussian profile plus noise now receives classification checks instead
  of being described as having nothing to classify.

## Upgrading from CRAN 0.9.21

Refit a model made with 0.9.21 when one of these applies; each is described
in full further down.

* **Gaussian indicators recorded in small units** (a proportion, metres, a
  rate). 0.9.21 floored variances at an absolute `1e-6` and stopped EM
  relative to `|log L|`: the same data multiplied by 0.001 and by 1e-6 gave
  log likelihoods 654 and 37,700 below the fit in the original units, with
  only a boundary warning. Both rules are now unit-free ("Fits that do not
  depend on the units of the data"), and so are SQUAREM and the
  quasi-Newton finishes ("SQUAREM and replication in any units").
* **Occasion labels such as `T1` ... `T10`.** 0.9.21 ordered them
  alphabetically (`T1`, `T10`, `T2`, ...), so `lta()` and `multilpa(time =)`
  fitted the occasions in the wrong order. Such labels are now refused with
  `latents_bad_time`; give the order with `factor(levels = )` or zero-pad
  the numbers.
* **Sampling weights with zeros.** A unit of weight zero is now absent from
  every count, criterion, entropy and start ("Sampling weights").
* **Several group classes.** Start 1 of a multilevel model is now
  deterministic and can number the classes the other way round ("Class
  numbering of multilevel models").
* **ICL, AWE and CLC** charge the entropy of the complete classification
  ("Information criteria"), and EII, EEI and EEE default to
  `search = "emEM"` ("Default start search").
* **Mixture regressions with `common = "intercept"`.** The argument had no
  effect in 0.9.21; it is now refused.

## Latent GOLD: mixture regression and sampling weights (2026-10-28)

Mixture regressions (Gaussian, Poisson, logistic; model-based and robust
standard errors) and sampling weights in `mixture_regression()`, `lpa()`,
`lta()` and two-level `multilpa()` were compared with Latent GOLD 6.1: 172
estimates, standard errors, log likelihoods and parameter counts, all within
its printed rounding (`equivalence/latentgold-wider/`). The comparison found:

* `coef()` of an `lta()` fit run by the general engine (sampling weights,
  transition or initial covariates, occasion-varying transitions, second
  order) returned each continuous mean minus the indicator's grand mean, the
  engine's internal coordinates. It now returns the means in the data's
  units, as the profile table always did. `vcov()`, standard errors and every
  table were unaffected.
* `mixture_regression()` stored its class variances under the names of their
  positions before the classes were relabelled. No reported number used the
  names.

## Release review (2026-10-27)

A whole-package review confirmed 40 defects. The repairs below each have a
regression test; SQUAREM and the replication count are in the next section.
The enumeration's `criteria` table, which duplicates part of its grid, is
kept because 0.9.21 code reads it.

### Results that change

* A latent transition origin that no expected move leaves (an empty row)
  reported its starting values as estimates and counted them as parameters.
  Its probabilities are now equal, flagged `estimated = FALSE`, left out of
  `n_parameters`, AIC and BIC (by K - 1 per empty row), and warned with
  `latents_empty_transition_row`. Wald inference reports `NA` for them and
  keeps the rest of the table, which it used to refuse as singular.
* A three-step class carried by a single group (or, under `vcov_type =
  "independent"`, a single unit) reported a standard error of about 1e-16
  and a zero-width interval under BCH, proportional and modal assignment, and
  its pairwise contrasts p = 0. They are now `NA`, with
  `latents_small_classes`.
* `sensitivity()` of an `lta()` fit refits with the fit's own `tol` and
  `max_iter`; it used `lta()`'s defaults.
* Zero-weight units no longer choose the starts or the maximum of a mixture
  regression or growth mixture, and a level or category seen only on
  zero-weight rows adds no parameter (regressions, ordinal indicators and
  `local_dependence` pairs).
* An autoregressive growth residual held on its bound reports `NA` for its
  standard error and interval, not zero.
* `diagnostics()` of a single-level fit, or of a two-level fit with one
  group class, no longer has the trivial "groups" row in its classification
  and average-posterior tables, and `get_results(fit, "all")` follows it.

### Refusals and warnings

* A latent class model with more free parameters than response patterns
  allow warns `latents_unidentified` (Goodman, 1974); its estimates are
  unchanged.
* `confint()`, `get_results()` and `summary()` of regression, growth and
  additive fits refuse a level outside (0, 1), which used to give an
  inverted interval, and refuse arguments they do not take. `vcov()` warns
  `latents_unused_argument` for arguments it ignores, and refuses
  `method = "bootstrap"` on covariate fits instead of returning the Wald
  matrix. `summary()` of a fit, an enumeration or a bootstrap LRT warns for
  an argument it ignores (`summary(fit, digits = 3)` used to do nothing).
* `get_results()` refuses an abbreviated table name.
* More input mistakes raise classed conditions: three-step contracts,
  unsupported three-step families, and the option checks of `multilpa()`
  (`latents_bad_argument`).

### Other changes

* `summary()` of a general LTA fit is a result object: read its table with
  `get_results(summary(fit), "table")`. Its transition and second-order
  tables gain an `estimated` column, and its assignments take the shared
  layout.
* `coef()` and `confint()` methods for pooled imputations; the intervals are
  Rubin's t intervals with each parameter's degrees of freedom.
* The default plot of an all-categorical fit with text category labels draws
  the response probabilities instead of failing.

### Known issue

* A covariate fit stops on the change in log likelihood, which along a flat
  direction of the logit coefficients cannot place them closer than about
  1e-6 to the maximum, and tightening `tol` does not help. On the
  equivalence data the coefficients differ from a reference optimiser by
  2.6e-6, against standard errors of 0.16 to 0.82. A stop on the score of
  the coefficients would remove this.

## SQUAREM and replication in any units

* SQUAREM extrapolated on the raw parameters, whose means are in the
  indicators' units and variances in their squares, so its step length,
  and with it the maximum a start reached, depended on the units. With
  `test_anxiety` of `srl` multiplied by 1000 (VVV, four profiles, five
  starts, seed 3) the fit and its grand-scaled refit differed by 4.92 in
  log likelihood. SQUAREM now steps in standard-deviation units in
  `multilpa()` (and so `lpa()`, `lca()` and `multilca()`), both `lta()`
  engines and the additive family, and the two fits agree to 1e-12. Plain
  EM was already unit-free and is unchanged bit for bit. Accelerated fits
  can stop at a different iteration and, from some starts, at a different
  maximum: on the package's 192 recorded reference fits, 179 log
  likelihoods moved by less than 1e-6, 9 by up to 1e-4, three rose (by
  0.19, 0.84 and 2.29) and one, a short staged fit, fell by 0.0006.
* The quasi-Newton finish of the general `lta()` engine (covariates,
  occasion-specific transitions, weights, mover-stayer, second order) and of
  growth mixtures stopped relative to `|log L|`, which a change of units
  moves, and searched the means in the units of the data. With the
  indicators multiplied by 0.001, a mover-stayer fit of `course_engagement`
  stopped up to 0.028 below the maximum it reaches in the original units.
  Both now stop on the change per observed value, the general `lta()`
  engine searches the means in standard deviations, and the fit is the
  same in any units.
* Starts whose log likelihoods are within the fit's own precision, `tol`
  times the number of observed values (the change at which its EM stops),
  are now treated as tied, and the earliest converged one among them is
  reported, in every family. The window used to be a relative `1e-10` of
  `|log L|`, which a change of units moves and which at `|log L|` in the
  thousands spans 4e-7: a weakly separated transition model fitted with
  `tol = 1e-12` reported a start 3.2e-7 below another, a gap depmixS4's EM
  then recovered from the reported point. With the default `tol = 1e-8` the
  window is wider than before, so between starts that agree to the precision
  asked for the earlier one is now reported more often.
* The additive family's search for between-group variances at zero
  compared a derivative with respect to a variance, which has units, with a
  slack proportional to `|log L|`. It now compares the derivative in
  reference-variance units with a slack per observed rating, and accepts a
  boundary fit within `tol` per observed rating, as EM stops.
* Group classes of the additive family that the data cannot tell apart are
  now flagged as collapsed when setting them equal lowers the log
  likelihood by less than the replication tolerance, not only when their
  parameters agree to 1e-8. EM approaches such a pair along a flat ridge and
  rarely brings it that close: a dispersion fit whose classes differ only in
  means was flagged from one seed in twelve, and its inference then failed
  with a singular information instead of the collapsed-class refusal. It is
  now flagged from every seed tried (`latents_weak_class`).
* A start counts as reaching the best log likelihood when it is within
  `1e-6` times the number of observed values of it, no longer
  `1e-6 * (1 + |log L|)`, which moved 40-fold with the units for the same
  gaps between starts. For data in ordinary units the two are of the same
  size, so the counts printed for such fits barely change. Every family
  uses this rule and stores it as `replication_tolerance`, which
  `sensitivity()` reads; the general `lta()` engine, mixture regressions
  and growth mixtures did not store one before.

## Astra review repairs

* Zero-weight people no longer determine the covariate coordinates or boundary
  checks of an LTA fit. Ordinary categorical levels are discovered on positive
  support; excluded-only categories do not add parameters or change BIC.
  Unidentified transition coefficients have unavailable bootstrap covariance
  rows and columns as well as unavailable standard errors.
* A constant Gaussian response is accepted when its design or offset leaves
  residual variation. Argument forwarding uses distinct internal bindings and
  evaluates each supplied expression once, including through `multilca()`.
* Diagnostics and regression enumeration use the fitted family's contributing
  units, sampling weights and classification levels. Newly supported diagnostics
  families expose the documented entropy and classification tables. Unfinished
  sensitivity refits remain visible and are not labelled as replicated maxima.
* Enumeration examples use the current complete-data entropy convention and
  updated ICL/AWE winners. The full-output equivalence harness checks the new
  replication likelihood, fails partial output coverage, and runs all reference
  suites in CI with explicit failures for unavailable dependencies.

## Class numbering of multilevel models

* Start 1 now takes its group classes from Ward's clustering of the groups'
  profile compositions, as its profiles already came from Ward's clustering,
  so it no longer depends on the seed. A model with more than one group
  class can reach the same maximum as in earlier versions with its group
  classes, and sometimes its profiles, numbered the other way round: the
  log likelihood and estimates are unchanged, only the labels swap. Scripts
  that refer to a class by its number should check which class it now is.

## diagnostics() as a structured report

* `diagnostics()` reports a fit in five domains (estimation, class size,
  classification, local independence, information), one row per quantity
  with its value, a status and, where one defines it, its reference: how
  many starts reached the best log likelihood, starts that failed or did not
  finish, estimates at a bound, the smallest class, classes with no assigned
  unit, relative entropy, the lowest average posterior and odds of correct
  classification, and the bivariate residuals significant at 0.05 after
  Benjamini-Hochberg adjustment. Values are reported without cut-offs: the
  published ones (AvePP .70, OCC 5, a 5% smallest class) are rules of thumb,
  and the package does not impose them. A status marks only facts: a best
  start that did not converge or a class estimated from less than one unit
  (`fail`); failed or unfinished starts, estimates at a bound, a class with no
  assigned unit, or dependent indicator pairs (`check`). The print names
  those first. `get_results(diagnostics(fit), "checks")` is the table; every
  value is read from a table the fit already computes, and no estimate
  changes. `as.data.frame()` of the result is still the entropy table.
* `diagnostics()` now reports on additive, cross-level, mixture-regression
  and growth-mixture fits, which it used to refuse; a quantity a family
  cannot support is `n/a` with the reason. A one-class fit and a class with
  no assigned unit no longer print `NA`.
* Starts that a search screens (short runs by design) are reported apart from
  finished starts, never as unconverged.

## A readable result when it is printed or clicked

* A fitted latent class model prints its class sizes (and group-class sizes)
  and a wide table of response probabilities, one row per indicator category
  and one column per class, rounded to three decimals; it used to print a
  long table of probabilities and thresholds at full precision.
* `str()` of a fit, an enumeration or a `diagnostics()` result is two lines:
  the result's header and where its tables are. A host that shows `str()`
  before `print()` when an object is clicked (Beatrina, RStudio's
  Environment pane) used to show some ninety components of internals first.
  `str(unclass(x))` still shows the raw list.

## How often the best solution was reached

* `print()` of a fit and `diagnostics()` now say how many starts reached the
  best log likelihood ("best log likelihood reached by 4 of 10 starts"), and
  flag a best solution found by a single start. `print()` used to show the
  winning start's index as "best start: 5/10", which read as a count; the
  count itself was only in `summary()`. The numbers are in
  `get_results(fit, "replication")` and `get_results(diagnostics(fit),
  "replication")`.

## Clearer refusals

* Every common input mistake is now refused with a `latents_*` condition
  class, no call, and a message that names what is wrong and how to fix it.
  Of 56 mistakes across the exported verbs, 23 used to raise a bare
  `simpleError`; none does now (contract test
  `test-review-1027-error-contract.R`).
* A misspelt, missing or repeated column, and an `id` listed among `vars`,
  each have their own message, with the likely spelling: "`vars` names `Z`,
  which `data` does not have. Did you mean `z`?" (the four used to share
  "Supply at least two rows, unique existing indicators, and one distinct
  group column."). The same for `id`, `time`, a three-step `outcome` and
  `descriptives(vars = )`.
* Missing values name the indicators and their counts and point to
  `missing = "fiml"`; infinite and NaN values are reported separately,
  because FIML treats only `NA` as missing. A non-numeric indicator in a
  profile model names the column and its class and points to `categorical`
  or `lca()`.
* Passing a fitted model or an enumeration where the data belong now says so
  by the caller's own name: "`data` must be a data frame, but `grid` is an
  enumeration result ... Pass the data frame the model is to be fitted on."
  Passing data where a fitted model belongs raises the new class
  `latents_not_a_fit` (also a `latents_bad_argument`), saying what was
  passed; `parameter_inference()`, `get_tna()` and `get_group_tna()` gain
  default methods that do the same instead of R's "no applicable method".
* The plot verbs (`plot_diagnostics()`, `plot_profiles()` and the rest)
  name the object they were given and, for a data frame, say that the plot
  is drawn from the fitted model itself.
* `candidate_fit()` no longer evaluates a call headed by `get_candidate`'s
  body with the arguments by value, so its conditions carry no data.

## Conditions never carry the data (0.97.4)

* An error raised while a fitting verb validated its input could carry the
  whole data set in its call: `multilca()` forwarded to `multilpa()` with
  `do.call()`, so the recorded call held the data by value, and `multilpa()`
  put the setup function's body (375 lines) at the head of the call it
  evaluated. R's own error printing truncates the call, but a host that
  deparses `conditionCall()` in full streamed every value of a 400,000-row
  data set (a 10,000-row object printed as 3.3 MB; now 83 bytes). The
  fitting verbs now evaluate calls of symbols only (`multilpa()`,
  `multilca()`, `lta()`), a non data frame is refused with
  `latents_bad_data` naming the object's class, and every candidate error
  that is stored on a grid or raised again has its call bounded
  (`<data.frame[400000 x 9]>` in place of the values).

### Third pass: interactions (2026-10-27)

* A plain categorical item whose labels happen to contain `::` is no longer
  read as a half-observed joined item (the one-class log likelihood was 47.8
  too high); only a declared `local_dependence` pair carries that meaning.
* `compare_models()` again accepts fits named `n_profiles`, `n_group_classes`
  or `cutoff`, as 0.9.21 did; numeric values keep their filtering meaning.
* A positive sampling weight whose ratio underflows stays a present unit.
* Regression, growth and additive fits apply the sampling-weight convention
  to their entropy, average posteriors and classification tables; the
  two-level regression ICL charges the joint classification entropy.
* Covariate inference keeps half-observed pair rows; zero-weight outcomes no
  longer set a Gaussian floor or a growth coordinate; the regression floor
  stays relative below machine epsilon; a random-intercept growth model whose
  every person has one observation is refused as `latents_not_identified`;
  Brownian residual bounds and starts follow the time units.

### Second review of every formula (2026-10-26, round 2)

* **Standard errors of centred fits.** Under `centering = "grand"` or
  `scaling = "grand"` the Wald standard errors treated the centring constants
  as known and were 22-31% too small; they are now estimated jointly (the
  centring equations are stacked into the information) and calibrate to
  within 10% in simulation. Under person centring or scaling the observed
  information is too small by construction, so such fits now default to
  `vcov_type = "robust"` and refuse `"observed"` and `"opg"`; the robust
  sandwich, which used to be refused as rank deficient, is correct there.
* **Information criteria.** ICL, AWE and CLC charge the soft entropy of the
  complete-data classification (group classes and profiles) under both
  conventions; the conventions differ only in the sample size. Under the old
  rule the "groups" convention charged only the group-class entropy, which is
  zero for every single-level fit, so `icl_groups` equalled BIC and
  `clc_groups` always chose the largest model. Transition fits keep one
  entropy per level. Fits with several group classes saved by 0.97.1 report
  `NA` for these three until refitted. KIC is documented as
  −2L + 3(q + 1), the form tidyLPA computes.
* **Sampling weights.** Every classification estimate (entropy sum of the
  criteria, relative entropy, average posteriors, odds of correct
  classification, the classification-error matrix and the BCH weights) uses
  the sampling-weighted posteriors; modal counts stay counts. Units of weight
  zero are absent from every count, estimate, start, refusal and warning in
  every family, including the cross-level fits (whose restricted log
  likelihood was inflated by 19/15 with zero-weight groups), grand centring
  and scaling, and `drop_empty = FALSE`.
* **A half-observed categorical pair** (one item of a `local_dependence` pair
  missing) now contributes its observed item under `missing = "fiml"`; it
  used to be dropped, which biased the estimates under MAR. The
  `latents_partial_pair` warning is no longer raised. Every pair of
  indicators now has a row in `get_results(x, "residuals")`; ordinal
  indicators enter the categorical table, and mixed or count pairs report an
  association without a test. The categorical chi-square is documented as
  conservative.
* **Units and origins.** The VEI, EVE and VVE M-steps no longer floor a
  unitless shape at `min_variance`; `min_variance = NULL` (the new default)
  means 1e-6 times the smallest indicator variance in the fitted units; EM
  stops on a log-likelihood change of `tol` times the number of observed
  values, which does not depend on the units; the membership regression's
  Newton step and the separation flag no longer depend on a covariate's
  origin; mixture-regression standard errors are origin-free up to shifts of
  1e6; growth mixtures standardise their predictors as well as the outcome
  (time in days, hours or calendar years gives the same fit); `lta()`
  covariate designs are fitted on a centred, scaled chart.
* **Shared negative-binomial dispersion** (`count_dispersion = "equal"`) is
  no longer forced to its floor on a local condition (the fit through `lta()`
  was 23.5 log-likelihood units low).
* **Bootstrap likelihood-ratio tests** take their p-value over the valid
  replicates and report `n_invalid`, with a `latents_failed_replicates`
  warning; the p-value is `NA` only when fewer than half the replicates are
  valid. Under the null a spurious extra class is expected to separate or
  stall, so an invalid replicate is information, not failure.
* **Smaller:** start 1 is deterministic for any number of group classes;
  covariate and cross-level fits record `n_best_replicated`; nominal and
  ordinal fits of the same items can be compared; the `r3step()` by-group-class
  table reports `NA` for a non-positive variance; `sensitivity()` counts
  optima by an anchor rule; bootstrap resamples are drawn before any refit; a
  CSV round trip of the data is accepted by every inference verb; `lta()`
  refuses rank-deficient covariate designs and counts only identified logits;
  person centring together with person scaling refuses data whose every group
  has fewer than three observations (`latents_atomic_groups` otherwise);
  `mixture_regression()` counts rows of positive weight in the observation
  nesting and no longer refuses a two-level model with as many groups as row
  classes.

### Fits that do not depend on the units of the data

* `min_variance` now defaults to `NULL`: one millionth of the smallest
  indicator variance (observed values in rows of positive weight, on the
  scale fitted), in `multilpa()`, `lta()`, `fit_staged()` and the group-class
  and cross-level families. The old default was an absolute `1e-6`, so data
  recorded in small units met it in ordinary fits: iris in metres gave a VVV
  solution whose log likelihood was 1.00 below the centimetre fit, and at
  `x 0.001` VVI and EEE were 414 and 465 below, with only a boundary warning.
  A number keeps its meaning, and the bound used is stored on the fit.
  `mixture_regression()` and growth models already used a relative bound.
* EM now stops when a step changes the log likelihood by at most `tol` times
  the number of observed values, rather than `tol * (1 + |log L|)`, which a
  change of units moved. Plain EM on overlapping profiles now stops at the
  same step at units `x 1e-6`, `x 1` and `x 1e6` (0.048 standard errors
  from the maximum); before, it stopped 0.06 to 0.19 standard errors away
  depending on the units. The reported `iterations` can differ from earlier
  versions.
* The VEI M-step compared a shape sum, which has no units, with
  `min_variance`. With indicators whose scales differ by about 1e9 (a rate
  and a population), every start failed with "EM likelihood decreased",
  where mclust fits. The EVE and VVE M-steps had an absolute floor that broke
  them for data in units of about 1e-9.
* Under SQUAREM, `iterations` counts the EM steps kept: a cycle that falls
  back to plain EM is two steps, not three.
* `centering = "person"` with `scaling = "person"` turns a group with two
  observations of an indicator into exactly -1/sqrt(2) and 1/sqrt(2), on
  which a profile can collapse. Such groups now raise `latents_atomic_groups`,
  and data in which every group is like that are refused.
* The `scaling` documentation now lists correctly which structures grand
  scaling leaves unchanged (EEI, VEI, EVI, VVI, EEE, VEE, EVV, VVV) and which
  it changes (EII, VII, EEV, VEV, EVE, VVE).

### Cores

* `cores` larger than the machine's core count is capped at that count with a
  `latents_too_many_cores` warning. On macOS with Apple's Accelerate BLAS, and
  on Windows, each worker is a separate R process holding the package and a
  copy of the data, so `cores = 24` on a ten-core machine started 24 R
  processes and could exhaust memory.

### Parallel seed sensitivity

* `sensitivity()` gains a trailing `cores` argument. The seeds are refitted in
  parallel, each refit on one core and seeding itself, so the table is
  identical to `cores = 1`, for profile, transition and growth fits alike.
  Before, every refit ran on one core whatever the original fit used.

### Classification error and three-step estimators

* Rows with no observed indicator are no longer counted in the classification
  error matrix, the BCH weights or any three-step estimator (`three_step()`,
  `r3step()`), as they already were not in entropy, the average posteriors and
  the information criteria. Such a row's posterior is the class prior, so it
  was assigned to the larger class and charged that class with errors: on 120
  rows of which 3 were empty, P(assigned 1 | true 2) was 0.025 instead of
  5e-5. A group with no observed row is likewise left out at the group level.
  `drop_empty = FALSE` keeps every row, as before.

### Standard errors

* A probability on its upper bound (one minus the floor, because every other
  category of its row is on the floor) has an `NA` standard error, as the
  documentation states. In a binary item it was reported as exactly 0, with
  the interval [1, 1], when it was the last category, and as `NA` when it was
  the first.
* When an item's last category is on the floor (a category one class never
  uses), only that probability is held at its bound. The whole row used to be
  held, so the row's other probabilities had no standard error, and every
  other error, the class shares included, was computed as if they were known
  and came out too small: on poLCA's `gss82` with three classes the share
  errors were 0.0396, 0.0307 and 0.0348 instead of 0.0521, 0.0462 and 0.0348.
  The errors no longer depend on the order of the categories.

### Sampling weights

* A unit of sampling weight zero now changes nothing: the weights are scaled
  to the number of units with a positive weight, and such units are in no
  sample size, criterion, entropy or classification table, in every family
  (`multilpa()`, `lta()`, the additive family, `mixture_regression()` and
  growth mixtures). Adding four zero-weight groups to fifteen used to scale
  the log likelihood by 19/15 and move BIC from 622.99 to 784.74. `nobs()`
  counts units of positive weight, as `nobs.glm()` does; `n_observations` and
  `n_groups` still count the data, and every row keeps its posterior.

### Growth mixtures in any units

* A growth mixture is now fitted in standard coordinates and mapped back
  exactly: the outcome is centred and scaled, and so is every column of the
  class-specific, common and random designs except the intercept. Its
  estimates, log likelihood, standard errors, posteriors, predicted random
  effects and predictions no longer depend on the units or origin of the
  outcome or of the predictors. With the outcome multiplied by 1e-5 a fit
  used to stop 0.06 standard errors short; with time in days instead of
  years a random slope's variance was held at the variance floor, 11.3
  log-likelihood units low, and reported converged; time in hours stopped
  0.12 standard errors short; and a calendar-year random slope stopped with
  "EM likelihood decreased". Outcome scales from 1e-5 to 1e6 (with shifts up
  to 1e8), and time in days, hours, units of 1e5 or calendar years, now give
  the same fit to about 1e-11 standard errors.
* Every growth fit now ends with Newton steps on the analytic score and is
  judged converged when it is within 1e-4 standard errors of the root, as
  fits with correlated residuals already were; the quasi-Newton stop alone
  left fits up to 0.35 standard errors short.
* `min_variance` floors each random effect's conditional variance (given the
  ones before it) in the standard coordinates, in EM and in the final steps
  alike, and a random-effect covariance is judged degenerate there, also by
  the eigenvalues of its correlation matrix. A slope in small time units is
  no longer flagged degenerate.
* The class table's `assigned_persons` and `average_posterior` leave out
  persons of sampling weight zero, like every other count.
* A multilevel growth mixture's ICL charges the entropy of the joint
  classification of clusters and persons. It used to add the persons'
  marginal entropies to the clusters' entropy, which overstated it.
* `parameter_inference()` drops the interval of a residual standard
  deviation only for the class whose variance is held at its floor, not for
  every class.

### Standard errors of mixture regressions

* Standard errors no longer depend on a predictor's units or origin. The
  singularity check and the inverse of the information are taken on its
  correlation form, and its numerical derivatives step by each parameter's own
  standard-error scale. With a calendar-year predictor every standard error
  used to be `NA`, even with one class, and the slope's error moved by 0.2%
  under a shift of origin. Growth mixtures and the latent-transition
  three-step share the corrected inverse.

### Smaller repairs from the review

* `positive_rows` in `get_results(x, "gaussian_diagnostics")` is the number of
  rows with a Gaussian value and a positive weight, as documented. It used to
  count the posteriors that had not underflowed, which for a component with
  SD 0.01 reported 87 of 300 rows.
* Additive-family inference and `bootstrap_lrt()` accept the fitted data with
  its rows reordered within groups; a changed rating is still refused.
* `r3step()` reports `NA` with a `latents_singular_information` warning for a
  coefficient whose variance is not positive, instead of a standard error of 0.
* The "densities vanished" stop is classed `latents_vanished_density`.
* A trajectory three-step check names why a weight column cannot be read.

### Refusals

* `enumerate_classes()`, `enumerate_lpa()` and `enumerate_lca()` raise the
  data error that `multilpa()` raises (`latents_bad_data`,
  `latents_bad_categorical`) when every candidate fails with it: a column that
  does not exist, an indicator that is constant or entirely missing, a text
  indicator. They used to return a grid of failed cells without a condition.
  Failures that differ between candidates still appear in the grid.
* `mixture_regression()` refuses a Gaussian outcome that takes one value in
  every row with `latents_bad_data` and says so. It used to report that every
  start degenerated and suggest fewer classes or more starts.

### Review of 2026-10-26 (numerical and validation defects)

* Exactly constant indicators and covariates are recognised by their range,
  not by a computed standard deviation, which is about 1e-17 for a repeated
  0.1.
* Standard errors are kept for category labels containing `:`, `,`, `[` or `]`.
* A direct-effect covariate constant where its indicator is observed is
  refused; a 17-digit ordinal label no longer merges categories on redraw;
  unobserved joint categories keep their labels.
* EVI/VEI variances of a near-empty profile are no longer floored (a 12%
  bias); growth models in small units no longer stop 380 log-likelihood units
  short; ordinal regressions at large predictor scales converge; the
  negative-binomial dispersion is no longer forced to its floor.
* Homogeneous and occasion-specific `lta()` fits of the same data can be
  compared; person-centred fits accept their data re-read from CSV.

# latents 0.97.1

### Enumeration on several cores

* With `seed` given, `enumerate_classes()` prepares the candidates of a
  profile, latent class or transition grid on the workers. Previously all of
  them were prepared in one process before the workers started; for a large
  grid this took a long time on one core. Each candidate seeds its own starts,
  so results are identical to `cores = 1`. Unseeded grids are still prepared
  in grid order in the calling process.
* `enumerate_classes()` gains a trailing `verbose` argument. The default,
  `cores > 1`, reports the preparation, the fitting of the starts and the
  finished grid in at most four `latents_progress` messages, with the elapsed
  time. Group-class families are not reported.

### Clearer reports

* An enumeration's print separates candidates that failed to fit from fitted
  candidates that did not converge. Before, both were reported as "did not
  converge".
* A start whose k-means starting values stop at 100 iterations is recorded as
  a failed start, with that reason in `get_results(fit, "starts")`. R's own
  k-means warning, which read as though the model had not converged, is no
  longer passed on.
* Printed enumeration tables omit the `model` column when it is empty in every
  row, as in a latent class grid. `get_results()` keeps it.
* `odds_correct_classification` is `Inf` for a class classified perfectly
  (average posterior one to within rounding), instead of `NA` or a number
  such as 7e15.
* The knee documentation cites Satopaa et al. (2011).

### Review corrections

* Plot methods take their 0.9.21 positional arguments again (for example
  `plot(grid, "aic")` for family and transition grids, whose `plot(grid)`
  still draws every criterion as a panel), and `multilpa()` keeps its
  released argument positions, with `boundary` last. The
  enumeration summary and `get_results()` again provide `criteria` beside
  `best`. Released base-graphics style arguments to bootstrap-LRT and pooled
  plots are accepted with a `latents_unused_argument` warning; these plots are
  drawn with ggplot2 and returned.
* Joined local-dependence items keep their full grid of category
  combinations, each with an estimated probability, so unobserved
  combinations are part of the model.
* Direct-effect slopes are solved on the weighted design rather than its
  cross-product, which no longer appears singular when the covariates have
  full rank.
* The homogeneous transition recursion falls back to the log-space recursion
  when subnormal transition probabilities would overflow its scaled messages.
* For a model whose inference provides no per-unit scores, the
  Vuong-Lo-Mendell-Rubin test is reported as unavailable (a classed
  `latents_invalid_lmr` refusal, recorded in the `note` column) instead of
  failing with an unclassed error.
* Parallel workers use the parent's random number generator kinds, are
  stopped when a call ends abnormally, and retry a lost task at most once.
  `options(warn = 2)` is honoured on every worker backend.

# latents 0.97.0

* `enumerate_classes()` gains a trailing `type` argument to request `"lca"`,
  `"lpa"` or `"lta"` explicitly. `"lca"` treats numeric response codes as
  categorical; conflicting routing or measurement arguments are refused.
  `type = NULL` keeps existing fitting behavior. The resolved type is recorded
  on the enumeration and shown in its print header; tables are unchanged.
* A numeric indicator fitted as categorical that has non-integer values or
  more than ten categories now raises a `latents_suspect_categorical`
  warning, once per call: every distinct value is fitted as a category, and
  a continuous indicator may call for a model with Gaussian indicators.

### Faster latent transition models

* Profile and latent class grids in `enumerate_classes()` now share the
  worker pool across candidates and their starts, including both emEM phases.
  Ordered starting values, selection and finalisation stay in the parent,
  preserving seeded fits. Membership covariates and joined categorical
  local-dependence models retain serial starts as candidate tasks; group-class
  families retain candidate-level parallelism.

* Lost parallel worker results are retried once with fewer workers instead of
  discarding the grid. Starts or candidates still lost are recorded as failed,
  and one `latents_worker_lost` warning reports losses and recovery. Socket
  dispatch retains completed chunks; tasks within a broken chunk may be rerun.
* `lta()` gains `cores` for parallel random starts, including emEM
  screening. Transition grids in `enumerate_classes()` dynamically schedule
  individual starts in one worker pool, longest estimated first, without
  nested workers or waiting for candidate batches. Both emEM screening and
  finalist runs use the pool. Seeded results are identical to `cores = 1`.
* `lta()` gains `acceleration`, defaulting to `"squarem"` as in
  `multilpa()`. Fits usually need several times fewer EM steps and can
  change slightly, at the same or a higher maximum.
  `acceleration = "none"` restores the previous fits.
* The forward-backward recursion of homogeneous transition models runs as
  matrix products, several times faster per EM step with identical results.
* `enumerate_classes()` with `time` or a group-class `family` no longer
  fails every candidate on socket workers (Windows, and macOS with
  Accelerate BLAS).

### Default start search

* Models whose profiles share one covariance matrix (EII, EEI, EEE) now
  default to `search = "emEM"`; all others keep `"kmeans"`. Their results
  can change, usually to a higher maximum, and fitting takes about twice
  as long. `search = "kmeans"` restores the previous fits.

### Review corrections

* General latent transition fits detect covariance eigenvalues at the
  variance floor, including occasion-specific measurement. Such fits are
  marked as boundary solutions and refuse Wald inference.
* Zero-weight units no longer determine the internal Gaussian origin,
  imputation, clustering geometry or initial measurement parameters in
  profile, membership-covariate and weighted transition fits. Excluded
  outliers cannot consume an initial profile or erase supported-row
  precision. Observation and group mappings are preserved.
* Gaussian support warnings in ordinary and membership-covariate fits
  exclude categorical-only and retained empty rows. Cross-level models
  apply sampling weights and check both the individual and group Gaussian
  components.
* `get_results(fit, "gaussian_diagnostics")` exposes per-component Gaussian
  membership mass, positive contributing rows, weight effective sample size,
  standardized variances and covariance eigenvalues, and the absolute
  boundary status. The narrowest standardized direction is also printed.
  These are descriptive diagnostics, with no automatic spike cutoff or
  change to start selection. Normalized membership mass is not an
  independent sample count; thin-profile warnings are heuristics, not
  proofs of singularity or nonidentification.
* Bootstrap inference for `fit_staged()` repeats the measurement stage's
  own search controls before refitting the group-class stage, preserving
  an explicit emEM pool even when a single supplied start disables
  screening at stage two.
* `latents_thin_profiles` counts the free Gaussian measurement parameters:
  covariances held at zero and held means are excluded, and held variances
  do not trigger the warning.
* In latent transition models, thin-profile warnings use weighted
  memberships of rows with observed Gaussian indicators. Empty or
  categorical-only rows cannot hide a thin Gaussian profile, even under
  `drop_empty = FALSE`. Occasion-specific measurement is checked separately
  at each occasion, which the warning identifies. Fitted estimates and
  classification counts are unchanged.

### Start search

* `search = "emEM"` (new) finds starts by the emEM search
  (Biernacki, Celeux & Govaert, 2003): a pool of `n_screen` starts (five
  times `n_starts` by default) is screened by a short EM run each, and the
  `n_starts` finalists are run to convergence. A successful first start
  (Ward or supplied) reserves one finalist slot; the others are ranked by
  screening likelihood. With `n_starts = 1`, a valid Ward start is therefore
  always finalized; a single supplied start skips screening. This existing
  protection is now documented. The pool begins with the kmeans
  search's own starts. `search = "kmeans"` reproduces the fits of 0.96.1,
  apart from the corrections listed below.
  On the start-search benchmark (180 fits), emEM attained the finite
  above-cutoff reference in 136 fits and the kmeans search in 130. The
  benchmark's descriptive cutoff was 0.01 for the smallest standardized
  variance or covariance eigenvalue; values below it do not prove false
  classes. On the matched 80 heteroscedastic fits, both searches attained
  that reference in 40 fits, while emEM returned 25 below-cutoff solutions
  against 12 for the kmeans search. Retained runs took about twice as long
  by median time and 1.45 times as long in total.
* Either search reports the highest likelihood among the finished starts.
  The new `latents_thin_profiles` warning fires when a profile with its own
  variances has Gaussian membership mass below its number of free Gaussian
  parameters. It prompts inspection of covariance and support diagnostics.
  The warning is heuristic and the estimate is unchanged.
* Under emEM, `get_results(fit, "starts")` lists every start of the pool,
  with its `kind`, `stage` and likelihood after screening, and "Best
  likelihood replicated" counts the starts that were run to convergence.
* `bootstrap_lrt()` refits each model with the search, starts, stopping
  rule, selection and boundary policy it was fitted with, by default
  (`n_starts`, `max_iter` and `tol` now default to `NULL`). Staged and
  covariate refits keep the pool size too.
* `cores` distributes the emEM screening and the finished starts, with
  identical results.

### Corrections found by the feature equivalence suites

* `three_step(method = "ml")` runs its EM from a fixed set of starts (Latent
  GOLD's proportional start first, the modal means, quantile placements,
  tail starts and twenty draws from a fixed internal seed) and keeps the
  highest maximum. One start missed a higher regular maximum for skewed or
  heavy-tailed outcomes. The result does not depend on, and does not move,
  the caller's random seed.
* `three_step(method = "ml")` reports `NA` standard errors, with a
  `latents_small_classes` warning, for a class carried by a single unit
  (by a single group under `vcov_type = "cluster"`); its sandwich was
  exactly zero.
* Covariate fits record their `max_iter` and `tol`, so a bootstrap refit
  repeats them.
* `direct_effects` values that are `NA`, not character, or end in an empty
  term (`"y3 ~ z +"`) are refused as `latents_bad_argument`; a trailing `+`
  was silently dropped.
* Reported covariance matrices are exactly symmetric (they could differ by
  one unit in the last place between the triangles).
* Bootstrap inference drops, and counts, a resample that misses a category
  of a categorical or joined item, instead of stopping with an internal
  error.

# latents 0.96.1

### Inference and numerical corrections

* Entirely missing groups under `drop_empty = TRUE`, and zero-weight units,
  no longer satisfy robust-inference minimum-unit checks or inflate the
  number of independent units used for the CR1 correction. The correction
  is consistent across profile, covariate, transition, additive, regression
  and growth models where those units are supported.
* Direct-effect information is evaluated with centered covariates and its
  full covariance transformed back to the reported origin. Translating a
  predictor no longer causes a false singular-information error.
* ML three-step outcomes are fitted and differentiated in centered, scaled
  coordinates, preserving standard errors when the outcome has a large
  origin relative to its spread.
* VLMR calculations retain well-conditioned information coordinates and
  use an equivalent symmetric eigenvalue calculation. Imhof probabilities
  handle negative statistics and retain genuine small distribution weights.
  Comparisons whose reference weights cannot be resolved numerically report
  unavailable probabilities.
* Structured covariance updates handle active variance floors with a
  feasible solver that checks optimality, avoiding numerical failures near
  a bound while preserving the covariance structure and EM monotonicity.
* Socket-worker cleanup closes surviving and failed connections separately,
  preserving the typed parallel-failure error when a worker dies.

### Paired indicators and weighted trees

* Pair preparation rejects generated names that conflict with existing raw
  columns and applies native categorical validation before joining labels.
  Distinct numeric categories cannot silently merge during pair creation.
* Regression and growth trees use the sampling weights of the units their
  posteriors classify, including row, group and person adapters.
* The independent mover-stayer equivalence reference uses stable probability
  normalization, so large initial logits no longer overflow its path sum.

# latents 0.96.0

### Direct effects

* `multilpa(direct_effects = "y ~ z")` lets a covariate move a continuous
  indicator's mean by the same amount in every profile (differential item
  functioning), Latent GOLD's `y <- 1 + Cluster + z`. The effects are
  estimated with the profiles, by a conditional maximization step in EM on
  centred covariates, have standard errors in `parameter_inference()`, and are
  reported by `get_results(fit, "direct_effects")`. `predict()`, refits, the
  bivariate residuals and parametric draws use the covariates. Against Latent
  GOLD 6.1 the log likelihood, BIC, effect, intercepts, variances and their
  standard errors agree to every printed digit. Effects on categorical
  indicators and effects that differ by profile are not yet available.

### Compatibility

* Every call that worked in 0.95.12 gives the same result: the new arguments
  (`local_dependence`, `direct_effects`, `three_step(variance = )`) come last
  with defaults that change nothing, and `three_step(method = "m")` still
  means `"modal"`; `"ml"` is taken by its full name.

### Local dependence

* `local_dependence = "item1:item2"` in `lca()`, `multilca()`, `multilpa()`
  and the enumeration functions models a pair of categorical indicators
  jointly within each class, the class-specific association Latent GOLD fits
  with `item1 <-> item2 | Cluster`; the other indicators stay independent
  given the class. `get_results(fit, "local_dependence")` gives each class's
  joint probabilities beside their values under independence. On simulated
  data the log likelihood, parameter count, class sizes and item marginals
  agree with Latent GOLD 6.1 to every printed digit.
* A pair of continuous indicators, `local_dependence = "y1:y2"`, gets its own
  covariance within each class while every other covariance stays zero:
  between local independence (VVI) and a full covariance (VVV). The log
  likelihood, parameter count, covariances and their standard errors agree
  with Latent GOLD 6.1 (`y1 <-> y2 | Cluster`) to every printed digit.
* `compare_models()` compares a model with paired items and one without, which
  are fitted to the same data.

### Lo-Mendell-Rubin tests

* `get_results(models, "lmr")` reports the Vuong-Lo-Mendell-Rubin (VLMR) and
  adjusted Lo-Mendell-Rubin (aLMR) tests for every pair of converged
  candidates of a grid that differ by one profile or one group class, as
  Latent GOLD computes them: the likelihood-ratio statistic referred to a
  weighted chi-square distribution whose weights are the eigenvalues of
  Vuong's matrix, with p-values by Imhof's method. On two test cases the
  statistics, the mean and standard deviation of the reference distribution
  and the p-values agree with Latent GOLD 6.1 to every printed digit.

### Maximum-likelihood three-step distal outcomes

* `three_step(method = "ml")` estimates the class means of a distal outcome
  by maximum likelihood, with the step-one classification errors held fixed
  (Vermunt 2010; Bakk, Tekle and Vermunt 2013). `variance = "varying"`
  (default) or `"equal"` sets the outcome's variances; the table adds the
  estimated `variance` of each class. Standard errors are the clustered
  sandwich, and `contrast = "pairs"` compares the classes, as for the other
  methods. BCH remains the default, being more robust to non-normal
  outcomes (Bakk and Vermunt 2016).

# latents 0.95.12

### Quieter and stricter checks

* A covariate model no longer passes on the estimation warnings of the
  one-start base fit it uses to check its input (a variance at its bound, no
  convergence, failed starts): it checks its own estimate for each of them,
  so one boundary gave two warnings and now gives one, about the model
  returned. It now also warns, with `latents_small_classes`, when one of its
  own classes has effective membership below one.
* A named plot function asked for one view (`plot_profiles(fit, what =
  "heatmap", scale = "standardized")`) warns with `latents_unused_argument`
  about an argument that view does not use, instead of ignoring it without
  a word. `what = "all"` stays silent, as each view takes its share.

### Faster standard errors

* `parameter_inference()` computes the observed information of a profile,
  class, covariate or transition model once per R session for that fit and
  its data, and reuses it. The first call takes as long as before; later
  calls, including those made by `plot_profiles()`, `plot_responses()` and
  `get_results()` when they draw or tabulate intervals, skip the two score
  evaluations per parameter. On a two-level model of 20,000 observations
  with 94 parameters, a repeated call took 0.23 s instead of 13 s. Results
  are identical to a computation without reuse, and any change to the fit or
  its data computes the information anew.

# latents 0.95.11

### Named plot functions on every fitted model

* `plot_profiles()`, `plot_diagnostics()`, `plot_prevalence()` and
  `plot_transitions()` take every fitted model of the package and draw the
  view of that fit's own `plot()` that shows the same quantity:
  `plot_profiles()` draws a growth mixture's `"trajectories"` (also of a
  mixture regression with `class_level = "group"`), another mixture
  regression's `"fitted"` lines, an additive or dispersion fit's `"means"`,
  a cross-level fit's `"profiles"` and a pooled fit's `"measurement"`
  parameters; `plot_diagnostics(what = "posteriors")` draws the
  `"posteriors"` of growth mixtures and mixture regressions;
  `plot_prevalence()` draws a cross-level fit's `"composition"` and a pooled
  fit's `"profile"` parameters (each profile's probability within each group
  class);
  `plot_transitions()` draws a general `lta()` fit's transition matrix, as
  before. `what = "all"` draws these and names the views it leaves out.
* Any other combination (for example `plot_profiles(what = "bars")` on a
  growth mixture, or `plot_diagnostics()` on an additive fit) is refused with
  `latents_nothing_to_plot`, which is also of class `latents_bad_argument`
  (the class of the former refusal), and the message names the views the
  fit's `plot()` offers and the named functions that draw them. An argument
  the mapped view does not take (for example `intervals` for a trajectory)
  is refused by name with `latents_bad_argument`.
* For profile, class and transition models nothing changes, except that
  `plot_transitions()` on a model without transitions, and a named function
  whose view the fit's `plot()` lacks (for example `plot_prevalence()` on a
  transition fit), now raise `latents_nothing_to_plot` with
  `latents_bad_argument`, where they raised one of the two.

### One contract for every grid

* `plot_enumeration()`, `plot_selection()` and `plot_tree()` take every grid
  the package makes: profile and class grids, transition grids
  (`enumerate_classes(time = )`), group-class family grids
  (`enumerate_classes(family = )`) and regression grids
  (`enumerate_regressions()`). Each grid offers the views its candidates
  carry. Transition grids draw every view, the tree included. Family grids
  draw their criteria and classification quality against the number of group
  classes, and a tree of their group classes for each family; their
  candidates have no profiles, so the selection views are not drawn.
  Regression grids draw every view, the tree of their classes included.
  `get_results(x, "tree", ...)` returns the shared masses on every grid.
* A view a grid does not offer is refused with the new condition
  `latents_unsupported_view` (also of class `latents_bad_argument`), whose
  message gives the reason and the views the grid offers. An unknown
  criterion of a grid is now also of class `latents_bad_argument`.
* `plot_views(x)` of a grid lists every view of that grid, whether it is
  offered, and why not.
* `plot()` of a transition or family grid draws through the same views.
  Without `what` it draws the grid's summary, AIC, BIC and BIC over
  individuals as panels of one plot, as `plot_enumeration()` does (0.9.21
  drew BIC alone). The 0.9.21 calls that name a criterion still work:
  `plot(x, what = "aic")`, abbreviations such as `"bic_i"`, and
  `criterion =`. It also takes `"all"` and every other view of the grid, and
  passes `combine`, `labels`, `mark_minimum`, `cutoff` and `using` on to
  `plot_enumeration()`, which 0.9.21 ignored. The plots have the look of the
  other grid plots: a criterion in its own panel, the lowest value ringed.
* Regression grids get `plot()`, which draws AIC, BIC, SABIC and ICL across
  the classes, and `summary()`, which returns their table (as
  `as.data.frame()` does). Formerly `plot()` failed with a base R error and
  `summary()` gave the default list summary.
* The selection views of a class (LCA) grid name the singular "class"
  instead of "classe".

### Criterion panels

* A grid whose every candidate has one group class draws each criterion once,
  under its plain name: the group versions of AWE, ICL and CLC are left out
  when their individual versions are drawn, since without group classes ICL
  over groups equals BIC over groups. `what = "all"` on such a grid draws AIC,
  KIC, BIC, SABIC, CAIC, AWE, ICL and CLC.


# latents 0.95.10

### Tree plots

* `plot_tree()`, `plot(x, what = "tree")` and
  `plot_enumeration(x, what = "tree")` return one plot per tree, printed one
  at a time on an interactive screen device; `combine = TRUE` returns every
  tree as a panel of one plot, the figure of earlier versions. A class has
  the colour of the class it matches in the first tree, and the classes of
  one row have distinct colours. Every band that carries posterior mass is
  plotted.
* New arguments choose the trees and rows: `n_group_classes`, `n_profiles`
  and `model`. A single value of `n_profiles` compares the solutions with
  that number of profiles across numbers of group classes. `level = "groups"`
  plots how the group classes divide. `from = c(3, 3), to = c(4, 3)` (or
  `candidates` by name) compares two chosen solutions.
* `get_results(x, "tree", ...)` returns the shared posterior mass behind any
  tree, one row per pair of classes, with the same arguments.
* Trees use the rows and groups the fit counts (`drop_empty`) and their
  sampling weights at both levels, as the effective class counts do. A fit
  with a noise component is refused with `latents_unsupported_noise`, and a
  named candidate that did not converge with `latents_unconverged`. A grid
  with no tree raises `latents_no_tree`, also of class
  `latents_nothing_to_plot`.
* `plot_enumeration()` refuses an argument it does not use with
  `latents_bad_argument`.

### Documentation

* New guide, *Tree plot*.

# latents 0.95.9

### Documentation

* The package website links each author's name to their own website
  (https://saqr.me, https://sonsoles.me), as the sibling packages' sites do.

# latents 0.95.8

### Documentation

* The vignettes, the guides, the case studies and the README use the named
  plot functions: `plot_profiles()` (with `what = "bars"`, `"heatmap"`,
  `"raincloud"`), `plot_responses()`, `plot_prevalence()`,
  `plot_sequences()`, `plot_transitions()`, `plot_diagnostics()`,
  `plot_enumeration()`, `plot_selection()` and `plot_tree()`, in place of
  `plot(x, what = )`. The function tables of the enumeration guides list
  them. The model families without named functions (growth mixtures,
  trajectory classes, mixture regression, additive and cross-level models)
  keep `plot(x, what = )`.
* The LPA and enumeration vignettes mention `cores`.
* Every help page is revised: titles and descriptions now say what each
  function does and returns, without design rationale or history, and every
  argument states its default. The revision corrected statements that did
  not match the code, among them the refusals listed by `enumerate_lca()`,
  `enumerate_classes()`, `bootstrap_lrt()`, `compare_models()`,
  `enumerate_regressions()`, `prior_control()` and `fit_staged()`; the
  columns of `starting_values()` under full covariance; the tables offered by
  `get_results()` (starting values, ordinal and count tables, diagnostics);
  when `latents_no_group_classes` and `latents_unsupported_inference` are
  raised; and the default `vcov_type` of weighted fits.

# latents 0.95.7

### Bug fixes

* Line plots of profiles and response probabilities (`plot_profiles()`,
  `plot_responses()`) no longer run long indicator names together: names
  that cannot sit side by side are turned. The profile labels at the right
  end get a margin sized to the longest label, so they are no longer cut off
  when there are many indicators, and they are spaced on the axis as drawn,
  so lines that end close together keep their labels apart.

# latents 0.95.6

### New features

* The named plot functions now cover every view of a fit, grouped by
  subject so that typing `plot_` lists them all, and each group's `what`
  takes `"all"` to draw every view it has:
  * `plot_profiles(x, what = )` adds `"bars"`, `"heatmap"`, `"raincloud"`,
    `"parallel"` and `"pairs"` to the lines (`what` is its last argument, so
    existing calls are unchanged);
  * new `plot_sequences(x)`;
  * new `plot_selection(x, what = )` for an enumeration: `"deltas"`,
    `"heatmap"`, `"smallest"`, `"trust"`, `"diagnostics"`, or `"all"`
    (default), the same as `plot_enumeration(x, what = "selection")`.

  Each draws exactly what `plot(x, what = )` draws.

### Bug fixes

* `plot(x, what = "all")` and `report()` left out the heatmap of a model
  with only categorical indicators, which draws its response probabilities.

# latents 0.95.3

### New features

* `multilpa()`, and so `multilca()`, `lpa()` and `lca()`, take `cores` to run
  the random starts of one fit in parallel. The starting values are drawn
  first, in order, so the fit is identical to `cores = 1`. Models with
  membership covariates still run their starts on one core.

### Bug fixes

* `multilca()` and `enumerate_lca()` no longer refuse `categorical`, so a call written for
  `enumerate_classes()` or `multilpa()` (`categorical = vars`) works as is.
  Every indicator in `vars` stays categorical whatever `categorical` names;
  a name that is not in `vars` is still refused.

* `plot(what = "heatmap")` on polytomous categorical indicators now draws one
  panel per profile, with indicators as rows and response categories as
  columns. The former layout, a block of columns per indicator, ran to
  indicators x categories columns, so with many indicators the numbers
  overlapped and the indicator names were cut off. Binary indicators keep
  their one-column-per-indicator layout.

### Documentation

* The model-selection and class-enumeration guides, the enumeration and
  evaluation vignettes, the LPA workflow, the README and the help pages of
  `lpa()`, `enumerate_lpa()` and `enumerate_classes()` describe the
  degenerate one-observation solution and `boundary = "exclude"`.

# latents 0.95.2

Fixes found by Equivalence 2.0, the strengthened external validation.

### New features

* New argument `boundary = c("allow", "exclude")` in `lpa()`, `multilpa()`,
  `enumerate_lpa()` and the covariate models. `"exclude"` sets aside starts
  with a variance at `min_variance`: solutions in which a Gaussian profile
  shrinks around one observation, whose likelihood gain is set by
  `min_variance` rather than by the data. The more starts are run, the
  likelier one finds such a solution: in a simulation study (29 datasets), a
  4-profile dataset fitted with 500 starts gave a 5-profile solution with a
  one-observation profile that made BIC choose 5; with `boundary = "exclude"`
  BIC chose 4, as Mplus did. Opt-in: the default (`"allow"`) and every
  existing result are unchanged. When every start is at the bound, the best
  is reported with a `latents_all_starts_boundary` warning. Categorical
  indicators alone are unaffected.

### Bug fixes

* `lta()` beyond homogeneous transitions (`transitions = "occasion"`,
  covariates, second order, mover-stayer) stopped its quasi-Newton finish
  early near a transition estimated at zero: since 0.9.13 the finish's
  stopping rule was tied to `tol`, and at the default `tol = 1e-10` an
  occasion-specific model of LMest's reference was left 2.8e-6 below its
  maximum with transitions 1e-4 off. The finish again runs to L-BFGS-B's
  tight rule, as in 0.9.1, and agrees with LMest to 1.7e-6 or better. Fits
  move closer to their maximum; nothing else changes.
* `r3step()` could stop with `latents_no_converge` when a covariate nearly
  separates the classes: the likelihood with the classification error held
  fixed is not concave, and BFGS from zero settled on a separation ridge
  (Mplus UG example 7.3: -log L 81.7 on the ridge, 62.6 at the maximum). The
  search now also starts from the proportional-assignment estimate, which is
  concave and lies in the maximum's basin, and keeps the better converged
  end; it reproduces Mplus's R3STEP estimates.
* `enumerate_classes(cores = )` on macOS with R switched to Apple's
  Accelerate (vecLib) BLAS: forked workers died with a segmentation fault
  once the session had used linear algebra, and the enumeration then stopped
  with an internal error at printing. Accelerate is not safe in a forked
  process, so there the candidates are now fitted by a socket cluster, as on
  Windows; results are identical to `cores = 1`, as before. Elsewhere the
  fork is kept.
* A worker killed by a signal (a crash, out of memory) returned an empty
  candidate that was not recognized as lost; it now raises
  `latents_parallel_failure`, as a worker that fails with an error does, and
  so does a socket worker that dies.
* The condition catalogue lists `latents_parallel_failure`, raised by
  `enumerate_classes(cores = )` when a worker process fails.

### Equivalence

* The external validation now includes 160 simulated datasets with known
  parameters fitted by Mplus and Latent GOLD (likelihood, parameters,
  criteria, posteriors, entropy, AvePP, OCC, class counts, group level,
  `predict()`, enumeration grids, `compare_models()`, `cores`), 42 genuine
  Mplus runs kept in JStats, `scaling` against Mplus, sampling weights in a
  mixture against Mplus WEIGHT with robust standard errors, `bootstrap_lrt()`
  against TECH14, and the group-class families at Latent GOLD's finest
  quadrature. Comparisons that were scripts now run as tests, and the Latent
  GOLD kit is compared with the current package rather than with targets
  saved when it was built. See `equivalence/README.md`.
* Every documented claim of agreement with another program or a published
  result is registered in `equivalence/CLAIMS.csv` with the test that backs
  it, and a test fails when a claim or its test changes. Two claims gained
  tests (`bic_rows` equals flexmix's BIC; a noise observation is profile 0,
  as in mclust's classification); one was narrowed to what is tested (plain
  EM follows mclust's path from the same start).
* A sabotage battery (`equivalence/sabotage/run.R`) injects realistic
  defects into the loaded package and confirms that the tests detect them.

### Documentation

* The class-enumeration guide is rewritten in scientific register, with the
  methods cited.

# latents 0.95.0

### Parallel enumeration

* `enumerate_classes()`, and through it `enumerate_lpa()` and
  `enumerate_lca()`, gain `cores = 1L`, which fits the candidates on several
  processor cores: forked processes on macOS and Linux, a socket cluster on
  Windows. Every candidate sets and restores its own seed, so with `seed`
  given the result is identical to `cores = 1` for profile, transition
  (`time`) and group-class family grids; only the time differs (a 24-candidate
  grid of `course_engagement` took 11.7 s on one core, 4.2 s on four and
  3.3 s on eight). Candidates' warnings and messages are passed on in grid
  order, with their condition classes. `cores = 1`, the default, runs as
  before. `parallel`, part of R, is now in Imports.

# latents 0.9.41

### `drop_empty`

* `multilpa()` and `lta()` gain `drop_empty = TRUE`, at the end of their
  signatures; enumerations, `lpa()`, `lca()` and covariate fits pass it on
  through `...`. It governs what the reported statistics do with rows whose
  every indicator is missing and groups whose every row is such a row. The
  estimates and log likelihood are the same either way.
  * `TRUE` (default): such rows and groups are left out of the sample sizes
    of the information criteria, of entropy, AvePP and the entropy charge of
    ICL, AWE and CLC, of the classification tables and plots, and now also of
    the effective class counts, as Mplus leaves them out. Class counts now
    agree with Mplus's counts based on posterior probabilities.
  * `FALSE`: every number is what latents 0.9.37 reported, so results
    published with earlier versions can be reproduced (checked against
    0.9.37 for enumeration grids, covariate fits, both `lta()` engines and
    `compare_models()`).
* Fits with different `drop_empty` are refused by `compare_models()` and
  `bootstrap_lrt()` with `latents_incomparable_models`.

# latents 0.9.40

### Groups with nothing observed

* A group whose every row has every indicator missing is now left out of the
  group sample size (`bic_groups` and the other group-convention criteria)
  and of group-level entropy, AvePP, the classification and entropy tables and
  the group diagnostic plots, as Mplus leaves it out of its clusters. Its
  posterior is the group-class prior, not a classification. A group with at
  least one observed row is unaffected; group class sizes and assignments
  keep every group. Fits record `n_informative_groups`. In a single-level
  fit each row is its own group, so `bic_groups` now equals
  `bic_individual` when rows are empty. Fits without such groups report
  exactly what they did.

# latents 0.9.39

Fixes found by an independent review of 0.9.38, and new references for
missing indicators.

### Bug fixes

* `compare_models()` on a general `lta()` fit kept rows with nothing observed
  in the entropy charge of ICL, AWE and CLC; it now uses the rows the fit
  classified, as every other fit does.
* `predict()` failed on a fit with no continuous indicator when it was fitted
  with `scaling`, or with `centering = "grand"`; such a fit is neither centred
  nor scaled, and now predicts.
* `scaling = "person"` judged whether a unit's standard deviation was zero
  against the largest value of every unit, so a real but small spread beside
  much larger values was refused, as was a unit with no observed value of an
  indicator. Each unit is now judged on its own values, and a unit with
  nothing observed is left unscaled.

### Equivalence

* Five datasets with different missing-value patterns (missing completely at
  random, rows with nothing observed, missing at random on an observed
  indicator, one indicator mostly missing, and a two-level set with empty
  groups), fitted by genuine Mplus 9 and Latent GOLD 6.1. Log likelihood,
  sample size, BIC, entropy, AvePP and posteriors agree; the programs'
  differing treatment of rows and groups with nothing observed is recorded.

### Documentation

* The model-selection guide describes the choices made before fitting
  (random starts, centring and scaling, missing values), why entropy and AvePP
  can disagree, the two levels of classification, and refitting with more
  starts. The class-enumeration guide points to them.

# latents 0.9.38

Scaling of the indicators, group-level diagnostic plots, and classification
statistics that leave out rows with nothing observed. Every addition is a new
argument at the end of a signature or a new value of `what`.

### Classification quality with missing indicators (bug fix)

* Under `missing = "fiml"`, a row with every indicator missing was counted in
  relative entropy, AvePP, the classification and entropy tables, the
  entropy charge of ICL, AWE and CLC, and the case-level diagnostic plots.
  Such a row contributes nothing to the likelihood, so its posterior is the
  prior its group gives it, not a classification; counting it pulled every
  one of those statistics down (in one example relative entropy from 0.78 to
  0.63). They are now taken over the rows with at least one observed
  indicator, the same rows the individual-level BIC already counted
  (`n_informative`). Estimates, log likelihoods and class sizes are
  unchanged; fits without such rows report exactly what they did. The rows
  used are kept on the fit as `informative`.
* The same rows are now left out of the individual-level BIC of `lta()` and
  of covariate fits, which counted every row; `multilpa()` already did. A
  covariate fit's `bic_individual` therefore agrees with its own
  `get_results(x, "information_criteria")`, which it did not when rows were
  empty. The additive and cross-level families require complete data and
  are unaffected.
* The classification-quality and criteria plots of an enumeration labelled
  each line at its end in every panel, which ran off the edge beside a
  second panel and overlapped where values were close. With several panels
  they now share one legend, keyed by colour and shape; a single panel keeps
  its direct labels, spaced further apart.

### Scaling the indicators (`scaling`)

* `multilpa()` gains `scaling = c("none", "person", "grand")`, at the end of
  its signature, beside `centering`. `"grand"` divides each continuous
  indicator by its overall standard deviation and `"person"` by each unit's
  own, so `centering = "grand", scaling = "grand"` fits grand z-scores and
  `centering = "person", scaling = "person"` within-person z-scores. It
  reaches `lpa()`, the covariate fits and every enumeration through `...`.
  The scales are kept on the fit: `get_results(x, "data")` returns the
  columns supplied, `predict()` scales new data the same way, and the
  bootstrap and staged fits treat a scaled fit as they treat a centred one.
  With structures that give every indicator its own variance, grand scaling
  changes only the units (the log likelihood moves by `n * sum(log(sd))`,
  memberships are identical); with `EII` and `VII` it changes the solution.
  An undefined or zero standard deviation is refused with `latents_bad_data`.
* The print of a grand-centred fit no longer calls its profiles "profiles of
  change"; that describes person centring only.

### Group-level classification plots

* `plot_diagnostics()` draws the group classes as well as the profiles. The
  new views `"group_sizes"`, `"group_posteriors"`, `"group_entropy"` and
  `"group_avepp"` are the four individual views drawn for the groups, and
  `what = "all"` includes them when a fit has two or more group classes. The
  individual views are unchanged. A fit with one group class refuses the
  group views with `latents_no_group_classes`.


### Documentation

* A new guide, *Model selection: the evidence for the number of classes*
  (`vignettes/articles/model-selection.Rmd`), describes every piece of the
  enumeration evidence: the criteria and their definitions, the two
  sample-size conventions, deltas, the knee rule and its stability notes,
  classification quality, the estimation checks, comparisons and the
  bootstrap test.
* `enumerate_classes()` documents the arguments its `...` passes to every
  candidate, leading with `n_starts`.
* The function tables of the enumeration vignette and guide name the minima
  table `"best"`.

# latents 0.9.37

Growth mixtures and latent transition models; nothing in the 0.9.21 interface
changes, and every addition is a new argument at the end of a signature, a new
`type`, or a verb that used to refuse.

### Enumeration grids

* The table of each criterion's minimum is now `get_results(x, "best")`, the
  name the transition and family grids already use. `"criteria"` read as
  "every candidate's criteria", which is `"candidates"`. `"criteria"` still
  returns the same table, so existing scripts are unaffected.

### Growth mixture models (`mixture_regression(random = )`)

* `predict()` works on growth mixtures (it used to refuse). The persons in
  `newdata` are new persons: `type = "response"` and `"class_response"` give
  class trajectories, which need no outcome; given the outcomes,
  `"posterior"` classifies each person, `"class_conditional"` and
  `"conditional"` give each person's own curve with their empirical Bayes
  random effects, and `"random_effects"` gives those effects by class. At
  lcmm's own estimates they equal `lcmm::predictY()` (marginal and
  subject-specific), `predictClass()` and `predictRE()` (three
  configurations; `equivalence/test-growth-predict-lcmm.R`). A multilevel
  growth mixture classifies the clusters in `newdata` as new clusters.
* `parameter_inference()` has a growth-mixture method: one row per natural
  parameter (coefficients, membership logits, residual and random-effect
  standard deviations, variances and correlations, class shares) by Wald, or
  by a person bootstrap (`method = "bootstrap"`) whose replicates are matched
  to the fitted classes by their posterior memberships of the fitted persons.
  The Wald rows reproduce the fit's own tables.
* `r3step()` works on trajectory models (latent class growth and growth
  mixtures; it used to refuse): person covariates predicting class
  membership, with the classification errors held fixed.
* `sensitivity()` refits a growth mixture under each seed and reports how many
  persons keep their class.
* Correlated residuals (new arguments `residual_correlation = "ar"` or
  `"brownian"` and `residual_time`): a continuous-time autoregressive process
  or Brownian motion, shared by the classes, on top of the random effects and
  the residual variance, as `lcmm::hlme(cor = AR(time))` and `cor = BM(time)`.
  At hlme's estimates the likelihood equals hlme's to 1e-8 (three
  configurations; `equivalence/test-growth-correlated-lcmm.R`). A process
  the data cannot separate from the residual variance or the random effects
  is left on its bound with a `latents_random_boundary` warning, and its
  rows have no standard error; the rest of the fit is inferred with those
  coordinates held.

### Latent transition models (`lta()`)

* `r3step(fit, data, initial_covariates =, transition_covariates =)` on an
  `lta()` fit is the three-step latent Markov model (Di Mari, Oberski and
  Vermunt 2016): covariates of the first profile and of the transitions, with
  the classification errors held fixed. It agrees with Latent GOLD 6.1's
  `step3 modal ml` latent Markov model to its printed precision (four
  decimals) on three designs: log likelihoods, coefficients and standard
  errors (`equivalence/latentgold-step3-markov/`).
* Wald standard errors for the extended `lta()` model (covariate,
  occasion-varying or second-order transitions) cover every covariance
  structure (`model =`) and `covariance_model = "full"`; they were limited to
  EEI and VVI. The analytic scores agree with numerical gradients of the
  likelihood to about 2e-7.
* `sensitivity()` takes `lta()` fits (with `data`) and reports how many
  occasions keep their profile under each seed.
* `plot(fit, what = "bars")`, `"profiles"` and `"responses"` of a transition
  fit draw the 95% intervals of `parameter_inference()` (new argument
  `intervals`, default `TRUE`); they used to draw point estimates only.

# latents 0.9.36

* `plot_enumeration()` of a grid over group classes at one number of
  profiles puts the group classes on the x-axis; it used to stack every
  candidate in a single column with overlapping labels. Grids over profiles
  plot as before.

# latents 0.9.35

* A posterior probability that rounds a hair past 1 no longer gives a
  negative `uncertainty` in `get_results(x, "assignments")`, and no longer
  drops that case, with a warning, from the classification-certainty and
  entropy plots of `diagnostics()` and `plot()`. The estimates are unchanged.
* Base R's "empty cluster: try a better set of initial centers" warning from
  the k-means of a random start is no longer passed on; the fit already
  reports failed starts in its own warning and in `get_results(x, "starts")`.

# latents 0.9.34

* latents is on CRAN (0.9.21, published 2026-10-16). The README carries the CRAN
  badge and opens with a list of every vignette, guide and case study.

# latents 0.9.33

* `plot()` of an occasion-varying `lta()` fit orders its panels by occasion
  number (with ten or more occasions, "occasion 10" used to come before
  "occasion 2") and puts the occasion on a second line of the panel title,
  which narrow panels used to cut off.
* Top legends with more than four keys wrap into rows instead of running off
  the right edge (for example four group classes and the shared key of a
  group-class variances plot).

# latents 0.9.32

A function-by-function review of the package and its documentation. The
maths it re-derived (analytic scores against numerical gradients for every
model family, delta-method Jacobians, information criteria, entropy, AvePP,
OCC, the knee) held; the corrections below are to edge cases, input
validation, consistency between tables, plots and warnings, and documentation.

## Changed results

* Group-class families (`family = "additive"`, `"dispersion"`,
  `"additive_dispersion"`) use SQUAREM-accelerated EM. Converged fits are
  unchanged to rounding (8 Latent GOLD cases to 7.5e-16 relative); fits that
  plain EM left crawling now converge within `max_iter`, so bootstrap
  replicates that used to be withheld are kept.
* The membership-covariate separation check, the one-step membership
  regression, `r3step()` (pooled and `by_group_class = TRUE`) and covariate-fit
  `parameter_inference()` now work in orthonormal design coordinates, so
  results no longer depend on where a covariate's zero lies: a covariate with
  mean 100 and sd 1 is no longer mistaken for separation, and a covariate
  shifted by 1e6 gives the same estimates and standard errors as the
  centred one. `r3step()` also converges more tightly (about 1e-4 closer to
  its maximum).
* With sampling weights, `small_classes` and the `latents_small_classes`
  warning use the weighted effective counts the fit reports; group-class
  families report weighted class counts and compute effective groups (and
  plot intervals) from the robust covariance they use for inference.
* In weighted general `lta()` fits, the model-implied `"transitions"`,
  `"initial"` and `"second_order_transitions"` tables average groups with
  their sampling weights.
* `lta()` and its `parameter_inference()` refits default to `max_iter = 5000`
  (was `1000`). The tighter `tol = 1e-10` of 0.9.31 needs more updates:
  three- and four-profile models of two `course_engagement` indicators took
  1,000 to 1,500 and stopped unconverged. Fits that converged before are
  unchanged.
* Transition and family grids (`enumerate_classes(time = )`,
  `enumerate_classes(family = )`) pass their candidates' warnings on, as
  profile grids do; they used to keep them only in the `warnings` column.
* `get_results(grid, "best")` on transition and family grids considers only
  converged candidates; when none converged, `"best"` raises
  `latents_failed_candidate`.
* Response-probability whiskers in plots are the logit-scale intervals of
  `parameter_inference()`, so plot and table agree.

## Refusals instead of silent wrong answers

* `lta()` refuses (`latents_bad_data`) transition covariates when
  `occasions = "grid"` leaves a position skipped inside a group's span;
  the covariate used to be taken as 0 there.
* `multilpa()`, `lta()` and covariate fits refuse (`latents_bad_time`) a
  `time` column of numbers stored as text whose alphabetical order differs
  from their numeric order ("10" before "2"). Zero-padded numbers and ISO
  dates are still accepted.
* `fit_staged(measurement = )` refuses (`latents_bad_stage`) a measurement
  fitted with weights, centring or ordinal/count indicators instead of
  holding it under a different likelihood.
* `mixture_regression(family = "ordinal")` refuses (`latents_not_identified`)
  a single-level mixture whose predictors take too few distinct values to
  identify it.
* An ordinal indicator with fewer than two observed categories is refused
  (`latents_bad_data`) even when its factor declares more levels; a declared
  but never-observed level warns (`latents_unused_category`).
* Categorical columns whose distinct values print alike are refused instead of
  getting duplicate labels. `prior_control()` refuses a `scale` that is not a
  square, symmetric, positive-definite matrix.
* `three_step()` refuses `level = "groups"` for trajectory models instead of
  ignoring it. `vcov(fit, scale = "unconstrained", method = "bootstrap")`
  raises `latents_unsupported_inference` instead of returning `NULL`.
* `predict()` on membership-covariate and latent transition fits raises
  `latents_unsupported_prediction`.

## Fixes

* `lpa()` and `lca()` record the call as written, so `update()` refits from
  the caller's data.
* `centering = "grand"` or `"person"` (with or without covariates) no longer
  refuses indicators recorded in very small units.
* `mixture_regression()`: factor predictors can be shared across classes
  through `common`. Regression and growth fits no longer store their fitting
  frame, so saved fits are several times smaller; enumerated growth
  candidates record a call that names the data.
* Growth mixtures: `plot(what = "random")` draws when a class's random-effect
  covariance is degenerate and is identical on every draw; trajectories work
  for intercept-only class formulas; `confint()` honours `parm`; failed starts
  warn (`latents_failed_starts`); `plot(what = "classification")` names each
  class on the bars instead of relying on colour alone.
* `get_results(fit, "recovery")` reports a missing known class as `NA`, and on
  a membership-covariate fit compares group-level truths with group classes.
* `diagnostics()` and `report()` work on fits with a noise component.
  `descriptives()` of a single-level fit no longer reports `n_groups`/`icc`
  for the internal row numbering.
* `print(summary())` names single-level and latent class fits correctly;
  `rows = 0` prints headings only; empty `"profiles"`/`"responses"` tables keep
  their standard-error columns.
* General `lta()` fits gain a `"responses"` table for categorical indicators,
  flag an emptied stayer class as a small class, and no longer warn spuriously
  for a dotted covariate name. Cross-level fits store `profile_prevalences`,
  reported as a `prevalence` column of `get_results(fit, "profiles")`.
* Plots: `plot_profiles()`, `plot_responses()`, `plot_prevalence()`,
  `plot_diagnostics()`, `plot_enumeration()` and `plot_tree()` refuse objects
  they cannot draw with `latents_bad_argument` and name the right plot; a
  noise fit's prevalence plot gives the noise share; profile views name the
  indicators they leave out; `plot()` of a general `lta()` fit uses the
  package's probability scale; the group-class variances plot is keyed by
  class, panelled by level.
* Capped VEI covariance M-steps are reported once per fit with a count.
* Methods that ignored extra arguments (`get_results()` on general transition
  fits and transition or family grids, `vcov()` on transition fits) now warn
  (`latents_unused_argument`); `vcov(fit, vcov_type = "robust")` on a general
  transition fit used to return the observed covariance silently.
* Classed conditions replace unclassed errors for invalid seeds, controls,
  `parm`, `start`, `n_group_classes`, `lca(n_profiles = )` and covariate-model
  inputs; `sensitivity()` and `enumerate_classes()` validate their controls
  and seed before fitting.
* Robust standard errors (CR1 since 0.9.31) are exactly `sqrt(G / (G - 1))`
  times those of Mplus's MLR, which applies no small-sample factor: 0.8%
  larger with 60 groups. `?parameter_inference` says so, and the external
  Mplus comparison now checks that relationship to 1e-6.
* `compare_models()` no longer prints "<- lowest BIC" and its plot no longer
  circles each criterion's lowest value: `delta_bic` and `bic_weight` carry
  the evidence and the choice stays with the reader.
* Plots keep the Okabe-Ito grey for categories that are not classes (a noise
  component, terms shared by every class); classes cycle through the other
  colours, so the two never coincide. Noise tiles were `grey75`.
* Documentation: `?latents-conditions` lists every condition the code raises
  (checked by a test); `?lta`, `?multilpa`, `?parameter_inference`,
  `?bootstrap_lrt`, `?get_results`, `?three_step`, `?fit_staged`, the plot
  help and `course_engagement` are corrected where they disagreed with the
  code.

# latents 0.9.31

* Cluster-robust standard errors now carry the CR1 small-sample factor
  `G / (G - 1)` for `G` independent units, the default of
  `sandwich::vcovCL()`. It applies to every `vcov_type = "robust"` result
  (profile, covariate, transition, family, mixture-regression and growth fits)
  and to `three_step()` (both `vcov_type`s, every `method`) and `r3step()`.
  Robust standard errors rise by `sqrt(G / (G - 1))`: about 2.6% with 20
  groups and 0.5% with 100. A one-class
  robust regression now reproduces `sandwich::vcovCL()` exactly.
* Wald inference is refused (`latents_boundary_fit`) for a membership-covariate
  fit whose logit coefficients are heading to infinity because a covariate
  separates the classes. Previously the refusal held only when the observed
  information happened to be numerically singular; otherwise a coefficient
  such as -504 was reported with a standard error of 669 and p = 0.45.
* The separation check (and the `latents_extreme_coefficients` warning)
  measures each coefficient per root-mean-square unit of its covariate, so it
  no longer depends on the units a covariate is recorded in: a real effect of
  a covariate in tiny units is no longer flagged, and a separating covariate
  in large units no longer escapes.
* `lta()` defaults to `tol = 1e-10` (was `1e-8`), as do the bootstrap refits of
  its `parameter_inference()`. The group-class shares and the mean group
  posteriors, one EM step apart, now agree to about `2e-5` instead of `2e-4`,
  for roughly 1.5 times the iterations.

# latents 0.9.30

* Released function and method names and argument spellings remain available.
  Tables retain their existing columns (new columns are added; in an
  enumeration grid's tables `candidate` now comes first). Corrections to
  statistical results and input validation are described below.
* Audit corrections: bootstrap inference verifies supplied data, ordering and
  sampling weights before resampling; membership-covariate models apply and
  retain requested centering; constrained covariance updates preserve their
  structure at an active variance floor; growth likelihoods retain residual
  information at extreme random-effect scales; R3STEP inference is stable
  under predictor-unit changes; pooling structural constants reports zero
  missing information; descriptives refuse matrix-valued columns.
* Weighted residual tables use weighted descriptive moments and withhold
  unsupported design-based test statistics and p-values. Weighted profile
  trees use the fitted weights for shares and flows. Noise-model entropy
  includes noise mass; named noise diagnostics retain the supported sizes
  view and refuse unsupported entropy, posterior and AvePP views.
* Comparisons support homogeneous/general transition and membership-covariate
  fits both directly and through enumeration; classification follows the
  retained comparison rows. Data identity includes raw indicators, category
  labels, centering, grouping, occasions and weights. Covariate enumeration
  tolerates unavailable start-replication counts. `get_results(x, "all")`
  returns the available tables for regression comparisons and family/transition
  enumerations. Smallest-profile curves distinguish covariance structures.
* Aggregated-binomial prediction consistently reports expected successes.
  New-data trial counts are inferred from the response columns or supplied
  with the named-only `trials` argument (a number, vector or column name).
  When neither is available, a classed data error requests trial counts;
  `trials = 1` explicitly requests per-trial probabilities.
* Every candidate of an enumeration grid has a name, and the one name
  reaches it everywhere. `classes3_groups4` is 3 classes and 4 group classes
  (`profiles` when an indicator is continuous; no `_groups` part on a
  single-level grid; `_VVI` etc. when the grid crosses structures;
  `<model>_groups2`, such as `additive_equal_groups2`, on a group-class
  family grid). The name is the first
  column, `candidate`, of the grid's table and names each pick in
  `get_results(grid, "criteria")`; `grid$classes3_groups4` is that model, and
  an editor completes the names after `grid$`. Grids saved by earlier versions
  are named the same way.
* New `get_candidate()`, beside `get_results()`, takes one model out of a
  grid: `get_candidate(grid, "classes3_groups4")`, `get_candidate(grid, 3, 4)`
  or `get_candidate(grid, c(3, 4))`. `candidate_fit()`, its 0.9.21 name,
  still works.
* `compare_models()` compares candidates of a grid: by name,
  `compare_models(grid, "classes3_groups4", "classes3_groups5")`; by class
  counts, `compare_models(grid, c(3, 4), c(3, 5))`; or by selection,
  `compare_models(grid, n_profiles = 3)`, or every candidate,
  `compare_models(grid)` (candidates that failed are left out, with a
  message). It gives one row per candidate with its information criteria
  (the candidates table's numbers), entropies, AvePP summaries, smallest
  share, convergence, `delta_bic` and `bic_weight`, and takes a `cutoff` for
  the AvePP count. It also compares fits from `multilpa()`, `lpa()`, `lca()`
  or `get_candidate()` directly, refusing fits to different data. The
  printout marks the row with the lowest BIC as "lowest BIC" (it read
  "best", also for mixture-regression comparisons, whose columns are
  unchanged).
* A comparison carries the models it compared:
  `get_candidate(comparison, "classes3_groups4")` returns the model behind a
  row, and `get_results(comparison, "classification")` gives every class's
  AvePP and sizes for compared profile, class and transition models.
* Enumeration grids are much smaller in memory and on disk. Each candidate
  used to record, as its `call`, the whole fitting function and a full copy of
  the data, so a grid held one copy of the data per candidate. A candidate's
  call is now a short call that names the data, for example
  `multilpa(data = d, vars = ..., n_profiles = 3, n_group_classes = 8, ...)`.
  This applies to `enumerate_classes()` (profile, transition and group-class
  family grids) and `enumerate_regressions()`.
* Model-selection evidence for an enumeration grid, none of which selects a
  model, each reading the criterion in `using` (default `"bic_groups"`):
  `get_results(grid, "deltas")` gives, for every candidate, how much the
  criterion fell from one profile fewer and from one group class fewer;
  `get_results(grid, "knees")` gives a maximum perpendicular chord-distance
  knee, a simplified Kneedle-inspired rule without smoothing or sensitivity
  thresholds, along the profiles, along the group classes and
  over every candidate by its number of parameters, with the range it was
  found on and whether it holds without the last candidate.
  `plot(grid, what = "deltas")`, `"heatmap"`, `"smallest"` and `"trust"`
  draw the evidence; `what = "selection"` draws them all with the entropy
  and AvePP views.
* `plot(grid, what = "entropy")`, `"avepp"` and `"diagnostics"` (both) draw
  classification quality across the candidates of an enumeration grid, in
  the words of a single model's plot: relative entropy and the lowest AvePP
  of each candidate against the number of classes, one panel per level, the
  AvePP `cutoff` (0.8) dashed. `plot_enumeration()` takes them too.
* An enumeration grid compares classification quality across candidates. Its
  candidates table gains, beside the two entropies, the smallest, largest and
  mean average posterior probability (AvePP) of each candidate's profiles and
  group classes and the number of classes below a cutoff:
  `avepp_min`, `avepp_max`, `avepp_mean`, `avepp_below` and the same four
  prefixed `group_` for the group classes.
  The cutoff defaults to 0.8 and is set with
  `get_results(grid, "candidates", cutoff = )`. The columns are computed from
  the stored fits when the table is read, so grids fitted by earlier versions
  report them without refitting. `get_results(grid, "classification")` gives
  every class's AvePP, modal and estimated sizes, odds of correct
  classification and relative entropy for every candidate in one table.
* Named plot functions for the main views, each with its own listed arguments
  so an editor completes them: `plot_profiles()`, `plot_responses()`,
  `plot_prevalence()`, `plot_diagnostics()`, `plot_transitions()` and
  `plot_tree()`, alongside `plot_enumeration()`. Each draws exactly what the
  matching `plot(x, what = )` draws, and an argument the fit's family does not
  take is refused by name. `plot()` on a `diagnostics()` result now draws
  through `plot_diagnostics()`.
* Enumeration plots choose their criteria with `what`, as every other plot
  does: `plot(grid, what = "bic_groups")`, `plot(grid, what = "all")`,
  `plot(grid, what = "tree")`, and `plot_enumeration(grid, what = )`. The
  0.9.21 spelling, `criterion =` (alone or with `what = "enumeration"`),
  draws the same plot; `criterion` with any other `what` is an error of class
  `latents_bad_argument` (with `what = "tree"`, a warning, as the tree uses no
  criterion). This applies to transition and group-class family grids too.
* Every plot now returns a ggplot object. The last base-graphics plots are
  ported: `plot()` on a pooled imputation fit, on a `bootstrap_lrt()` result
  and on a mixture regression (`"fitted"`, `"coefficients"`, `"posteriors"`).
  Style them with ggplot2. The bootstrap plot still takes `palette`; its
  base-graphics `style` constants have nothing to set on a ggplot, so passing
  them draws the plot with a warning of class `latents_unused_argument`. The
  posterior strip of a mixture regression spreads its points deterministically
  and no longer draws from the random number stream.
  These formerly base plots draw when called, including inside loops, and
  return the ggplot invisibly. Code that used their previous returned fit
  object should keep its original fit separately; the returned object is now
  the reusable plot.
* `plot()` on a pooled fit chooses its parameters with `what` (`level`, which
  elsewhere means a confidence level, still works) and takes `subtitle`.
* No `plot()` method drops an argument it does not use without a word. Those
  that ignored one silently in 0.9.21 (comparison, growth mixture, mixture
  regression, pooled, additive, cross-level, general transition fits,
  bootstrap tests, transition and group-class family grids) now draw with a
  warning of class `latents_unused_argument`; the others refuse it with an
  error of class `latents_bad_argument`, as they did.
* `graphics` is no longer imported.
* The reference manual consolidates the former 74 help pages into 38 topics.
  Released help aliases are retained, and plot/prediction method usages list
  their actual arguments and defaults. `print()`, `summary()` and
  `as.data.frame()` are described once on
  `?get_results`, which also lists every family's tables (the eight
  `get_results.<class>` pages are folded in). All plotting, the fourteen
  `plot.<class>` pages and the `plot_*()` functions, is on one page,
  `?"latents-plots"`. `predict()` and `simulate()` share `?"latents-predict"`.
* The README's likelihood equation renders on GitHub: `\left\{` was stripped
  to `\left{` by GitHub's markdown, so the braces are now `\lbrace`/`\rbrace`.

# latents 0.9.23

* New `plot_enumeration()` draws the information criteria of an enumeration
  grid with its arguments listed on its own, so an editor completes them; it
  is the view `plot()` draws for an enumeration.
* `plot()` and `plot_enumeration()` take `criterion = "all"` for every
  information criterion in the grid. Pairs that coincide in a single-level
  grid are drawn once.

# latents 0.9.22

* Numeric categorical models now support `plot(what = "profiles")`, with
  `statistic = "mean"`, `"median"`, or `"mode"` on the original category
  scale. These summaries use fitted probabilities, support FIML and membership
  covariates, and do not compute standard errors. Mixed-model profile plots
  continue to show continuous indicators.

# latents 0.9.21

* The version submitted to CRAN; the same code as 0.9.20.

# latents 0.9.20

* Tests fit each model from a single start unless the test is about starts,
  which halves the time the tests take on CRAN. No package code changed.

# latents 0.9.19

* Every help page has an example: the `plot()` methods for model comparisons,
  family and transition enumerations, pooled imputations, cross-level and
  general transition fits, and `predict()` for growth mixtures, which shows
  the refusal and where the trajectories are.
* The transition-enumeration example in `?get_results.latents_transition_enumeration`
  fits the first 50 students, so it runs in a few seconds on Windows as well
  (it took 12 s there).
* The remaining tests that compare latents with another implementation move to
  `tests/equivalence/`: the depmixS4 and LMest fits, exact path sums, dense
  likelihoods, enumerated EM steps and independent Hessians, together with
  their reference helpers and fixtures. The installed tests no longer contain
  any such comparison.

# latents 0.9.18

* The growth mixture guide moves from the vignettes to the package website
  (same address), so the package ships three vignettes: `lpa`, `lca` and
  `lta`. Its model fits made it most of the vignette build time.

# latents 0.9.17

* The logo reads "latents" (it still showed the earlier name, "multilatent"),
  and the website icons are regenerated from it (they showed "multilpa").
* The reference manual is shorter: the `print()`, `summary()`,
  `as.data.frame()`, `coef()`, `vcov()`, `confint()`, `logLik()` and `nobs()`
  methods that had a page each are documented together on four pages
  (`?"latents-print"`, `?"latents-summary"`, `?"latents-as-data-frame"`,
  `?"latents-model-methods"`). The methods are unchanged.
* Fewer tests run on CRAN: tests that verify the estimators against their
  definitions, tests that only check a plot draws, and the regression tests
  from past audits are skipped there. All of them still run locally and in CI.

# latents 0.9.16

* Every help page now has an example that runs in the automatic checks. The
  growth mixture pages (`get_results()`, `plot()`, `simulate()`) and
  `compare_models()` open with a fit on the first 60 students of
  `growth_scores`, which runs in about a second; the full-data examples stay
  in `\donttest{}`. The `growth_scores`, `growth_schools` and `student_esm`
  pages open with a summary of the data. The `compare_models()` example no
  longer fits a random intercept that `course_engagement` cannot support.
* The package ships four vignettes: `lpa`, `lca`, `lta` and `growth-mixture`.
  The workflow guides (`workflow-lpa`, `workflow-lca`) and the guides to
  evaluation, covariates, mixture regression, latent class growth analysis and
  group-class models move to articles on the package website, with unchanged
  addresses. This keeps the installed package small and the checks quick.
* Tests that compare latents with other implementations (`numDeriv`,
  `mclust`, `MASS`, `glm()`, Latent GOLD and Mplus results, reference
  implementations, the numbers in the workflow guides) move to
  `tests/equivalence/`, which is not part of the package and runs in its own
  CI job. `numDeriv` is no longer a suggested package. Slow tests that refit
  many models are skipped on CRAN; every model family keeps tests that run
  there.
* The README's enumeration example names its covariance structure
  (`model = "VVI"`), as `candidate_fit()` needs once several structures are
  crossed, and its text matches the regenerated output.

# latents 0.9.15

* Vignettes and worked examples save and restore the user's display options
  (`digits` and `width`) after use, including the scripts installed in
  `inst/doc`.

* Corrects regression and growth inference: robust covariance checks the number
  and rank of independent units, multilevel growth BIC uses clusters, enumerated
  models retain the appropriate inference, and boundary parameters no longer
  erase unrelated delta-method standard errors.
* Stabilizes Poisson and negative-binomial regression densities, NB dispersion
  derivatives near the Poisson limit, and growth residual calculations for
  outcomes with large locations. Single-class binary regressions are accepted;
  nonfinite outcomes and evaluated predictors are refused explicitly.
* Trajectory tables include offsets and report response-scale standard errors.
  Alternate time columns work in tables and plots; random-effect spread holds
  other predictors fixed, and observed means respect sampling weights.
* Withholds regression bootstrap likelihood-ratio p-values when refits fail
  validation. Model comparisons also check trials, weights and cluster layout.
* Handles unassigned regression classes, weighted transition class counts and
  matrix-valued prediction inputs correctly. Noise cases retain their mass in
  profile percentages and appear as class zero in sequence tables and plots.

# latents 0.9.14

* Declares `MASS` and `withr`, used by the tests, in Suggests (the 0.9.13
  check warned about them as undeclared).

# latents 0.9.13

## One engine

* Every estimator now runs on one shared kernel (see
  `validation/ENGINE_DESIGN.md`): one EM driver (with SQUAREM for the profile
  models), one start runner, one quasi-Newton finish, one measurement block
  (Gaussian, categorical, ordinal and count indicators) shared by the
  profile, covariate, cross-level and both transition engines, one
  forward-backward pass for both transition engines, shared latent
  structures for the regression and growth models, and shared inference,
  alignment and simulation services. Results are unchanged: 251 recorded
  fits of every engine and option reproduce bit for bit
  (`equivalence/engine-golden/`), except where a defect below was fixed.

## New models

* Multilevel growth mixture models: `mixture_regression(cluster = )` nests
  the persons (`id`) in clusters (schools, clinics) with `n_group_classes`
  group classes of clusters that shift how probable each trajectory class is,
  and `group_membership` covariates of the group class. The clusters are the
  independent units of the likelihood, standard errors and BIC. New tables
  `"group_classes"` and `"clusters"`; `enumerate_regressions()` compares the
  number of group classes; `simulate()` draws the two levels. New dataset
  `growth_schools`; the growth vignette has a section on it.
* Negative-binomial regression mixtures: `mixture_regression(family =
  "negative_binomial")`, NB2 with a dispersion per class or shared
  (`variance = "equal"`). One class reproduces `MASS::glm.nb()`. A
  dispersion at the Poisson limit warns (`latents_boundary`) and the other
  standard errors are conditional on it.
* Ordinal regression mixtures: `mixture_regression(family = "ordinal")`, a
  proportional-odds (cumulative logit) regression within each class with
  class-specific ordered thresholds in place of the intercept. The outcome
  is an ordered factor, a factor or whole-number categories. Thresholds are
  rows `"threshold:<lower>|<upper>"` of the coefficient table; class means
  are expected category scores; `predict(type = "probabilities")` gives the
  category probabilities and `simulate()` draws categories in the outcome's
  own type. One class reproduces `MASS::polr()`. Mixtures that one outcome
  per class assignment cannot identify (two categories, or no predictors)
  are refused with `latents_not_identified`.

## Documentation

* New vignette "Latent class growth analysis: classes of trajectories"
  (`vignette("trajectory-classes")`): trajectory classes without random
  effects, every plot view, curved trajectories, pass/fail and ordered-band
  trajectories, distal outcomes, and when to prefer the growth mixture model.
  The growth mixture vignette now focuses on random effects and the
  multilevel model, with its plots.
* Coefficient forest plots of trajectory fits label every estimate with the
  same number of decimals per panel (two, or three for small values).
* The recovery table of a regression mixture prints in the report layout.

## Bug fixes

* Growth mixture models gain a start from the persons' own least-squares
  trajectories. Every previous start came from the model without random
  effects, which can split persons by their level; the growth fit then
  stopped at a lower maximum (by 9 to 11 log-likelihood units on test data
  where classes differ in slope).
* A fit with a noise component: `summary()` and `get_results(x, "all")` no
  longer fail; `"classification_errors"` and `"bch_weights"` (which were
  wrong for noise units) refuse with `latents_unsupported_noise`.
  `sensitivity()` no longer reports `NA` agreement for noise fits.
* `parameter_inference()` of a general `lta()` fit now uses the caller's
  `step` in its table, not only in the attached covariance.
* `get_tna()` and `get_group_tna()` on a general `lta()` fit refuse with
  `latents_unsupported_tna` instead of a bare "no applicable method".
* An additive-family information matrix that is negative definite is
  refused with `latents_singular_information` instead of an unclassed
  Cholesky error.
* A group covariate that varies within a group, and all starts failing in
  the covariate and homogeneous transition engines, now raise classed
  conditions (`latents_bad_data`, `latents_all_starts_failed`).
* Two-level regression summaries label a class logit's intercept within a
  group class "Intercept, group class 1" rather than
  "Intercept x group_class_1".

# latents 0.9.12

## Growth mixture models

* `mixture_regression()` gains `random`, `random_covariance` and
  `random_diagonal`: random effects within classes, so persons scatter around
  their class trajectory (the Gaussian growth mixture model). The likelihood
  is exact; estimation is EM with the random effects as missing data,
  finished by quasi-Newton with analytic scores. Validated against an
  independent dense likelihood, numerical derivatives and `lcmm::hlme()`
  (six configurations; the same likelihood at `hlme`'s estimates and at
  least its maximum).
* Variable names may replace one-sided formulas: `random = "wave"`,
  `random = "intercept"`, `common = "x"`, `membership = "age"`.
* A degenerate random-effect covariance (a variance at zero, or perfectly
  correlated random effects) warns (`latents_random_boundary`) and reports
  no interval for it.
* New `compare_models()`: fits to the same data side by side, with BIC
  differences and Schwarz weights; `plot()` draws the comparison.
* Trajectory models (`id` with `class_level = "group"`, with or without
  `random`) share new ggplot2 views: class trajectories with confidence
  bands, the spread of persons around them and observed class means;
  individual trajectories; a coefficient forest; predicted random effects;
  and a classification structure plot.
* Result tables of these fits print in a report layout (readable labels,
  intervals as `[low, high]`, APA p-values, fit statistics as a card) while
  remaining plain numeric data frames.
* Growth fits work with `enumerate_regressions()`, including its bootstrap
  likelihood-ratio test, have a `simulate()` method, and support distal
  outcomes in `three_step()` (person-level, bias-corrected); `three_step()`
  results print in the report layout. `predict()` and `r3step()` refuse
  growth fits with classed conditions.
* New dataset `growth_scores` and vignette "Growth mixture models".

# latents 0.9.11

* Transition fits no longer report false non-convergence. The quasi-Newton
  finish counted only L-BFGS-B's code 0 as converged, but its line search can
  stop at the maximum (code 52, "ABNORMAL_TERMINATION_IN_LNSRCH") when rounding
  hides any further ascent, which is platform-dependent. It now restarts once
  from where it stopped and accepts the fit when the log likelihood cannot be
  raised by more than `tol`. In 120 negative-binomial test datasets 5 were
  flagged unconverged at points with scores of about 1e-6 that a restart
  improved by less than 1e-12; none are now. This, not the iteration cap, is
  what failed the 0.9.8-0.9.10 CI checks on macOS and R 4.1 (the 0.9.10
  test change is reverted).

# latents 0.9.10

* Test-only fix: a negative-binomial transition test required convergence
  within `max_iter = 5`, which held on some platforms only (macOS CI and R 4.1
  needed more iterations, and the package correctly warned). The test now
  allows 500 iterations; its checks of dispersion floors, monotone likelihood
  and stationary scores are unchanged.

# latents 0.9.9

## Transition models: full-covariance simulation and bootstrap inference

* Transition fits with full covariance structures (EEE, VVV and the other
  non-diagonal structures, homogeneous or extended) can now be simulated, so
  `bootstrap_lrt()` compares them. Each profile's rows are drawn through the
  Cholesky factor of its covariance; diagonal fits keep their random stream.
* `parameter_inference(method = "bootstrap", data = )` works on transition
  fits. Persons are resampled with replacement and refitted; each replicate's
  profiles are matched to the original's on the measurement and its group
  classes on what the matched profiles make of them, with initial logits
  rebased when the reference profile moves. Standard errors are the
  replicates' standard deviation and intervals their percentiles. This gives
  inference for extended models whose Wald inference is unavailable (full
  covariance structures) and supports sampling weights.
* Fixed a 0.9.8 regression: transition fits refused `boundary = "fix"`, which
  `get_results(fit, "responses")` and `summary()` request, so categorical
  transition response tables lost their standard errors. Transition inference
  holds nothing at a bound, so `"fix"` now agrees with `"error"`; a fit with an
  active bound is still refused with `latents_boundary_fit`.

# latents 0.9.8

## Model review fixes

* Ordinal/count indicators now work through `lca()`, `multilca()` and
  `enumerate_lca()`, as documented. Bootstrap and imputation label alignment
  move their measurement parameters together, restore the ordinal reference
  profile, and match profiles using the extra indicators.
* Negative-binomial bootstrap comparisons preserve the fitted count model and
  refuse incompatible indicator specifications. Supplied fitting data are
  checked against all measured indicators and model-specific sufficient
  statistics before inference or comparison.
* Poisson log densities and negative-binomial dispersion derivatives use
  numerically stable calculations. Invalid extra-indicator column shapes and
  nonfinite ordinal values are refused. Single-level covariate inference
  accepts the original data without demanding an internal identifier.
* Second-order occasion transitions omit coefficient coordinates that never
  apply. Stayer tables report identity transitions; second-order inference
  includes the second-order coefficient block. Count boundaries retain their
  constraints during transition optimization and refuse Wald inference.
* Extended transition `coef()` now carries complete measurement covariance
  coordinates for all fourteen structures, with likelihood-preserving decoding.
  The existing restrictions on extended-model Wald inference remain explicit.
* Transition fits with only categorical, ordinal or count indicators return all
  result tables successfully, including an empty continuous-profile table.
  Transition data validation checks original category labels and occasion values.
* One-profile NB fits with varying dispersion retain the correct profile label
  and report their dispersion standard error.
* Group-class bootstrap comparisons verify group sizes and scatter and refuse
  equivalent one-class family specifications. Weighted restricted cross-level
  fitting avoids collisions with indicator names used for internal weights.
* Sampling weights normalize safely at extreme common scales. Regression
  mixture class counts, shares and composition respect sampling weights;
  invalid integer controls and seeds are refused before conversion.
* Confidence levels and numerical differentiation controls receive explicit
  validation on the reviewed model surfaces.

## Follow-up review fixes

* Negative-binomial fits near the Poisson limit no longer stop short of the
  maximum. The dispersion objective is convex there in log dispersion, and
  the old Newton step could stall while still reporting convergence (fitted
  dispersion 4.2e-7 against a profile maximum of 1.4e-5, with singular Wald
  inference). The M-step now uses modified Newton with reflected eigenvalues.
* `bootstrap_lrt()` refuses a full-covariance transition null up front with
  `latents_unsupported_inference`, rather than failing every replicate and
  returning `p_value = NA` with a generic warning.
* `parameter_inference()` on transition fits refuses `method = "bootstrap"`
  and `boundary = "fix"` with `latents_unsupported_inference` instead of
  silently returning the Wald table. On general `lta()` fits, `adjust` now
  produces a `p_adjusted` column.
* `bootstrap_lrt()` accepts the caller's own data frame for single-level fits.
* Under sampling weights, the `"classification"` table's `estimated_n` and
  `estimated_proportion` (and so the odds of correct classification and
  `diagnostics()`' smallest class) use weighted posterior totals. The same
  applies to restricted cross-level `composition` and `group_classes` counts.

# latents 0.9.7

## Ordinal and count indicators everywhere; negative-binomial counts

* Membership covariates (`profile_covariates`, `group_covariates`) and
  `lta()` (every extension, invariant or occasion-specific measurement) now
  take `ordinal` and `count` indicators, with Wald inference, weights and
  simulation (`bootstrap_lrt()` for transition fits).
* `count_model = "negative_binomial"` models counts more variable than a
  Poisson within a profile (NB2: variance `mu + alpha mu^2`), with
  `count_dispersion = "varying"` (per profile) or `"equal"` (shared).
  `get_results(fit, "count_means")` adds the dispersion and its standard
  error. A dispersion at zero is the Poisson limit: the fit warns
  (`latents_boundary`) and Wald inference is refused there.
* Checked against Latent GOLD 6.1 (`poisson overdispersed`), class-varying
  and shared dispersion: log likelihoods agree to 4e-5 and parameters to its
  printed precision; seven reference cases in all, stored as a test fixture.
* The Newton line searches of the membership logits, ordinal and count
  M-steps accept a step that changes the objective by rounding only, instead
  of halving it up to forty times at the maximum (30x faster on a negative-
  binomial fit).

# latents 0.9.6

## Bug fix

* `multilpa()` with `profile_covariates` or `group_covariates` silently fitted
  any covariance structure other than EEI, VVI, EEE and VVV as one of those
  four (for example `model = "VEI"` gave VVI) and did not record the
  structure. It now refuses them with `latents_unsupported_structure`; fit the
  covariate-free model with the structure and use `three_step()` or
  `r3step()`.
* `?"latents-conditions"` documents the condition classes added in 0.9.3–0.9.6.

# latents 0.9.5

## Ordinal and count indicators

* `multilpa()`, `lpa()`, `lca()` and `multilca()` take `ordinal` and `count`,
  naming indicators in `vars`:
  * `ordinal`: an adjacent-category logit with category intercepts shared by
    every profile and one location per profile (Latent GOLD's default
    ordinal model): `(K - 1) + (C - 1)` parameters per indicator instead of
    `C (K - 1)` as `categorical`. Ordered factors or whole numbers.
  * `count`: Poisson, one mean per profile.
* New tables: `get_results(fit, "ordinal")` (category probabilities and
  locations, with standard errors) and `get_results(fit, "count_means")`.
* Mixed freely with continuous and categorical indicators, in single- and
  two-level fits, with `missing = "fiml"` and `weights`. Wald (analytic
  scores) and bootstrap inference, `bootstrap_lrt()`, `predict()`,
  `sensitivity()` and `enumerate_classes()` handle them.
* Checked against Latent GOLD 6.1 on five datasets (mixed, ordinal only,
  count only, three classes, two-level): log likelihoods agree to 4e-5 and
  every parameter to 5e-5, the precision Latent GOLD prints
  (`equivalence/latentgold-ordinal/`); the results are a test fixture.
* Not yet available with membership covariates, `start`, `fixed`, `prior`,
  `noise`, `lta()` or the group-class families (`latents_unsupported_indicator`).

# latents 0.9.4

* The membership M-step of covariate fits (`profile_covariates`,
  `group_covariates`) now solves its weighted multinomial logits by
  Newton-Raphson with the exact Hessian instead of BFGS. BFGS stalled above
  its score tolerance on about half of the steps of an ordinary fit, which
  left fits flagged unconverged (and refused standard errors) on some
  platforms but not others. Estimates are unchanged; fits converge in fewer
  EM iterations.
* Continuous integration passes again: the R 4.1 floor job no longer tries to
  install `mice` (its current dependency chain needs R >= 4.4), one lint is
  fixed, and the checkout action is updated.

# latents 0.9.3

## Sampling weights

* `multilpa()`, `lpa()`, `lca()`, `multilca()`, `lta()` and
  `mixture_regression()` take `weights`, the name of a column holding one
  sampling weight per independent unit (each `id` group, or each row of a
  single-level fit). The fit maximizes the pseudo log likelihood
  `sum_j w_j log L_j`, with the weights scaled to sum to the number of units
  (Mplus's convention), so AIC and BIC stay on the sample's scale.
* Covered: every covariance structure, categorical indicators,
  `missing = "fiml"`, membership covariates, every `family` (additive,
  dispersion, additive-dispersion, restricted and full cross-level), every
  extension of `lta()`, and all three nestings of `mixture_regression()`.
* Integer weights reproduce the fit to the data with each unit repeated that
  many times, to 1e-9 or better in the log likelihood, for every engine
  (tested). A one-profile weighted fit reproduces Mplus 9's weighted
  regression (`TYPE = COMPLEX`, `WEIGHT`) and closed-form weighted least
  squares (`equivalence/weights/`).
* Standard errors are the sandwich: `parameter_inference()`, `vcov()`,
  `confint()` and `get_results()` default to `vcov_type = "robust"` for a
  weighted fit and refuse `"observed"` and `"opg"`
  (`latents_unsupported_weights`). `method = "bootstrap"` resamples units
  with their weights.
* A weighted fit prints its weight column and Kish's effective sample size;
  `get_results(fit, "model")` gains a `weights` column. Effective counts and
  proportions are weighted; classification and posterior tables report each
  unit's own posterior.
* Refused with `latents_unsupported_weights`: `bootstrap_lrt()` and
  `enumerate_regressions(bootstrap = )` (the parametric bootstrap ignores the
  design), `three_step()`, `r3step()`, `prior` and `noise`. Weights that vary
  within a unit, are negative or missing raise `latents_bad_weights`.

# latents 0.9.2

## Choosing and testing transition models

* `enumerate_classes(..., time = )` enumerates latent transition models over
  `n_profiles` and `n_group_classes` (and covariance structures named in
  `model`), with every other `lta()` argument held fixed; read with
  `get_results()`, `plot()` and `candidate_fit(grid, n_profiles = ,
  n_group_classes = )`.
* `bootstrap_lrt()` compares nested transition fits (`data` required): more
  profiles or group classes, occasion-varying against homogeneous
  transitions, added covariates, occasion-specific measurement, second
  order, a stayer class. The null model is simulated with the data's groups,
  occasions, covariates and missing values kept.
* `lta(mover_stayer = TRUE)` adds a class of stayers who never change
  profile (Goodman's mover-stayer model), with its own initial distribution.
* The extended transition model now takes `missing = "fiml"` and covariance
  structures (`model =`, `covariance_model = "full"`); standard errors with
  diagonal covariances (including FIML).
* A transition or initial probability estimated at zero (a move never
  observed, common with occasion-varying transitions on sparse late
  occasions) is reported as a boundary fit (`latents_boundary`, standard
  errors withheld) instead of as non-convergence.

# latents 0.9.1

## Latent transition analysis: beyond homogeneous transitions

* `lta()` gains arguments that relax each default assumption of the model:
  * `transition_covariates`: covariates (fixed or changing over occasions)
    shift the transition probabilities through a multinomial logit per
    origin profile, the log odds of moving to each other profile rather
    than staying;
  * `initial_covariates`: covariates shift the starting profile;
  * `transitions = "occasion"`: a separate transition matrix for each move;
  * `measurement = "occasion"`: profile means and variances (or response
    probabilities) per occasion;
  * `order = 2`: second-order transitions, the next profile depending on the
    two previous ones;
  * `model =`: the fourteen covariance structures (`"VEI"`, `"EEV"`, ...)
    for the homogeneous model, as in `multilpa()`.
* The extensions return a `multilpa_lta` fit with tidy tables
  (`get_results(fit, "transition_coefficients")`, `"transitions"` per
  occasion, `"second_order_transitions"`, `"initial_coefficients"`,
  `"profiles"`, ...), Wald standard errors from analytic scores (observed,
  robust or OPG), `coef()`, `vcov()`, `logLik()`, `nobs()`, `BIC()` and
  `plot()`. EM is finished by a quasi-Newton search on the exact likelihood.
* External agreement, 30 of 30 quantities (`equivalence/lta-extensions/`):
  Mplus User's Guide examples 8.13 and 8.14 and an Mplus model with
  occasion-specific thresholds (log-likelihoods within Mplus's printed
  precision, estimates within 2e-3), depmixS4 covariate transitions
  (log-likelihood within 1e-7) and LMest time-heterogeneous transitions
  (within 4e-8). Second-order transitions, which no external program
  fits, are checked against an exact sum over every profile path.

# latents 0.9.0

## Experimental: the Houle et al. (2026) multilevel families

* `multilpa(family = )` fits the multilevel latent profile families of
  Houle, Morin & Harvey (2026) beside the existing profile model
  (`family = "profiles"`, their dispersion-heterogeneity model):
  * group-class families with no individual profiles, each group carrying a
    Gaussian intercept per indicator: `"additive"` (classes differ in means;
    within-group variance shared), `"dispersion"` (classes differ in
    within-group variances) and `"additive_dispersion"` (both);
    `between_variance = "varying"` or `"equal"`. The likelihood is exact (no
    quadrature); zero between-group variances are reached by a boundary
    maximization with a Karush-Kuhn-Tucker check; interior fits are finished
    by Newton steps.
  * cross-level families with individual profiles and group classes from the
    same indicators, following the published manifest-aggregation
    specification: `"restricted_cross_level"` and `"full_cross_level"`.
    Their likelihood is a working likelihood (group means reuse the
    ratings), so it compares only cross-level fits of the same data; no
    standard errors.
* Tidy tables for every family through `get_results()` (`"parameters"`,
  `"group_classes"`, `"groups"`, `"intercepts"`, `"recovery"`, ... for the
  group-class families; `"profiles"`, `"composition"`, ... for the
  cross-level ones), `summary()`, `plot()`, `coef()`, `vcov()`, `confint()`,
  `logLik()`, `nobs()` and `parameter_inference()` (observed, cluster-robust
  or OPG covariance from analytic group scores).
* `enumerate_classes(family = c("additive", "dispersion",
  "additive_dispersion"), n_group_classes = )` crosses families,
  between-variance restrictions and class counts; `candidate_fit(grid,
  n_group_classes = , model = )` picks one. `bootstrap_lrt()` tests nested
  group-class fits by parametric bootstrap.
* `latents_weak_class` warning with `effective_groups` in the
  `"group_classes"` table: a class supported by fewer than 50 effective
  groups (small, poorly separated, or both) has intervals that can be
  miscalibrated. Threshold set by a predeclared rule and evaluated on
  held-out simulations (flags 100% of a rare, weakly separated condition and
  0.8% of well-separated ones). Also flags designs with fewer than about 50
  groups per class.
* Evidence (`validation/ADDITIVE_SIMULATION.md`,
  `equivalence/latentgold-families/`): predeclared simulations with 1000
  datasets per condition (additive: 31 of 34 reference-condition parameters
  meet every gate; dispersion and additive-dispersion: 45 of 46; every miss
  is the ML small-sample downward bias of a between-group variance, coverage
  0.925-0.947). Latent GOLD 6.1 agrees on parameter counts in all eight
  external cases, and its likelihoods and posteriors converge to latents'
  as its quadrature is refined (largest remaining gaps 1.6e-4 and 2.0e-4).
  Robust standard errors are needed for heavy-tailed ratings.
* New vignette: `vignette("additive")`.
* Not yet available for these families: missing data, covariates, and an
  external Mplus comparison (Mplus is not available here).

# latents 0.8.8

## plot() returns ggplot objects (breaking)

* `plot()` on a fit (`multilpa`, covariate and transition fits), on an
  enumeration and on a `diagnostics()` result now returns ggplot objects
  instead of drawing with base graphics. Print one to draw it, save it with
  `ggplot2::ggsave()`, or restyle it with `+ ggplot2::theme()`. ggplot2 stays
  in Suggests: a plot without it raises `latents_missing_package`, and
  `diagnostics()` and `report()` print a message and carry on.
* `what = "all"`, `plot()` on a diagnostics result and an enumeration plotted
  with `combine = FALSE` return a `latents_plots` list, named by view, that
  draws every plot when printed.
* The styling arguments `palette`, `symbols`, `linetypes`, `style` and the
  style constants passed through `...` are gone; an unknown argument now
  raises `latents_bad_argument` instead of being ignored.
* Every view shares one profile order (largest first), one profile share
  (the posterior share every table reports) and Okabe-Ito colours paired with
  shapes. Bars start at zero, the heatmaps have colour keys, and the case
  diagnostics are strips of every case with each profile's mean marked.
* New views: `"parallel"` (every case as a line, one panel per profile),
  `"pairs"` (scatter-plot matrix with each profile's 95% covariance ellipse)
  and, for enumerations, `what = "tree"` (an icicle of how profiles split as
  more are added). `plot_views()` lists them.
* Model comparison switches from direct labels to a legend beyond five
  series; failed candidates are marked with a cross on the panel floor.
* `plot(diagnostics(fit))` on a one-profile fit returns the sizes and average
  posterior views instead of refusing.

# latents 0.8.7

## Messages

* The single-level notice from `multilpa()`, `multilca()` and
  `enumerate_classes()` with `id = NULL` is now one line ("Fitted a
  single-level latent profile model."). The error for a missing `id` is
  shorter too.

## ggplot2 views (in development)

* Internal ggplot2 builders for every clustering view: profiles, bars (from
  zero), heatmap (mixture-standardized, with a colour key), raincloud, sizes,
  entropy and posteriors (per-case strips), average posterior probability and
  model comparison. New views: parallel coordinates, a scatter-plot matrix
  with 95% covariance ellipses and uncertainty-sized points, and an icicle of
  how profiles split across numbers of profiles (Zappia and Oshlack, 2018),
  with bands carrying posterior mass and colours following each profile's
  lineage.
* Every view uses one profile share (the posterior share), one display order
  (largest first) and Okabe-Ito colours paired with shapes.
* ggplot2 is in Suggests. `plot()` still draws the base-graphics views; the
  switch comes in a later release.

# latents 0.8.6

## Revision review fixes

* Mixture-regression prediction retains the training membership factor levels,
  contrasts and transformed terms. Posterior predictions rebuild group IDs
  for new rows; fitted plots retain grouping and membership columns and accept
  transformed regression predictors.
* A group-membership formula with one slope and no intercept is optimized as
  a slope, rather than mistaken for an intercept-only model.
* `r3step(by_group_class = TRUE)` now handles the one-group-class limit and
  agrees with the pooled regression for both observed and robust covariance.
* `enumerate_lpa()` preserves an omitted `model`, allowing the documented
  covariance switches and all-categorical inputs through the wrapper.
* `bootstrap_lrt()` refuses prior-fitted models whose likelihoods are evaluated
  at posterior modes. Its fixed-mask missingness assumptions are clarified,
  and multiple-imputation guidance distinguishes algebraic equivalence from
  evidence about coverage and imputation-model validity.

# latents 0.8.5

## Single-level verbs, covariance-structure grids and workflow vignettes

* New `lpa()` and `lca()` fit single-level latent profile and latent class
  models. They are `multilpa(id = NULL)` and `multilca(id = NULL)` under
  their ordinary names, return the same object, and refuse `id` and
  `n_group_classes`.
* `multilpa()` and `lpa()` take `model`, the three-letter covariance code
  (`"EEE"`, `"VVI"`, ...), as an alternative to `variance_model` /
  `covariance_model` or `volume` / `shape` / `orientation`; giving both is an
  error.
* New `enumerate_lpa()` and `enumerate_lca()` compare single-level latent
  profile and latent class models, as `lpa()` and `lca()` fit one:
  `enumerate_lpa(data, vars, n_profiles = 1:6)` crosses the counts with the
  covariance structures, and `enumerate_lca(data, vars, n_classes = 1:6)`
  treats every indicator as categorical. `enumerate_classes()` remains the
  verb for nested data and needs `id`.
* The two-level verbs (`multilpa()`, `multilca()`, `enumerate_classes()`)
  with `id = NULL` say which analysis was estimated -- a single-level latent
  profile or latent class analysis -- and name the single-level verbs; called
  without `id` they are refused with the same pointer.
* `enumerate_classes(model = )` now defaults to `"basic"`: the four
  structures that combine equal or varying variances with covariances absent
  or present (`EEI`, `VVI`, `EEE`, `VVV`). `"all"` fits the 14 structures;
  codes can still be named. With one continuous indicator the structures that
  coincide are fitted once; with only categorical indicators no structure is
  crossed; a structure set through the other arguments replaces the default,
  and `model = NULL` restores the previous grid. Because a grid now holds
  several structures per class count, `candidate_fit()` needs `model` to
  pick one of them.
* New data set `srl`: five self-regulated learning scales for 300
  respondents simulated by a large language model.
* New vignettes `vignette("workflow-lpa")` and `vignette("workflow-lca")`:
  a latent profile and a latent class analysis from the data to a reported
  model.
* `get_results(fit, "profiles")` and `get_results(fit, "responses")` carry
  `mean_standard_error`, `variance_standard_error` and
  `probability_standard_error` without the data being passed; they are `NA`,
  with a message, when the information matrix cannot be formed.
* `diagnostics()` draws its classification plots by default
  (`plots = TRUE`): profile sizes, posterior probabilities, entropy
  contributions and the average posterior probability matrix.
* `plot(what = "raincloud")`: for each profile and indicator, the density,
  a quartile box and the observations.
* Wald intervals for probabilities are formed on the logit scale and for
  variances on the log scale, so they stay inside the parameter's range;
  `confint()` returns exactly the intervals `parameter_inference()` reports.
* A probability at or below 1e-6 counts as on its bound for inference: EM
  approaches zero slowly and stopped at 1e-8, leaving a singular information
  matrix that `boundary = "fix"` could not hold.
* `r3step()` and `three_step()` name their adjusted p-value `p_adjusted` and
  their classes `profile_1`, `group_class_1`, as `parameter_inference()` does
  (were `p_value_adjusted` and `class_1`).
* `id = NULL` now raises its `latents_single_level` notice as a message, not
  a warning: a single-level model is a legitimate choice.
* `diagnostics()` prints only the individual level for a fit with one group
  class, and a single-level fit's assignments carry no `group_class` column.
* `plot(what = "profiles")` draws 95% intervals when the fit has standard
  errors (`intervals = TRUE`, profiles dodged apart);
  `plot(what = "responses")` draws intervals too.
* `plot(what = "bars")` and the profile plot no longer fail on a bound-active
  or unconverged fit: they draw without whiskers and the subtitle says why.
* `plot(what = "heatmap")` on an all-categorical fit draws the response
  probabilities: one row per class, one column per category, values printed
  where they fit.
* Covariate fits draw every measurement and classification view
  (`"bars"`, `"heatmap"`, `"raincloud"`, `"sizes"`, `"avepp"` are new for
  them).
* The enumeration plot uses display names (AIC, BIC, ICL), draws one BIC
  panel when the group-level and individual-level values coincide (a
  single-level grid), labels series only by what distinguishes them, and
  spaces labels by the text height. Gridlines and zero lines stay inside the
  panel. Okabe-Ito yellow now comes seventh, after the six colours that read
  well on the light panel.
* `summary()` of an enumeration prints a compact candidates table and names
  `get_results()` for the rest.
* An all-categorical fit prints as a latent class analysis rather than with
  a residual covariance, and a single-level fit's header reads "BIC".
* A covariate fit with one group class names its membership intercept
  `(Intercept)`, not `group_class_1`.
* Help-page examples use the bundled data sets and verbs, with no `$` or
  bracket indexing; a test keeps it that way.

## Bootstrap likelihood-ratio tests with missing data and covariates

* `bootstrap_lrt()` now accepts fits made with `missing = "fiml"`: every
  simulated replicate is given the observed missing cells before it is
  refitted, so the reference distribution loses the same information as the
  observed statistic. The pattern is treated as fixed (independent of the
  profiles).
* `bootstrap_lrt()` now accepts membership-covariate fits: the covariates are
  held at their observed values and memberships are drawn from the fitted
  logits. A one-group-class null may be compared with an alternative that adds
  group covariates or slopes by group class; otherwise both models must use the
  same membership regressions (`latents_bad_nesting`). A covariate replicate
  whose likelihood has converged counts as valid even if a membership logit is
  still drifting (an over-fitted class), since the statistic reads only the
  maximized likelihood; the replicate table's new `logits_settled` column
  records it.
* Models with different missing-data handling are now refused as
  incomparable (`latents_incomparable_models`).

## Standard errors for latent transition models

* `parameter_inference()`, `vcov()` and `confint()` now work on `lta()` fits
  instead of refusing them: observed-information, robust (clustered on
  sequences) and OPG standard errors for the measurement model, the
  group-class shares, the initial profile probabilities and the transition
  probabilities. The scores follow the Fisher identity (expected counts from
  the forward-backward pass minus the counts the fitted probabilities imply)
  and match numerical differentiation of the likelihood; the standard errors
  match a Richardson-extrapolated Hessian to 1e-6 and depmixS4's numerical
  ones to within its own finite-difference error. Over 300 simulated panels,
  95% intervals covered 0.947-0.960. The `latents_no_inference` condition is
  retired.
* `coef()` on an `lta()` fit now uses the names and order of `vcov()`: under
  `variance_model = "equal"` it reports the shared variances once, under
  `covariance_model = "full"` the covariance matrices, and the group-class
  shares come before the initial and transition probabilities.

## Multiple imputation for missing covariates

* New `pool_imputations()`: fits `multilpa()` to every completed data set
  (a list of data frames, or a `mids` object from mice), aligns each fit's
  profile and group-class labels to the first, and pools with Rubin's rules.
  The table adds `df`, `within`, `between`, `riv` and `fmi` to the columns of
  `parameter_inference()`; `get_results(x, "imputations")` and
  `get_results(x, "fits")` show each imputation and the relabelling it needed,
  and `plot()` shows the imputations' estimates beside the pooled interval.
  This is the route for missing covariates, which `missing = "fiml"` cannot
  integrate out. Aligning a covariate fit reparameterizes its membership
  logits, whose reference category moves with the labels. Pooled values agree
  with `mice::pool.scalar()` to machine precision. A failing imputation raises
  `latents_pooling_failed` rather than being dropped.

## First-stage uncertainty in staged fits

* `parameter_inference(method = "bootstrap")` now accepts a `fit_staged()`
  result instead of refusing it. Every resample of the groups refits both
  stages, so the standard errors and percentile intervals carry the
  measurement's sampling variability into the group-class estimates; the
  Wald default still conditions on the measurement. The table reports every
  parameter, the measurement included. On a simulated two-class design the
  within-class profile probabilities' errors were 14% and 32% wider than the
  conditional ones, and a bootstrap that held the first stage reproduced the
  conditional Wald errors to within 3%.

## Profile-covariate slopes by group class

* `multilpa(profile_covariates = , profile_slopes = "group_class")` lets each
  group class have its own profile-covariate slopes, so a covariate can predict
  profile membership differently in different kinds of group (a cross-level
  interaction). The coefficients are reported with terms such as
  `z:group_class_1`, and inference, robust and OPG covariances cover them.
  The default, `"shared"`, is the previous model. With one group class the two
  are the same model. Naming it without `profile_covariates` raises
  `latents_bad_argument`.

## Missing indicators in membership-covariate models

* `multilpa(profile_covariates = , group_covariates = , missing = "fiml")`
  now fits: missing indicators are integrated out of the measurement density
  exactly as in the covariate-free model, for diagonal and full residual
  covariance and for categorical and mixed indicators. Previously this
  combination was refused with `latents_bad_argument`. Standard errors
  (observed, robust and OPG) cover it through the Fisher identity: a missing
  value enters the scores through its conditional mean and covariance given
  the row's observed values. Covariates themselves must still be complete; a
  missing covariate now raises the classed `latents_bad_data`.

## Faster estimation

* `multilpa(acceleration = "squarem")`, the new default, accelerates EM with
  SQUAREM (Varadhan & Roland 2008). A cycle extrapolates along the last two
  EM steps and finishes with an ordinary EM step, falling back to plain EM
  whenever the extrapolation would lower the likelihood, so the path stays
  monotone; convergence is judged over the whole cycle, so it is never looser
  than plain EM's. `acceleration = "none"` reproduces plain EM's path (use it
  to follow mclust or Mplus step by step). Prior fits, membership-covariate
  models always use plain EM. Results can
  differ slightly from earlier versions: the accelerated fit usually ends
  closer to the maximum, and on a flat, weakly identified likelihood it may
  stop at a different point of the ridge.
* The first start of `multilpa()` is now Ward's hierarchical clustering, with
  each cluster's own means and spread; the rest remain random k-means. Over
  14 structures on 3 datasets it never lowered the best-of-10 likelihood and
  raised it in 3 of 36 cases, for example on iris VVV from -186.57 to
  mclust's -180.19. Fits may therefore land on a better maximum, or on the
  same maximum with profiles numbered differently.
* EVE and VVE take a fixed number of warm-started orientation steps per EM
  iteration (a generalized EM step) instead of solving the orientation to
  convergence inside every M-step, which took a median of 220 inner steps per
  iteration. The same maxima are reached 3.5 to 8 times faster, and SQUAREM
  now applies to them.
* `enumerate_classes(id = NULL)` enumerates single-level models, the search
  `mclust::mclustBIC()` performs, with one single-level notice for the grid.
* New `predict()` method for `multilpa()` and `multilca()` fits: the modal
  profile (`"class"`), every profile's probability (`"posterior"`), or each
  row's log density under the fitted mixture (`"density"`) for new rows,
  prepared as the fit prepared its own (categorical levels, centering,
  missing data, noise). Groups in the new data are classified as new groups.
  On the training rows it reproduces the fit's posteriors, and it matches
  `mclust::estep()` at the same parameters.
* `descriptives(by = )` lists strata in a fixed order (a factor's levels,
  otherwise sorted) rather than in order of first appearance.
* The E-step no longer computes a missingness pattern per row for complete
  data and takes row maxima without a per-row `apply()`: about 14 times faster
  per iteration, with identical results.
* On 5,000 rows, 4 indicators and 3 profiles, 13 of the 14 covariance
  structures now fit faster than `mclust::Mclust()` from one start (VVV is
  the exception), and 13 of 14 reach a higher likelihood.

## Mixture regression

* New `mixture_regression()` fits finite mixtures of regressions (clusterwise, latent
  class or regression-mixture models) for `"gaussian"`, `"binomial"` (0/1,
  factor or `cbind(successes, failures)`) and `"poisson"` (with `offset()`)
  outcomes. Classes can belong to each row, to a whole group
  (`class_level = "group"`, flexmix's `y ~ x | id`), or to each row with a
  second-level group class that shifts the class shares within groups
  (`n_group_classes`, Vermunt 2003). Membership covariates at both levels
  (`membership`, `group_membership`), coefficients shared across classes
  (`common`), and equal residual variances (`variance = "equal"`).
* Analytic scores (Fisher's identity) give observed-information, sandwich
  (clustered on the independent unit, or on `id` for a single-level fit) and
  outer-product standard errors; validated against a numerical Hessian of the
  likelihood in every nesting. A Monte Carlo study
  (`validation/mixture-regression-recovery.R`, 200 replications) found 91.5--97% coverage
  of nominal 95% intervals for the two-level model.
* `get_results()` tables: `coefficients` (with Benjamini-Hochberg-adjusted
  slopes and odds/rate ratios), `classes`, `membership`, `group_classes`,
  `fit` (AIC, BIC, SABIC, ICL, entropy), `assignments`, `groups`, `fitted`,
  `starts`, `classification` and `recovery` (against a known classification).
  `predict()`, `simulate()`, `plot()`, `summary()`, `coef()`, `vcov()`,
  `confint()`, `logLik()` and `nobs()` methods.
* New `enumerate_regressions()` compares class counts in one table, with an
  optional parametric bootstrap likelihood-ratio test.
* A binary outcome with one trial per class assignment is refused
  (`latents_not_identified`), since that mixture is not identified.
* Single-level and group-level fits reproduce flexmix's likelihood at its
  estimates and its binomial and Poisson estimates
  (`equivalence/test-mixture-regression-flexmix.R`).
* New dataset `study_hours` and `vignette("mixture-regression")`.

## Fixes

* `multilca(id = NULL)` now fits one group class by default, as `multilpa()`
  does, instead of refusing the default `n_group_classes = 2`.
* `equivalence/run.R` evaluated tests in the retired `multilpa` package's
  namespace when it was still installed; it now uses `latents`.

## mclust parity

* `parameter_inference()`, `vcov()` and `confint()` give Wald standard errors
  for all fourteen covariance structures, not only EEI, VVI, EEE and VVV. The
  ten that constrain the volume, the shape or the orientation across profiles
  are differentiated in their own free coordinates --- log volumes, log shapes
  with the determinant-one constraint built in, and Cayley-transform
  orientations (log-Cholesky factors for VEE and EVV) --- with exactly as many
  dimensions as the structure has parameters, and carried to the reported
  variances and covariances by the delta method. Robust (`"robust"`, `"opg"`)
  errors work for them too. Validated against an independent numerical Hessian
  of the observed-data likelihood in a different chart.
* `multilpa(prior = prior_control())` gives maximum a posteriori estimates
  under the conjugate prior of Fraley and Raftery (2007), mclust's
  `priorControl()`: the same default hyperparameters and the same M-steps as
  `mclust::me()`, for the ten structures mclust defines a prior for.
  `log_likelihood`, and the criteria built from it, are the unpenalized
  likelihood at the posterior mode, as in mclust. Wald inference is refused for
  such a fit; the bootstrap refits with the prior. New condition class
  `latents_unsupported_prior`.
* `multilpa(noise = TRUE)` adds mclust's uniform noise component over the
  data's hypervolume (`hypvol()`), for one group class and continuous, complete
  indicators, with any covariance structure and with or without `prior`. The
  noise is profile `0` in `get_results(fit, "assignments")` (with a
  `posterior_noise` column), `"posteriors"`, `"profile_probabilities"` and the
  classification diagnostics. The component counts two parameters, as in
  mclust, so BIC equals mclust's. New condition class
  `latents_unsupported_noise`, raised also by the verbs that do not yet account
  for a noise component.

## New controls

* `parameter_inference()` and `vcov()` give standard errors for a one-step
  covariate model with categorical or mixed indicators
  (`multilca(profile_covariates = )`, `multilpa(categorical = , profile_covariates = )`).
  Each category's probability is reported with a delta-method standard error
  and no Wald test, as for the covariate-free model.
* `parameter_inference(vcov_type = "opg")` inverts the outer product of the
  group scores (the BHHH estimate). It is the estimator `glca` reports, and
  reproduces its standard errors; `"observed"` remains the default.
* `parameter_inference(boundary = "fix")` holds categorical response
  probabilities that sit on `min_probability` at that bound and reports the
  other parameters conditionally on them. The default, `"error"`, refuses as
  before, now with a message that names the alternative.
* `multilpa()`, `multilca()` and `lta()` gain `select_start`. `"converged"`
  reports the best converged start whenever one converged, instead of an
  unconverged start that edges it out by a negligible likelihood.
* `r3step()` gains `by_group_class`. `TRUE` gives every group class its own
  profile intercepts, with shared slopes and one latent group class per group,
  the two-level form of the one-step model. The pooled regression attenuates
  the slopes when group classes differ in their profile mix.

## Changed defaults

* `r3step()` now defaults to `vcov_type = "robust"`, standard errors clustered
  on groups, matching `three_step()`. Pass `vcov_type = "observed"` for the
  previous behaviour.

## Bug fixes

* A covariate fit now stores `min_probability`, so inference recognises
  response probabilities on their bound and refuses with
  `latents_boundary_fit` instead of a numerically singular information.

# latents 0.8.4

* The website URL in DESCRIPTION ends in a slash, and the README gives the
  CRAN installation command.
* The package check is about a third faster: vignettes use fewer random
  starts, the slowest examples use smaller data, and six exhaustive test
  files are skipped on CRAN (they run on GitHub Actions).

# latents 0.8.3

* The README opens with the package description, and every vignette and
  article names its authors.

# latents 0.8.2

* The package description is rewritten.

# latents 0.8.1

* A documentation website, built with pkgdown, is published at
  <https://pak.dynasite.org/latents/>.

# latents 0.8.0

## Package renamed

* The package is renamed from multilpa to **latents**, a name that covers
  latent profile, latent class and latent transition models and the planned
  extensions. The model functions keep their names (`multilpa()`,
  `multilca()`, `lta()`), as do the result classes. Condition classes now use
  the `latents_` prefix (for example `latents_bad_argument`, previously
  `multilpa_bad_argument`), `multilpa_plot_types()` is `plot_views()`,
  and the conditions catalogue is `?"latents-conditions"`. The vignettes are
  `vignette("lpa", package = "latents")`, `"evaluation"`, `"covariates"`,
  `"lca"` and `"lta"`.

## New features

* `multilca()` fits a two-level latent class model: `multilpa()` with every
  indicator categorical, so the items are named once. Mixed models remain
  `multilpa(categorical = )`.

# multilpa 0.7.0

## Behaviour changes

* `enumerate_classes(structure =)` and `candidate_fit(structure =)` are now
  `model =`, and the candidate grid's `structure` column is `model`. An
  argument `multilpa()` does not take, including the old `structure`, is
  refused with `multilpa_bad_argument` before any candidate is fitted, instead
  of failing inside every candidate.
* Printing an enumeration grid shows AIC, BIC under both conventions, ICL,
  entropy at both levels, and the diagnostics `boundary` and
  `n_best_replicated`.
* `plot()` on an enumeration grid draws several criteria as line plots, one
  panel per criterion and one line per covariance model and group-class
  count; the default shows AIC, both BICs and ICL. `combine = FALSE` draws
  each criterion as a separate figure.

## Documentation

* `?course_engagement` and the case studies now say that the indicators are
  simulated on the `log1p` scale, not transformed from counts.
* New case study on two-level latent class analysis of `student_esm`; the two
  latent transition case studies are merged into one. Case studies and
  vignettes print whole tables instead of filtering them.

## Tests

* The bivariate-residual test for a covariate fit no longer asserts whether
  a separated fixture converges, which differed on R-devel.

# multilpa 0.6.1

## Tests

* Two score tests compared the analytic and numerical gradients at the fitted
  maximum, where both are near zero; they failed on Linux and Windows. They now
  compare at a point away from the maximum, where a wrong score is detected.
* The test suite runs on small data with few starts and iterations: no test
  fits the full bundled dataset, shared fixtures are fitted once, and tests
  whose inequality depended on random starts reaching the global maximum now
  warm-start the wider model. The full suite takes about a fifth of its
  previous time, and slow tests no longer need to be skipped on CRAN.

# multilpa 0.6.0

## New data

* `student_esm`: experience-sampling data on the leisure activities of 100
  university students at 2,582 non-study prompts, with four affect ratings,
  from openESM dataset 0062 (Neubauer & Schmiedek, 2024; CC-BY 4.0). Eight
  yes/no activity items make it the categorical counterpart of
  `course_engagement`, used for two-level latent class and mixed-measurement
  examples.

## Behaviour changes

* One seed rule for every verb: `seed` (and each of `sensitivity()`'s `seeds`)
  is any whole number that `set.seed()` accepts, negative ones included.
  `multilpa()`, `fit_staged()`, `lta()`, the covariate fit, `bootstrap_lrt()`
  and `sensitivity()` previously refused negative seeds, and the bootstrap in
  `parameter_inference()` accepted a fraction that `set.seed()` truncates, so
  seeds 1.2 and 1.7 drew the same stream. Every refusal is now
  `multilpa_bad_argument`.
* Consistent condition classes for the three-step verbs:
  - A measurement indicator offered as an external variable raises
    `multilpa_indicator_reused` from both verbs, alongside
    `multilpa_bad_outcome` from `three_step()` and the new
    `multilpa_bad_covariate` from `r3step()` (previously
    `multilpa_bad_argument`).
  - `r3step()` raises `multilpa_bad_covariate` for a predictor that varies
    within a group at `level = "groups"` (previously `multilpa_bad_outcome`)
    and for a rank-deficient design (previously `multilpa_bad_inference_data`).
  - `data` with the wrong number of rows raises `multilpa_bad_inference_data`
    from `three_step()`, `r3step()` and the bivariate residuals, as it already
    did from the assignments table and `sensitivity()`. These were unclassed.

## Internal

* The fallback row for a transition state with no outgoing moves in
  `get_tna()` is computed for every state in one matrix product rather than a
  loop; the result is bit-identical.

# multilpa 0.5.0

## Behaviour changes

* `get_results(x, "assignments", truth = )` counts each group once when a
  truth column is compared against group classes, rather than once per
  observation, so larger groups no longer weigh more in recovery. Unused factor
  levels are omitted.
* The `multilpa_unverified_alignment` warning now also fires when the supplied
  `data` shares only columns that cannot establish row order, such as a
  repeated group identifier on its own.
* `three_step()` and `r3step()` refuse a membership-covariate fit
  (`multilpa_unsupported_three_step`) and refuse measurement indicators as
  outcomes or predictors. Both now check `data` alignment against the fit.
* `bootstrap_lrt()` requires matching covariance structure and centering,
  repeats grand-mean centering in every refit, and refuses person-centred fits
  (`multilpa_unsupported_bootstrap`).
* `get_tna()` and `get_group_tna()` no longer turn a state with no outgoing
  moves into a certain self-transition. The fit's own row is kept and
  `multilpa_empty_transition_row` is warned.
* `fit_staged(measurement = )` refuses a first stage with more than one group
  class (`multilpa_bad_stage`).
* `sensitivity()` works on `fit_staged()` results, refuses directly fixed fits,
  checks `data` alignment, and requires non-negative whole-number seeds.

## Bug fixes

* `sensitivity()` now refits with the reference fit's own `max_iter` and `tol`;
  `multilpa()` did not store them, so the defaults were always used.
* Relabelling group classes during bootstrap alignment now also reorders
  `group_posteriors`, `group_classes` and `effective_group_counts`.
* Covariate fits no longer store categorical indicators in `indicator_data`,
  which coerced the matrix to character.
* Categorical indicators round-trip with their original type (numeric, factor)
  through `get_results(x, "data")` and bootstrap simulation.
* BCH `three_step()` no longer reports a spurious significant contrast for a
  constant outcome; zero standard errors give `NA` statistics and p-values in
  `three_step()` and `r3step()`.
* `descriptives()` excludes rows with a missing group ID from `n_groups` and
  `icc`, and describes repeated `vars` once.
* Enumeration summaries and plots distinguish covariance structures; the plot
  marks the minimum among converged candidates only.

# multilpa 0.4.1

`vignettes/figure`, a leftover `knitr` figure directory from a preview render,
is added to `.Rbuildignore`. It was entering the tarball and drawing a
"looks like a leftover from 'knitr'" NOTE. The tarball's `vignettes/` now holds
the four `.Rmd` sources and nothing else.

# multilpa 0.4.0

## `get_data()` is renamed `get_results()`

**Breaking, and deliberately without an alias.** The package has never been on
CRAN, so there is no installed base to migrate.

The old name was ambiguous in a way that showed up inside the verb itself:
`get_data(x, what = "data")` used the same word for the accessor and for one of
its own values, and "data" in an R package normally means the input rather than
what a model produced. `get_results()` says what comes back, and keeps the
`get_*` family consistent with `get_tna()` and `get_group_tna()`.

All thirteen S3 methods, the four vignettes and the README move with it. The
`what =` values are unchanged, so `get_results(x, what = "profiles")` returns
exactly what `get_data(x, what = "profiles")` did. Earlier entries in this file
still say `get_data()`, because that is what those versions shipped.

# multilpa 0.3.0

## Two features moved to `future/`

The package now exports 17 verbs rather than 19. The test applied was whether
an export makes a claim the package can back; two did not, and they are kept in
`future/`, which `.Rbuildignore` excludes, with a note on what would bring each
one back.

* **`fit_random_intercept()` is removed.** Its independent validation reached
  only the one-profile limit -- multiple-profile parameter recovery had never
  been checked against anything -- and it refused `parameter_inference()`,
  `vcov()` and `confint()` outright, so nothing in the package could say when
  it was wrong. It is also a different model family: a continuous group effect
  rather than the discrete group classes the package is named for. It appeared
  in no vignette.
* **`lmr_lrt()` is removed.** It returned a `p_value` column that was always
  `NA`, because the Vuong-Lo-Mendell-Rubin reference distribution is not
  reproduced. `bootstrap_lrt()` gives a calibrated p-value for the same
  comparison. The Lo-Mendell-Rubin adjustment factor itself had been verified
  against two genuine Mplus `TECH11` runs; that evidence is recorded in
  `future/README.md` against the day it is restored.

Fifteen S3 methods and twelve help pages go with them. `multilpa_plot_types()`
no longer lists `"random_intercepts"`, and `get_data()` no longer offers the
`"random_intercepts"` table.

Two conditions changed status as a consequence, both of which had the random
intercept as their only trigger:

* `multilpa_quadrature_check` is removed from `?"multilpa-conditions"`. Nothing
  raises it any more.
* `multilpa_no_group_classes` is still raised by guards in `R/diagnostics.R`
  and `R/get-data.R`, but every model family this version fits has discrete
  group classes, so no call can reach it. The catalogue says so rather than
  describing a refusal a caller cannot produce.

## `sensitivity()` asks whether the solution survives a different seed

EM converges to a local maximum from the starts it was given. The package could
report how many starts within one seed's stream reached the best likelihood, and
nothing else: a user could not vary the seed on their own data at all. That left
the package unable to meet, on a user's data, the standard it applies to its own.

`sensitivity(fit, seeds = )` refits under each seed and returns one tidy row per
seed: `log_likelihood`, `converged`, `iterations`, `optimum` (which distinct
maximum it reached, `1` being the best), `best`, and `agreement` (the proportion
of observations assigned as the reference fit assigned them).

`agreement` is computed after each refit's arbitrary profile labels have been
matched to the reference fit's, so two identical solutions that happened to
number their profiles differently report agreement of one rather than zero.

Writing it exposed a defect in `.multilpa_permute_profiles()`, which has been
fixed. It permuted `means`, `variances`, `covariances`,
`profile_probabilities` and `response_probabilities`, but left
`subject_posteriors`, `subject_profiles` and `effective_profile_counts` behind.
The bootstrap it was written for reads only parameters, so a permuted fit that
described one labelling and assigned another went unnoticed. `subject_profiles`
holds labels rather than positions and takes the inverse permutation.

`multilpa()` fits only. `lta()` and covariate fits raise
`multilpa_unsupported_sensitivity`, and a seed whose refit fails contributes an
`NA` row with a `multilpa_sensitivity_dropped` warning naming how many.

## `fit_transitions()` is renamed `lta()`

**Breaking, and deliberately without an alias.** The package has never been on
CRAN, so there is no installed base to migrate and a compatibility alias would
be permanent debt paid for nobody.

The old name described the mechanism; every other fitting verb here is named
for its method. `multilpa()` is latent profile analysis, so latent transition
analysis is `lta()`. Under the old name the capability was undiscoverable by the
name the literature uses: searching the index for "lta" found nothing.

`multilta()` was considered and rejected. It is one character from `multilpa()`,
and the two take the same leading arguments, so a typo would fit a different
model and still return a result.

Nothing else moves. The fitted object is still `multilpa_transitions`, the table
is still `get_data(x, "transitions")`, and the help page is now `?lta`.

## A latent transition fit draws, instead of refusing

`plot()` on an `lta()` fit raised `multilpa_no_plot` for every view. That was a
blanket refusal, not a limitation: the measurement model is the one `multilpa()`
fits, and the fit already carried `means`, `variances`, `subject_posteriors`,
`subject_profiles` and `sequence_lengths`. Six of the seven plot helpers worked
against one unchanged; the seventh was called with the wrong signature.

`plot.multilpa_transitions()` now offers ten views: `"transitions"`,
`"profiles"`, `"bars"`, `"heatmap"`, `"responses"`, `"sequences"`, `"sizes"`,
`"entropy"`, `"posteriors"`, `"avepp"`, and `"all"`.

`what = "transitions"` is new and is this family's own view: the estimated
transition matrix as a heatmap, one panel per group class, row the current
profile and column the next, so the diagonal is persistence. Every cell prints
its probability, and the fill is the white-to-blue ramp the average-posterior
matrix already uses, so darker is a higher probability anywhere in the package.
A row with no data support has its label parenthesised: such a row is uniform by
construction rather than estimated.

`"bars"` draws point estimates with no whiskers here. The interval it draws on a
`multilpa()` fit is a Wald interval, and this family has no standard errors, so
there is none to draw; the subtitle says so rather than the plot implying an
uncertainty it does not have.

`multilpa_plot_types()` lists thirteen views. `multilpa_no_plot` keeps its entry
in `?"multilpa-conditions"`, reworded to what still raises it: `what = "all"` on
an object whose method names no views.

## `get_tna()` and `get_group_tna()` are documented

Both were added in 0.11.7 and appeared nowhere in the README. The README now
carries a section and a feature-matrix row for them, and every one of the 17
exports appears there again.

## Build

`^future$`, `^docs$`, `^vignettes/\.beatrina` and `^vignettes/.*\.tex$` are
added to `.Rbuildignore`. The hidden `vignettes/.beatrina-*` directories were in
`.gitignore` but not in `.Rbuildignore`, so they stayed out of version control
and shipped in the tarball anyway, along with a leftover `vignettes/multilpa.tex`
from a preview render.

`R CMD check --as-cran` on this version reports 1 NOTE, 0 errors and 0 warnings.
The NOTE is "New submission".

# multilpa 0.11.9

## Vignette figures reconciled with the standardized data

Every number in the prose of all four vignettes was checked against the value
its own chunk actually prints, and the stale ones corrected. Seventy-one figures
changed. Three passages needed rewriting rather than renumbering, because the
conclusion moved with the data:

* In the transition guide, the student class with the *highest* probability of
  staying in lower activity is now class 3 at 0.947, where the text had class 3
  as the least persistent at 0.584. The near-boundary transitions are now near
  zero rather than near one, so the caution about sparse cells was rewritten
  around the estimates that are actually at a bound.
* In the evaluation guide, the full-covariance candidate whose best solution
  occurs only once is the three-profile, two-class model, not the three-by-three.
* In the covariate guide, the one-step fit now reports `converged FALSE`: its
  membership logits still carry a non-negligible score at `max_iter`. The text
  says so, and reads the coefficients as the direction of the association rather
  than as a converged maximum.

## Authors and installation

Sonsoles López-Pernas is recorded as an author, and both authors carry their
ORCIDs. The README's installation section now gives
`remotes::install_github("mohsaqr/multilpa")` before the local
`R CMD INSTALL .`, which was the only instruction it offered.

# multilpa 0.11.8

## The introductory vignette, rewritten

`vignette("multilpa")` now opens with what the model found. The order is the
profile means, then how big each profile and class is, then the profile
probabilities, then everything else -- estimation diagnostics, the plots, the
sequences. Printing the fit gives the means one row per profile with each
profile's size, so the guide starts from a result rather than from plumbing.

Two passages are gone. One introduced the notation `pi_{k|m}` and `omega_m`,
said the probabilities lie between zero and one and sum to one, and never used
the symbols again; what it was reaching for is now said where the table is read.
The other listed the model's assumptions in a block, of which the clause that
does any work -- the default diagonal covariance treats the indicators as
independent within a profile, which a shared association between activity
measures can violate -- has moved to Limitations, where a reader can act on it.
A sentence restating the ICC formula's symbols in words went with them.

Every figure in the text was recomputed against the standardized indicators.

# multilpa 0.11.7

## Transition networks, with `get_tna()` and `get_group_tna()`

A fitted transition model hands itself to the tna package:

```r
moves <- fit_transitions(course_engagement, vars = activity, id = "student",
                         time = "sequence", n_profiles = 3, n_group_classes = 3)

get_tna(moves)        # one network for the whole sample
get_group_tna(moves)  # one network per latent group class
```

`get_tna()` returns a `tna` model and `get_group_tna()` a `group_tna`, so every
verb of that package applies: `centralities()`, `communities()`, `cliques()`,
`compare()`, `plot()`. `tna` is in Suggests and is required only by these two.

What is handed over is the **estimates**. tna's own method for a fitted mixture
model, `group_model.mhmm()`, takes the model's grouping and counts transitions
between modal assignments, discarding the estimated matrices; these verbs pass
the matrices themselves, maximised over the full posteriors, so an observation
split 0.6/0.4 between two profiles contributes to both rather than entirely to
one. The two answer different questions and will differ wherever classification
is uncertain.

The whole-sample network is the **marginal** transition matrix: the expected
counts summed across classes and then normalised, not the class-probability
average of the class matrices. A class holding a tenth of the units but a fifth
of the transitions counts for its transitions, which is what a marginal
probability means. A state no observation ever leaves is reported as a
self-transition rather than as `NaN`.

# multilpa 0.11.6

## Single-level fits, with `id = NULL`

`multilpa()` accepts `id = NULL`, which treats the observations as independent
-- each row is its own unit, `n_group_classes` becomes one -- and fits the
ordinary Gaussian or latent-class mixture the two-level model reduces to:

```r
multilpa(faithful, vars = c("eruptions", "waiting"), id = NULL, n_profiles = 3)
```

`id` has **no default**, and omitting it still raises, now with
`multilpa_bad_argument` naming both ways forward rather than R's bare
"argument \"id\" is missing". A single-level model is not the default in a
package for multilevel models, and a forgotten grouping must not quietly become
a different model: on data with 106 students in 1,422 enrolments it would report
1,422 independent observations and change every standard error. Every
`id = NULL` fit raises a `multilpa_single_level` warning saying what it fitted.

Nothing about the likelihood changes, and the fit is identical to the one you
get by numbering the rows yourself and asking for one group class; that
equivalence is asserted in the tests. Every verb of the package works on it,
including `parameter_inference()`, whose cluster bootstrap degenerates to the
ordinary nonparametric bootstrap when each unit holds one row -- which is the
right bootstrap for independent observations -- and `vcov_type = "robust"`,
whose per-group scores become per-row scores.

The unit column is fabricated internally as `.observation`. It is not handed
back by `get_data(x, "data")` or `get_data(x, "assignments")`, and data that
already carry a column of that name raise `multilpa_bad_data` rather than having
it silently overwritten. Asking for more than one group class without an `id`
raises `multilpa_bad_argument`: one observation per unit leaves no composition
for a second-level class to differ in.

A single-level fit prints as what it is --- `Latent profile analysis: 2
profiles` over `160 observations` --- rather than reporting one group class and
160 groups of one.

`fit_staged()` and `fit_transitions()` still require `id`, and say so in their
own documentation rather than inheriting `multilpa()`'s: both have a second
level by construction.

# multilpa 0.11.5

## Every covariance structure now reports uncertainty

`parameter_inference()` gains `method = "bootstrap"`. It resamples *groups*
with replacement, refits inside the same covariance family, undoes the
relabelling each refit comes back with, and reports the percentile interval and
the standard deviation of the replicates:

```r
fit <- multilpa(course_engagement, activity, id = "student", n_profiles = 2,
                n_group_classes = 1, volume = "varying", shape = "equal",
                orientation = "axis")
parameter_inference(fit, method = "bootstrap", iter = 999, seed = 1)
```

This closes the gap the fourteen structures opened. Ten of them --- EII, VII,
VEI, EVI, VEE, EVE, VVE, EEV, VEV and EVV --- constrain the volume, the shape
or the orientation, which the Wald coordinates cannot express, so the Wald path
refused them. `enumerate_classes(structure = )` could therefore select a model
the package would not then do inference on: on this package's own
`course_engagement` data the best of the fourteen is `VEI`. Every structure the
grid can select now reports an interval.

`vcov()` and `confint()` reach the same path through `...`. `confint()` returns
the percentile interval the table reports rather than rebuilding
`estimate +/- z * se` from the bootstrap standard error, so one fit does not
give two different intervals.

Groups are the resampling unit, not rows, so the interval carries the same
independence assumption as `vcov_type = "robust"` and not the stronger one
`"observed"` makes. Against the cluster-robust sandwich on `course_engagement`,
where both are available, the bootstrap standard errors agree to within 11%
across every estimated parameter. `statistic`, `p_value` and `p_adjusted` are
`NA` on this path: the interval is the inference, and putting a normal
approximation back on top of the replicates it was read from is the assumption
the path exists to avoid.

A fit that holds a measurement block is refused with
`multilpa_unsupported_inference`, because the held values came from another fit
and resampling these data does not resample them.

## Fixes

* `fit_staged()` accepted a measurement fitted with a constrained covariance
  structure and silently fitted the second stage as the nearest structure
  `variance_model` and `covariance_model` can name: a `VEI` measurement came
  back recorded as `VVI`. The held means and variances were the right numbers,
  but under the wrong model, and `n_parameters_with_measurement` counted the
  spread as `VVI`'s. `multilpa()` already refuses to hold the variances of such
  a structure, so a staged fit of one was never supported; it is now refused up
  front with `multilpa_bad_stage`, naming the structure. Refit the measurement
  as EEI, VVI, EEE or VVV, or fit both levels together with `multilpa()`, which
  estimates all fourteen directly.

* `bootstrap_lrt()` refitted each replicate from `variance_model` and
  `covariance_model` alone, which cannot express a constrained structure: a
  `VEI` model came back as `VVI`, two parameters wider, so the reference
  distribution belonged to a different pair of models than the statistic
  compared against it. Replicates are now refitted from the structure the model
  records.

## Corrections to the 0.11.3 notes below

The 0.11.3 section was written in two passes and the earlier one was left
standing. Two of its statements are wrong, and are corrected here rather than
edited away:

* "Six of `mclust`'s fourteen models --- the ones with a non-axis-parallel
  orientation --- remain unavailable" was superseded within that same release
  by "All fourteen mclust covariance structures". All fourteen ship, and their
  parameter counts match `mclust`'s own.
* The refusal list "`parameter_inference()`, `vcov()` and `confint()` refuse
  EII, VII, VEI and EVI" named four of the ten that were actually refused. The
  rule, rather than the list: Wald inference was available exactly for the four
  structures reachable by `variance_model` and `covariance_model` alone --- EEI,
  VVI, EEE and VVV. As of this release the bootstrap covers all fourteen.

# multilpa 0.11.4

## Equivalence tests no longer ship with the package

Every test that compares multilpa against other software or published output
now lives in `equivalence/` in the source repository, and that folder is not
part of the built package. This covers Mplus, mclust, tidySEM, depmixS4 and
the in-house enumeration. The shipped `tests/testthat/` keeps unit and
invariant tests only. Run the comparisons with `Rscript equivalence/run.R`.

`mclust` and `tidySEM` are no longer in `Suggests`; nothing shipped uses them.

## Latent transition references

`fit_transitions()` is now checked against the external references pinned by
the JStats project: depmixS4 Gaussian and binary panel models, two Mplus LTA
runs, and the plain-LTA row of Table 5 of Muthen & Asparouhov (2022,
*Psychological Methods*). Log-likelihoods, criteria and parameters agree to
the reference's precision. On one binary depmixS4 fixture (K = 3), the
recorded depmixS4 solution is a local maximum: multilpa reaches a likelihood
0.016 higher. A brute-force path-sum likelihood confirms both values.

# multilpa 0.11.3

## Within-person profiles, without doing the transform yourself

`multilpa(centering = "person")` subtracts each group's own mean from its rows
before the measurement model sees them, so the profiles become profiles of
*change* rather than of level. This is the person-mean-centred, within-person
design (Quintana, 2021; Voelkle, Brose, Schmiedek, & Lindenberger, 2014), and
until now it meant transforming the data frame first --- which left the fit
holding numbers that no longer matched the frame it was given, so every verb
that checks row alignment was comparing centred values against raw ones.

```r
multilpa(course_engagement,
         c("browse", "lectures", "forum_read", "forum_post", "attendance"),
         id = "student", n_profiles = 3, n_group_classes = 1,
         centering = "person")
```

The offsets stay on the fit, so nothing downstream has to guess which scale it
is looking at: `get_data(x, "data")` returns the columns you supplied,
alignment is still checked against them, and the bivariate residuals are
computed on the scale the model was estimated on. `"grand"` subtracts one mean
per indicator instead. `"none"` is the default and changes nothing.

Centring removes exactly the between-unit variation, so `"person"` refuses with
`multilpa_bad_data` when it leaves an indicator constant --- which is what
happens when a unit has one observation of it. With `"person"` the group
classes become types of *change pattern*, not types of unit.

## All fourteen mclust covariance structures

`volume`, `shape` and `orientation` are the three pieces a covariance
decomposes into --- `Sigma_k = lambda_k * D_k * A_k * D_k'` --- and naming
them reaches every model `mclust` fits, where `variance_model` and
`covariance_model` reached four.

| `volume` | `shape` | `orientation` | model |
|---|---|---|---|
| equal / varying | spherical | --- | EII, VII |
| equal / varying | equal / varying | axis | EEI, VEI, EVI, VVI |
| equal / varying | equal / varying | equal | EEE, VEE, EVE, VVE |
| equal / varying | equal / varying | varying | EEV, VEV, EVV, VVV |

All three default to `NULL`, which follows `variance_model` and
`covariance_model`, so a call that names none of them fits exactly what it did
before. An ellipsoidal structure carries a full covariance array whatever
`covariance_model` said, because an orientation cannot live in a diagonal
block.

The estimates are those of Celeux and Govaert (1995). Ten have a closed form
or a scalar fixed point; EVE and VVE share one orientation across profiles
while letting the shapes differ, which has no closed form, and use the
minorize-maximize step of Browne and McNicholas (2014) --- majorize each trace
term linearly, then solve the resulting orthogonal Procrustes problem by
singular value decomposition.

**Verified against `mclust`.** Parameter counts match its own `df` for all
fourteen. Against `mclust::mstep()` on the same responsibilities, twelve agree
to 1e-11 or better. EVE and VVE differ by at most 5.7e-4 --- and attain a
*lower* value of the objective the M-step minimizes, by 2.1e-7 and 2.8e-5, so
the difference is `mclust` stopping short rather than a disagreement about the
estimate.

## Selecting over structures as well as class counts

`enumerate_classes(structure = )` crosses the covariance structures with the
class counts, which is the grid model-based clustering is usually selected
over, and adds a `structure` column to the candidate table:

```r
enumerate_classes(course_engagement, activity, id = "student",
                  n_profiles = 2:4, n_group_classes = 1,
                  structure = c("EEI", "VEI", "EVI", "VVI"),
                  centering = "person")
```

`candidate_fit()` gains `structure` to name one candidate out of a grid where
the class counts alone name several; asking without it raises
`multilpa_unknown_candidate` listing the structures it could have meant.

## Also new

The two structures whose M-step iterates over orientations now warm start from
the previous iteration's answer, which is worth about 1.2 to 1.3 times on them.
An M-step that stops at its iteration cap is **counted and reported once per
fit, with the count**, rather than warning from inside every M-step: one step
in three hundred is a different thing from all of them, and a warning raised
hundreds of times says neither.

`get_data(x, "assignments")` gains an `uncertainty` column: one minus the
posterior of the profile each row was assigned to, which is exactly what modal
assignment discards.

## The constrained diagonal covariance family

`variance_model` ties a profile's covariance volume to its shape: `"equal"`
constrains both and `"varying"` frees both. `volume` and `shape` free them
separately, which makes the six axis-parallel `mclust` models reachable
instead of two:

| `volume` | `shape` | model | spread parameters |
|---|---|---|---|
| equal | spherical | EII | 1 |
| varying | spherical | VII | K |
| equal | equal | EEI | d |
| varying | equal | VEI | K + d − 1 |
| equal | varying | **EVI** | 1 + K(d − 1) |
| varying | varying | VVI | Kd |

Each profile's covariance decomposes as `Sigma_k = lambda_k * A_k`, a volume
`lambda_k = |Sigma_k|^(1/d)` and a shape `A_k` with determinant one; the
estimates are those of Celeux and Govaert (1995). `volume`/`shape` default to
`NULL`, which follows `variance_model`, so every existing fit is unchanged.

Verified against `mclust`: the parameter counts match its own for all six
models, and at `mclust`'s maximum-likelihood estimate the log likelihood agrees
with one written from the mixture definition to 1e-10. The maximized
likelihoods are monotone along the whole nesting lattice.

Three consequences worth knowing:

* `parameter_inference()`, `vcov()` and `confint()` refuse EII, VII, VEI and
  EVI with `multilpa_unsupported_inference`. The free coordinates this package
  differentiates are log variances, one per profile and indicator, which is the
  wrong chart for a constrained volume or shape --- it has more coordinates than
  the model has parameters. EEI, VVI, EEE and VVV are unaffected.
* A constrained structure is maximized across every profile at once, so it
  cannot be combined with a held `variances` block; `fixed = "means"` is fine.
* A `start` taken from a wider structure is projected onto the requested family
  before the first iteration, because EM only promises a likelihood that does
  not decrease *within* the family it is maximizing over.

`covariance_model = "full"` constrains neither volume nor shape, so naming
either argument with it raises `multilpa_bad_argument`. Six of `mclust`'s
fourteen models --- the ones with a non-axis-parallel orientation --- remain
unavailable.

# multilpa 0.11.2

## One verb for every table

Eleven table verbs and the `what =` argument of seven `as.data.frame()`
methods are **replaced by one generic**, `get_data(x, what = )`. A reader no
longer has to know which verb owns which table before they can ask for it:

```r
get_data(fit)                            # the primary table, the measurement model
get_data(fit, "entropy")                 # was entropy_table(fit)
get_data(fit, "transitions")             # was transitions(fit)
get_data(fit, "profile_probabilities")   # was as.data.frame(fit, what = "...")
names(get_data(fit, "all"))              # every table this fit can produce
```

`what = "all"` returns a named list of every table, built from the same
definitions a single `what` uses, so the two cannot disagree. A table this
particular fit cannot produce -- sequences for a fit made without `time`,
bivariate residuals for a family with no discrete group classes -- is left out
rather than erroring, so the names of the list say what was available.

`summary()` now carries every table and prints all of them, truncated to
`rows = 10` each with the `get_data()` call that returns the rest.
`as.data.frame()` is plain coercion to the primary table; it no longer takes
`what`, and an argument it cannot use raises `multilpa_bad_argument` naming it
rather than quietly returning a different table.

### Migration

| Was | Now |
|---|---|
| `assignments(fit, data)` | `get_data(fit, "assignments", data = data)` |
| `average_posteriors(fit)` | `get_data(fit, "average_posteriors")` |
| `bch_weights(fit)` | `get_data(fit, "bch_weights")` |
| `bivariate_residuals(fit, data)` | `get_data(fit, "residuals", data = data)` |
| `classification_errors(fit)` | `get_data(fit, "classification_errors")` |
| `classification_table(fit)` | `get_data(fit, "classification")` |
| `entropy_table(fit)` | `get_data(fit, "entropy")` |
| `information_criteria(fit)` | `get_data(fit, "information_criteria")` |
| `sequences(fit)` | `get_data(fit, "sequences")` |
| `sequence_summary(fit)` | `get_data(fit, "sequence_summary")` |
| `transitions(fit)` | `get_data(fit, "transitions")` |
| `as.data.frame(fit, what = "x")` | `get_data(fit, "x")` |
| `as.data.frame(summary(fit), what = "fit")` | `get_data(fit, "model")` |
| `as.data.frame(diagnostics(fit), what = "posteriors")` | `get_data(diagnostics(fit), "average_posteriors")` |
| `fit_covariates(...)` | `multilpa(..., profile_covariates = )` |

Three tables that could only be reached through a `summary()` object are now
on the fit itself: `"counts"` (effective class memberships), `"covariances"`
(the within-profile residual covariance matrices) and `"model"` (the one-row
fit summary, previously `what = "fit"`). An enumeration grid gains
`"criteria"`, the per-criterion minima its summary already printed.

### Two defaults changed

`"classification"`, `"average_posteriors"`, `"classification_errors"` and
`"bch_weights"` default to **both levels** on a fit that has discrete group
classes, where the old verbs defaulted to `"individuals"`. This is the rule
`diagnostics()` already followed. Pass `level = "individuals"` for the old
behaviour. `"classification_errors"` and `"bch_weights"` accept
`level = "both"` for the first time; weight a regression with one level, not
with both.

`as.data.frame()` on a fitted transition model returns the **measurement
model**, not the transition matrix, because the primary table is now the same
kind of thing for every fitted family. `get_data(fit, "transitions")` is the
transition matrix.

## Covariates are arguments to `multilpa()`, not a separate verb

`fit_covariates()` is **removed**. Membership covariates are part of the model
`multilpa()` fits, so they are arguments to it:

```r
multilpa(course_engagement, activity, id = "student", n_profiles = 2,
         n_group_classes = 2, profile_covariates = "previous_grade")
```

The result is a `multilpa_covariates` object exactly as before, and it records
the `multilpa()` call the caller wrote. The covariate likelihood has no
observed-data form, no starting-value contract of the shape the covariate-free
EM uses and no held-measurement machinery, so `start`, `missing = "fiml"` and
`fixed` are refused by name with `multilpa_bad_argument` rather than accepted
and ignored.

## Every plot a fit supports, in one call

`plot(x, what = "all")` draws every view the fit has the ingredients for,
re-issuing the caller's own call once per view so style arguments carry into
each panel. A view that refuses is named at the end rather than stopping the
sequence.

## Recovery is one call

`get_data(x, "assignments", truth = )` cross-tabulates the model's labels
against known ones, so checking a fit against a generating truth no longer
needs `xtabs()` on the assignment frame:

```r
get_data(fit, "assignments", data = course_engagement,
         truth = c("engagement", "student_type"))
```

Each truth column is compared against the level it describes, read off the
data rather than guessed: a column taking one value within every group is a
property of the group and goes against `group_class`, and one that varies
inside any group goes against `profile`. The `assignment` column records which
comparison was made.

## A fitted model prints its own estimates

`print()` on any fitted model now shows the measurement model under the header,
so `fit` alone is the whole result. `fit` followed by `as.data.frame(fit)` was
two calls showing one thing.

```r
fit <- multilpa(course_engagement,
                c("browse", "lectures", "forum_read", "forum_post", "attendance"),
                id = "student", n_profiles = 2, n_group_classes = 2)
fit
#> Two-level latent profile analysis: 2 profiles, 2 group classes
#> ...
#>  profile  indicator  mean variance standard_deviation
#>        1     browse 2.854   0.3521             0.5934
#>        ...
#> Every other table: get_data(x, what = ), or get_data(x, "all").
```

An all-categorical fit prints its response probabilities instead, rather than
an empty block under a heading promising means. `rows` bounds the printout;
`get_data()` returns any table whole.

## Also

* `get_data(x, "model")` renames the plain fit's `bic` column to `bic_groups`,
  which is what the covariate and random-intercept families already called it,
  so one verb no longer spells the same convention two ways. The table is
  family-shaped rather than a padded union: a random-intercept fit reports its
  integration diagnostics and a covariate fit its predictor counts.
  `"information_criteria"` is the table with one column set for every family,
  and is the one to compare fits across families with.
* `report()` gains `rows`, forwarded to `print(summary(x))`, because a first
  look at a fit with thousands of observations would otherwise be mostly
  posteriors.
* A covariate model with **no** covariates can no longer be requested: naming
  no covariate is a request for the covariate-free model. The two differ only
  in parameterisation, the intercept-only multinomial logits standing in for
  the mixing probabilities.
* `report()` and `diagnostics()` are unchanged in purpose. `diagnostics()`
  reports its tables through `get_data()` and renames its average-posterior
  table from `"posteriors"` to `"average_posteriors"`, which no longer
  collides with the individual posteriors.
* `multilpa_removed_argument` is **removed** from the condition catalogue. It
  could only be raised by `classification_table(detail = )`, and that verb no
  longer exists, so nothing can raise it. `?"multilpa-conditions"` no longer
  claims otherwise.
* A three-step method called on a random-intercept fit now refuses with
  `multilpa_no_group_classes` instead of an unclassed error, so
  `get_data(x, "all")` can tell a table that does not apply from a defect.

# multilpa 0.11.1

## One bundled dataset, covering every model family

`school_engagement` and `engagement_panel` are **removed** and replaced by
`course_engagement`: 1422 enrolments of 106 students across 32 courses, ordered
within each student. The same rows now carry every model the package fits.

```r
activity <- c("browse", "lectures", "forum_read", "forum_post", "attendance")

multilpa(course_engagement, activity, id = "student",
         n_profiles = 2, n_group_classes = 2)          # two-level
fit_transitions(course_engagement, activity, id = "student",
                time = "sequence", n_profiles = 2)      # the same students, in order
fit_covariates(course_engagement, activity, id = "student",
               n_profiles = 2, n_group_classes = 1,
               profile_covariates = "previous_grade")   # and what predicts it
```

The old pair could not do this. `school_engagement`'s `term` column was not
time -- it indexed twelve different students within a school -- so anything
about order needed the second dataset and a different set of people, and
neither dataset carried a covariate at all.

Two new columns make that possible. `previous_grade` is the standardised grade
the student earned in the course *before* this one, so it predicts engagement
rather than summarising it. `attendance` is the count of days the student was
active, following the same engagement as the click measures but recorded one
tick per day rather than one per click, so it carries more noise than they do.

The generating truth still ships beside the observations -- `engagement` for the
enrolment pattern and `student_type` for the kind of student -- so a fitted
model can still be checked against what produced it:

```r
xtabs(~ profile + engagement, data = assignments(fit, data = course_engagement))
```

The truth column is called `engagement` rather than `profile` because
`assignments()` adds a column named `profile`, and a dataset already carrying
that name would make the verb refuse to join its own output onto it.

The data are simulated. The simulation is calibrated to the shape of a real
learning-analytics export, using only aggregate constants from it: marginal
means and spreads on the log scale, the separation between an engaged and a
disengaged pattern, and how persistent that pattern is across a student's
courses. No row, identifier or value of any real student is present, and the
source is not distributed. See `data-raw/course-engagement.R`.

The activity indicators are natural-log counts, `log1p(events)`. Raw
learning-analytics counts are strongly right-skewed and on an unlogged draft the
skew alone looked like an extra class, exactly as a floor effect does.

This release is a correctness release. An independent review of 0.10.0 found
that several advertised paths returned plausible numbers that were wrong, rather
than failing. Every one of them is fixed below, each with a regression test that
reproduces the original defect. Nothing here adds a model family.

## Silently wrong results, now correct

### A covariate's units no longer change the fit

`fit_covariates()` ran its membership M-step on the coefficients as supplied.
BFGS stops on a *relative* tolerance measured in the coordinates it is given, so
a predictor in small units stalled the search at a point where the score was
still large, and `optim()` reported success. Rescaling one covariate by `1e-7`
moved the log likelihood from `-620.1281` to `-730.5748` and drove the slope to
zero, with no warning. The search now runs on root-mean-square-scaled
coordinates and convergence is judged by a unit-free score, not by the
optimiser's return code. Estimates are returned on your scale, unchanged:

```r
# identical to machine precision at every scale
fit_covariates(data, "y", "g", 2, profile_covariates = "z")
```

### Staged categorical measurement attaches to the right item

A held measurement was passed between stages as an unlabelled, positional list,
while the compatibility check compared only the *sorted set* of indicator names
-- which a permutation passes. Reordering `categorical = c("a", "b")` to
`c("b", "a")` swapped whole response distributions (`-587.5909` to `-671.5937`);
reversing a factor's levels swapped categories within an item (`-1047.9689`).
`starting_values()` now keeps the indicator names and category labels, and
`multilpa()` aligns on them. An encoding that cannot be aligned raises
`multilpa_bad_start` or `multilpa_bad_stage` instead of returning a number.
An unlabelled start is still read positionally, so stored starts keep working.

### Bootstrap comparisons keep the constraints they were given

`bootstrap_lrt()` accepted fixed-measurement models but dropped `fixed` and the
held values from every refit, so the observed statistic compared constrained
fits while the null distribution came from unconstrained ones -- a 1-versus-3
parameter statistic referred to a 9-versus-11 parameter null. Refits now carry
the constraint, and a pair whose constrained nesting cannot be established is
refused with `multilpa_bad_nesting`.

### Transition counts no longer overflow

Expected transition counts exponentiated a scaled forward factor and an unscaled
backward factor separately, so on a long sequence one overflowed to `Inf` where
the other underflowed to `0`. A converged 600-occasion fit reported `NaN` for
every count, and fitting 1,800-occasion sequences failed outright. The complete
pair log probability is now formed before exponentiating, which cannot overflow.
Counts agree with explicit path enumeration to `1.8e-15`.

### `descriptives()` no longer invents or replaces values

An `NA` stratifier produced `NA` row selectors, and `x[NA, ]` fabricates a row
rather than dropping one: every stratum gained a phantom missing observation
while the real row vanished. Missing strata are now a labelled stratum, and the
per-stratum counts are asserted to account for every input row. Separately, the
profile assignment was written into the same frame as the indicators, so an
indicator named `profile` or `group_class` was overwritten by its own class
labels. The assignment is now kept beside the frame, never in it.

### Cluster-robust standard errors are refused when they cannot exist

`three_step()` and `r3step()` summed estimating-equation contributions within
groups without requiring enough independent groups. Those contributions sum to
zero at the estimate -- that sum *is* the stationarity condition -- so a single
group gave a standard error of `5.9e-17` and a confidence interval of zero
width, and `r3step()`'s robust path reported p-values of exactly zero. Both now
require more independent groups than the quantities reported, raising
`multilpa_too_few_groups`. `three_step(vcov_type = "independent")` offers an
explicit, labelled unclustered alternative.

### `diagnostics()` no longer ignores the argument you passed

`diagnostics()` and `report()` documented that `...` reached the verbs they
gather, but discarded it. Because `by` is a real argument of
`bivariate_residuals()`, `diagnostics(fit, by = "overall")` silently returned
per-profile residuals under the name you asked to pool. Both verbs now take `by`
explicitly and refuse any argument they cannot forward.

### A bootstrap comparison stops withholding its p-value

`bootstrap_lrt()` judged a replicate invalid if its statistic fell below a fixed
`-1e-5`. But the statistic is `2 * (ll_alt - ll_null)` on log likelihoods of
order 10^3 to 10^4, while EM stops on a *relative* tolerance, and under the null
the alternative converges to the null solution -- so every replicate sits at the
boundary with noise of about `2 * tol * |log likelihood|`. On
`course_engagement` that is 1.6e-04, sixteen times the window. Replicates that
had all converged were reported as "Nonconvergence or reversed likelihood" and
the p-value withheld. The window now scales with `tol` and the likelihood, so
the same comparison returns every replicate valid and a p-value.

### An ordinary converged fit stops being told it has not converged

`parameter_inference()` cut a fixed `0.01` on the scaled score. On 720
observations at the default `tol = 1e-8` that fires on a fit which converged
with all 20 starts replicating. The criterion is now the displacement the score
implies, expressed in the estimate's own standard errors: a score `g` against
information `I` moves the estimate by `g / I`, and the standard error is
`sqrt(1 / I)`, so `g * SE` is that displacement and is dimensionless. The
default-tolerance fit sits 0.17% of a standard error from stationarity and is
now quiet; a fit at `tol = 1e-4` sits 5% away and still warns, now saying so in
those words.

### Misaligned data is refused, not scored

`assignments()` documented that it checked row alignment but compared only the
row count, and `bivariate_residuals()` trusted row order. A reversed frame was
paired with the wrong posteriors and returned as a result. Both now compare
every column the fit recognises and raise `multilpa_bad_inference_data` on a
mismatch; a frame sharing no column with the fit warns
`multilpa_unverified_alignment` rather than passing in silence.

## Now supported

* `parameter_inference()`, `vcov()` and `confint()` work on fixed and staged
  fits, reporting standard errors conditional on the held measurement -- which
  `?fit_staged` already promised. They previously failed with
  `length(theta) == object$n_parameters is not TRUE`. Held blocks enter the
  likelihood at their held values and contribute no row to the information
  matrix; the conditional likelihood matches an independent enumeration to
  `5.7e-14` and the standard errors match an independent numerical information
  matrix to `1.3e-07` relative.
* `assignments()` and `descriptives(by = "profile")` work on a random-intercept
  fit. Both were documented to, but the fit stored no modal assignment, so the
  first failed with an unclassed `cbind()` error.
* `multilpa_plot_types()` lists `"random_intercepts"`, which
  `plot.multilpa_random_intercept()` has always accepted.
* `fit_transitions(max_iter = 0)` returns a usable fit instead of failing with
  "attempt to set an attribute on NULL".
* `plot()` draws `what = "entropy"` and `what = "posteriors"` for covariate and
  random-intercept fits, so `plot(diagnostics(fit))` works for every family.

## Behaviour changes

* `max_iter = 0` with a supplied `start` now evaluates that start alone and
  ignores `n_starts`, which is what `?multilpa` promised. It previously scored
  all random starts and returned the best of them.
* Restart selection now breaks ties deterministically rather than by
  `which.max()`, which was choosing between equally valid optima on
  floating-point noise.
* A separated covariate fit is now reported as unconverged rather than
  converged, so `parameter_inference()` refuses it instead of returning Wald
  intervals around an unbounded coefficient.
* Every warning the package raises now carries a documented condition class, so
  a simulation loop can muffle the qualification it expects and let the rest
  through. New classes: `multilpa_extreme_coefficients`,
  `multilpa_quadrature_check`, `multilpa_empty_transition_row`,
  `multilpa_unverified_alignment`, `multilpa_no_free_parameters`,
  `multilpa_held_parameter`.
* `three_step()` gains `vcov_type`; `diagnostics()` and `report()` gain `by`;
  `plot.multilpa_covariates()` and `plot.multilpa_random_intercept()` gain
  `"entropy"` and `"posteriors"`. All are additive and keep positional calls.
* `multilpa_bootstrap_lrt` objects carry `fixed`, and
  `as.data.frame(what = "test")` gains a `fixed` column.
* `parameter_inference()` on a fixed fit reports no row for a held block and
  carries a `fixed` attribute. The table is on the natural scale, where every
  probability in a set is reported, so it has `n_parameters` rows plus one for
  each reference category the estimation scale drops. `vcov()` is
  `n_parameters` square on the unconstrained scale; on the natural scale it is
  correspondingly larger and singular.
* `descriptives(by = )` can return one additional row: the labelled missing
  stratum.
* `profile_prevalence` for a transition fit is read from the final expectation
  rather than the last M-step. These agree at convergence; the new value is the
  correct one for an unconverged or zero-iteration fit.
* `starting_values()` keeps names and column labels on categorical response
  blocks, and `as.data.frame(what = "responses")` now reports those labels
  rather than positions, matching what the fitted object's own response table
  reports for the same quantity. A start carrying no labels is still read
  positionally and its tidy view still reports positions.
* `assignments()` raises `multilpa_bad_inference_data` rather than
  `multilpa_bad_nesting` for a row-count mismatch, which is a data contract
  failure, not a model-nesting one.
* Thirteen further user-reachable guards in `bootstrap_lrt()` and
  `fit_transitions()` now raise their catalogued class instead of a bare
  `simpleError`, so `?"multilpa-conditions"` is true of them.
* A `multilpa(fixed = )` fit records `n_parameters_with_measurement`, which only
  `fit_staged()` set before. `print()` consequently reported the same count
  twice on such a fit; it now reports the two counts it names.
* `parameter_inference()` results carry `score_displacement` beside
  `scaled_score`, the quantity the convergence warning is judged on.
* `bivariate_residuals()` gains `adjust`, and its table a `p_adjusted` column.
  Every pair of indicators is a separate test, so a five-indicator fit asks ten
  questions at once; the family size is now the verb's to declare rather than
  the reader's to reconstruct. The default is `"none"`, so `p_adjusted` equals
  `p_value` unless asked otherwise.

## Infrastructure

* GitHub Actions check the package on macOS, Windows and Linux across release,
  devel and oldrel, plus an explicit R 4.1 job for the declared `Depends` floor;
  a second workflow runs `lintr`.
* `DESCRIPTION` gains `URL` and `BugReports`.

# multilpa 0.10.0

## New features

### A fit carries the data it was built from

Post-fit verbs no longer ask for data the fit is already holding. The model
frame is rebuilt from the uncentred indicators, the categorical codes and their
levels, the identifier and the occasion — all of which were stored already — so
this costs no extra memory, and the round trip is exact for every model family.

```r
parameter_inference(fit)          # was parameter_inference(fit, students)
bivariate_residuals(fit)
vcov(fit); confint(fit)
bootstrap_lrt(smaller, larger)
as.data.frame(fit, what = "data") # the columns the model was fitted to
```

Supplying `data` still works and is still the stronger check, since a frame that
does not reproduce the fitted likelihood is rejected as before.
`three_step()` and `r3step()` keep requiring `data`, because an outcome or a
covariate is a column the model never saw. The `multilpa_data_required`
condition is retired: there is now something sensible to default to.

### Three convenience verbs

* `descriptives(x, vars, id)` — one row per variable with n, missing, mean, sd,
  range, and the **intraclass correlation**: the share of variance lying between
  groups. An ICC near zero says the groups do not differ, so a model built to
  tell them apart has nothing to find. It accepts a data frame before fitting or
  a fit afterwards, and `by = "profile"` splits it by assigned profile.
* `diagnostics(x)` — entropy, classification quality, average posteriors and
  bivariate residuals in one call, as a classed object with `print()`,
  `as.data.frame(x, what = )` and `plot()`. `plots = TRUE` draws them too. A
  model family that cannot supply a table says so rather than failing.
* `report(x)` — summary, descriptives, diagnostics and every plot the fit
  supports, in one call. It skips views a given fit has no ingredients for
  rather than stopping at the first refusal.

### `assignments()`

Every observation with the profile it was assigned to, the group class, and its
posteriors -- optionally with a frame of your own columns kept alongside:

```r
labelled <- assignments(fit, data = course_engagement)
xtabs(~ profile + engagement, data = labelled)
```

Comparing an assignment with anything else means putting them in the same row,
and doing that with two separate objects assumes they share an order. This verb
owns the alignment: a frame with the wrong number of rows is refused, and a
column the assignments would overwrite is an error rather than a replacement.

It is a verb rather than another `as.data.frame(what = )` value because it
computes a join with a frame the caller supplies. The line: `as.data.frame()`
represents what a fit already holds, and anything bringing in columns the model
never saw gets its own verb with its own documented arguments.

### `plot(fit, what = "bars")` draws its intervals without being handed data

The fit carries the columns it was built from, so the 95% intervals are drawn by
default; `data` is now only an override. A model family whose standard errors
are not implemented gets bars without whiskers rather than an error.

`report()` likewise asks each `plot()` method what views it accepts instead of
assuming the full catalogue, and skips a view the method refuses -- naming what
it could not draw rather than dropping it silently.

## Bug fixes

* The intraclass correlation is now `1` rather than `NA` when a variable has no
  within-group variance. That is a well-defined case — every observation equals
  its group's own value — not a degenerate one.

# multilpa 0.9.0

## Breaking changes

### `information_criteria()` returns the shape people report

The default is now one row, with one column per criterion — the shape a
model-comparison table is published in, and the same column names
`enumerate_classes()` uses, so a chosen model's row sits under the grid it was
chosen from without renaming anything.

```r
information_criteria(fit)
#>   log_likelihood n_parameters      aic      kic bic_groups bic_individual ...
#> 1      -3778.259           15 7586.518 7604.518   7617.933       7655.207 ...

information_criteria(fit, format = "long")   # the previous shape
```

* The `penalty` column has been removed. No convention reports it — not Mplus,
  not `tidyLPA::get_fit()`, not mclust — and it is the difference between two
  numbers the table already carries.
* `definitions = TRUE` describes one criterion per row, so it now requires
  `format = "long"` and raises `multilpa_bad_argument` otherwise.
* `as.data.frame(fit, what = "information_criteria")` follows the verb, and its
  `format` argument now actually reaches it. Previously `format` was captured by
  the accessor's own signature and never forwarded.
* `convention` remains `NA_character_` for `deviance`, `aic` and `kic`, now
  documented as meaning *the question does not arise* rather than *a value is
  missing*.

## New features

* Two example datasets ship with the package, so the first example runs as
  written: `school_engagement` (720 students in 60 schools) and
  `engagement_panel` (120 students over four waves). Both are simulated and both
  carry the kind each row was generated from — `engaged` and `state` — so a
  fitted model can be checked against what produced it.
* `vignette("multilpa")` is rebuilt on those datasets and gains a section
  cross-tabulating the fitted profiles against the truth column.

# multilpa 0.8.0

## Breaking changes

### Argument names now follow the sibling package `Nestimate`

No deprecation shim: the old spellings are gone and raise
`unused argument`.

| was | is | where |
|---|---|---|
| `object` | `x` | first formal of this package's own verbs |
| `indicators` | `vars` | `multilpa()`, `fit_covariates()`, `fit_staged()`, `fit_transitions()`, `fit_random_intercept()`, `enumerate_classes()` |
| `group` | `id` | the same six |
| `level_ci` | `ci_level` | `three_step()`, `r3step()` |
| `p_adjust` | `adjust` | `parameter_inference()`, `three_step()`, `r3step()` |
| `n_boot` | `iter` | `bootstrap_lrt()` |
| `profiles`, `group_classes` | `n_profiles`, `n_group_classes` | `enumerate_classes()` |

```r
# before
multilpa(students, indicators = vars, group = "school", n_profiles = 2)
# after
multilpa(students, vars = vars, id = "school", n_profiles = 2)
```

The S3 methods for the base generics `summary()`, `coef()`, `vcov()`,
`confint()`, `logLik()` and `nobs()` keep `object`, because those generics own
that name.

**Tidy output columns did not change.** Every posterior, sequence and
classification table keeps its `group` column, `level` keeps the value
`"group"`, and fitted objects keep `group_classes` and `n_group_classes`.

## New features

* `multilpa_plot_types()` lists every view the `plot()` methods accept, with the
  group each belongs to and what it answers. It carried `@export` in 0.7.0 but
  was missing from `NAMESPACE`, so it could not be called.
* `plot()` gains four views: `"bars"` (profile means as grouped bars, with 95%
  intervals when `data` is supplied), `"heatmap"` (means as deviations from each
  indicator's grand mean), `"entropy"` and `"posteriors"` (per-case
  classification quality as ridges).

## Bug fixes

* A refused plot no longer leaves the graphics device unusable. Restoring `mfg`
  as part of `par(no.readonly = TRUE)` switches `new` on as a documented side
  effect, so a view that validated and raised before drawing left every
  subsequent plot ready to overlay onto a blank panel.
* Direct labels are no longer clipped at the right edge on narrow panels. The
  gap was measured in user units while the margin reserved for it was measured
  in inches.
* Point area in `plot(fit)` encodes prevalence absolutely rather than being
  min-max normalised within the plot, which had mapped any spread onto the full
  size range and so carried no information.

# multilpa 0.7.0

* Every public verb returns a tidy `data.frame`. See `CHANGES.md` in the source
  repository for the full sweep.
