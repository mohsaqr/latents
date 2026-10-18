# Changelog

## latents 0.95.9

#### Documentation

- The package website links each author’s name to their own website
  (<https://saqr.me>, <https://sonsoles.me>), as the sibling packages’
  sites do.

## latents 0.95.8

#### Documentation

- The vignettes, the guides, the case studies and the README use the
  named plot functions:
  [`plot_profiles()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  (with `what = "bars"`, `"heatmap"`, `"raincloud"`),
  [`plot_responses()`](https://pak.dynasite.org/latents/reference/latents-plots.md),
  [`plot_prevalence()`](https://pak.dynasite.org/latents/reference/latents-plots.md),
  [`plot_sequences()`](https://pak.dynasite.org/latents/reference/latents-plots.md),
  [`plot_transitions()`](https://pak.dynasite.org/latents/reference/latents-plots.md),
  [`plot_diagnostics()`](https://pak.dynasite.org/latents/reference/latents-plots.md),
  [`plot_enumeration()`](https://pak.dynasite.org/latents/reference/latents-plots.md),
  [`plot_selection()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  and
  [`plot_tree()`](https://pak.dynasite.org/latents/reference/latents-plots.md),
  in place of `plot(x, what = )`. The function tables of the enumeration
  guides list them. The model families without named functions (growth
  mixtures, trajectory classes, mixture regression, additive and
  cross-level models) keep `plot(x, what = )`.
- The LPA and enumeration vignettes mention `cores`.
- Every help page is revised: titles and descriptions now say what each
  function does and returns, without design rationale or history, and
  every argument states its default. The revision corrected statements
  that did not match the code, among them the refusals listed by
  [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md),
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md),
  [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md),
  [`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md),
  [`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md),
  [`prior_control()`](https://pak.dynasite.org/latents/reference/prior_control.md)
  and
  [`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md);
  the columns of
  [`starting_values()`](https://pak.dynasite.org/latents/reference/starting_values.md)
  under full covariance; the tables offered by
  [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
  (starting values, ordinal and count tables, diagnostics); when
  `latents_no_group_classes` and `latents_unsupported_inference` are
  raised; and the default `vcov_type` of weighted fits.

## latents 0.95.7

#### Bug fixes

- Line plots of profiles and response probabilities
  ([`plot_profiles()`](https://pak.dynasite.org/latents/reference/latents-plots.md),
  [`plot_responses()`](https://pak.dynasite.org/latents/reference/latents-plots.md))
  no longer run long indicator names together: names that cannot sit
  side by side are turned. The profile labels at the right end get a
  margin sized to the longest label, so they are no longer cut off when
  there are many indicators, and they are spaced on the axis as drawn,
  so lines that end close together keep their labels apart.

## latents 0.95.6

#### New features

- The named plot functions now cover every view of a fit, grouped by
  subject so that typing `plot_` lists them all, and each group’s `what`
  takes `"all"` to draw every view it has:
  - `plot_profiles(x, what = )` adds `"bars"`, `"heatmap"`,
    `"raincloud"`, `"parallel"` and `"pairs"` to the lines (`what` is
    its last argument, so existing calls are unchanged);
  - new `plot_sequences(x)`;
  - new `plot_selection(x, what = )` for an enumeration: `"deltas"`,
    `"heatmap"`, `"smallest"`, `"trust"`, `"diagnostics"`, or `"all"`
    (default), the same as `plot_enumeration(x, what = "selection")`.

  Each draws exactly what `plot(x, what = )` draws.

#### Bug fixes

- `plot(x, what = "all")` and
  [`report()`](https://pak.dynasite.org/latents/reference/report.md)
  left out the heatmap of a model with only categorical indicators,
  which draws its response probabilities.

## latents 0.95.3

#### New features

- [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  and so
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md),
  [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md) and
  [`lca()`](https://pak.dynasite.org/latents/reference/lca.md), take
  `cores` to run the random starts of one fit in parallel. The starting
  values are drawn first, in order, so the fit is identical to
  `cores = 1`. Models with membership covariates still run their starts
  on one core.

#### Bug fixes

- [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
  and
  [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md)
  no longer refuse `categorical`, so a call written for
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  or
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  (`categorical = vars`) works as is. Every indicator in `vars` stays
  categorical whatever `categorical` names; a name that is not in `vars`
  is still refused.

- `plot(what = "heatmap")` on polytomous categorical indicators now
  draws one panel per profile, with indicators as rows and response
  categories as columns. The former layout, a block of columns per
  indicator, ran to indicators x categories columns, so with many
  indicators the numbers overlapped and the indicator names were cut
  off. Binary indicators keep their one-column-per-indicator layout.

#### Documentation

- The model-selection and class-enumeration guides, the enumeration and
  evaluation vignettes, the LPA workflow, the README and the help pages
  of [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
  [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md)
  and
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  describe the degenerate one-observation solution and
  `boundary = "exclude"`.

## latents 0.95.2

Fixes found by Equivalence 2.0, the strengthened external validation.

#### New features

- New argument `boundary = c("allow", "exclude")` in
  [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md)
  and the covariate models. `"exclude"` sets aside starts with a
  variance at `min_variance`: solutions in which a Gaussian profile
  shrinks around one observation, whose likelihood gain is set by
  `min_variance` rather than by the data. The more starts are run, the
  likelier one finds such a solution: in a simulation study (29
  datasets), a 4-profile dataset fitted with 500 starts gave a 5-profile
  solution with a one-observation profile that made BIC choose 5; with
  `boundary = "exclude"` BIC chose 4, as Mplus did. Opt-in: the default
  (`"allow"`) and every existing result are unchanged. When every start
  is at the bound, the best is reported with a
  `latents_all_starts_boundary` warning. Categorical indicators alone
  are unaffected.

#### Bug fixes

- [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) beyond
  homogeneous transitions (`transitions = "occasion"`, covariates,
  second order, mover-stayer) stopped its quasi-Newton finish early near
  a transition estimated at zero: since 0.9.13 the finish’s stopping
  rule was tied to `tol`, and at the default `tol = 1e-10` an
  occasion-specific model of LMest’s reference was left 2.8e-6 below its
  maximum with transitions 1e-4 off. The finish again runs to L-BFGS-B’s
  tight rule, as in 0.9.1, and agrees with LMest to 1.7e-6 or better.
  Fits move closer to their maximum; nothing else changes.
- [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md)
  could stop with `latents_no_converge` when a covariate nearly
  separates the classes: the likelihood with the classification error
  held fixed is not concave, and BFGS from zero settled on a separation
  ridge (Mplus UG example 7.3: -log L 81.7 on the ridge, 62.6 at the
  maximum). The search now also starts from the proportional-assignment
  estimate, which is concave and lies in the maximum’s basin, and keeps
  the better converged end; it reproduces Mplus’s R3STEP estimates.
- `enumerate_classes(cores = )` on macOS with R switched to Apple’s
  Accelerate (vecLib) BLAS: forked workers died with a segmentation
  fault once the session had used linear algebra, and the enumeration
  then stopped with an internal error at printing. Accelerate is not
  safe in a forked process, so there the candidates are now fitted by a
  socket cluster, as on Windows; results are identical to `cores = 1`,
  as before. Elsewhere the fork is kept.
- A worker killed by a signal (a crash, out of memory) returned an empty
  candidate that was not recognized as lost; it now raises
  `latents_parallel_failure`, as a worker that fails with an error does,
  and so does a socket worker that dies.
- The condition catalogue lists `latents_parallel_failure`, raised by
  `enumerate_classes(cores = )` when a worker process fails.

#### Equivalence

- The external validation now includes 160 simulated datasets with known
  parameters fitted by Mplus and Latent GOLD (likelihood, parameters,
  criteria, posteriors, entropy, AvePP, OCC, class counts, group level,
  [`predict()`](https://rdrr.io/r/stats/predict.html), enumeration
  grids,
  [`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md),
  `cores`), 42 genuine Mplus runs kept in JStats, `scaling` against
  Mplus, sampling weights in a mixture against Mplus WEIGHT with robust
  standard errors,
  [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
  against TECH14, and the group-class families at Latent GOLD’s finest
  quadrature. Comparisons that were scripts now run as tests, and the
  Latent GOLD kit is compared with the current package rather than with
  targets saved when it was built. See `equivalence/README.md`.
- Every documented claim of agreement with another program or a
  published result is registered in `equivalence/CLAIMS.csv` with the
  test that backs it, and a test fails when a claim or its test changes.
  Two claims gained tests (`bic_rows` equals flexmix’s BIC; a noise
  observation is profile 0, as in mclust’s classification); one was
  narrowed to what is tested (plain EM follows mclust’s path from the
  same start).
- A sabotage battery (`equivalence/sabotage/run.R`) injects realistic
  defects into the loaded package and confirms that the tests detect
  them.

#### Documentation

- The class-enumeration guide is rewritten in scientific register, with
  the methods cited.

## latents 0.95.0

#### Parallel enumeration

- [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md),
  and through it
  [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md)
  and
  [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md),
  gain `cores = 1L`, which fits the candidates on several processor
  cores: forked processes on macOS and Linux, a socket cluster on
  Windows. Every candidate sets and restores its own seed, so with
  `seed` given the result is identical to `cores = 1` for profile,
  transition (`time`) and group-class family grids; only the time
  differs (a 24-candidate grid of `course_engagement` took 11.7 s on one
  core, 4.2 s on four and 3.3 s on eight). Candidates’ warnings and
  messages are passed on in grid order, with their condition classes.
  `cores = 1`, the default, runs as before. `parallel`, part of R, is
  now in Imports.

## latents 0.9.41

#### `drop_empty`

- [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  and [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) gain
  `drop_empty = TRUE`, at the end of their signatures; enumerations,
  [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
  [`lca()`](https://pak.dynasite.org/latents/reference/lca.md) and
  covariate fits pass it on through `...`. It governs what the reported
  statistics do with rows whose every indicator is missing and groups
  whose every row is such a row. The estimates and log likelihood are
  the same either way.
  - `TRUE` (default): such rows and groups are left out of the sample
    sizes of the information criteria, of entropy, AvePP and the entropy
    charge of ICL, AWE and CLC, of the classification tables and plots,
    and now also of the effective class counts, as Mplus leaves them
    out. Class counts now agree with Mplus’s counts based on posterior
    probabilities.
  - `FALSE`: every number is what latents 0.9.37 reported, so results
    published with earlier versions can be reproduced (checked against
    0.9.37 for enumeration grids, covariate fits, both
    [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) engines
    and
    [`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)).
- Fits with different `drop_empty` are refused by
  [`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)
  and
  [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
  with `latents_incomparable_models`.

## latents 0.9.40

#### Groups with nothing observed

- A group whose every row has every indicator missing is now left out of
  the group sample size (`bic_groups` and the other group-convention
  criteria) and of group-level entropy, AvePP, the classification and
  entropy tables and the group diagnostic plots, as Mplus leaves it out
  of its clusters. Its posterior is the group-class prior, not a
  classification. A group with at least one observed row is unaffected;
  group class sizes and assignments keep every group. Fits record
  `n_informative_groups`. In a single-level fit each row is its own
  group, so `bic_groups` now equals `bic_individual` when rows are
  empty. Fits without such groups report exactly what they did.

## latents 0.9.39

Fixes found by an independent review of 0.9.38, and new references for
missing indicators.

#### Bug fixes

- [`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)
  on a general
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit kept
  rows with nothing observed in the entropy charge of ICL, AWE and CLC;
  it now uses the rows the fit classified, as every other fit does.
- [`predict()`](https://rdrr.io/r/stats/predict.html) failed on a fit
  with no continuous indicator when it was fitted with `scaling`, or
  with `centering = "grand"`; such a fit is neither centred nor scaled,
  and now predicts.
- `scaling = "person"` judged whether a unit’s standard deviation was
  zero against the largest value of every unit, so a real but small
  spread beside much larger values was refused, as was a unit with no
  observed value of an indicator. Each unit is now judged on its own
  values, and a unit with nothing observed is left unscaled.

#### Equivalence

- Five datasets with different missing-value patterns (missing
  completely at random, rows with nothing observed, missing at random on
  an observed indicator, one indicator mostly missing, and a two-level
  set with empty groups), fitted by genuine Mplus 9 and Latent GOLD 6.1.
  Log likelihood, sample size, BIC, entropy, AvePP and posteriors agree;
  the programs’ differing treatment of rows and groups with nothing
  observed is recorded.

#### Documentation

- The model-selection guide describes the choices made before fitting
  (random starts, centring and scaling, missing values), why entropy and
  AvePP can disagree, the two levels of classification, and refitting
  with more starts. The class-enumeration guide points to them.

## latents 0.9.38

Scaling of the indicators, group-level diagnostic plots, and
classification statistics that leave out rows with nothing observed.
Every addition is a new argument at the end of a signature or a new
value of `what`.

#### Classification quality with missing indicators (bug fix)

- Under `missing = "fiml"`, a row with every indicator missing was
  counted in relative entropy, AvePP, the classification and entropy
  tables, the entropy charge of ICL, AWE and CLC, and the case-level
  diagnostic plots. Such a row contributes nothing to the likelihood, so
  its posterior is the prior its group gives it, not a classification;
  counting it pulled every one of those statistics down (in one example
  relative entropy from 0.78 to 0.63). They are now taken over the rows
  with at least one observed indicator, the same rows the
  individual-level BIC already counted (`n_informative`). Estimates, log
  likelihoods and class sizes are unchanged; fits without such rows
  report exactly what they did. The rows used are kept on the fit as
  `informative`.
- The same rows are now left out of the individual-level BIC of
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) and of
  covariate fits, which counted every row;
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  already did. A covariate fit’s `bic_individual` therefore agrees with
  its own `get_results(x, "information_criteria")`, which it did not
  when rows were empty. The additive and cross-level families require
  complete data and are unaffected.
- The classification-quality and criteria plots of an enumeration
  labelled each line at its end in every panel, which ran off the edge
  beside a second panel and overlapped where values were close. With
  several panels they now share one legend, keyed by colour and shape; a
  single panel keeps its direct labels, spaced further apart.

#### Scaling the indicators (`scaling`)

- [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  gains `scaling = c("none", "person", "grand")`, at the end of its
  signature, beside `centering`. `"grand"` divides each continuous
  indicator by its overall standard deviation and `"person"` by each
  unit’s own, so `centering = "grand", scaling = "grand"` fits grand
  z-scores and `centering = "person", scaling = "person"` within-person
  z-scores. It reaches
  [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md), the
  covariate fits and every enumeration through `...`. The scales are
  kept on the fit: `get_results(x, "data")` returns the columns
  supplied, [`predict()`](https://rdrr.io/r/stats/predict.html) scales
  new data the same way, and the bootstrap and staged fits treat a
  scaled fit as they treat a centred one. With structures that give
  every indicator its own variance, grand scaling changes only the units
  (the log likelihood moves by `n * sum(log(sd))`, memberships are
  identical); with `EII` and `VII` it changes the solution. An undefined
  or zero standard deviation is refused with `latents_bad_data`.
- The print of a grand-centred fit no longer calls its profiles
  “profiles of change”; that describes person centring only.

#### Group-level classification plots

- [`plot_diagnostics()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  draws the group classes as well as the profiles. The new views
  `"group_sizes"`, `"group_posteriors"`, `"group_entropy"` and
  `"group_avepp"` are the four individual views drawn for the groups,
  and `what = "all"` includes them when a fit has two or more group
  classes. The individual views are unchanged. A fit with one group
  class refuses the group views with `latents_no_group_classes`.

#### Documentation

- A new guide, *Model selection: the evidence for the number of classes*
  (`vignettes/articles/model-selection.Rmd`), describes every piece of
  the enumeration evidence: the criteria and their definitions, the two
  sample-size conventions, deltas, the knee rule and its stability
  notes, classification quality, the estimation checks, comparisons and
  the bootstrap test.
- [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  documents the arguments its `...` passes to every candidate, leading
  with `n_starts`.
- The function tables of the enumeration vignette and guide name the
  minima table `"best"`.

## latents 0.9.37

Growth mixtures and latent transition models; nothing in the 0.9.21
interface changes, and every addition is a new argument at the end of a
signature, a new `type`, or a verb that used to refuse.

#### Enumeration grids

- The table of each criterion’s minimum is now `get_results(x, "best")`,
  the name the transition and family grids already use. `"criteria"`
  read as “every candidate’s criteria”, which is `"candidates"`.
  `"criteria"` still returns the same table, so existing scripts are
  unaffected.

#### Growth mixture models (`mixture_regression(random = )`)

- [`predict()`](https://rdrr.io/r/stats/predict.html) works on growth
  mixtures (it used to refuse). The persons in `newdata` are new
  persons: `type = "response"` and `"class_response"` give class
  trajectories, which need no outcome; given the outcomes, `"posterior"`
  classifies each person, `"class_conditional"` and `"conditional"` give
  each person’s own curve with their empirical Bayes random effects, and
  `"random_effects"` gives those effects by class. At lcmm’s own
  estimates they equal `lcmm::predictY()` (marginal and
  subject-specific), `predictClass()` and `predictRE()` (three
  configurations; `equivalence/test-growth-predict-lcmm.R`). A
  multilevel growth mixture classifies the clusters in `newdata` as new
  clusters.
- [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  has a growth-mixture method: one row per natural parameter
  (coefficients, membership logits, residual and random-effect standard
  deviations, variances and correlations, class shares) by Wald, or by a
  person bootstrap (`method = "bootstrap"`) whose replicates are matched
  to the fitted classes by their posterior memberships of the fitted
  persons. The Wald rows reproduce the fit’s own tables.
- [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md)
  works on trajectory models (latent class growth and growth mixtures;
  it used to refuse): person covariates predicting class membership,
  with the classification errors held fixed.
- [`sensitivity()`](https://pak.dynasite.org/latents/reference/sensitivity.md)
  refits a growth mixture under each seed and reports how many persons
  keep their class.
- Correlated residuals (new arguments `residual_correlation = "ar"` or
  `"brownian"` and `residual_time`): a continuous-time autoregressive
  process or Brownian motion, shared by the classes, on top of the
  random effects and the residual variance, as
  `lcmm::hlme(cor = AR(time))` and `cor = BM(time)`. At hlme’s estimates
  the likelihood equals hlme’s to 1e-8 (three configurations;
  `equivalence/test-growth-correlated-lcmm.R`). A process the data
  cannot separate from the residual variance or the random effects is
  left on its bound with a `latents_random_boundary` warning, and its
  rows have no standard error; the rest of the fit is inferred with
  those coordinates held.

#### Latent transition models (`lta()`)

- `r3step(fit, data, initial_covariates =, transition_covariates =)` on
  an [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit is
  the three-step latent Markov model (Di Mari, Oberski and Vermunt
  2016): covariates of the first profile and of the transitions, with
  the classification errors held fixed. It agrees with Latent GOLD 6.1’s
  `step3 modal ml` latent Markov model to its printed precision (four
  decimals) on three designs: log likelihoods, coefficients and standard
  errors (`equivalence/latentgold-step3-markov/`).
- Wald standard errors for the extended
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) model
  (covariate, occasion-varying or second-order transitions) cover every
  covariance structure (`model =`) and `covariance_model = "full"`; they
  were limited to EEI and VVI. The analytic scores agree with numerical
  gradients of the likelihood to about 2e-7.
- [`sensitivity()`](https://pak.dynasite.org/latents/reference/sensitivity.md)
  takes [`lta()`](https://pak.dynasite.org/latents/reference/lta.md)
  fits (with `data`) and reports how many occasions keep their profile
  under each seed.
- `plot(fit, what = "bars")`, `"profiles"` and `"responses"` of a
  transition fit draw the 95% intervals of
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  (new argument `intervals`, default `TRUE`); they used to draw point
  estimates only.

## latents 0.9.36

- [`plot_enumeration()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  of a grid over group classes at one number of profiles puts the group
  classes on the x-axis; it used to stack every candidate in a single
  column with overlapping labels. Grids over profiles plot as before.

## latents 0.9.35

- A posterior probability that rounds a hair past 1 no longer gives a
  negative `uncertainty` in `get_results(x, "assignments")`, and no
  longer drops that case, with a warning, from the
  classification-certainty and entropy plots of
  [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  and [`plot()`](https://rdrr.io/r/graphics/plot.default.html). The
  estimates are unchanged.
- Base R’s “empty cluster: try a better set of initial centers” warning
  from the k-means of a random start is no longer passed on; the fit
  already reports failed starts in its own warning and in
  `get_results(x, "starts")`.

## latents 0.9.34

- latents is on CRAN (0.9.21, published 2026-10-16). The README carries
  the CRAN badge and opens with a list of every vignette, guide and case
  study.

## latents 0.9.33

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of an
  occasion-varying
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit
  orders its panels by occasion number (with ten or more occasions,
  “occasion 10” used to come before “occasion 2”) and puts the occasion
  on a second line of the panel title, which narrow panels used to cut
  off.
- Top legends with more than four keys wrap into rows instead of running
  off the right edge (for example four group classes and the shared key
  of a group-class variances plot).

## latents 0.9.32

A function-by-function review of the package and its documentation. The
maths it re-derived (analytic scores against numerical gradients for
every model family, delta-method Jacobians, information criteria,
entropy, AvePP, OCC, the knee) held; the corrections below are to edge
cases, input validation, consistency between tables, plots and warnings,
and documentation.

### Changed results

- Group-class families (`family = "additive"`, `"dispersion"`,
  `"additive_dispersion"`) use SQUAREM-accelerated EM. Converged fits
  are unchanged to rounding (8 Latent GOLD cases to 7.5e-16 relative);
  fits that plain EM left crawling now converge within `max_iter`, so
  bootstrap replicates that used to be withheld are kept.
- The membership-covariate separation check, the one-step membership
  regression,
  [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md)
  (pooled and `by_group_class = TRUE`) and covariate-fit
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  now work in orthonormal design coordinates, so results no longer
  depend on where a covariate’s zero lies: a covariate with mean 100 and
  sd 1 is no longer mistaken for separation, and a covariate shifted by
  1e6 gives the same estimates and standard errors as the centred one.
  [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md)
  also converges more tightly (about 1e-4 closer to its maximum).
- With sampling weights, `small_classes` and the `latents_small_classes`
  warning use the weighted effective counts the fit reports; group-class
  families report weighted class counts and compute effective groups
  (and plot intervals) from the robust covariance they use for
  inference.
- In weighted general
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fits, the
  model-implied `"transitions"`, `"initial"` and
  `"second_order_transitions"` tables average groups with their sampling
  weights.
- [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) and its
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  refits default to `max_iter = 5000` (was `1000`). The tighter
  `tol = 1e-10` of 0.9.31 needs more updates: three- and four-profile
  models of two `course_engagement` indicators took 1,000 to 1,500 and
  stopped unconverged. Fits that converged before are unchanged.
- Transition and family grids (`enumerate_classes(time = )`,
  `enumerate_classes(family = )`) pass their candidates’ warnings on, as
  profile grids do; they used to keep them only in the `warnings`
  column.
- `get_results(grid, "best")` on transition and family grids considers
  only converged candidates; when none converged, `"best"` raises
  `latents_failed_candidate`.
- Response-probability whiskers in plots are the logit-scale intervals
  of
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md),
  so plot and table agree.

### Refusals instead of silent wrong answers

- [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) refuses
  (`latents_bad_data`) transition covariates when `occasions = "grid"`
  leaves a position skipped inside a group’s span; the covariate used to
  be taken as 0 there.
- [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) and
  covariate fits refuse (`latents_bad_time`) a `time` column of numbers
  stored as text whose alphabetical order differs from their numeric
  order (“10” before “2”). Zero-padded numbers and ISO dates are still
  accepted.
- `fit_staged(measurement = )` refuses (`latents_bad_stage`) a
  measurement fitted with weights, centring or ordinal/count indicators
  instead of holding it under a different likelihood.
- `mixture_regression(family = "ordinal")` refuses
  (`latents_not_identified`) a single-level mixture whose predictors
  take too few distinct values to identify it.
- An ordinal indicator with fewer than two observed categories is
  refused (`latents_bad_data`) even when its factor declares more
  levels; a declared but never-observed level warns
  (`latents_unused_category`).
- Categorical columns whose distinct values print alike are refused
  instead of getting duplicate labels.
  [`prior_control()`](https://pak.dynasite.org/latents/reference/prior_control.md)
  refuses a `scale` that is not a square, symmetric, positive-definite
  matrix.
- [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
  refuses `level = "groups"` for trajectory models instead of ignoring
  it. `vcov(fit, scale = "unconstrained", method = "bootstrap")` raises
  `latents_unsupported_inference` instead of returning `NULL`.
- [`predict()`](https://rdrr.io/r/stats/predict.html) on
  membership-covariate and latent transition fits raises
  `latents_unsupported_prediction`.

### Fixes

- [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md) and
  [`lca()`](https://pak.dynasite.org/latents/reference/lca.md) record
  the call as written, so
  [`update()`](https://rdrr.io/r/stats/update.html) refits from the
  caller’s data.
- `centering = "grand"` or `"person"` (with or without covariates) no
  longer refuses indicators recorded in very small units.
- [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md):
  factor predictors can be shared across classes through `common`.
  Regression and growth fits no longer store their fitting frame, so
  saved fits are several times smaller; enumerated growth candidates
  record a call that names the data.
- Growth mixtures: `plot(what = "random")` draws when a class’s
  random-effect covariance is degenerate and is identical on every draw;
  trajectories work for intercept-only class formulas;
  [`confint()`](https://rdrr.io/r/stats/confint.html) honours `parm`;
  failed starts warn (`latents_failed_starts`);
  `plot(what = "classification")` names each class on the bars instead
  of relying on colour alone.
- `get_results(fit, "recovery")` reports a missing known class as `NA`,
  and on a membership-covariate fit compares group-level truths with
  group classes.
- [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  and [`report()`](https://pak.dynasite.org/latents/reference/report.md)
  work on fits with a noise component.
  [`descriptives()`](https://pak.dynasite.org/latents/reference/descriptives.md)
  of a single-level fit no longer reports `n_groups`/`icc` for the
  internal row numbering.
- `print(summary())` names single-level and latent class fits correctly;
  `rows = 0` prints headings only; empty `"profiles"`/`"responses"`
  tables keep their standard-error columns.
- General [`lta()`](https://pak.dynasite.org/latents/reference/lta.md)
  fits gain a `"responses"` table for categorical indicators, flag an
  emptied stayer class as a small class, and no longer warn spuriously
  for a dotted covariate name. Cross-level fits store
  `profile_prevalences`, reported as a `prevalence` column of
  `get_results(fit, "profiles")`.
- Plots:
  [`plot_profiles()`](https://pak.dynasite.org/latents/reference/latents-plots.md),
  [`plot_responses()`](https://pak.dynasite.org/latents/reference/latents-plots.md),
  [`plot_prevalence()`](https://pak.dynasite.org/latents/reference/latents-plots.md),
  [`plot_diagnostics()`](https://pak.dynasite.org/latents/reference/latents-plots.md),
  [`plot_enumeration()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  and
  [`plot_tree()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  refuse objects they cannot draw with `latents_bad_argument` and name
  the right plot; a noise fit’s prevalence plot gives the noise share;
  profile views name the indicators they leave out;
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a general
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit uses
  the package’s probability scale; the group-class variances plot is
  keyed by class, panelled by level.
- Capped VEI covariance M-steps are reported once per fit with a count.
- Methods that ignored extra arguments
  ([`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
  on general transition fits and transition or family grids,
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) on transition fits) now
  warn (`latents_unused_argument`); `vcov(fit, vcov_type = "robust")` on
  a general transition fit used to return the observed covariance
  silently.
- Classed conditions replace unclassed errors for invalid seeds,
  controls, `parm`, `start`, `n_group_classes`, `lca(n_profiles = )` and
  covariate-model inputs;
  [`sensitivity()`](https://pak.dynasite.org/latents/reference/sensitivity.md)
  and
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  validate their controls and seed before fitting.
- Robust standard errors (CR1 since 0.9.31) are exactly
  `sqrt(G / (G - 1))` times those of Mplus’s MLR, which applies no
  small-sample factor: 0.8% larger with 60 groups.
  [`?parameter_inference`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  says so, and the external Mplus comparison now checks that
  relationship to 1e-6.
- [`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)
  no longer prints “\<- lowest BIC” and its plot no longer circles each
  criterion’s lowest value: `delta_bic` and `bic_weight` carry the
  evidence and the choice stays with the reader.
- Plots keep the Okabe-Ito grey for categories that are not classes (a
  noise component, terms shared by every class); classes cycle through
  the other colours, so the two never coincide. Noise tiles were
  `grey75`.
- Documentation: `?latents-conditions` lists every condition the code
  raises (checked by a test);
  [`?lta`](https://pak.dynasite.org/latents/reference/lta.md),
  [`?multilpa`](https://pak.dynasite.org/latents/reference/multilpa.md),
  [`?parameter_inference`](https://pak.dynasite.org/latents/reference/parameter_inference.md),
  [`?bootstrap_lrt`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md),
  [`?get_results`](https://pak.dynasite.org/latents/reference/get_results.md),
  [`?three_step`](https://pak.dynasite.org/latents/reference/three_step.md),
  [`?fit_staged`](https://pak.dynasite.org/latents/reference/fit_staged.md),
  the plot help and `course_engagement` are corrected where they
  disagreed with the code.

## latents 0.9.31

- Cluster-robust standard errors now carry the CR1 small-sample factor
  `G / (G - 1)` for `G` independent units, the default of
  [`sandwich::vcovCL()`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html).
  It applies to every `vcov_type = "robust"` result (profile, covariate,
  transition, family, mixture-regression and growth fits) and to
  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
  (both `vcov_type`s, every `method`) and
  [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md).
  Robust standard errors rise by `sqrt(G / (G - 1))`: about 2.6% with 20
  groups and 0.5% with 100. A one-class robust regression now reproduces
  [`sandwich::vcovCL()`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html)
  exactly.
- Wald inference is refused (`latents_boundary_fit`) for a
  membership-covariate fit whose logit coefficients are heading to
  infinity because a covariate separates the classes. Previously the
  refusal held only when the observed information happened to be
  numerically singular; otherwise a coefficient such as -504 was
  reported with a standard error of 669 and p = 0.45.
- The separation check (and the `latents_extreme_coefficients` warning)
  measures each coefficient per root-mean-square unit of its covariate,
  so it no longer depends on the units a covariate is recorded in: a
  real effect of a covariate in tiny units is no longer flagged, and a
  separating covariate in large units no longer escapes.
- [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) defaults
  to `tol = 1e-10` (was `1e-8`), as do the bootstrap refits of its
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md).
  The group-class shares and the mean group posteriors, one EM step
  apart, now agree to about `2e-5` instead of `2e-4`, for roughly 1.5
  times the iterations.

## latents 0.9.30

- Released function and method names and argument spellings remain
  available. Tables retain their existing columns (new columns are
  added; in an enumeration grid’s tables `candidate` now comes first).
  Corrections to statistical results and input validation are described
  below.
- Audit corrections: bootstrap inference verifies supplied data,
  ordering and sampling weights before resampling; membership-covariate
  models apply and retain requested centering; constrained covariance
  updates preserve their structure at an active variance floor; growth
  likelihoods retain residual information at extreme random-effect
  scales; R3STEP inference is stable under predictor-unit changes;
  pooling structural constants reports zero missing information;
  descriptives refuse matrix-valued columns.
- Weighted residual tables use weighted descriptive moments and withhold
  unsupported design-based test statistics and p-values. Weighted
  profile trees use the fitted weights for shares and flows. Noise-model
  entropy includes noise mass; named noise diagnostics retain the
  supported sizes view and refuse unsupported entropy, posterior and
  AvePP views.
- Comparisons support homogeneous/general transition and
  membership-covariate fits both directly and through enumeration;
  classification follows the retained comparison rows. Data identity
  includes raw indicators, category labels, centering, grouping,
  occasions and weights. Covariate enumeration tolerates unavailable
  start-replication counts. `get_results(x, "all")` returns the
  available tables for regression comparisons and family/transition
  enumerations. Smallest-profile curves distinguish covariance
  structures.
- Aggregated-binomial prediction consistently reports expected
  successes. New-data trial counts are inferred from the response
  columns or supplied with the named-only `trials` argument (a number,
  vector or column name). When neither is available, a classed data
  error requests trial counts; `trials = 1` explicitly requests
  per-trial probabilities.
- Every candidate of an enumeration grid has a name, and the one name
  reaches it everywhere. `classes3_groups4` is 3 classes and 4 group
  classes (`profiles` when an indicator is continuous; no `_groups` part
  on a single-level grid; `_VVI` etc. when the grid crosses structures;
  `<model>_groups2`, such as `additive_equal_groups2`, on a group-class
  family grid). The name is the first column, `candidate`, of the grid’s
  table and names each pick in `get_results(grid, "criteria")`;
  `grid$classes3_groups4` is that model, and an editor completes the
  names after `grid$`. Grids saved by earlier versions are named the
  same way.
- New
  [`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md),
  beside
  [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md),
  takes one model out of a grid:
  `get_candidate(grid, "classes3_groups4")`, `get_candidate(grid, 3, 4)`
  or `get_candidate(grid, c(3, 4))`.
  [`candidate_fit()`](https://pak.dynasite.org/latents/reference/get_candidate.md),
  its 0.9.21 name, still works.
- [`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)
  compares candidates of a grid: by name,
  `compare_models(grid, "classes3_groups4", "classes3_groups5")`; by
  class counts, `compare_models(grid, c(3, 4), c(3, 5))`; or by
  selection, `compare_models(grid, n_profiles = 3)`, or every candidate,
  `compare_models(grid)` (candidates that failed are left out, with a
  message). It gives one row per candidate with its information criteria
  (the candidates table’s numbers), entropies, AvePP summaries, smallest
  share, convergence, `delta_bic` and `bic_weight`, and takes a `cutoff`
  for the AvePP count. It also compares fits from
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
  [`lca()`](https://pak.dynasite.org/latents/reference/lca.md) or
  [`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md)
  directly, refusing fits to different data. The printout marks the row
  with the lowest BIC as “lowest BIC” (it read “best”, also for
  mixture-regression comparisons, whose columns are unchanged).
- A comparison carries the models it compared:
  `get_candidate(comparison, "classes3_groups4")` returns the model
  behind a row, and `get_results(comparison, "classification")` gives
  every class’s AvePP and sizes for compared profile, class and
  transition models.
- Enumeration grids are much smaller in memory and on disk. Each
  candidate used to record, as its `call`, the whole fitting function
  and a full copy of the data, so a grid held one copy of the data per
  candidate. A candidate’s call is now a short call that names the data,
  for example
  `multilpa(data = d, vars = ..., n_profiles = 3, n_group_classes = 8, ...)`.
  This applies to
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  (profile, transition and group-class family grids) and
  [`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md).
- Model-selection evidence for an enumeration grid, none of which
  selects a model, each reading the criterion in `using` (default
  `"bic_groups"`): `get_results(grid, "deltas")` gives, for every
  candidate, how much the criterion fell from one profile fewer and from
  one group class fewer; `get_results(grid, "knees")` gives a maximum
  perpendicular chord-distance knee, a simplified Kneedle-inspired rule
  without smoothing or sensitivity thresholds, along the profiles, along
  the group classes and over every candidate by its number of
  parameters, with the range it was found on and whether it holds
  without the last candidate. `plot(grid, what = "deltas")`,
  `"heatmap"`, `"smallest"` and `"trust"` draw the evidence;
  `what = "selection"` draws them all with the entropy and AvePP views.
- `plot(grid, what = "entropy")`, `"avepp"` and `"diagnostics"` (both)
  draw classification quality across the candidates of an enumeration
  grid, in the words of a single model’s plot: relative entropy and the
  lowest AvePP of each candidate against the number of classes, one
  panel per level, the AvePP `cutoff` (0.8) dashed.
  [`plot_enumeration()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  takes them too.
- An enumeration grid compares classification quality across candidates.
  Its candidates table gains, beside the two entropies, the smallest,
  largest and mean average posterior probability (AvePP) of each
  candidate’s profiles and group classes and the number of classes below
  a cutoff: `avepp_min`, `avepp_max`, `avepp_mean`, `avepp_below` and
  the same four prefixed `group_` for the group classes. The cutoff
  defaults to 0.8 and is set with
  `get_results(grid, "candidates", cutoff = )`. The columns are computed
  from the stored fits when the table is read, so grids fitted by
  earlier versions report them without refitting.
  `get_results(grid, "classification")` gives every class’s AvePP, modal
  and estimated sizes, odds of correct classification and relative
  entropy for every candidate in one table.
- Named plot functions for the main views, each with its own listed
  arguments so an editor completes them:
  [`plot_profiles()`](https://pak.dynasite.org/latents/reference/latents-plots.md),
  [`plot_responses()`](https://pak.dynasite.org/latents/reference/latents-plots.md),
  [`plot_prevalence()`](https://pak.dynasite.org/latents/reference/latents-plots.md),
  [`plot_diagnostics()`](https://pak.dynasite.org/latents/reference/latents-plots.md),
  [`plot_transitions()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  and
  [`plot_tree()`](https://pak.dynasite.org/latents/reference/latents-plots.md),
  alongside
  [`plot_enumeration()`](https://pak.dynasite.org/latents/reference/latents-plots.md).
  Each draws exactly what the matching `plot(x, what = )` draws, and an
  argument the fit’s family does not take is refused by name.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
  [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  result now draws through
  [`plot_diagnostics()`](https://pak.dynasite.org/latents/reference/latents-plots.md).
- Enumeration plots choose their criteria with `what`, as every other
  plot does: `plot(grid, what = "bic_groups")`,
  `plot(grid, what = "all")`, `plot(grid, what = "tree")`, and
  `plot_enumeration(grid, what = )`. The 0.9.21 spelling, `criterion =`
  (alone or with `what = "enumeration"`), draws the same plot;
  `criterion` with any other `what` is an error of class
  `latents_bad_argument` (with `what = "tree"`, a warning, as the tree
  uses no criterion). This applies to transition and group-class family
  grids too.
- Every plot now returns a ggplot object. The last base-graphics plots
  are ported: [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  on a pooled imputation fit, on a
  [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
  result and on a mixture regression (`"fitted"`, `"coefficients"`,
  `"posteriors"`). Style them with ggplot2. The bootstrap plot still
  takes `palette`; its base-graphics `style` constants have nothing to
  set on a ggplot, so passing them draws the plot with a warning of
  class `latents_unused_argument`. The posterior strip of a mixture
  regression spreads its points deterministically and no longer draws
  from the random number stream. These formerly base plots draw when
  called, including inside loops, and return the ggplot invisibly. Code
  that used their previous returned fit object should keep its original
  fit separately; the returned object is now the reusable plot.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a pooled
  fit chooses its parameters with `what` (`level`, which elsewhere means
  a confidence level, still works) and takes `subtitle`.
- No [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
  drops an argument it does not use without a word. Those that ignored
  one silently in 0.9.21 (comparison, growth mixture, mixture
  regression, pooled, additive, cross-level, general transition fits,
  bootstrap tests, transition and group-class family grids) now draw
  with a warning of class `latents_unused_argument`; the others refuse
  it with an error of class `latents_bad_argument`, as they did.
- `graphics` is no longer imported.
- The reference manual consolidates the former 74 help pages into 38
  topics. Released help aliases are retained, and plot/prediction method
  usages list their actual arguments and defaults.
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) are
  described once on
  [`?get_results`](https://pak.dynasite.org/latents/reference/get_results.md),
  which also lists every family’s tables (the eight
  `get_results.<class>` pages are folded in). All plotting, the fourteen
  `plot.<class>` pages and the `plot_*()` functions, is on one page,
  [`?"latents-plots"`](https://pak.dynasite.org/latents/reference/latents-plots.md).
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) share
  [`?"latents-predict"`](https://pak.dynasite.org/latents/reference/latents-predict.md).
- The README’s likelihood equation renders on GitHub: `\left\{` was
  stripped to `\left{` by GitHub’s markdown, so the braces are now
  `\lbrace`/`\rbrace`.

## latents 0.9.23

- New
  [`plot_enumeration()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  draws the information criteria of an enumeration grid with its
  arguments listed on its own, so an editor completes them; it is the
  view [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws
  for an enumeration.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`plot_enumeration()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  take `criterion = "all"` for every information criterion in the grid.
  Pairs that coincide in a single-level grid are drawn once.

## latents 0.9.22

- Numeric categorical models now support `plot(what = "profiles")`, with
  `statistic = "mean"`, `"median"`, or `"mode"` on the original category
  scale. These summaries use fitted probabilities, support FIML and
  membership covariates, and do not compute standard errors. Mixed-model
  profile plots continue to show continuous indicators.

## latents 0.9.21

CRAN release: 2026-10-16

- The version submitted to CRAN; the same code as 0.9.20.

## latents 0.9.20

- Tests fit each model from a single start unless the test is about
  starts, which halves the time the tests take on CRAN. No package code
  changed.

## latents 0.9.19

- Every help page has an example: the
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods for
  model comparisons, family and transition enumerations, pooled
  imputations, cross-level and general transition fits, and
  [`predict()`](https://rdrr.io/r/stats/predict.html) for growth
  mixtures, which shows the refusal and where the trajectories are.
- The transition-enumeration example in
  [`?get_results.latents_transition_enumeration`](https://pak.dynasite.org/latents/reference/get_results.md)
  fits the first 50 students, so it runs in a few seconds on Windows as
  well (it took 12 s there).
- The remaining tests that compare latents with another implementation
  move to `tests/equivalence/`: the depmixS4 and LMest fits, exact path
  sums, dense likelihoods, enumerated EM steps and independent Hessians,
  together with their reference helpers and fixtures. The installed
  tests no longer contain any such comparison.

## latents 0.9.18

- The growth mixture guide moves from the vignettes to the package
  website (same address), so the package ships three vignettes: `lpa`,
  `lca` and `lta`. Its model fits made it most of the vignette build
  time.

## latents 0.9.17

- The logo reads “latents” (it still showed the earlier name,
  “multilatent”), and the website icons are regenerated from it (they
  showed “multilpa”).
- The reference manual is shorter: the
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) and
  [`nobs()`](https://rdrr.io/r/stats/nobs.html) methods that had a page
  each are documented together on four pages
  ([`?"latents-print"`](https://pak.dynasite.org/latents/reference/get_results.md),
  [`?"latents-summary"`](https://pak.dynasite.org/latents/reference/get_results.md),
  [`?"latents-as-data-frame"`](https://pak.dynasite.org/latents/reference/get_results.md),
  [`?"latents-model-methods"`](https://pak.dynasite.org/latents/reference/latents-model-methods.md)).
  The methods are unchanged.
- Fewer tests run on CRAN: tests that verify the estimators against
  their definitions, tests that only check a plot draws, and the
  regression tests from past audits are skipped there. All of them still
  run locally and in CI.

## latents 0.9.16

- Every help page now has an example that runs in the automatic checks.
  The growth mixture pages
  ([`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`simulate()`](https://rdrr.io/r/stats/simulate.html)) and
  [`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)
  open with a fit on the first 60 students of `growth_scores`, which
  runs in about a second; the full-data examples stay in `\donttest{}`.
  The `growth_scores`, `growth_schools` and `student_esm` pages open
  with a summary of the data. The
  [`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)
  example no longer fits a random intercept that `course_engagement`
  cannot support.
- The package ships four vignettes: `lpa`, `lca`, `lta` and
  `growth-mixture`. The workflow guides (`workflow-lpa`, `workflow-lca`)
  and the guides to evaluation, covariates, mixture regression, latent
  class growth analysis and group-class models move to articles on the
  package website, with unchanged addresses. This keeps the installed
  package small and the checks quick.
- Tests that compare latents with other implementations (`numDeriv`,
  `mclust`, `MASS`, [`glm()`](https://rdrr.io/r/stats/glm.html), Latent
  GOLD and Mplus results, reference implementations, the numbers in the
  workflow guides) move to `tests/equivalence/`, which is not part of
  the package and runs in its own CI job. `numDeriv` is no longer a
  suggested package. Slow tests that refit many models are skipped on
  CRAN; every model family keeps tests that run there.
- The README’s enumeration example names its covariance structure
  (`model = "VVI"`), as
  [`candidate_fit()`](https://pak.dynasite.org/latents/reference/get_candidate.md)
  needs once several structures are crossed, and its text matches the
  regenerated output.

## latents 0.9.15

- Vignettes and worked examples save and restore the user’s display
  options (`digits` and `width`) after use, including the scripts
  installed in `inst/doc`.

- Corrects regression and growth inference: robust covariance checks the
  number and rank of independent units, multilevel growth BIC uses
  clusters, enumerated models retain the appropriate inference, and
  boundary parameters no longer erase unrelated delta-method standard
  errors.

- Stabilizes Poisson and negative-binomial regression densities, NB
  dispersion derivatives near the Poisson limit, and growth residual
  calculations for outcomes with large locations. Single-class binary
  regressions are accepted; nonfinite outcomes and evaluated predictors
  are refused explicitly.

- Trajectory tables include offsets and report response-scale standard
  errors. Alternate time columns work in tables and plots; random-effect
  spread holds other predictors fixed, and observed means respect
  sampling weights.

- Withholds regression bootstrap likelihood-ratio p-values when refits
  fail validation. Model comparisons also check trials, weights and
  cluster layout.

- Handles unassigned regression classes, weighted transition class
  counts and matrix-valued prediction inputs correctly. Noise cases
  retain their mass in profile percentages and appear as class zero in
  sequence tables and plots.

## latents 0.9.14

- Declares `MASS` and `withr`, used by the tests, in Suggests (the
  0.9.13 check warned about them as undeclared).

## latents 0.9.13

### One engine

- Every estimator now runs on one shared kernel (see
  `validation/ENGINE_DESIGN.md`): one EM driver (with SQUAREM for the
  profile models), one start runner, one quasi-Newton finish, one
  measurement block (Gaussian, categorical, ordinal and count
  indicators) shared by the profile, covariate, cross-level and both
  transition engines, one forward-backward pass for both transition
  engines, shared latent structures for the regression and growth
  models, and shared inference, alignment and simulation services.
  Results are unchanged: 251 recorded fits of every engine and option
  reproduce bit for bit (`equivalence/engine-golden/`), except where a
  defect below was fixed.

### New models

- Multilevel growth mixture models: `mixture_regression(cluster = )`
  nests the persons (`id`) in clusters (schools, clinics) with
  `n_group_classes` group classes of clusters that shift how probable
  each trajectory class is, and `group_membership` covariates of the
  group class. The clusters are the independent units of the likelihood,
  standard errors and BIC. New tables `"group_classes"` and
  `"clusters"`;
  [`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md)
  compares the number of group classes;
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) draws the two
  levels. New dataset `growth_schools`; the growth vignette has a
  section on it.
- Negative-binomial regression mixtures:
  `mixture_regression(family = "negative_binomial")`, NB2 with a
  dispersion per class or shared (`variance = "equal"`). One class
  reproduces
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html). A
  dispersion at the Poisson limit warns (`latents_boundary`) and the
  other standard errors are conditional on it.
- Ordinal regression mixtures: `mixture_regression(family = "ordinal")`,
  a proportional-odds (cumulative logit) regression within each class
  with class-specific ordered thresholds in place of the intercept. The
  outcome is an ordered factor, a factor or whole-number categories.
  Thresholds are rows `"threshold:<lower>|<upper>"` of the coefficient
  table; class means are expected category scores;
  `predict(type = "probabilities")` gives the category probabilities and
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) draws categories
  in the outcome’s own type. One class reproduces
  [`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html). Mixtures
  that one outcome per class assignment cannot identify (two categories,
  or no predictors) are refused with `latents_not_identified`.

### Documentation

- New vignette “Latent class growth analysis: classes of trajectories”
  (`vignette("trajectory-classes")`): trajectory classes without random
  effects, every plot view, curved trajectories, pass/fail and
  ordered-band trajectories, distal outcomes, and when to prefer the
  growth mixture model. The growth mixture vignette now focuses on
  random effects and the multilevel model, with its plots.
- Coefficient forest plots of trajectory fits label every estimate with
  the same number of decimals per panel (two, or three for small
  values).
- The recovery table of a regression mixture prints in the report
  layout.

### Bug fixes

- Growth mixture models gain a start from the persons’ own least-squares
  trajectories. Every previous start came from the model without random
  effects, which can split persons by their level; the growth fit then
  stopped at a lower maximum (by 9 to 11 log-likelihood units on test
  data where classes differ in slope).
- A fit with a noise component:
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  `get_results(x, "all")` no longer fail; `"classification_errors"` and
  `"bch_weights"` (which were wrong for noise units) refuse with
  `latents_unsupported_noise`.
  [`sensitivity()`](https://pak.dynasite.org/latents/reference/sensitivity.md)
  no longer reports `NA` agreement for noise fits.
- [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  of a general
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit now
  uses the caller’s `step` in its table, not only in the attached
  covariance.
- [`get_tna()`](https://pak.dynasite.org/latents/reference/get_tna.md)
  and
  [`get_group_tna()`](https://pak.dynasite.org/latents/reference/get_group_tna.md)
  on a general
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit
  refuse with `latents_unsupported_tna` instead of a bare “no applicable
  method”.
- An additive-family information matrix that is negative definite is
  refused with `latents_singular_information` instead of an unclassed
  Cholesky error.
- A group covariate that varies within a group, and all starts failing
  in the covariate and homogeneous transition engines, now raise classed
  conditions (`latents_bad_data`, `latents_all_starts_failed`).
- Two-level regression summaries label a class logit’s intercept within
  a group class “Intercept, group class 1” rather than “Intercept x
  group_class_1”.

## latents 0.9.12

### Growth mixture models

- [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
  gains `random`, `random_covariance` and `random_diagonal`: random
  effects within classes, so persons scatter around their class
  trajectory (the Gaussian growth mixture model). The likelihood is
  exact; estimation is EM with the random effects as missing data,
  finished by quasi-Newton with analytic scores. Validated against an
  independent dense likelihood, numerical derivatives and `lcmm::hlme()`
  (six configurations; the same likelihood at `hlme`’s estimates and at
  least its maximum).
- Variable names may replace one-sided formulas: `random = "wave"`,
  `random = "intercept"`, `common = "x"`, `membership = "age"`.
- A degenerate random-effect covariance (a variance at zero, or
  perfectly correlated random effects) warns (`latents_random_boundary`)
  and reports no interval for it.
- New
  [`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md):
  fits to the same data side by side, with BIC differences and Schwarz
  weights; [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  draws the comparison.
- Trajectory models (`id` with `class_level = "group"`, with or without
  `random`) share new ggplot2 views: class trajectories with confidence
  bands, the spread of persons around them and observed class means;
  individual trajectories; a coefficient forest; predicted random
  effects; and a classification structure plot.
- Result tables of these fits print in a report layout (readable labels,
  intervals as `[low, high]`, APA p-values, fit statistics as a card)
  while remaining plain numeric data frames.
- Growth fits work with
  [`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md),
  including its bootstrap likelihood-ratio test, have a
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) method, and
  support distal outcomes in
  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
  (person-level, bias-corrected);
  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
  results print in the report layout.
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md)
  refuse growth fits with classed conditions.
- New dataset `growth_scores` and vignette “Growth mixture models”.

## latents 0.9.11

- Transition fits no longer report false non-convergence. The
  quasi-Newton finish counted only L-BFGS-B’s code 0 as converged, but
  its line search can stop at the maximum (code 52,
  “ABNORMAL_TERMINATION_IN_LNSRCH”) when rounding hides any further
  ascent, which is platform-dependent. It now restarts once from where
  it stopped and accepts the fit when the log likelihood cannot be
  raised by more than `tol`. In 120 negative-binomial test datasets 5
  were flagged unconverged at points with scores of about 1e-6 that a
  restart improved by less than 1e-12; none are now. This, not the
  iteration cap, is what failed the 0.9.8-0.9.10 CI checks on macOS and
  R 4.1 (the 0.9.10 test change is reverted).

## latents 0.9.10

- Test-only fix: a negative-binomial transition test required
  convergence within `max_iter = 5`, which held on some platforms only
  (macOS CI and R 4.1 needed more iterations, and the package correctly
  warned). The test now allows 500 iterations; its checks of dispersion
  floors, monotone likelihood and stationary scores are unchanged.

## latents 0.9.9

### Transition models: full-covariance simulation and bootstrap inference

- Transition fits with full covariance structures (EEE, VVV and the
  other non-diagonal structures, homogeneous or extended) can now be
  simulated, so
  [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
  compares them. Each profile’s rows are drawn through the Cholesky
  factor of its covariance; diagonal fits keep their random stream.
- `parameter_inference(method = "bootstrap", data = )` works on
  transition fits. Persons are resampled with replacement and refitted;
  each replicate’s profiles are matched to the original’s on the
  measurement and its group classes on what the matched profiles make of
  them, with initial logits rebased when the reference profile moves.
  Standard errors are the replicates’ standard deviation and intervals
  their percentiles. This gives inference for extended models whose Wald
  inference is unavailable (full covariance structures) and supports
  sampling weights.
- Fixed a 0.9.8 regression: transition fits refused `boundary = "fix"`,
  which `get_results(fit, "responses")` and
  [`summary()`](https://rdrr.io/r/base/summary.html) request, so
  categorical transition response tables lost their standard errors.
  Transition inference holds nothing at a bound, so `"fix"` now agrees
  with `"error"`; a fit with an active bound is still refused with
  `latents_boundary_fit`.

## latents 0.9.8

### Model review fixes

- Ordinal/count indicators now work through
  [`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
  and
  [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md),
  as documented. Bootstrap and imputation label alignment move their
  measurement parameters together, restore the ordinal reference
  profile, and match profiles using the extra indicators.
- Negative-binomial bootstrap comparisons preserve the fitted count
  model and refuse incompatible indicator specifications. Supplied
  fitting data are checked against all measured indicators and
  model-specific sufficient statistics before inference or comparison.
- Poisson log densities and negative-binomial dispersion derivatives use
  numerically stable calculations. Invalid extra-indicator column shapes
  and nonfinite ordinal values are refused. Single-level covariate
  inference accepts the original data without demanding an internal
  identifier.
- Second-order occasion transitions omit coefficient coordinates that
  never apply. Stayer tables report identity transitions; second-order
  inference includes the second-order coefficient block. Count
  boundaries retain their constraints during transition optimization and
  refuse Wald inference.
- Extended transition [`coef()`](https://rdrr.io/r/stats/coef.html) now
  carries complete measurement covariance coordinates for all fourteen
  structures, with likelihood-preserving decoding. The existing
  restrictions on extended-model Wald inference remain explicit.
- Transition fits with only categorical, ordinal or count indicators
  return all result tables successfully, including an empty
  continuous-profile table. Transition data validation checks original
  category labels and occasion values.
- One-profile NB fits with varying dispersion retain the correct profile
  label and report their dispersion standard error.
- Group-class bootstrap comparisons verify group sizes and scatter and
  refuse equivalent one-class family specifications. Weighted restricted
  cross-level fitting avoids collisions with indicator names used for
  internal weights.
- Sampling weights normalize safely at extreme common scales. Regression
  mixture class counts, shares and composition respect sampling weights;
  invalid integer controls and seeds are refused before conversion.
- Confidence levels and numerical differentiation controls receive
  explicit validation on the reviewed model surfaces.

### Follow-up review fixes

- Negative-binomial fits near the Poisson limit no longer stop short of
  the maximum. The dispersion objective is convex there in log
  dispersion, and the old Newton step could stall while still reporting
  convergence (fitted dispersion 4.2e-7 against a profile maximum of
  1.4e-5, with singular Wald inference). The M-step now uses modified
  Newton with reflected eigenvalues.
- [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
  refuses a full-covariance transition null up front with
  `latents_unsupported_inference`, rather than failing every replicate
  and returning `p_value = NA` with a generic warning.
- [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  on transition fits refuses `method = "bootstrap"` and
  `boundary = "fix"` with `latents_unsupported_inference` instead of
  silently returning the Wald table. On general
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fits,
  `adjust` now produces a `p_adjusted` column.
- [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
  accepts the caller’s own data frame for single-level fits.
- Under sampling weights, the `"classification"` table’s `estimated_n`
  and `estimated_proportion` (and so the odds of correct classification
  and
  [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)’
  smallest class) use weighted posterior totals. The same applies to
  restricted cross-level `composition` and `group_classes` counts.

## latents 0.9.7

### Ordinal and count indicators everywhere; negative-binomial counts

- Membership covariates (`profile_covariates`, `group_covariates`) and
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) (every
  extension, invariant or occasion-specific measurement) now take
  `ordinal` and `count` indicators, with Wald inference, weights and
  simulation
  ([`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
  for transition fits).
- `count_model = "negative_binomial"` models counts more variable than a
  Poisson within a profile (NB2: variance `mu + alpha mu^2`), with
  `count_dispersion = "varying"` (per profile) or `"equal"` (shared).
  `get_results(fit, "count_means")` adds the dispersion and its standard
  error. A dispersion at zero is the Poisson limit: the fit warns
  (`latents_boundary`) and Wald inference is refused there.
- Checked against Latent GOLD 6.1 (`poisson overdispersed`),
  class-varying and shared dispersion: log likelihoods agree to 4e-5 and
  parameters to its printed precision; seven reference cases in all,
  stored as a test fixture.
- The Newton line searches of the membership logits, ordinal and count
  M-steps accept a step that changes the objective by rounding only,
  instead of halving it up to forty times at the maximum (30x faster on
  a negative- binomial fit).

## latents 0.9.6

### Bug fix

- [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  with `profile_covariates` or `group_covariates` silently fitted any
  covariance structure other than EEI, VVI, EEE and VVV as one of those
  four (for example `model = "VEI"` gave VVI) and did not record the
  structure. It now refuses them with `latents_unsupported_structure`;
  fit the covariate-free model with the structure and use
  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
  or [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md).
- [`?"latents-conditions"`](https://pak.dynasite.org/latents/reference/latents-conditions.md)
  documents the condition classes added in 0.9.3–0.9.6.

## latents 0.9.5

### Ordinal and count indicators

- [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
  [`lca()`](https://pak.dynasite.org/latents/reference/lca.md) and
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
  take `ordinal` and `count`, naming indicators in `vars`:
  - `ordinal`: an adjacent-category logit with category intercepts
    shared by every profile and one location per profile (Latent GOLD’s
    default ordinal model): `(K - 1) + (C - 1)` parameters per indicator
    instead of `C (K - 1)` as `categorical`. Ordered factors or whole
    numbers.
  - `count`: Poisson, one mean per profile.
- New tables: `get_results(fit, "ordinal")` (category probabilities and
  locations, with standard errors) and
  `get_results(fit, "count_means")`.
- Mixed freely with continuous and categorical indicators, in single-
  and two-level fits, with `missing = "fiml"` and `weights`. Wald
  (analytic scores) and bootstrap inference,
  [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md),
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`sensitivity()`](https://pak.dynasite.org/latents/reference/sensitivity.md)
  and
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  handle them.
- Checked against Latent GOLD 6.1 on five datasets (mixed, ordinal only,
  count only, three classes, two-level): log likelihoods agree to 4e-5
  and every parameter to 5e-5, the precision Latent GOLD prints
  (`equivalence/latentgold-ordinal/`); the results are a test fixture.
- Not yet available with membership covariates, `start`, `fixed`,
  `prior`, `noise`,
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) or the
  group-class families (`latents_unsupported_indicator`).

## latents 0.9.4

- The membership M-step of covariate fits (`profile_covariates`,
  `group_covariates`) now solves its weighted multinomial logits by
  Newton-Raphson with the exact Hessian instead of BFGS. BFGS stalled
  above its score tolerance on about half of the steps of an ordinary
  fit, which left fits flagged unconverged (and refused standard errors)
  on some platforms but not others. Estimates are unchanged; fits
  converge in fewer EM iterations.
- Continuous integration passes again: the R 4.1 floor job no longer
  tries to install `mice` (its current dependency chain needs R \>=
  4.4), one lint is fixed, and the checkout action is updated.

## latents 0.9.3

### Sampling weights

- [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
  [`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md),
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) and
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
  take `weights`, the name of a column holding one sampling weight per
  independent unit (each `id` group, or each row of a single-level fit).
  The fit maximizes the pseudo log likelihood `sum_j w_j log L_j`, with
  the weights scaled to sum to the number of units (Mplus’s convention),
  so AIC and BIC stay on the sample’s scale.
- Covered: every covariance structure, categorical indicators,
  `missing = "fiml"`, membership covariates, every `family` (additive,
  dispersion, additive-dispersion, restricted and full cross-level),
  every extension of
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md), and all
  three nestings of
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md).
- Integer weights reproduce the fit to the data with each unit repeated
  that many times, to 1e-9 or better in the log likelihood, for every
  engine (tested). A one-profile weighted fit reproduces Mplus 9’s
  weighted regression (`TYPE = COMPLEX`, `WEIGHT`) and closed-form
  weighted least squares (`equivalence/weights/`).
- Standard errors are the sandwich:
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html) and
  [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
  default to `vcov_type = "robust"` for a weighted fit and refuse
  `"observed"` and `"opg"` (`latents_unsupported_weights`).
  `method = "bootstrap"` resamples units with their weights.
- A weighted fit prints its weight column and Kish’s effective sample
  size; `get_results(fit, "model")` gains a `weights` column. Effective
  counts and proportions are weighted; classification and posterior
  tables report each unit’s own posterior.
- Refused with `latents_unsupported_weights`:
  [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
  and `enumerate_regressions(bootstrap = )` (the parametric bootstrap
  ignores the design),
  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md),
  [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md),
  `prior` and `noise`. Weights that vary within a unit, are negative or
  missing raise `latents_bad_weights`.

## latents 0.9.2

### Choosing and testing transition models

- `enumerate_classes(..., time = )` enumerates latent transition models
  over `n_profiles` and `n_group_classes` (and covariance structures
  named in `model`), with every other
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) argument
  held fixed; read with
  [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  `candidate_fit(grid, n_profiles = , n_group_classes = )`.
- [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
  compares nested transition fits (`data` required): more profiles or
  group classes, occasion-varying against homogeneous transitions, added
  covariates, occasion-specific measurement, second order, a stayer
  class. The null model is simulated with the data’s groups, occasions,
  covariates and missing values kept.
- `lta(mover_stayer = TRUE)` adds a class of stayers who never change
  profile (Goodman’s mover-stayer model), with its own initial
  distribution.
- The extended transition model now takes `missing = "fiml"` and
  covariance structures (`model =`, `covariance_model = "full"`);
  standard errors with diagonal covariances (including FIML).
- A transition or initial probability estimated at zero (a move never
  observed, common with occasion-varying transitions on sparse late
  occasions) is reported as a boundary fit (`latents_boundary`, standard
  errors withheld) instead of as non-convergence.

## latents 0.9.1

### Latent transition analysis: beyond homogeneous transitions

- [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) gains
  arguments that relax each default assumption of the model:
  - `transition_covariates`: covariates (fixed or changing over
    occasions) shift the transition probabilities through a multinomial
    logit per origin profile, the log odds of moving to each other
    profile rather than staying;
  - `initial_covariates`: covariates shift the starting profile;
  - `transitions = "occasion"`: a separate transition matrix for each
    move;
  - `measurement = "occasion"`: profile means and variances (or response
    probabilities) per occasion;
  - `order = 2`: second-order transitions, the next profile depending on
    the two previous ones;
  - `model =`: the fourteen covariance structures (`"VEI"`, `"EEV"`, …)
    for the homogeneous model, as in
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).
- The extensions return a `multilpa_lta` fit with tidy tables
  (`get_results(fit, "transition_coefficients")`, `"transitions"` per
  occasion, `"second_order_transitions"`, `"initial_coefficients"`,
  `"profiles"`, …), Wald standard errors from analytic scores (observed,
  robust or OPG), [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html),
  [`BIC()`](https://rdrr.io/r/stats/AIC.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html). EM is
  finished by a quasi-Newton search on the exact likelihood.
- External agreement, 30 of 30 quantities
  (`equivalence/lta-extensions/`): Mplus User’s Guide examples 8.13 and
  8.14 and an Mplus model with occasion-specific thresholds
  (log-likelihoods within Mplus’s printed precision, estimates within
  2e-3), depmixS4 covariate transitions (log-likelihood within 1e-7) and
  LMest time-heterogeneous transitions (within 4e-8). Second-order
  transitions, which no external program fits, are checked against an
  exact sum over every profile path.

## latents 0.9.0

### Experimental: the Houle et al. (2026) multilevel families

- `multilpa(family = )` fits the multilevel latent profile families of
  Houle, Morin & Harvey (2026) beside the existing profile model
  (`family = "profiles"`, their dispersion-heterogeneity model):
  - group-class families with no individual profiles, each group
    carrying a Gaussian intercept per indicator: `"additive"` (classes
    differ in means; within-group variance shared), `"dispersion"`
    (classes differ in within-group variances) and
    `"additive_dispersion"` (both); `between_variance = "varying"` or
    `"equal"`. The likelihood is exact (no quadrature); zero
    between-group variances are reached by a boundary maximization with
    a Karush-Kuhn-Tucker check; interior fits are finished by Newton
    steps.
  - cross-level families with individual profiles and group classes from
    the same indicators, following the published manifest-aggregation
    specification: `"restricted_cross_level"` and `"full_cross_level"`.
    Their likelihood is a working likelihood (group means reuse the
    ratings), so it compares only cross-level fits of the same data; no
    standard errors.
- Tidy tables for every family through
  [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
  (`"parameters"`, `"group_classes"`, `"groups"`, `"intercepts"`,
  `"recovery"`, … for the group-class families; `"profiles"`,
  `"composition"`, … for the cross-level ones),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html) and
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  (observed, cluster-robust or OPG covariance from analytic group
  scores).
- `enumerate_classes(family = c("additive", "dispersion", "additive_dispersion"), n_group_classes = )`
  crosses families, between-variance restrictions and class counts;
  `candidate_fit(grid, n_group_classes = , model = )` picks one.
  [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
  tests nested group-class fits by parametric bootstrap.
- `latents_weak_class` warning with `effective_groups` in the
  `"group_classes"` table: a class supported by fewer than 50 effective
  groups (small, poorly separated, or both) has intervals that can be
  miscalibrated. Threshold set by a predeclared rule and evaluated on
  held-out simulations (flags 100% of a rare, weakly separated condition
  and 0.8% of well-separated ones). Also flags designs with fewer than
  about 50 groups per class.
- Evidence (`validation/ADDITIVE_SIMULATION.md`,
  `equivalence/latentgold-families/`): predeclared simulations with 1000
  datasets per condition (additive: 31 of 34 reference-condition
  parameters meet every gate; dispersion and additive-dispersion: 45 of
  46; every miss is the ML small-sample downward bias of a between-group
  variance, coverage 0.925-0.947). Latent GOLD 6.1 agrees on parameter
  counts in all eight external cases, and its likelihoods and posteriors
  converge to latents’ as its quadrature is refined (largest remaining
  gaps 1.6e-4 and 2.0e-4). Robust standard errors are needed for
  heavy-tailed ratings.
- New vignette: `vignette("additive")`.
- Not yet available for these families: missing data, covariates, and an
  external Mplus comparison (Mplus is not available here).

## latents 0.8.8

### plot() returns ggplot objects (breaking)

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a fit
  (`multilpa`, covariate and transition fits), on an enumeration and on
  a
  [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  result now returns ggplot objects instead of drawing with base
  graphics. Print one to draw it, save it with
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html),
  or restyle it with `+ ggplot2::theme()`. ggplot2 stays in Suggests: a
  plot without it raises `latents_missing_package`, and
  [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  and [`report()`](https://pak.dynasite.org/latents/reference/report.md)
  print a message and carry on.
- `what = "all"`,
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
  diagnostics result and an enumeration plotted with `combine = FALSE`
  return a `latents_plots` list, named by view, that draws every plot
  when printed.
- The styling arguments `palette`, `symbols`, `linetypes`, `style` and
  the style constants passed through `...` are gone; an unknown argument
  now raises `latents_bad_argument` instead of being ignored.
- Every view shares one profile order (largest first), one profile share
  (the posterior share every table reports) and Okabe-Ito colours paired
  with shapes. Bars start at zero, the heatmaps have colour keys, and
  the case diagnostics are strips of every case with each profile’s mean
  marked.
- New views: `"parallel"` (every case as a line, one panel per profile),
  `"pairs"` (scatter-plot matrix with each profile’s 95% covariance
  ellipse) and, for enumerations, `what = "tree"` (an icicle of how
  profiles split as more are added).
  [`plot_views()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  lists them.
- Model comparison switches from direct labels to a legend beyond five
  series; failed candidates are marked with a cross on the panel floor.
- `plot(diagnostics(fit))` on a one-profile fit returns the sizes and
  average posterior views instead of refusing.

## latents 0.8.7

### Messages

- The single-level notice from
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
  and
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  with `id = NULL` is now one line (“Fitted a single-level latent
  profile model.”). The error for a missing `id` is shorter too.

### ggplot2 views (in development)

- Internal ggplot2 builders for every clustering view: profiles, bars
  (from zero), heatmap (mixture-standardized, with a colour key),
  raincloud, sizes, entropy and posteriors (per-case strips), average
  posterior probability and model comparison. New views: parallel
  coordinates, a scatter-plot matrix with 95% covariance ellipses and
  uncertainty-sized points, and an icicle of how profiles split across
  numbers of profiles (Zappia and Oshlack, 2018), with bands carrying
  posterior mass and colours following each profile’s lineage.
- Every view uses one profile share (the posterior share), one display
  order (largest first) and Okabe-Ito colours paired with shapes.
- ggplot2 is in Suggests.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) still draws
  the base-graphics views; the switch comes in a later release.

## latents 0.8.6

### Revision review fixes

- Mixture-regression prediction retains the training membership factor
  levels, contrasts and transformed terms. Posterior predictions rebuild
  group IDs for new rows; fitted plots retain grouping and membership
  columns and accept transformed regression predictors.
- A group-membership formula with one slope and no intercept is
  optimized as a slope, rather than mistaken for an intercept-only
  model.
- `r3step(by_group_class = TRUE)` now handles the one-group-class limit
  and agrees with the pooled regression for both observed and robust
  covariance.
- [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md)
  preserves an omitted `model`, allowing the documented covariance
  switches and all-categorical inputs through the wrapper.
- [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
  refuses prior-fitted models whose likelihoods are evaluated at
  posterior modes. Its fixed-mask missingness assumptions are clarified,
  and multiple-imputation guidance distinguishes algebraic equivalence
  from evidence about coverage and imputation-model validity.

## latents 0.8.5

### Single-level verbs, covariance-structure grids and workflow vignettes

- New [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md) and
  [`lca()`](https://pak.dynasite.org/latents/reference/lca.md) fit
  single-level latent profile and latent class models. They are
  `multilpa(id = NULL)` and `multilca(id = NULL)` under their ordinary
  names, return the same object, and refuse `id` and `n_group_classes`.
- [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  and [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md) take
  `model`, the three-letter covariance code (`"EEE"`, `"VVI"`, …), as an
  alternative to `variance_model` / `covariance_model` or `volume` /
  `shape` / `orientation`; giving both is an error.
- New
  [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md)
  and
  [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md)
  compare single-level latent profile and latent class models, as
  [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md) and
  [`lca()`](https://pak.dynasite.org/latents/reference/lca.md) fit one:
  `enumerate_lpa(data, vars, n_profiles = 1:6)` crosses the counts with
  the covariance structures, and
  `enumerate_lca(data, vars, n_classes = 1:6)` treats every indicator as
  categorical.
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  remains the verb for nested data and needs `id`.
- The two-level verbs
  ([`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md),
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md))
  with `id = NULL` say which analysis was estimated – a single-level
  latent profile or latent class analysis – and name the single-level
  verbs; called without `id` they are refused with the same pointer.
- `enumerate_classes(model = )` now defaults to `"basic"`: the four
  structures that combine equal or varying variances with covariances
  absent or present (`EEI`, `VVI`, `EEE`, `VVV`). `"all"` fits the 14
  structures; codes can still be named. With one continuous indicator
  the structures that coincide are fitted once; with only categorical
  indicators no structure is crossed; a structure set through the other
  arguments replaces the default, and `model = NULL` restores the
  previous grid. Because a grid now holds several structures per class
  count,
  [`candidate_fit()`](https://pak.dynasite.org/latents/reference/get_candidate.md)
  needs `model` to pick one of them.
- New data set `srl`: five self-regulated learning scales for 300
  respondents simulated by a large language model.
- New vignettes `vignette("workflow-lpa")` and
  `vignette("workflow-lca")`: a latent profile and a latent class
  analysis from the data to a reported model.
- `get_results(fit, "profiles")` and `get_results(fit, "responses")`
  carry `mean_standard_error`, `variance_standard_error` and
  `probability_standard_error` without the data being passed; they are
  `NA`, with a message, when the information matrix cannot be formed.
- [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  draws its classification plots by default (`plots = TRUE`): profile
  sizes, posterior probabilities, entropy contributions and the average
  posterior probability matrix.
- `plot(what = "raincloud")`: for each profile and indicator, the
  density, a quartile box and the observations.
- Wald intervals for probabilities are formed on the logit scale and for
  variances on the log scale, so they stay inside the parameter’s range;
  [`confint()`](https://rdrr.io/r/stats/confint.html) returns exactly
  the intervals
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  reports.
- A probability at or below 1e-6 counts as on its bound for inference:
  EM approaches zero slowly and stopped at 1e-8, leaving a singular
  information matrix that `boundary = "fix"` could not hold.
- [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md) and
  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
  name their adjusted p-value `p_adjusted` and their classes
  `profile_1`, `group_class_1`, as
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  does (were `p_value_adjusted` and `class_1`).
- `id = NULL` now raises its `latents_single_level` notice as a message,
  not a warning: a single-level model is a legitimate choice.
- [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  prints only the individual level for a fit with one group class, and a
  single-level fit’s assignments carry no `group_class` column.
- `plot(what = "profiles")` draws 95% intervals when the fit has
  standard errors (`intervals = TRUE`, profiles dodged apart);
  `plot(what = "responses")` draws intervals too.
- `plot(what = "bars")` and the profile plot no longer fail on a
  bound-active or unconverged fit: they draw without whiskers and the
  subtitle says why.
- `plot(what = "heatmap")` on an all-categorical fit draws the response
  probabilities: one row per class, one column per category, values
  printed where they fit.
- Covariate fits draw every measurement and classification view
  (`"bars"`, `"heatmap"`, `"raincloud"`, `"sizes"`, `"avepp"` are new
  for them).
- The enumeration plot uses display names (AIC, BIC, ICL), draws one BIC
  panel when the group-level and individual-level values coincide (a
  single-level grid), labels series only by what distinguishes them, and
  spaces labels by the text height. Gridlines and zero lines stay inside
  the panel. Okabe-Ito yellow now comes seventh, after the six colours
  that read well on the light panel.
- [`summary()`](https://rdrr.io/r/base/summary.html) of an enumeration
  prints a compact candidates table and names
  [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
  for the rest.
- An all-categorical fit prints as a latent class analysis rather than
  with a residual covariance, and a single-level fit’s header reads
  “BIC”.
- A covariate fit with one group class names its membership intercept
  `(Intercept)`, not `group_class_1`.
- Help-page examples use the bundled data sets and verbs, with no `$` or
  bracket indexing; a test keeps it that way.

### Bootstrap likelihood-ratio tests with missing data and covariates

- [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
  now accepts fits made with `missing = "fiml"`: every simulated
  replicate is given the observed missing cells before it is refitted,
  so the reference distribution loses the same information as the
  observed statistic. The pattern is treated as fixed (independent of
  the profiles).
- [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
  now accepts membership-covariate fits: the covariates are held at
  their observed values and memberships are drawn from the fitted
  logits. A one-group-class null may be compared with an alternative
  that adds group covariates or slopes by group class; otherwise both
  models must use the same membership regressions
  (`latents_bad_nesting`). A covariate replicate whose likelihood has
  converged counts as valid even if a membership logit is still drifting
  (an over-fitted class), since the statistic reads only the maximized
  likelihood; the replicate table’s new `logits_settled` column records
  it.
- Models with different missing-data handling are now refused as
  incomparable (`latents_incomparable_models`).

### Standard errors for latent transition models

- [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) and
  [`confint()`](https://rdrr.io/r/stats/confint.html) now work on
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fits
  instead of refusing them: observed-information, robust (clustered on
  sequences) and OPG standard errors for the measurement model, the
  group-class shares, the initial profile probabilities and the
  transition probabilities. The scores follow the Fisher identity
  (expected counts from the forward-backward pass minus the counts the
  fitted probabilities imply) and match numerical differentiation of the
  likelihood; the standard errors match a Richardson-extrapolated
  Hessian to 1e-6 and depmixS4’s numerical ones to within its own
  finite-difference error. Over 300 simulated panels, 95% intervals
  covered 0.947-0.960. The `latents_no_inference` condition is retired.
- [`coef()`](https://rdrr.io/r/stats/coef.html) on an
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit now
  uses the names and order of
  [`vcov()`](https://rdrr.io/r/stats/vcov.html): under
  `variance_model = "equal"` it reports the shared variances once, under
  `covariance_model = "full"` the covariance matrices, and the
  group-class shares come before the initial and transition
  probabilities.

### Multiple imputation for missing covariates

- New
  [`pool_imputations()`](https://pak.dynasite.org/latents/reference/pool_imputations.md):
  fits
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  to every completed data set (a list of data frames, or a `mids` object
  from mice), aligns each fit’s profile and group-class labels to the
  first, and pools with Rubin’s rules. The table adds `df`, `within`,
  `between`, `riv` and `fmi` to the columns of
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md);
  `get_results(x, "imputations")` and `get_results(x, "fits")` show each
  imputation and the relabelling it needed, and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) shows the
  imputations’ estimates beside the pooled interval. This is the route
  for missing covariates, which `missing = "fiml"` cannot integrate out.
  Aligning a covariate fit reparameterizes its membership logits, whose
  reference category moves with the labels. Pooled values agree with
  [`mice::pool.scalar()`](https://amices.org/mice/reference/pool.scalar.html)
  to machine precision. A failing imputation raises
  `latents_pooling_failed` rather than being dropped.

### First-stage uncertainty in staged fits

- `parameter_inference(method = "bootstrap")` now accepts a
  [`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
  result instead of refusing it. Every resample of the groups refits
  both stages, so the standard errors and percentile intervals carry the
  measurement’s sampling variability into the group-class estimates; the
  Wald default still conditions on the measurement. The table reports
  every parameter, the measurement included. On a simulated two-class
  design the within-class profile probabilities’ errors were 14% and 32%
  wider than the conditional ones, and a bootstrap that held the first
  stage reproduced the conditional Wald errors to within 3%.

### Profile-covariate slopes by group class

- `multilpa(profile_covariates = , profile_slopes = "group_class")` lets
  each group class have its own profile-covariate slopes, so a covariate
  can predict profile membership differently in different kinds of group
  (a cross-level interaction). The coefficients are reported with terms
  such as `z:group_class_1`, and inference, robust and OPG covariances
  cover them. The default, `"shared"`, is the previous model. With one
  group class the two are the same model. Naming it without
  `profile_covariates` raises `latents_bad_argument`.

### Missing indicators in membership-covariate models

- `multilpa(profile_covariates = , group_covariates = , missing = "fiml")`
  now fits: missing indicators are integrated out of the measurement
  density exactly as in the covariate-free model, for diagonal and full
  residual covariance and for categorical and mixed indicators.
  Previously this combination was refused with `latents_bad_argument`.
  Standard errors (observed, robust and OPG) cover it through the Fisher
  identity: a missing value enters the scores through its conditional
  mean and covariance given the row’s observed values. Covariates
  themselves must still be complete; a missing covariate now raises the
  classed `latents_bad_data`.

### Faster estimation

- `multilpa(acceleration = "squarem")`, the new default, accelerates EM
  with SQUAREM (Varadhan & Roland 2008). A cycle extrapolates along the
  last two EM steps and finishes with an ordinary EM step, falling back
  to plain EM whenever the extrapolation would lower the likelihood, so
  the path stays monotone; convergence is judged over the whole cycle,
  so it is never looser than plain EM’s. `acceleration = "none"`
  reproduces plain EM’s path (use it to follow mclust or Mplus step by
  step). Prior fits, membership-covariate models always use plain EM.
  Results can differ slightly from earlier versions: the accelerated fit
  usually ends closer to the maximum, and on a flat, weakly identified
  likelihood it may stop at a different point of the ridge.
- The first start of
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  is now Ward’s hierarchical clustering, with each cluster’s own means
  and spread; the rest remain random k-means. Over 14 structures on 3
  datasets it never lowered the best-of-10 likelihood and raised it in 3
  of 36 cases, for example on iris VVV from -186.57 to mclust’s -180.19.
  Fits may therefore land on a better maximum, or on the same maximum
  with profiles numbered differently.
- EVE and VVE take a fixed number of warm-started orientation steps per
  EM iteration (a generalized EM step) instead of solving the
  orientation to convergence inside every M-step, which took a median of
  220 inner steps per iteration. The same maxima are reached 3.5 to 8
  times faster, and SQUAREM now applies to them.
- `enumerate_classes(id = NULL)` enumerates single-level models, the
  search
  [`mclust::mclustBIC()`](https://mclust-org.github.io/mclust/reference/mclustBIC.html)
  performs, with one single-level notice for the grid.
- New [`predict()`](https://rdrr.io/r/stats/predict.html) method for
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  and
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
  fits: the modal profile (`"class"`), every profile’s probability
  (`"posterior"`), or each row’s log density under the fitted mixture
  (`"density"`) for new rows, prepared as the fit prepared its own
  (categorical levels, centering, missing data, noise). Groups in the
  new data are classified as new groups. On the training rows it
  reproduces the fit’s posteriors, and it matches
  [`mclust::estep()`](https://mclust-org.github.io/mclust/reference/estep.html)
  at the same parameters.
- `descriptives(by = )` lists strata in a fixed order (a factor’s
  levels, otherwise sorted) rather than in order of first appearance.
- The E-step no longer computes a missingness pattern per row for
  complete data and takes row maxima without a per-row
  [`apply()`](https://rdrr.io/r/base/apply.html): about 14 times faster
  per iteration, with identical results.
- On 5,000 rows, 4 indicators and 3 profiles, 13 of the 14 covariance
  structures now fit faster than
  [`mclust::Mclust()`](https://mclust-org.github.io/mclust/reference/Mclust.html)
  from one start (VVV is the exception), and 13 of 14 reach a higher
  likelihood.

### Mixture regression

- New
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
  fits finite mixtures of regressions (clusterwise, latent class or
  regression-mixture models) for `"gaussian"`, `"binomial"` (0/1, factor
  or `cbind(successes, failures)`) and `"poisson"` (with
  [`offset()`](https://rdrr.io/r/stats/offset.html)) outcomes. Classes
  can belong to each row, to a whole group (`class_level = "group"`,
  flexmix’s `y ~ x | id`), or to each row with a second-level group
  class that shifts the class shares within groups (`n_group_classes`,
  Vermunt 2003). Membership covariates at both levels (`membership`,
  `group_membership`), coefficients shared across classes (`common`),
  and equal residual variances (`variance = "equal"`).
- Analytic scores (Fisher’s identity) give observed-information,
  sandwich (clustered on the independent unit, or on `id` for a
  single-level fit) and outer-product standard errors; validated against
  a numerical Hessian of the likelihood in every nesting. A Monte Carlo
  study (`validation/mixture-regression-recovery.R`, 200 replications)
  found 91.5–97% coverage of nominal 95% intervals for the two-level
  model.
- [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
  tables: `coefficients` (with Benjamini-Hochberg-adjusted slopes and
  odds/rate ratios), `classes`, `membership`, `group_classes`, `fit`
  (AIC, BIC, SABIC, ICL, entropy), `assignments`, `groups`, `fitted`,
  `starts`, `classification` and `recovery` (against a known
  classification). [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`simulate()`](https://rdrr.io/r/stats/simulate.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) and
  [`nobs()`](https://rdrr.io/r/stats/nobs.html) methods.
- New
  [`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md)
  compares class counts in one table, with an optional parametric
  bootstrap likelihood-ratio test.
- A binary outcome with one trial per class assignment is refused
  (`latents_not_identified`), since that mixture is not identified.
- Single-level and group-level fits reproduce flexmix’s likelihood at
  its estimates and its binomial and Poisson estimates
  (`equivalence/test-mixture-regression-flexmix.R`).
- New dataset `study_hours` and `vignette("mixture-regression")`.

### Fixes

- `multilca(id = NULL)` now fits one group class by default, as
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  does, instead of refusing the default `n_group_classes = 2`.
- `equivalence/run.R` evaluated tests in the retired `multilpa`
  package’s namespace when it was still installed; it now uses
  `latents`.

### mclust parity

- [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) and
  [`confint()`](https://rdrr.io/r/stats/confint.html) give Wald standard
  errors for all fourteen covariance structures, not only EEI, VVI, EEE
  and VVV. The ten that constrain the volume, the shape or the
  orientation across profiles are differentiated in their own free
  coordinates — log volumes, log shapes with the determinant-one
  constraint built in, and Cayley-transform orientations (log-Cholesky
  factors for VEE and EVV) — with exactly as many dimensions as the
  structure has parameters, and carried to the reported variances and
  covariances by the delta method. Robust (`"robust"`, `"opg"`) errors
  work for them too. Validated against an independent numerical Hessian
  of the observed-data likelihood in a different chart.
- `multilpa(prior = prior_control())` gives maximum a posteriori
  estimates under the conjugate prior of Fraley and Raftery (2007),
  mclust’s `priorControl()`: the same default hyperparameters and the
  same M-steps as
  [`mclust::me()`](https://mclust-org.github.io/mclust/reference/me.html),
  for the ten structures mclust defines a prior for. `log_likelihood`,
  and the criteria built from it, are the unpenalized likelihood at the
  posterior mode, as in mclust. Wald inference is refused for such a
  fit; the bootstrap refits with the prior. New condition class
  `latents_unsupported_prior`.
- `multilpa(noise = TRUE)` adds mclust’s uniform noise component over
  the data’s hypervolume (`hypvol()`), for one group class and
  continuous, complete indicators, with any covariance structure and
  with or without `prior`. The noise is profile `0` in
  `get_results(fit, "assignments")` (with a `posterior_noise` column),
  `"posteriors"`, `"profile_probabilities"` and the classification
  diagnostics. The component counts two parameters, as in mclust, so BIC
  equals mclust’s. New condition class `latents_unsupported_noise`,
  raised also by the verbs that do not yet account for a noise
  component.

### New controls

- [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  and [`vcov()`](https://rdrr.io/r/stats/vcov.html) give standard errors
  for a one-step covariate model with categorical or mixed indicators
  (`multilca(profile_covariates = )`,
  `multilpa(categorical = , profile_covariates = )`). Each category’s
  probability is reported with a delta-method standard error and no Wald
  test, as for the covariate-free model.
- `parameter_inference(vcov_type = "opg")` inverts the outer product of
  the group scores (the BHHH estimate). It is the estimator `glca`
  reports, and reproduces its standard errors; `"observed"` remains the
  default.
- `parameter_inference(boundary = "fix")` holds categorical response
  probabilities that sit on `min_probability` at that bound and reports
  the other parameters conditionally on them. The default, `"error"`,
  refuses as before, now with a message that names the alternative.
- [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
  and [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) gain
  `select_start`. `"converged"` reports the best converged start
  whenever one converged, instead of an unconverged start that edges it
  out by a negligible likelihood.
- [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md)
  gains `by_group_class`. `TRUE` gives every group class its own profile
  intercepts, with shared slopes and one latent group class per group,
  the two-level form of the one-step model. The pooled regression
  attenuates the slopes when group classes differ in their profile mix.

### Changed defaults

- [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md) now
  defaults to `vcov_type = "robust"`, standard errors clustered on
  groups, matching
  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md).
  Pass `vcov_type = "observed"` for the previous behaviour.

### Bug fixes

- A covariate fit now stores `min_probability`, so inference recognises
  response probabilities on their bound and refuses with
  `latents_boundary_fit` instead of a numerically singular information.

## latents 0.8.4

- The website URL in DESCRIPTION ends in a slash, and the README gives
  the CRAN installation command.
- The package check is about a third faster: vignettes use fewer random
  starts, the slowest examples use smaller data, and six exhaustive test
  files are skipped on CRAN (they run on GitHub Actions).

## latents 0.8.3

- The README opens with the package description, and every vignette and
  article names its authors.

## latents 0.8.2

- The package description is rewritten.

## latents 0.8.1

- A documentation website, built with pkgdown, is published at
  <https://pak.dynasite.org/latents/>.

## latents 0.8.0

### Package renamed

- The package is renamed from multilpa to **latents**, a name that
  covers latent profile, latent class and latent transition models and
  the planned extensions. The model functions keep their names
  ([`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md),
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md)), as do
  the result classes. Condition classes now use the `latents_` prefix
  (for example `latents_bad_argument`, previously
  `multilpa_bad_argument`), `multilpa_plot_types()` is
  [`plot_views()`](https://pak.dynasite.org/latents/reference/latents-plots.md),
  and the conditions catalogue is
  [`?"latents-conditions"`](https://pak.dynasite.org/latents/reference/latents-conditions.md).
  The vignettes are
  [`vignette("lpa", package = "latents")`](https://pak.dynasite.org/latents/articles/lpa.md),
  `"evaluation"`, `"covariates"`, `"lca"` and `"lta"`.

### New features

- [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
  fits a two-level latent class model:
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  with every indicator categorical, so the items are named once. Mixed
  models remain `multilpa(categorical = )`.
