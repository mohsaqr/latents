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
