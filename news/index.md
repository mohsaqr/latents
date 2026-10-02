# Changelog

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
- New vignette:
  [`vignette("additive")`](https://pak.dynasite.org/latents/articles/additive.md).
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
  [`plot_views()`](https://pak.dynasite.org/latents/reference/plot_views.md)
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
  [`candidate_fit()`](https://pak.dynasite.org/latents/reference/candidate_fit.md)
  needs `model` to pick one of them.
- New data set `srl`: five self-regulated learning scales for 300
  respondents simulated by a large language model.
- New vignettes
  [`vignette("workflow-lpa")`](https://pak.dynasite.org/latents/articles/workflow-lpa.md)
  and
  [`vignette("workflow-lca")`](https://pak.dynasite.org/latents/articles/workflow-lca.md):
  a latent profile and a latent class analysis from the data to a
  reported model.
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
- New dataset `study_hours` and
  [`vignette("mixture-regression")`](https://pak.dynasite.org/latents/articles/mixture-regression.md).

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
  [`plot_views()`](https://pak.dynasite.org/latents/reference/plot_views.md),
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
