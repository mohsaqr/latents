# latents 1.4.2

Changes since the CRAN release 1.1.9.

## Refit a model made with 1.1.9 if

* it is a two-step model with covariates and categorical indicators
  (`fixed = "measurement"`). Its standard errors now hold the response
  probabilities, as documented, and can be substantially smaller;
* you compared a staged fit with a jointly estimated one in
  `compare_models()`. Criteria now count the staged fit's measurement
  parameters, which can reverse the ranking;
* you used the local-independence check of `diagnostics()`. It used
  unadjusted p-values and flagged a pair in up to half of correctly
  specified models; it now applies the Benjamini-Hochberg adjustment it
  states;
* you filtered on `p_adjusted` from `get_results(fit, "residuals")`. It was
  the raw p-value and is now Benjamini-Hochberg, as in `get_residuals()`;
  `adjust = "none"` gives the raw values;
* it uses grand scaling. The standard deviation is now the population one,
  weighted by any sampling weights, so likelihoods and scales change
  slightly;
* it is an `lta()` fit with homogeneous transitions and
  `drop_empty = FALSE`. Its BIC now counts every row with a positive
  weight, as the general transition engine's does.

## Warnings

* A cross-level fit with group classes warns that its likelihood counts the
  ratings twice, so its criteria favour group classes and the full form even
  on data with no group structure.
* Enumerations name the candidates that reached their best log likelihood
  from one start only; refit those with more `n_starts` before comparing
  their criteria.
* A transition row with fewer than one expected outgoing move warns and is
  marked unsupported; its Wald tests are refused.
* `descriptives()` warns that it describes the recorded rows without
  weights.

## New features

* One function per plot or table, taking any object the package returns:
  `plot_sizes()`, `plot_entropy()`, `plot_posteriors()`, `plot_avepp()`,
  `plot_heatmap()`, `plot_bars()`, `plot_raincloud()`, `plot_parallel()`,
  `plot_pairs()`, `plot_trajectories()`, `plot_individuals()`,
  `plot_classification()`, `plot_coefficients()`, `plot_criteria()`,
  `plot_deltas()`, `plot_trust()` and `plot_all()`; `get_checks()`,
  `get_entropy()`, `get_classification()`, `get_avepp()`, `get_residuals()`
  and `get_replication()`. `overview()` prints and draws everything about an
  object.
* Enumeration grids have a verb per table (`get_candidates()`,
  `get_best()`, `get_deltas()`, `get_knees()`, `get_tree()`) and print a
  one-screen summary that warns about edge minima, single-start maxima and
  low average posteriors.
* `align_labels(fit, reference)` numbers one fit's profiles and group
  classes like another's, for comparing samples or waves.
* `profile_covariates` and `group_covariates` work with every covariance
  structure and take `start` and `fixed`; holding the measurement model
  gives the two-step estimator.
* `three_step()` takes a categorical outcome, tests that the classes share
  one mean or distribution (`contrast = "overall"`), propagates the
  measurement model's uncertainty (`step_one = "estimated"`) and can leave
  units with a missing outcome out (`missing = "exclude"`).
* `parameter_inference(boundary = "fix")`, `vcov()` and `confint()` work for
  `lta()` fits without covariates, occasion-specific measurement, second
  order or stayers, and `get_results()` then reports the standard errors of
  the means and variances.
* `plot()` of a bootstrap from `parameter_inference(method = "bootstrap")`
  draws each estimate with its interval; `as.data.frame(boot, what =
  "replicates")` returns the replicates.
* General transition fits have the `"posteriors"`, `"entropy"` and
  `"information_criteria"` tables, and `"covariances"` with
  `covariance_model = "full"`.
* More views of `lta()` fits: `plot_sequences()` and `plot_profiles()` for
  ordinal, count, weighted and mover-stayer fits; the transition network
  (`get_tna()`, `get_group_tna()`) when transitions are homogeneous and
  first order; `plot_raincloud()`, `plot_parallel()` and `plot_pairs()` with
  homogeneous transitions, each answered occasion a case.
* Fits whose indicators are all ordinal draw `plot_profiles()` (each
  profile's expected category), `plot_responses()` and `plot_heatmap()`.
* `plot_sequences()` sorts the units within each class by their most
  frequent profile.
* `descriptives()` with `id` reports the within- and between-person
  standard deviations.
* `plot_individuals()` takes `persons`: how many persons, or which.
* `plot_pairs()` marks each profile's density with its marker, so the
  curves are told apart by shape as well as colour.
* Every function has its own help page with a verb-led title; each page
  lists every argument's values and default, the fits it accepts and
  refuses, and the columns it returns.

## Bug fixes

* `compare_models(test = "lrt")` accepted pairs that are not nested; they
  are now refused. Comparing shared with group-class slopes failed.
* `lca(count_dispersion = )` fitted the count indicator as Gaussian.
* `confint()` failed on a covariate fit with `fixed = "measurement"`.
* `get_results(fit)` and `as.data.frame(fit)` returned an empty profile
  table for a fit without continuous indicators, `lca()` fits included;
  they now return its response table.
* `three_step()`: BCH contrasts involving a probability outside [0, 1] were
  tested; an outcome category seen only on excluded rows removed every
  categorical ML standard error.
* `mixture_regression()` under separation stopped with an error from
  `eigen()`; it now returns the fit with `NA` standard errors and a warning.
* A bootstrap table whose rows were reordered attached replicates to the
  wrong parameters; such tables are now refused.
* Person centring of an indicator observed once in a group fitted
  artificial zero deviations; it is now refused.
* Unused factor levels became membership or random-effect parameters.
* A ragged transition bootstrap crashed; incompatible resamples now count
  as failures.
* Weighted fits: some plot subtitles and the posterior-strip means were
  unweighted, and cross-level prevalence and sample counts included
  zero-weight groups.
* `diagnostics()` reported "2 of 2 starts" for restricted cross-level fits
  and false thin-profile warnings for constrained covariate fits.
* `predict()` ignored a misnamed argument; `get_results(pooled, "all")` was
  refused; `vcov(pooled, type = )` ignored `type`.
* Count-model arguments that do not apply were dropped in silence; they now
  warn.
* An invalid option value raises an error that names the argument and its
  accepted values.
* Assignment tables lost an identifier column named `class`; small
  coefficients lost precision beside large intercepts; regression
  intercepts were shown as odds or rate ratios.
* The growth individuals view dropped classes beyond the fourth.
* `plot_sequences()` failed on an occasion no group answered.
* `plot_parallel()` bridged a missing middle indicator and did not draw a
  case with one answered value between gaps.
* `plot_all()` of an ordinal `lta()` fit failed.
* Plot legends left out the key of a status no candidate had.
* A platform that does not report its processor count failed in `fork()`;
  it now runs on one core with a warning.

# latents 1.1.9

Changes since the CRAN release 0.9.21.

## Refit a model made with 0.9.21 if

* its Gaussian indicators are recorded in small units (proportions, metres,
  rates). The variance floor and the EM stopping rule now do not depend on
  the units, and fits in small units can reach a much higher likelihood;
* its occasions are labelled `T1` ... `T10`. They were ordered
  alphabetically; such labels are now refused, so give the order with
  `factor(levels = )`;
* it uses sampling weights with zeros. Units of weight zero now count for
  nothing;
* it has several group classes. The classes may now be numbered the other
  way round, with the same estimates;
* it is an EII, EEI or EEE model, which now defaults to `search = "emEM"`
  and can reach a higher maximum (`search = "kmeans"` gives the old fit);
* you compared models on ICL, AWE or CLC, which now use the entropy of the
  complete classification.

## New features

* `compare_models(..., test = "lrt")` tests each nested model against the
  one before it.
* `get_results(grid, "lmr")`: Vuong-Lo-Mendell-Rubin and adjusted
  Lo-Mendell-Rubin tests.
* `local_dependence = "item1:item2"` models a pair of indicators jointly
  within each class.
* `direct_effects = "y ~ z"` lets a covariate shift an indicator's mean in
  every profile.
* `three_step(method = "ml")` for distal outcomes, and `r3step()` on `lta()`
  fits.
* `predict()` for growth mixtures.
* `diagnostics()` reports estimation, class sizes, classification, local
  independence and information, one row per quantity.
* `print()` says how many starts reached the best log likelihood.
* `cores` runs the starts of one fit, an enumeration or `sensitivity()` in
  parallel, with the same result as one core.
* New arguments: `scaling` (grand or within-person z-scores),
  `boundary = "exclude"` (set aside starts with a variance at its bound) and
  `drop_empty` (rows with no observed indicator).
* New and extended plots: `plot_sequences()`, `plot_selection()`, more views
  in `plot_profiles()`, group-level views in `plot_diagnostics()`, and one
  plot per tree in `plot_tree()` (`combine = TRUE` for one panel).

## Bug fixes

* `coef()` of some `lta()` fits returned centred means.
* `sensitivity(fit, data = )` dropped every seed for single-level fits.
* Rows with no observed indicator were assigned a class and counted in the
  classification error and three-step estimates.
* `mixture_regression(common = "intercept")` had no effect; it is now
  refused.

# latents 0.9.21

* First release on CRAN.
