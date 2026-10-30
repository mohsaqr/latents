# Extract model estimates and analysis results

`get_results()` provides access to the estimates and summaries produced
by models and analyses in `latents`. It supports examination of the
estimated latent classes, their membership probabilities, and model-fit
statistics, as well as the preparation of results for reporting and
further analysis.

The results available depend on the supplied object. Fitted models
provide measurement estimates, class proportions, posterior membership
probabilities, and other model-specific quantities, such as transition
probabilities. Enumeration and comparison objects provide summaries of
candidate models, while diagnostic and bootstrap objects provide their
corresponding checks and test results.

The `what` argument selects the requested table. Each table is returned
as a standard R data frame, allowing its contents to be inspected,
filtered, combined, or exported. Omitting `what` returns the object's
primary table; `what = "all"` returns a named collection of the tables
included by default for that object.

## Usage

``` r
get_results(x, what = NULL, ...)
```

## Arguments

- x:

  A fitted model of this package, an enumeration grid, a model
  comparison, a bootstrap comparison, a starting-value set, a
  classification-diagnostics object, or the summary of any of them.

- what:

  Name of the table to return, or `"all"` for the collection of tables
  included by default for the object. `NULL` (default) returns the
  primary table: the measurement model for a fit, the candidate grid for
  an enumeration, the test row for a bootstrap comparison and the
  entropy table for a diagnostics object. A name the object does not
  offer raises `latents_bad_argument`, listing the names it does offer.

- ...:

  Arguments for the requested table. Each table accepts only the
  arguments listed for it under *Tables*; any other argument raises
  `latents_bad_argument`.

## Value

A base `data.frame` for a single `what`, whose columns are given under
*Tables*, or a named list of such data frames for `what = "all"`.

## Details

Some tables calculate additional quantities when requested. For example,
measurement tables may obtain standard errors through
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md).
Their availability therefore depends on the fitted model and the
requirements of the corresponding calculation.

The table catalogue below specifies the contents and accepted arguments
of each table, including results that require an explicit request.

Every object also has [`print()`](https://rdrr.io/r/base/print.html),
which shows a compact overview,
[`summary()`](https://rdrr.io/r/base/summary.html), which collects its
main tables (taking `level` and `vcov_type` where the model has standard
errors), and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), which
returns its primary table, the same as `get_results(x)`.

Before combining supplied `data` with model results, the function checks
row alignment using columns retained by the fitted object, such as
identifiers, occasions, and indicators. Numeric values serialized by
[`utils::write.csv()`](https://rdrr.io/r/utils/write.table.html) are
also accepted when serialization preserves distinct values. Columns
absent from the fitted object cannot establish alignment. Inference and
model comparison apply additional checks to verify that supplied data
reproduce the estimation data.

## Tables

The following tables apply to profile and class models and, where
supported, membership-covariate and transition models. Tables for other
model families are described in *Tables of the other model families*.

- `"profiles"`:

  The Gaussian measurement model, one row per profile and continuous
  indicator: `profile`, `indicator`, `mean`, `variance`,
  `standard_deviation`, `mean_standard_error`,
  `variance_standard_error`. Zero rows, with the same columns, when
  every indicator is categorical. Takes `scale` and `data`.

- `"responses"`:

  The categorical measurement model, one row per profile, categorical
  indicator and category: `profile`, `indicator`, `category`,
  `probability`, `threshold`, `probability_standard_error`. `threshold`
  is `qlogis(P(y <= category))` and is `NA_real_` for each indicator's
  final category; a probability on its bound has an `NA_real_` standard
  error. Zero rows, with the same columns, when no indicator is
  categorical. Takes `data`.

- `"covariances"`:

  The within-profile residual covariance matrices, one row per profile
  and ordered pair of continuous indicators: `profile`, `indicator`,
  `indicator_2`, `covariance`. Under the diagonal parameterization the
  off-diagonal covariances are reported as exact zeros.

- `"ordinal"`, `"count_means"`:

  Fits with ordinal or count indicators only: the ordinal and count
  measurement models. Take `data`.

- `"counts"`:

  Effective class sizes, one row per class at each level of the fit:
  `level`, `class`, `effective_count` (the summed posteriors) and
  `effective_proportion`. Effective counts use posterior probabilities;
  modal counts in `"classification"` use discrete class assignments.

- `"model"`:

  One row describing the fit: its dimensions, likelihood, parameter
  counts, `aic`, `bic_groups`, `bic_individual` and its convergence and
  boundary diagnostics. The columns differ by family (a covariate fit
  adds its predictor counts); `"information_criteria"` has the same
  columns for every family.

- `"posteriors"`:

  One row per individual and profile: `row`, `group`, `profile`,
  `posterior`, `modal`. Takes `format`. Excluded observations retain
  their rows, with `posterior` and `modal` set to `NA`. In wide format
  their assigned `profile` and posterior columns are `NA`.

- `"assignments"`:

  One row per observation: the columns the fit was built from (or those
  of `data`), then `profile`, `group_class` where the model has one,
  `uncertainty` (one minus the posterior of the assigned profile) and
  one `posterior_profile_*` column per profile, plus `posterior_noise`
  for a fit with a noise component (whose cases are profile 0). Takes
  `data` and `truth`; with `truth` it returns a recovery table instead
  (see *Recovery*). Row order is preserved. Under the default
  `drop_empty = TRUE`, observations with no observed indicators have
  `NA` for `profile`, `uncertainty` and the posterior columns. Group
  assignments are missing only when the entire group is excluded; an
  empty observation can therefore retain a group assignment supported by
  other observations. Zero-weight units are unclassified under either
  `drop_empty` setting.

- `"classification"`:

  Classification quality, one row per level and class: `level`, `class`,
  `n_modal`, `proportion_modal`, `estimated_n`, `estimated_proportion`,
  `average_posterior`, `odds_correct_classification`. Takes `level`.
  Under sampling weights every column that estimates the population is
  weighted: `estimated_n` and `estimated_proportion` are weighted
  posterior totals and `average_posterior` a weighted mean, so the odds
  of correct classification compare the two on the same population;
  `n_modal` and `proportion_modal` stay counts of the units supplied.
  The odds are `Inf` for a class classified perfectly (average posterior
  one to within rounding) and `NA` where they are undefined: a class
  with no modal members, or a single class, whose proportion is one.

- `"average_posteriors"`:

  For the units modally assigned to each class, their mean posterior for
  every class: `level`, `assigned_class`, `class`, `n_assigned`,
  `average_posterior`. Rows sum to one within each assigned class. Under
  sampling weights `average_posterior` is the weighted mean and
  `n_assigned` stays a count. This conditions on the assigned class,
  whereas `"classification_errors"` conditions on the true class. Takes
  `level`.

- `"classification_errors"`:

  The probability that a unit in one true class is assigned to another:
  `level`, `true_class`, `assigned_class`, `probability`, summing to one
  within each `true_class`. Takes `level`. It is
  `sum_i w_i p_it 1(W_i = s) / sum_i w_i p_it` for true class `t` and
  assigned class `s` (Vermunt 2010), with `w_i` the unit's sampling
  weight under sampling weights and 1 otherwise, so a weighted fit
  reports the population's error rates. On a membership-covariate fit
  this is a marginal summary over the fitted units, not conditional on
  the covariates. A fit with a noise component raises
  `latents_unsupported_noise`, as
  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
  does, and `"all"` and
  [`summary()`](https://rdrr.io/r/base/summary.html) omit the table.

- `"bch_weights"`:

  The BCH inverse-error weights, one row per unit and class: `level`,
  `unit`, `assigned_class`, `class`, `weight`. Takes `level`; a weighted
  regression should use one level, not `"both"`. They invert the
  `"classification_errors"` matrix, weighted as it is. On a
  membership-covariate fit the weights invert the marginal error matrix
  and are not a conditional three-step correction. Not available for a
  fit with a noise component, as for `"classification_errors"`.

- `"entropy"`:

  One row per level: `level`, `n_classes`, `n_units`, `entropy_sum`,
  `relative_entropy`. Under sampling weights `entropy_sum` is
  `sum_i w_i e_i` and `relative_entropy`
  `1 - sum_i w_i e_i / (sum_i w_i log K)`; `n_units` stays a count.

- `"replication"`:

  One row: `n_starts` (starts run to the end), `n_best_replicated`
  (those whose log likelihood is within `tolerance` of the best),
  `share`, `tolerance`, `best_start` (the selected start's index) and
  `log_likelihood`. A best solution reached by one start of several may
  not be the global maximum: refit with more starts.

- `"residuals"`:

  Bivariate residuals for every pair of indicators, checking
  within-profile independence: `profile`, `indicator_1`, `indicator_2`,
  `kind` (`"gaussian"`, `"categorical"`, `"mixed"` or `"count"`),
  `observed`, `expected`, `residual`, `effective_n`, `statistic`, `df`,
  `p_value`, `p_adjusted`. Takes `data`, `by` and `adjust`. See
  *Bivariate residuals* for what each kind computes.

- `"information_criteria"`:

  Every information criterion the package computes (see *Information
  criteria*). Takes `format` and `definitions`.

- `"sequences"`:

  The assignments in occasion order, one row per classified observation
  (a row with nothing observed is left out, or kept under
  `drop_empty = FALSE`): `group`, `group_class`, `time`, `profile`.
  Takes `format`. Raises `latents_no_time` for a fit made without
  `time`.

- `"sequence_summary"`:

  One row per group class: `group_class`, `groups`, `observations`,
  `mean_length`, `median_length`, `shortest`, `longest`, `complete`,
  `gaps`. Summarises classified groups only; lengths describe their
  recorded occasions, including occasions with missing indicators.
  Raises `latents_no_time` for a fit made without `time`.

- `"starts"`:

  One row per EM start: `start`, `log_likelihood`, `converged`,
  `iterations`, `error` and, where the family records it, `boundary`.
  Under `search = "emEM"`, one row per start of the pool, with its
  `kind`, its `stage` (`"final"` when run to convergence, `"screened"`
  otherwise) and its `screened_log_likelihood`.

- `"data"`:

  The data the model was fitted to, one row per observation in input
  order: the identifier (except for a single-level fit), the occasion
  where there is one, the sampling weights of a weighted fit (as
  rescaled for fitting) and the indicators under their original names.
  Columns the model did not use, including membership covariates, are
  not included;
  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
  and [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md)
  therefore take `data`.

Fits with discrete group classes also offer:

- `"group_posteriors"`:

  One row per observed group and group class: `group`, `group_size`,
  `log_likelihood`, `group_class`, `posterior`, `modal`. Takes `format`.
  Excluded groups retain their rows, with missing posterior
  probabilities and modal assignments, as for `"posteriors"`.

Tables specific to one family or object:

- `"gaussian_diagnostics"`:

  Gaussian profile, membership-covariate, transition and cross-level
  fits that record these diagnostics. One row per Gaussian component and
  measurement context: `context`, `profile`, `gaussian_mass`,
  `positive_rows`, `weight_effective_n`, `min_variance_ratio`,
  `min_standardized_eigenvalue`, `boundary`. `gaussian_mass` sums
  posterior times normalized sampling weight over rows with at least one
  observed Gaussian indicator; `positive_rows` counts those rows (the
  same for every component, since a Gaussian density is positive
  everywhere). `weight_effective_n` is the concentration measure
  `(sum(a)^2 / sum(a^2))` for those contributions, not an identification
  guarantee; Gaussian mass is not this effective sample size.
  `min_variance_ratio` is the smallest component variance divided by its
  observed sampling-weighted marginal variance on the fitted scale,
  including between-component variation. `min_standardized_eigenvalue`
  is the smallest eigenvalue of the component covariance after
  standardizing every indicator by its observed sampling-weighted
  marginal standard deviation. Missing reference variation gives `NA`
  ratios. `boundary` concerns the absolute variance or eigenvalue floor.
  These descriptive measures use no automatic spike cutoff and do not
  alter start selection.

- `"profile_probabilities"`, `"stages"`:

  A `multilpa` fit. `"profile_probabilities"` has one row per group
  class and profile: `group_class`, `profile`, `probability`,
  `group_class_probability`. `"stages"` has one row per estimation
  stage: `stage`, `group_classes`, `fixed`, `log_likelihood`,
  `parameters`, `parameters_with_measurement`, `converged`. An ordinary
  fit has one `"joint"` row and a
  [`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
  result has two; `fixed` is `NA_character_` for a stage that estimated
  every block.

- `"coefficients"`:

  A covariate fit. One row per estimated membership coefficient: `level`
  (`"profile"` or `"group"`), `outcome`, `term`, `parameter` (always
  `"logit"`), `estimate`. Standard errors are given by
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md).

- `"transitions"`, `"initial"`, `"sequence_lengths"`:

  A transition fit. `"transitions"` has one row per group class and
  ordered pair of profiles: `group_class`, `from`, `to`, `probability`,
  `expected_count`, `stable`, `estimated`, `group_class_probability`; it
  takes `estimated` and `stable`. `"initial"` has one row per group
  class and profile: `group_class`, `profile`, `probability`,
  `prevalence`, `group_class_probability`. `"sequence_lengths"` has one
  row per group: `group`, `group_class`, `occasions`, `observations`,
  `complete`.

- `"candidates"`, `"best"`, `"classification"`:

  An enumeration grid. `"candidates"` has one row per candidate model:
  its name (`candidate`, first), its class counts, covariance `model`,
  log likelihood, parameter count, every criterion under both
  sample-size conventions, both entropies, the average posterior
  probability (AvePP) summaries of each level (`avepp_min`, `avepp_max`,
  `avepp_mean`, `avepp_below` and the same four prefixed `group_`) and
  the convergence, boundary, replication, warning and error diagnostics,
  including `min_standardized_eigenvalue` and `min_gaussian_mass` where
  Gaussian diagnostics are available. These are descriptive summaries,
  without an automatic cutoff or change to candidate selection. Failed
  candidates are kept with `NA` estimates and their error text. It takes
  `cutoff` (default `0.8`): `avepp_below` and `group_avepp_below` count
  the classes whose AvePP is below it. A class with no modally assigned
  unit has no AvePP; it is excluded from the minimum, maximum and mean
  and counted as below the cutoff. A level with one class has `NA`
  summaries. `"best"` has one row per information criterion and
  sample-size convention: `criterion`, `convention`, `candidate` (the
  converged candidate with the lowest value), `n_profiles`,
  `n_group_classes`, `model` and `value`; `what = "criteria"` is
  accepted as a deprecated name for it. `"classification"` has one row
  per fitted candidate, level and class: `candidate`, `n_profiles`,
  `n_group_classes`, `model`, then the columns of a fit's
  `"classification"` table and that level's `relative_entropy`. Failed
  candidates and single-class levels contribute no rows.

- `"deltas"`, `"knees"`:

  Model-selection evidence of an enumeration grid, each taking `using`
  (default `"bic_groups"`). Neither selects a model. `"deltas"` has one
  row per candidate: `candidate`, its class counts and `model`,
  `criterion`, `value`, `delta_profile` (the decrease in the criterion
  from the converged candidate with one profile fewer, the same group
  classes and the same structure; positive means this candidate is
  lower), `from_profile` naming that candidate, and `delta_group_class`
  with `from_group_class` for one group class fewer. `value` and both
  deltas are `NA` for a candidate that did not converge, and a delta is
  `NA` when the smaller neighbour is missing or did not converge.
  `"knees"` has one row per line of the grid: `direction` (`"profiles"`,
  along the profiles at fixed group classes; `"group_classes"`, along
  the group classes at fixed profiles; `"parameters"`, every candidate
  by its number of parameters, one line per structure), `line`,
  `criterion`, `from`, `to`, `n_points`, `knee_candidate`, `distance`,
  `knee_without_largest` (the knee when the last candidate of the line
  is left out), `stable` and `note`. The knee is the candidate furthest
  below the straight line from the first to the last candidate, measured
  perpendicular to it with both axes rescaled to 0-1. This is a
  simplified rule inspired by Kneedle (Satopaa et al., 2011), without
  its smoothing or sensitivity thresholds, so a very small positive
  distance can still name a knee. A line needs three candidates for a
  knee and four for `knee_without_largest` and `stable`; otherwise they
  are `NA` and `note` gives the reason.

- `"local_dependence"`:

  For a fit with `local_dependence`: `pair`, `kind` (`"categorical"` or
  `"continuous"`), `profile`, `item_1`, `category_1`, `item_2`,
  `category_2`, `probability`, `independent`, `ratio`, `covariance` and
  `correlation`. A categorical pair has one row per profile and category
  combination, including unobserved combinations subject to the
  `min_probability` floor and included in the parameter count: the
  fitted joint `probability`, the product of the two marginal
  probabilities (`independent`, the joint probability under independence
  within the profile) and their `ratio`; its covariance columns are
  `NA`. A continuous pair has one row per profile, with the pair's
  `covariance` and the `correlation` it implies; its category and
  probability columns are `NA`. A fit without pairs raises
  `latents_no_local_dependence`, and `what = "all"` leaves the table
  out.

- `"direct_effects"`:

  For a fit with `direct_effects`, one row per effect: `indicator`,
  `covariate`, `estimate` (the change in the indicator's mean per unit
  of the covariate, the same in every profile), `standard_error`,
  `statistic`, `p_value`, `conf_low` and `conf_high` (Wald, from
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)).
  A fit without direct effects raises `latents_no_direct_effects`, and
  `what = "all"` leaves the table out.

- `"lmr"`:

  The Vuong-Lo-Mendell-Rubin (VLMR) and adjusted Lo-Mendell-Rubin (aLMR)
  tests of an enumeration grid (Vuong 1989; Lo, Mendell and Rubin 2001;
  Vermunt 2024), one row per pair of converged candidates that differ by
  one profile (`direction = "profiles"`) or one group class
  (`"group_classes"`) at the same structure: `direction`,
  `null_candidate`, `alternative_candidate`, `n_profiles`,
  `n_group_classes` and `model` of the alternative, `n_units` (the
  independent units scored: groups, or observations without groups),
  `parameters_added`, `lr` (twice the log-likelihood difference),
  `vlmr_p`, `adjusted_lr`
  (`lr / (1 + 1 / (parameters_added * log(n_units)))`), `almr_p`,
  `reference_mean` and `reference_sd` of the weighted chi-square
  reference distribution, whose weights are the eigenvalues of Vuong's
  matrix of the two models' scores and Hessians, and `note`. The
  p-values are computed by Imhof's method to about `1e-9`. A pair whose
  information matrix is singular, or whose candidate is on a boundary,
  keeps its row with `NA` statistics and the reason in `note`; warnings
  raised while scoring a candidate are recorded there too. Every
  candidate's information is computed, so the table is built only when
  asked for, and `what = "all"` leaves it out. A grid with no such pair
  raises `latents_no_lmr`.

- `"tree"`:

  The flows of
  [`plot_tree()`](https://pak.dynasite.org/latents/reference/latents-plots.md),
  on every grid, taking its arguments `n_group_classes`, `n_profiles`
  (`n_classes` on a regression grid), `model`, `candidates`, `from`,
  `to` and `level` (default `"individuals"`). One row per pair of
  classes in consecutive rows of each tree: `tree` (the tree's name, as
  in the plot list), `level`, `from_candidate`, `from_class`,
  `from_share`, `to_candidate`, `to_class`, `to_share`, `mass` (the
  share of all units the two classes share, the sum over units of the
  product of their posterior probabilities, weighted by sampling weights
  at the individual level), `share_of_from` and `share_of_to` (`mass` as
  a part of each class). Every pair is listed, including those that
  share no mass. Units are the rows and groups the fit counts (see
  `drop_empty`), weighted by their sampling weights at both levels. A
  grid with no tree raises `latents_no_tree`, a fit with a noise
  component `latents_unsupported_noise`, and `what = "all"` leaves the
  table out. The summary of a grid carries the default table only,
  without these arguments.

- `"comparison"`, `"classification"`:

  A
  [`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)
  result. `"comparison"` is its table, one row per model.
  `"classification"`, offered when every model is a profile, class,
  membership-covariate or transition model, has one row per model, level
  and class: `model`, then the columns of a fit's `"classification"`
  table and that level's `relative_entropy`; single-class levels are
  left out. The models themselves are returned by
  [`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md).

- `"test"`, `"replicates"`:

  A bootstrap comparison. `"test"` is the one-row result; `"replicates"`
  has one row per simulated dataset.

- `"profiles"`, `"responses"`, `"profile_probabilities"`,
  `"covariances"`:

  A starting-value set: its Gaussian means and variances, response
  probabilities, profile probabilities and, one row per profile and
  ordered pair of continuous indicators, its covariances (`profile`,
  `indicator`, `indicator_2`, `covariance`).

- `"entropy"`, `"classification"`, `"average_posteriors"`,
  `"residuals"`, `"replication"`:

  A
  [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  object: the tables of the same names described above.

## Tables of the other model families

Each family's first table is its default `what`; `"all"` returns a named
list of every table (without `"recovery"`, which needs `data` and
`truth`). Tables with intervals take `level`, and `vcov_type`
(`"observed"`, `"robust"` or `"opg"`; `NULL` reuses the type stored with
a regression fit). These objects also have
[`coef()`](https://rdrr.io/r/stats/coef.html), `vcov(type = )`,
`confint(parm = , level = )`,
[`logLik()`](https://rdrr.io/r/stats/logLik.html) and
[`nobs()`](https://rdrr.io/r/stats/nobs.html) methods.
[`confint()`](https://rdrr.io/r/stats/confint.html) of a mixture
regression or growth mixture uses the covariance stored with the fit and
refuses any other argument with `latents_bad_argument` (an interval on
another covariance is `get_results(x, "coefficients", vcov_type = )`);
[`vcov()`](https://rdrr.io/r/stats/vcov.html) warns
`latents_unused_argument` for an argument it does not take, such as
`vcov_type =` in place of `type =`.

- Mixture regressions
  ([`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)):

  `"coefficients"` (one row per class and term, with Wald statistics,
  `p_adjusted` by Benjamini-Hochberg over non-intercepts and
  `exp_estimate` for the binomial, Poisson and ordinal families; ordinal
  thresholds lead each class as terms `"threshold:<lower>|<upper>"`),
  `"classes"`, `"membership"` (multinomial logits and odds ratios),
  `"group_classes"` (rows for two-level fits only), `"fit"` (`bic`
  penalised by the number of independent units, `bic_rows` by the number
  of rows, as flexmix reports it), `"assignments"`, `"groups"`,
  `"fitted"`, `"starts"` (`stage` `"screened"` or `"completed"`),
  `"classification"` and `"recovery"` (`by = "class"` or
  `"group_class"`).

- Growth mixtures (`mixture_regression(random = )`):

  `"coefficients"`, `"classes"` (each class's share with its interval,
  `expected_persons`, `assigned_persons`, `average_posterior` and
  `residual_sd`), `"membership"`, `"fit"`, `"assignments"` (one row per
  person: `class`, `probability` and one `posterior_*` column per
  class), `"starts"` and `"recovery"`, as for a mixture regression, plus
  `"random"` (random-effect variances, standard deviations and
  correlations, with intervals on the log and Fisher-z scales),
  `"random_effects"` (one row per person), `"trajectories"` (each
  class's mean trajectory with a band) and `"individual"` (each person's
  predicted curve); the last two take `time`. They have no `"groups"`,
  `"fitted"` or `"classification"` table. A multilevel growth mixture
  adds `"group_classes"` and `"clusters"`, and its `"fit"` uses the
  number of clusters in `bic` (persons in `bic_persons`).

- Additive fits (`multilpa(family = "additive")`):

  `"parameters"` (class means, between- and within-group variances and
  class weights; variance intervals on the log scale, weights on the
  logit scale; `adjust` sets the p-value adjustment, `"none"` by
  default), `"group_classes"`, `"groups"`, `"intercepts"` and
  `"conditional_intercepts"` (each group's posterior intercept, overall
  or given each class), `"assignments"`, `"classification"`, `"fit"`,
  `"recovery"` and `"starts"`. `effective_groups` in `"group_classes"`
  is `w^2 (1 - w) / Var(w)`, the number of perfectly classified groups
  that would estimate the class weight as precisely; below 50 the
  class's intervals carry a `latents_weak_class` warning.

- Cross-level fits (`multilpa(family = "full_cross_level")` or
  `"restricted_cross_level"`):

  `"profiles"` (individual profiles, then the group classes' means of
  group means), `"group_classes"`, `"composition"` (profile shares
  within each class: model-implied for the full family, a
  posterior-weighted description for the restricted one), `"groups"`,
  `"assignments"`, `"fit"` and `"starts"`. The likelihood is a working
  likelihood, comparable only between cross-level fits of the same data.

- General transition fits
  ([`lta()`](https://pak.dynasite.org/latents/reference/lta.md) with
  `transitions = "occasion"`, transition or initial covariates,
  `measurement = "occasion"`, `order = 2`, `mover_stayer = TRUE`,
  `weights`, or `ordinal`/`count` indicators):

  `"transitions"` (model-implied, averaged over the observed groups,
  with their sampling weights in a weighted fit; `"initial"` and
  `"second_order_transitions"` likewise), `"transition_coefficients"`
  (log odds of each move against staying, with odds ratios),
  `"second_order_transitions"` and `"second_order_coefficients"`
  (`order = 2`), `"initial"` and `"initial_coefficients"`, `"profiles"`,
  `"responses"` (categorical indicators: `occasion`, `profile`,
  `indicator`, `category`, `probability`), `"ordinal"` and
  `"count_means"` (ordinal and count indicators), `"group_classes"`,
  `"assignments"`, `"fit"` and `"starts"`.

- Pooled imputations
  ([`pool_imputations()`](https://pak.dynasite.org/latents/reference/pool_imputations.md)):

  `"estimates"` (the pooled table), `"imputations"` (each imputation's
  aligned estimates) and `"fits"` (one row per imputation).

- Transition and group-class family grids (`enumerate_classes(time = )`,
  `enumerate_classes(family = )`):

  `"candidates"` (one row per candidate, its name first, with its
  likelihood, criteria, entropy and diagnostics) or `"best"`, the row
  with the lowest `criterion` (`"bic"` by groups, the default, `"aic"`
  or `"bic_individual"`). `"best"` reports the lowest value of a
  criterion; it does not select a model. With one group class every
  family is the same model, so those rows share a likelihood. Candidates
  are extracted by name with
  [`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md);
  transition grids can also be passed to
  [`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md).

## What `"all"` returns

A named list of tables included by default for the object, in catalogue
order. Some results require an explicit request, including the `"lmr"`
table for an enumeration object. Each included table is built with its
default arguments; supplying any further argument raises
`latents_bad_argument`. A table this object cannot produce (for example
a sequence table for a fit made without `time`, or bivariate residuals
for a family without discrete group classes) is omitted, so the names of
the list show what was available. Only the conditions that mean the
table does not exist for this object are skipped; any other error is
raised.

## Table arguments

These are passed through `...`; each table takes only those listed for
it under *Tables*.

- `data`:

  A data frame with one row per observation of the fit, in the order it
  was fitted in. `NULL` uses the columns the fit carries. The standard
  errors of `"profiles"` and `"responses"` come from
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  on these data. Without `data`, a fit that cannot supply them
  (bound-active, unconverged, singular information, a noise component)
  gives `NA` with a `latents_no_standard_errors` message; with `data`,
  the underlying condition is raised.

- `scale`:

  For `"profiles"`: `"raw"` (default) reports estimates in input units;
  `"standardized"` divides each indicator's deviation from its grand
  mean by that indicator's observed standard deviation, as
  `plot(x, scale = "standardized")` draws. The standard errors are
  divided by the same constant.

- `format`:

  For the posterior tables and `"sequences"`, `"long"` (default) gives
  one row per unit and class and `"wide"` one row per unit with one
  column per class. For `"information_criteria"`, `"wide"` (default)
  gives one row and `"long"` one row per criterion and sample-size
  convention.

- `level`:

  The classification level to report: `"individuals"`, `"groups"` or
  `"both"`. `NULL` (default) gives `"both"` for a fit with discrete
  group classes and `"individuals"` otherwise, as in
  [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md).

- `truth`:

  For `"assignments"`, one or more column names of `data` holding known
  labels to compare the assignments with. See *Recovery*.

- `by`:

  For `"residuals"`, `"profile"` (default) assesses each profile
  separately and `"overall"` pools them.

- `adjust`:

  For `"residuals"`, the multiplicity correction applied across the
  indicator pairs: `"none"` (default) or one of the
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) methods.

- `definitions`:

  For `"information_criteria"` in long format, `TRUE` adds each
  criterion's formula and reference. Default `FALSE`.

- `estimated`, `stable`:

  For `"transitions"`, `TRUE` or `FALSE` restricts the table to the rows
  estimated (or not) from observed occupancy, or to the diagonal (or
  off-diagonal) cells. `NULL` (default) keeps every row.

## Recovery

`what = "assignments"` with `truth` returns a cross-tabulation of the
model's labels against known ones instead of the per-observation table:
one row per class and truth value, with columns `assignment` (the model
label compared against), `class`, `truth` (the column name), `value`,
`n` and `proportion` (the share of units with that truth value assigned
to that class). Several truth columns are stacked. Excluded observations
or groups appear under `class = NA`, so the table distinguishes
unclassified units from assignments to an estimated class.

A truth column that is constant within every group is compared against
`group_class`; one that varies within any group is compared against
`profile`. The `assignment` column records which comparison was made.
Group-level recovery counts each group once. Unused factor levels are
omitted. For a fit without discrete group classes every column is
compared against `profile`.

Class labels are arbitrary, so a truth value assigned mostly to a class
with a different label indicates relabelling, not misclassification.

## Information criteria

Let `q` be the number of free parameters, `n` the chosen sample size,
and `EN` the entropy of the complete-data classification (defined
below). The criteria are `deviance = -2L`, `aic = -2L + 2q`,
`bic = -2L + q log(n)`, `sabic = -2L + q log((n + 2) / 24)`,
`caic = -2L + q (log(n) + 1)`, `awe = -2(L - EN) + 2q (1.5 + log(n))`,
`icl = -2L + q log(n) + 2 EN`, `kic = -2L + 3(q + 1)`, and
`clc = -2L + 2 EN`. `clc`, like `icl` and `awe`, uses the entropy sum;
`tidyLPA` reports a `CLC` based on relative entropy, and the two values
are not comparable. The `kic` form is Cavanaugh's (1999) criterion as
`tidyLPA` computes it (Rosenberg et al. 2018), with `q` the number of
free parameters; Cavanaugh's general form `-2L + 3q` is 3 smaller for
every model, so differences between models and the model with the
smallest value are the same under both.

`EN` is the soft entropy of the complete-data classification, the
entropy term of the classification log likelihood `log L_c = log L - EN`
of ICL (Biernacki, Celeux and Govaert 2000), AWE (Banfield and Raftery
1993) and CLC (Biernacki and Govaert 1997). The latent classification of
a two-level model is the group class of every group and the profile of
every row, so
`EN = sum_j H(C_j | y_j) + sum_j sum_h P(C_j = h | y_j) sum_i H(X_ij | C_j = h, y_ij)`:
the group-class entropy plus each row's profile entropy within its group
class. With one group class, and in every single-level model, it is the
entropy of the profile posteriors. It is the soft form, the posterior
expectation of `-log P(classes | y)` (McLachlan and Peel 2000, section
6.10.3), not the entropy of the modal (hard) assignment that Biernacki,
Celeux and Govaert (2000) also discuss. A fit with several group classes
saved by latents 0.97.1 or earlier lacks the rows' within-class
entropies, and its `awe`, `icl` and `clc` are `NA`; refit it to obtain
them.

Every criterion that depends on a sample size is reported under both
multilevel conventions: `"groups"` with `n` the number of groups and
`"individuals"` with `n` the number of rows. `EN` is the same under
both, so the conventions differ only in `n`, and in a single-level model
they coincide. Use one convention consistently across compared
candidates. `deviance`, `aic` and `kic` carry
`convention = NA_character_`; `clc`, which uses no sample size, is the
same in both rows. A transition model
([`lta()`](https://pak.dynasite.org/latents/reference/lta.md)) keeps one
entropy per level instead: its profiles form a Markov chain within a
person, whose path entropy is not the sum of the occasions' entropies.
Its `"groups"` criteria use group-class entropy, whereas its
`"individuals"` criteria use occasion-level profile entropy. These
definitions are reported by
`get_results(x, "information_criteria", format = "long", definitions = TRUE)`.
Individual sample sizes exclude rows with no observed indicators, and
group sample sizes exclude groups in which every row has none
(`n_informative_groups`); the entropy and the classification statistics
exclude the same rows and groups, as in Mplus. Units of sampling weight
zero are in no sample size and no entropy, under either `drop_empty`.
For a fit with sampling weights, `L` is the pseudo log likelihood and
`EN` the sampling-weighted entropy, each unit's terms multiplied by its
weight, so that integer weights give the entropy of the data with each
unit repeated that many times. The `"entropy"` table,
`relative_entropy`, the average posteriors, the odds of correct
classification, the classification errors and the BCH weights are
weighted the same way; modal counts stay counts.

## Classification quality

The odds of correct classification for class `k` is
`(p / (1 - p)) / (r / (1 - r))`, where `p` is the average posterior in
the assigned class and `r` is the model-estimated class proportion.
Values near one indicate classification no better than the class
proportions alone. A value of five is a conventional descriptive
benchmark, rather than a formal test of class separation. It is
`NA_real_` when `average_posterior` is missing or zero, or when
`estimated_proportion` is zero or one. A perfect `average_posterior` of
one gives `Inf` at an interior estimated proportion; a single-class
solution always reports `NA_real_`. `average_posterior` is `NA_real_`
for a class with no modal members. `n_modal` and `estimated_n` use modal
assignments and posterior probabilities, respectively, and need not
coincide.

`relative_entropy` is `1 - EN_level / (n log K)` (Ramaswamy et al.
1993), between zero and one, higher meaning sharper classification, with
`EN_level` the level's own entropy sum, reported as `entropy_sum`, lower
meaning sharper. With one group class the individual level's
`entropy_sum` is the `EN` of the information criteria; with several,
that `EN` adds the group-class entropy to the profile entropy within
each group class (see *Information criteria*). `relative_entropy` is
`NA_real_` for a level with a single class.

A fit with individual profiles but no discrete group classes reports
individuals only for `level = "both"` and raises
`latents_no_group_classes` for `level = "groups"`.

## Bivariate residuals

`"residuals"` checks, pair by pair, the assumption that indicators are
independent within a profile, weighting each observation by its
posterior (times its sampling weight). Every pair of indicators has a
row (per profile under `by = "profile"`, once under `by = "overall"`).
`data` must be in the order the model was fitted in; the alignment is
checked. Large residuals suggest associations that the measurement model
does not adequately represent. Calculations for each `kind` are
described below:

- `"gaussian"`, two Gaussian indicators: `observed` is the
  posterior-weighted correlation of the residuals from the profile's
  means (pooled over the profiles under `"overall"`), `expected` the
  correlation the fitted covariance implies (zero for a diagonal one).
  `statistic` is the difference of their Fisher z transforms times
  `sqrt(effective_n - 3)`, referred to the standard normal (Fisher,
  1921), with a two-sided p-value; `df` is `NA`. This is an approximate
  reference distribution for a posterior-weighted diagnostic.

- `"categorical"`, two categorical or ordinal indicators (an ordinal
  indicator's categories are its levels): `statistic` is the Pearson X^2
  of the posterior-weighted two-way table against the table the model
  implies, which is Latent GOLD's bivariate residual times its `df`, and
  `p_value` refers it to chi-square with `df = (R - 1)(C - 1)`.
  `observed` and `expected` are the two tables' Cramer's V (Cramer,
  1946), zero for the model's table within a profile. The chi-square
  reference ignores that the model's table is estimated from the same
  data and can therefore be conservative (Oberski, van Kollenburg &
  Vermunt, 2013). Interpret `p_value` together with `residual` and
  `statistic`; a non-significant result does not establish conditional
  independence.

- `"mixed"` (a Gaussian and a categorical or ordinal indicator) and
  `"count"` (any pair with a count indicator): `observed` is the
  posterior-weighted correlation of the two indicators' values within
  the profile, centred on the profile's weighted means (pooled over the
  profiles under `"overall"`); a categorical or ordinal indicator enters
  as its category number in category order, which for a nominal
  indicator with more than two categories depends on that order.
  `expected` is zero, the correlation of independent indicators, and
  `statistic`, `df` and `p_value` are `NA`: no reference distribution is
  computed for these pairs.

A half-observed row of a joined `local_dependence` pair enters that
pair's tables only where both items are observed. `p_adjusted` applies
the correction named by `adjust` across the pairs that have a p-value. A
fit with sampling weights uses them in the descriptive associations and
`effective_n`; design-based inference is not implemented, so
`statistic`, `df`, `p_value` and `p_adjusted` are `NA` for weighted
fits.

## References

Fisher, R. A. (1921). On the "probable error" of a coefficient of
correlation deduced from a small sample. *Metron*, 1, 3–32.

Cramer, H. (1946). *Mathematical Methods of Statistics*. Princeton
University Press.

Schwarz, G. (1978). Estimating the dimension of a model. *Annals of
Statistics*, 6, 461–464.

Sclove, S. L. (1987). Application of model-selection criteria to some
problems in multivariate analysis. *Psychometrika*, 52, 333–343.

Bozdogan, H. (1987). Model selection and Akaike's information criterion.
*Psychometrika*, 52, 345–370.

Banfield, J. D., & Raftery, A. E. (1993). Model-based Gaussian and
non-Gaussian clustering. *Biometrics*, 49, 803–821.

Biernacki, C., & Govaert, G. (1997). Using the classification likelihood
to choose the number of clusters. *Computing Science and Statistics*,
29, 451–457.

Cavanaugh, J. E. (1999). A large-sample model selection criterion based
on Kullback's symmetric divergence. *Statistics and Probability
Letters*, 42, 333–343.

Rosenberg, J. M., Beymer, P. N., Anderson, D. J., Van Lissa, C. J., &
Schmidt, J. A. (2018). tidyLPA: An R package to easily carry out latent
profile analysis (LPA) using open-source or commercial software.
*Journal of Open Source Software*, 3(30), 978.

McLachlan, G., & Peel, D. (2000). *Finite Mixture Models*. Wiley.

Ramaswamy, V., DeSarbo, W. S., Reibstein, D. J., & Robinson, W. T.
(1993). An empirical pooling approach for estimating marketing mix
elasticities with PIMS data. *Marketing Science*, 12, 103–124.

Biernacki, C., Celeux, G., & Govaert, G. (2000). Assessing a mixture
model for clustering with the integrated completed likelihood. *IEEE
Transactions on Pattern Analysis and Machine Intelligence*, 22, 719–725.

Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent
structure models with categorical variables. *Political Analysis*, 12,
3–27.

Nagin, D. S. (2005). *Group-Based Modeling of Development*. Harvard
University Press.

Vermunt, J. K. (2010). Latent class modeling with covariates: two
improved three-step approaches. *Political Analysis*, 18, 450–469.

Satopaa, V., Albrecht, J., Irwin, D., & Raghavan, B. (2011). Finding a
"kneedle" in a haystack: Detecting knee points in system behavior.
*Proceedings of the 31st International Conference on Distributed
Computing Systems Workshops*, 166–171.

Oberski, D. L., van Kollenburg, G. H., & Vermunt, J. K. (2013). A Monte
Carlo evaluation of three methods to detect local dependence in binary
data latent class models. *Advances in Data Analysis and
Classification*, 7, 267–279.

## See also

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), which
returns the primary table alone,
[`summary()`](https://rdrr.io/r/base/summary.html), which prints every
table, and
[`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
for the classification tables with a diagnostic summary.

## Examples

``` r
fit <- multilpa(
  course_engagement,
  vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
  id = "student", n_profiles = 2, n_group_classes = 2, n_starts = 4,
  seed = 1
)
get_results(fit)
#>    profile  indicator       mean  variance standard_deviation
#> 1        1     browse  0.5395774 0.5898715          0.7680309
#> 2        1   lectures  0.4274031 0.8446923          0.9190714
#> 3        1 forum_read  0.6198249 0.4814743          0.6938835
#> 4        1 forum_post  0.5304888 0.6887895          0.8299335
#> 5        1 attendance  0.6490498 0.3841263          0.6197793
#> 6        2     browse -0.7655998 0.5284684          0.7269583
#> 7        2   lectures -0.6064852 0.5383733          0.7337392
#> 8        2 forum_read -0.8792137 0.3631130          0.6025886
#> 9        2 forum_post -0.7526883 0.4211624          0.6489703
#> 10       2 attendance -0.9206599 0.3736293          0.6112523
#>    mean_standard_error variance_standard_error
#> 1           0.02710160              0.02948401
#> 2           0.03248649              0.04210054
#> 3           0.02473030              0.02439626
#> 4           0.02953027              0.03481645
#> 5           0.02226482              0.01955858
#> 6           0.03121144              0.03221700
#> 7           0.03093059              0.03209025
#> 8           0.02621510              0.02279295
#> 9           0.02768122              0.02607600
#> 10          0.02642702              0.02278996
get_results(fit, what = "profile_probabilities")
#>   group_class profile probability group_class_probability
#> 1           1       1   0.1704508                0.325297
#> 2           1       2   0.8295492                0.325297
#> 3           2       1   0.7859317                0.674703
#> 4           2       2   0.2140683                0.674703
get_results(fit, what = "entropy")
#>         level n_classes n_units entropy_sum relative_entropy
#> 1 individuals         2    1422   60.939978        0.9381731
#> 2      groups         2     106    4.434661        0.9396428
names(get_results(fit, what = "all"))
#>  [1] "profiles"              "responses"             "covariances"          
#>  [4] "gaussian_diagnostics"  "profile_probabilities" "counts"               
#>  [7] "posteriors"            "group_posteriors"      "assignments"          
#> [10] "classification"        "average_posteriors"    "classification_errors"
#> [13] "bch_weights"           "entropy"               "replication"          
#> [16] "residuals"             "information_criteria"  "model"                
#> [19] "stages"                "starts"                "data"                 

# Compare estimated assignments with the generating classes at both levels.
get_results(fit, what = "assignments", data = course_engagement,
         truth = c("engagement", "student_type"))
#>    assignment class        truth      value   n proportion
#> 1     profile     1   engagement disengaged  16 0.02735043
#> 2     profile     2   engagement disengaged 569 0.97264957
#> 3     profile     1   engagement    engaged 820 0.97968937
#> 4     profile     2   engagement    engaged  17 0.02031063
#> 5 group_class     1 student_type  committed   1 0.01612903
#> 6 group_class     2 student_type  committed  61 0.98387097
#> 7 group_class     1 student_type   wavering  34 0.77272727
#> 8 group_class     2 student_type   wavering  10 0.22727273
```
