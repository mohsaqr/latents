# Get any table of a fitted model, grid, comparison, test or diagnostics result as a data frame

`get_results()` returns the results of an object of this package as data
frames. Every fitted model, enumeration grid, model comparison,
bootstrap test, starting-value set and diagnostics result holds a set of
named tables, and `what` chooses among them. The sections below list
every table name for each kind of object, with what one row is, its
columns and the arguments it takes. Many tables also have a verb of
their own that takes any object, such as
[`get_entropy()`](https://pak.dynasite.org/latents/reference/get_entropy.md),
[`get_best()`](https://pak.dynasite.org/latents/reference/get_best.md)
or
[`get_tree()`](https://pak.dynasite.org/latents/reference/get_tree.md).

## Usage

``` r
get_results(x, what = NULL, ...)
```

## Arguments

- x:

  A fitted model, an enumeration grid, a model comparison, a bootstrap
  test, a starting-value set, a diagnostics result or pooled
  imputations; or the [`summary()`](https://rdrr.io/r/base/summary.html)
  of a fitted model, of a profile or class grid or of a bootstrap test,
  which holds the same tables.

- what:

  Which table to return:

  - `NULL` (default): the object's primary table, the same table as
    `as.data.frame(x)` and the first one listed for the object below. A
    profile or class fit without continuous indicators returns the first
    of `"responses"`, `"ordinal"` and `"count_means"` that it has; a
    general transition fit (class `multilpa_lta`) returns
    `"transitions"`;

  - one table name, matched in full and not by its first letters;

  - `"all"`: every table the object builds by default, as a named list.

- ...:

  Arguments of the requested table, named in its entry below and
  described there or under *Table arguments*.

## Value

For one table, a data frame: the table's entry below gives what one row
is and the columns. For `what = "all"`, a named list of data frames in
the order of the entries, each built with its default arguments. The
list leaves out the tables this object cannot produce (the sequence
tables of a fit made without `time`, for example) and those built only
when asked for (`"lmr"`, `"recovery"`, and `"group_posteriors"` of a fit
with fewer than two group classes).

Conditions raised:

- `latents_bad_argument`: `what` is not a table of the object (the
  message lists the names it has), `...` holds an argument the table
  does not take, or `x` is not an object of this package;

- `latents_unused_argument` (a warning): such an argument, given to a
  table of a cross-level fit or of a transition or family grid;

- `latents_no_time`: a sequence table of a fit made without `time`;

- `latents_no_local_dependence`, `latents_no_direct_effects`: the table
  of those terms, for a fit without them;

- `latents_unsupported_noise`: the classification errors, BCH weights or
  residuals of a fit with a noise component.

## Tables of a profile, class or transition fit

The tables of a fit from
[`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
[`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
or
[`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md),
with or without membership covariates, and of an
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit of
class `multilpa_transitions` (the basic transition model). Any other
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit has
class `multilpa_lta` and the tables listed under *Tables of the other
model families*.

- `"profiles"`:

  The Gaussian measurement model, one row per profile and continuous
  indicator: `profile`, `indicator`, `mean`, `variance`,
  `standard_deviation`, `mean_standard_error`,
  `variance_standard_error`. Takes `scale` and `data`.

- `"responses"`:

  The categorical measurement model, one row per profile, indicator and
  category: `profile`, `indicator`, `category`, `probability`,
  `threshold` (`qlogis(P(y <= category))`, `NA` for the last category),
  `probability_standard_error`. Takes `data`.

- `"covariances"`:

  One row per profile and ordered pair of continuous indicators:
  `profile`, `indicator`, `indicator_2`, `covariance` (exact zeros off
  the diagonal under a diagonal structure).

- `"ordinal"`, `"count_means"`:

  Only for fits with ordinal or count indicators. `"ordinal"` has one
  row per profile, indicator and category: `profile`, `indicator`,
  `category`, `probability`, `location` (the profile's location on the
  indicator's latent scale) and `location_standard_error`.
  `"count_means"` has one row per profile and indicator: `profile`,
  `indicator`, `mean`, `standard_error`. Take `data`.

- `"local_dependence"`, `"direct_effects"`:

  For fits with those arguments. Local dependence has `pair`, `kind`,
  `profile`, `item_1`, `category_1`, `item_2`, `category_2`,
  `probability`, `independent`, `ratio`, `covariance`, `correlation`: a
  categorical pair has one row per profile and category combination with
  the joint `probability`, the product of the marginals (`independent`)
  and their `ratio`; a continuous pair one row per profile with
  `covariance` and `correlation`. Direct effects have one row per
  effect: `indicator`, `covariate`, `estimate`, `standard_error`,
  `statistic`, `p_value`, `conf_low`, `conf_high`.

- `"counts"`:

  Effective class sizes (summed posteriors), one row per level and
  class: `level`, `class`, `effective_count`, `effective_proportion`.

- `"model"`, `"information_criteria"`:

  One row describing the fit: dimensions, likelihood, parameter counts,
  criteria, convergence and boundary status; `"information_criteria"`
  has the same columns for every family and takes `format` and
  `definitions` (see *Information criteria*).

- `"posteriors"`:

  One row per observation and profile: `row`, `group`, `profile`,
  `posterior`, `modal`. Takes `format`. Excluded observations keep their
  rows with `NA`.

- `"group_posteriors"`:

  One row per group and group class: `group`, `group_size`,
  `log_likelihood`, `group_class`, `posterior`, `modal`. Takes `format`.

- `"assignments"`:

  One row per observation: the columns the fit was built from (or those
  of `data`), `profile`, `group_class` where there is one, `uncertainty`
  (one minus the posterior of the assigned profile), one
  `posterior_profile_*` column per profile, and `posterior_noise` for a
  noise fit (profile 0). Takes `data` and `truth`; with `truth` it
  returns a recovery table (*Table arguments*). Rows with no observed
  indicator have `NA` assignments under `drop_empty = TRUE`, but keep a
  group class supported by other rows.

- `"classification"`:

  One row per level and class: `level`, `class`, `n_modal`,
  `proportion_modal`, `estimated_n`, `estimated_proportion`,
  `average_posterior`, `odds_correct_classification` (see
  *Classification quality*). Takes `level`.

- `"average_posteriors"`:

  For the units assigned to each class, their mean posterior for every
  class: `level`, `assigned_class`, `class`, `n_assigned`,
  `average_posterior`; rows sum to one within an assigned class. Takes
  `level`.

- `"classification_errors"`, `"bch_weights"`:

  The probability that a unit of one true class is assigned to another
  (`level`, `true_class`, `assigned_class`, `probability`), and the BCH
  weights that invert that matrix (`level`, `unit`, `assigned_class`,
  `class`, `weight`). Take `level`. Marginal over the units on a
  covariate fit.

- `"entropy"`:

  One row per level: `level`, `n_classes`, `n_units`, `entropy_sum`,
  `relative_entropy`.

- `"replication"`:

  One row: the numbers of starts run, screened, converged and failed,
  and `n_best_replicated`, the starts within `tolerance` of the best log
  likelihood. See
  [`get_replication()`](https://pak.dynasite.org/latents/reference/get_replication.md).

- `"starts"`:

  One row per start: `start`, `log_likelihood`, `converged`,
  `iterations`, `error`, `boundary` where recorded, and under
  `search = "emEM"` `kind`, `stage` (`"final"` or `"screened"`) and
  `screened_log_likelihood`.

- `"residuals"`:

  Bivariate residuals of every indicator pair: `profile`, `indicator_1`,
  `indicator_2`, `kind`, `observed`, `expected`, `residual`,
  `effective_n`, `statistic`, `df`, `p_value`, `p_adjusted`. One row per
  pair and profile, or per pair under `by = "overall"`. Takes `data`,
  `by` and `adjust`. See *Bivariate residuals*.

- `"gaussian_diagnostics"`:

  Only for fits with continuous indicators. One row per Gaussian
  component: `context`, `profile`, `gaussian_mass`, `positive_rows`,
  `weight_effective_n`, `min_variance_ratio`,
  `min_standardized_eigenvalue`, `boundary`: how much posterior mass a
  component carries, how concentrated it is, and how small its variances
  are relative to the data. Descriptive, with no cutoff and no effect on
  start selection.

- `"sequences"`, `"sequence_summary"`:

  Fits with `time`: the assignments in occasion order (`group`,
  `group_class`, `time`, `profile`; takes `format`) and one row per
  group class (`group_class`, `groups`, `observations`, `mean_length`,
  `median_length`, `shortest`, `longest`, `complete`, `gaps`).

- `"data"`:

  The data the model was fitted to, in input order: the identifier, the
  occasion, the rescaled weights and the indicators under their original
  names; not the covariates, so
  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
  and [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md)
  take `data`.

- `"profile_probabilities"`, `"stages"`:

  A fit without membership covariates or transitions: one row per group
  class and profile (`group_class`, `profile`, `probability`,
  `group_class_probability`), and one row per estimation stage (`stage`,
  `group_classes`, `fixed`, `log_likelihood`, `parameters`,
  `parameters_with_measurement`, `converged`; two rows for a
  [`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
  result).

- `"coefficients"`:

  A fit with membership covariates (`profile_covariates` or
  `group_covariates`): one row per membership coefficient (`level`
  `"profile"` or `"group"`, `outcome`, `term`, `parameter`, always
  `"logit"`, `estimate`); standard errors come from
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md).

- `"transitions"`, `"initial"`, `"sequence_lengths"`:

  An [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit:
  one row per group class and ordered pair of profiles (`group_class`,
  `from`, `to`, `probability`, `expected_count`,
  `origin_expected_count`, `supported`, `stable`, `estimated`,
  `group_class_probability`; takes `estimated` and `stable`); one row
  per group class and profile (`group_class`, `profile`, `probability`,
  `prevalence`, `group_class_probability`); one row per group (`group`,
  `group_class`, `occasions`, `observations`, `complete`).

Under sampling weights every population quantity, here and in the grid
tables, is weighted: effective and estimated counts and proportions,
average posteriors, the odds of correct classification, entropy,
classification errors, BCH weights, the criteria and the tree's `mass`.
Modal counts stay unweighted counts of the units of positive weight; a
unit of weight zero keeps its row and posteriors but enters no
statistic.

## Tables of an enumeration grid

A profile or class grid from
[`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md),
[`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md)
or
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
has the tables below.
[`get_candidates()`](https://pak.dynasite.org/latents/reference/get_candidates.md),
[`get_best()`](https://pak.dynasite.org/latents/reference/get_best.md),
[`get_deltas()`](https://pak.dynasite.org/latents/reference/get_deltas.md),
[`get_knees()`](https://pak.dynasite.org/latents/reference/get_knees.md)
and
[`get_tree()`](https://pak.dynasite.org/latents/reference/get_tree.md)
return them too, and their pages describe the columns one by one.

- `"candidates"`:

  One row per candidate, failed ones included: `candidate`, its class
  counts and `model`, log likelihood, parameter count, every criterion
  under both sample-size conventions, both entropies, the AvePP
  summaries of each level, convergence, boundary, replication, the
  Gaussian diagnostics, `warnings` and `error`. Takes `cutoff` (default
  0.8) for `avepp_below`; a class with no assigned unit counts as below
  it.
  [`get_candidates()`](https://pak.dynasite.org/latents/reference/get_candidates.md)
  drops the columns empty in every row.

- `"best"`:

  One row per criterion and sample-size convention, with the converged
  candidate that has the lowest value. `"criteria"` is its older name,
  still in `"all"`.

- `"classification"`:

  One row per fitted candidate, level and class: `candidate`, its class
  counts and `model`, the columns of a fit's `"classification"` table
  and the level's `relative_entropy`.

- `"deltas"`, `"knees"`:

  Model-selection evidence, taking `using`, one criterion of the grid
  (default `"bic_groups"`); neither selects a model. `"deltas"` has one
  row per candidate, with the decrease in the criterion from the
  converged candidate with one profile fewer (`delta_profile`) and with
  one group class fewer (`delta_group_class`): positive when this
  candidate is lower, `NA` where the neighbour is missing or
  unconverged. `"knees"` has one row per line of the grid (along the
  profiles, the group classes or the parameters), with the candidate
  furthest below the chord of the line (`knee_candidate`), a simplified
  Kneedle rule; a line needs three candidates for a knee and four for
  `stable`.

- `"tree"`:

  The flows of
  [`plot_tree()`](https://pak.dynasite.org/latents/reference/plot_tree.md),
  one row per pair of classes in consecutive rows of a tree, with
  `mass`, the share of all units the two classes share (summed over the
  products of their posteriors). Takes the arguments of
  [`plot_tree()`](https://pak.dynasite.org/latents/reference/plot_tree.md):
  `n_group_classes`, `n_profiles` (`n_classes` on a regression grid),
  `model`, `candidates`, `from`, `to` and `level` (default
  `"individuals"`).

- `"lmr"`:

  The Vuong-Lo-Mendell-Rubin and adjusted Lo-Mendell-Rubin tests, one
  row per pair of converged candidates differing by one profile or one
  group class at the same structure: `direction`, `null_candidate`,
  `alternative_candidate`, the alternative's class counts and `model`,
  `n_units`, `parameters_added`, `lr`, `vlmr_p`, `adjusted_lr`
  (`lr / (1 + 1 / (parameters_added * log(n_units)))`), `almr_p`,
  `reference_mean`, `reference_sd` and `note`. The reference is a
  weighted chi-square with the eigenvalues of Vuong's matrix as weights,
  evaluated by Imhof's method; a pair with a singular information matrix
  or a boundary candidate keeps its row with `NA` and the reason in
  `note`. Built only when asked for.

A transition grid (`enumerate_classes(time = )`) or a family grid
(`enumerate_classes(family = )`) has three tables:

- `"candidates"`:

  One row per candidate, with `model`, `warnings`, `error` and the
  columns
  [`get_candidates()`](https://pak.dynasite.org/latents/reference/get_candidates.md)
  lists for the grid.

- `"best"`:

  One row: the row of `"candidates"` with the lowest value of one
  criterion among the converged candidates. Takes `criterion`: `"bic"`
  (default), `"aic"` or `"bic_individual"`.

- `"tree"`:

  As for a profile or class grid. Not in `"all"`.

A regression grid from
[`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md)
has:

- `"fit"`:

  One row per candidate, with the columns
  [`get_candidates()`](https://pak.dynasite.org/latents/reference/get_candidates.md)
  lists for a regression grid and `note`; `best_bic` is `TRUE` for the
  converged candidate with the lowest BIC.

- `"tree"`:

  As for a profile or class grid, with `n_classes` in place of
  `n_profiles`.

- `"model"`:

  Not a table: the fitted model of one candidate, chosen with
  `n_classes` and `n_group_classes` (default 1). Not in `"all"`.

The fitted models of the other grids and of a comparison are extracted
with
[`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md).

## Tables of the other objects

- A
  [`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)
  result:

  `"comparison"`: one row per model, with the columns
  [`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)
  lists. `"classification"`: one row per model, level and class, with
  the columns of a fit's `"classification"` table and
  `relative_entropy`; only when every model compared is a profile,
  class, covariate or transition model.

- A
  [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
  result:

  `"test"`: one row, the test. `"replicates"`: one row per simulated
  data set.

- A starting-value set
  ([`starting_values()`](https://pak.dynasite.org/latents/reference/starting_values.md)):

  `"profiles"`, `"responses"`, `"profile_probabilities"` and
  `"covariances"`: the rows of a fit's tables of those names, with the
  estimates only (no standard errors or thresholds).

- A
  [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  result:

  `"entropy"`, `"classification"`, `"average_posteriors"`, `"residuals"`
  and `"replication"`, each as the table of that name of a fit, and
  `"checks"`, one row per check (columns in
  [`get_checks()`](https://pak.dynasite.org/latents/reference/get_checks.md)).
  `"residuals"` exists only for a profile, class or
  `multilpa_transitions` fit.

## Tables of the other model families

Each family has its own tables; the first is the primary one. Tables
with intervals take `level` (the confidence level, default 0.95) and
`vcov_type` (`"observed"`, `"robust"` or `"opg"`; `NULL` reuses the type
stored with a regression fit).

A mixture regression
([`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)):

- `"coefficients"`:

  One row per class and term: `class`, `term`, `estimate`, `std_error`,
  `statistic`, `p_value`, `conf_low`, `conf_high`, `p_adjusted`
  (Benjamini-Hochberg over the non-intercepts) and, where the link makes
  it an odds or rate ratio, `exp_estimate` (`NA` for intercepts).

- `"classes"`:

  One row per class: `class`, `share`, `count`, `n_assigned`,
  `mean_posterior`, `prior`, `prior_std_error`, and for a Gaussian model
  the residual standard deviation `sigma` with `sigma_std_error`.

- `"membership"`:

  One row per coefficient of the class membership model: the logit
  `estimate` against a `reference` class, with its standard error, test,
  interval and `odds_ratio`.

- `"group_classes"`:

  Two-level fits (empty otherwise): one row per group class and class,
  with the `probability` of the class within the group class and
  `group_share`.

- `"fit"`:

  One row: the family, nesting, class counts, numbers of observations
  and groups, log likelihood, number of parameters, `aic`, `bic` (on the
  rows under observation nesting, on the groups under group or two-level
  nesting), `bic_rows`, `sabic`, `icl`, `entropy`, `group_entropy`,
  `smallest_share`, convergence and replication.

- `"assignments"`:

  One row per observation: `row`, the identifier, `class`, `posterior`
  and one `probability_class_*` column per class.

- `"groups"`:

  Fits with `id` (empty otherwise): one row per group, with its number
  of rows, its class or group class and the posterior probabilities.

- `"fitted"`:

  One row per observation: `row`, `observed`, `fitted`, `fitted_modal`,
  `residual`, `class`. An ordinal outcome is scored `1,...,C`.

- `"starts"`:

  One row per start: `start`, `kind`, `log_likelihood`, `converged`,
  `iterations`, `stage` (`"screened"` or `"completed"`), `degenerate`,
  `selected`.

- `"classification"`:

  For the units assigned to each class, their mean posterior of every
  class: `assigned`, `class`, `mean_posterior`, `n_assigned`.

- `"recovery"`:

  The assigned classes against known ones: one row per assigned class
  and truth value, with `n` and `share`. Needs `data` and `truth` (one
  column of `data`), and takes `by = "class"` (default) or
  `"group_class"`. Not in `"all"`.

Identifier columns that collide with output names use an `id_` or
`cluster_` prefix, with a unique suffix if needed.

A growth mixture
([`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
with `random`):

- `"coefficients"`:

  One row per class and term of the growth curve, with standard error,
  test and interval.

- `"classes"`:

  One row per class: `share` with its standard error and interval,
  `expected_persons`, `assigned_persons`, `average_posterior`,
  `residual_sd`.

- `"random"`:

  Variances, standard deviations and correlations of the random effects,
  one row per class and parameter, with intervals on the log and Fisher
  scales.

- `"membership"`:

  The coefficients of the class membership model.

- `"fit"`:

  One row: class count, persons, observations, the random-effect
  structure, number of parameters, log likelihood, `aic`, `bic`,
  `bic_rows`, `sabic`, `icl`, `entropy`, `smallest_share`, convergence
  and replication. In a multilevel growth mixture `bic` is on the
  clusters and `bic_persons` on the persons.

- `"assignments"`:

  One row per person: the identifier, `class`, `probability` and one
  `posterior_*` column per class.

- `"random_effects"`:

  One row per person: the identifier, `class` and the predicted random
  effects.

- `"trajectories"`, `"individual"`:

  Each class's mean curve with its band (`class`, the time variable,
  `estimate`, `std_error`, `conf_low`, `conf_high`), and each person's
  curve (the identifier, the time variable, `class`, `observed`,
  `class_mean`, `predicted`). Both take `time`, the name of the numeric
  time variable.

- `"starts"`:

  One row per start, as for a mixture regression.

- `"group_classes"`, `"clusters"`:

  Only in a multilevel growth mixture: its group classes and its
  clusters.

- `"recovery"`:

  As for a mixture regression, without `by`. Not in `"all"`.

An additive fit
([`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
with `family = "additive"`, `"dispersion"` or `"additive_dispersion"`):

- `"parameters"`:

  Class means, between- and within-group variances and class weights:
  `level`, `group_class`, `indicator`, `parameter`, `estimate`,
  `standard_error`, `statistic`, `p_value`, `p_adjusted`, `conf_low`,
  `conf_high`. Variance intervals are formed on the log scale and weight
  intervals on the logit scale. Takes `adjust` (`"none"` by default).

- `"group_classes"`:

  One row per group class: `group_class`, `weight`, `count`,
  `n_assigned`, `mean_posterior` and `effective_groups`, the number of
  perfectly classified groups that would estimate the weight as
  precisely (below 50, intervals warn).

- `"groups"`:

  One row per group: the identifier, `n`, `group_class`, `posterior` and
  one `probability_group_class_*` column per group class.

- `"intercepts"`, `"conditional_intercepts"`:

  The predicted intercept of each group on each indicator
  (`intercept_mean`, `intercept_variance`): one row per group and
  indicator, and one row per group, indicator and group class.

- `"assignments"`:

  One row per observation, with the group class of its group.

- `"classification"`:

  For the groups assigned to each group class, their mean posterior of
  every group class.

- `"fit"`, `"starts"`:

  One row describing the fit; one row per start.

- `"recovery"`:

  The assigned group classes against known ones. Needs `data` and
  `truth`. Not in `"all"`.

A cross-level fit (`multilpa(family = "restricted_cross_level")` or
`"full_cross_level"`):

- `"profiles"`:

  The profiles, then the group classes' means of the group means:
  `level`, `profile`, `indicator`, `mean`, `variance`, `prevalence`.

- `"group_classes"`:

  One row per group class: `group_class`, `weight`, `count`,
  `n_assigned`.

- `"composition"`:

  The share of each profile within each group class (`group_class`,
  `profile`, `share`): model-implied for the full family, a
  posterior-weighted description for the restricted one.

- `"groups"`, `"assignments"`:

  One row per group, and one row per observation, with their classes and
  posterior probabilities.

- `"fit"`, `"starts"`, `"gaussian_diagnostics"`:

  One row describing the fit; the record of the starts (for the
  restricted family, one row per part of the search); one row per
  Gaussian component.

The criteria of a cross-level fit come from a working likelihood and do
not choose the number of group classes or the form; see
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).

A general transition fit, the
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit of
class `multilpa_lta` that the extensions of the basic transition model
listed in [`lta()`](https://pak.dynasite.org/latents/reference/lta.md)
return:

- `"transitions"`, `"initial"`:

  The transition and initial probabilities the model implies, averaged
  over the groups and the occasions they apply to. `"transitions"` has
  one row per group class, occasion (`NA` when the transitions do not
  differ by occasion) and ordered pair of profiles: `group_class`,
  `occasion`, `from`, `to`, `probability`, `origin_expected_count`,
  `supported`, `estimated`. `"initial"` has `group_class`, `profile`,
  `probability`.

- `"transition_coefficients"`, `"initial_coefficients"`:

  The same parameters as log odds, one row per coefficient, with
  standard error, test, interval and `odds_ratio`.

- `"second_order_transitions"`, `"second_order_coefficients"`:

  Only for `order = 2`.

- `"profiles"`, `"covariances"`, `"responses"`, `"ordinal"`,
  `"count_means"`:

  The measurement model, each with an `occasion` column: profile means
  and variances, covariances (under `covariance_model = "full"`), and
  the categorical, ordinal and count models of fits with such
  indicators.

- `"group_classes"`:

  One row per group class: `group_class`, `weight`, `count`.

- `"assignments"`, `"posteriors"`, `"entropy"`,
  `"information_criteria"`:

  As for the other transition fits, with the same arguments.

- `"fit"`, `"starts"`, `"gaussian_diagnostics"`:

  One row describing the fit; one row per start; and, for continuous
  indicators, one row per Gaussian component.

`profile` is `"profile_1"`, ... in `"profiles"`, `"responses"` and the
transition tables of a general transition fit, and the profile's number
in its other tables.

Pooled imputations
([`pool_imputations()`](https://pak.dynasite.org/latents/reference/pool_imputations.md)):

- `"estimates"`:

  One row per pooled parameter: `level`, `outcome`, `term`, `parameter`,
  `estimate`, `standard_error`, `statistic`, `p_value`, `p_adjusted`,
  `conf_low`, `conf_high`, `df`, and the `within` and `between`
  variances with `riv` and `fmi`.

- `"imputations"`:

  Each imputation's aligned estimates, one row per imputation and
  parameter.

- `"fits"`:

  One row per imputation: `imputation`, `log_likelihood`,
  `n_parameters`, `converged` and the label orders used to align it.

## Table arguments

- `data`:

  One row per observation of the fit, in the fitted order; `NULL` uses
  the columns the fit carries. The standard errors of `"profiles"` and
  `"responses"` come from
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  on these data; a fit that cannot supply them gives `NA` with a message
  when `data` is omitted, and raises the condition when it is given.

- `scale`:

  For `"profiles"`: `"raw"` (default) or `"standardized"`, each
  indicator's deviation from its grand mean in units of the observed
  population standard deviation (a weighted mean and spread, without
  zero-weight rows, for a weighted fit).

- `format`:

  For the posterior tables and `"sequences"`, `"long"` (default, one row
  per unit and class) or `"wide"` (one row per unit). For
  `"information_criteria"`, `"wide"` (default, one row) or `"long"` (one
  row per criterion and convention).

- `level`:

  For the classification tables of a profile, class or transition fit:
  `"individuals"`, `"groups"` or `"both"`; `NULL` (default) is `"both"`
  for a fit with two or more group classes and `"individuals"`
  otherwise. For the other model families `level` is the confidence
  level of the intervals (default 0.95).

- `truth`:

  For `"assignments"`, column names of `data` holding known labels. The
  table then cross-tabulates the model's labels against them, one row
  per class and truth value: `assignment` (the label compared:
  `group_class` for a truth constant within every group of a fit with
  group classes, `profile` otherwise), `class`, `truth`, `value`, `n`
  and `proportion` (the share of units with that truth value assigned to
  that class). Excluded units appear under `class = NA`. Labels are
  arbitrary, so a truth value assigned mostly to a class with another
  label is relabelling, not misclassification.

- `by`, `adjust`:

  For `"residuals"`: `by` is `"profile"` (default, each profile
  separately) or `"overall"` (pooled over the profiles), and `adjust` is
  the multiplicity correction across the pairs: `"BH"` (default),
  `"none"` or another
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) method.

- `definitions`:

  For `"information_criteria"` in long format, `TRUE` adds each
  criterion's formula and reference; `FALSE` (default) leaves them out.

- `estimated`, `stable`:

  For `"transitions"` of a `multilpa_transitions` fit; `NULL` (default)
  keeps every row. `estimated = TRUE` keeps the rows estimated from
  observed occupancy and `FALSE` the others; `stable = TRUE` keeps the
  diagonal cells (staying) and `FALSE` the off-diagonal ones (moving).

`cutoff`, `using`, `criterion` and the tree arguments belong to the grid
tables and are described with them.

## Information criteria

With `q` free parameters, `n` the sample size and `EN` the entropy of
the classification: `deviance = -2L`, `aic = -2L + 2q`,
`bic = -2L + q log(n)`, `sabic = -2L + q log((n + 2) / 24)`,
`caic = -2L + q (log(n) + 1)`, `awe = -2(L - EN) + 2q (1.5 + log(n))`,
`icl = -2L + q log(n) + 2 EN`, `kic = -2L + 3(q + 1)` (as tidyLPA
computes it) and `clc = -2L + 2 EN` (on the entropy sum, so not
comparable with tidyLPA's CLC on relative entropy). `EN` is the soft
entropy of the complete-data classification, as in ICL, AWE and CLC; in
a two-level model it is the group-class entropy plus each row's profile
entropy within its group class.

Every criterion with a sample size is reported under two conventions,
`"groups"` (`n` the number of groups) and `"individuals"` (`n` the
number of rows), which coincide in a single-level model; use one
consistently across candidates. A transition fit keeps one entropy per
level instead, since its profiles form a chain within a person. Rows
with no observed indicator, groups made only of such rows, and units of
weight zero are in no sample size and no entropy.

## Classification quality

The odds of correct classification of a class is
`(p / (1 - p)) / (r / (1 - r))`, with `p` the average posterior in the
assigned class and `r` the estimated class proportion: near one, the
classification is no better than the class proportions; five is a
conventional benchmark, not a test. It is `NA` when `p` is missing or
zero or `r` is zero or one (so always for a single class), and `Inf`
when `p` is one at an interior `r`. `relative_entropy` is
`1 - entropy_sum / (n log K)`, with `entropy_sum` the level's own
entropy, between zero and one, higher meaning sharper classification;
`NA` for a level with one class. A fit without group classes reports
individuals only.

## Bivariate residuals

`"residuals"` checks, pair by pair, that the indicators are independent
within a profile, weighting each observation by its posterior. Two
Gaussian indicators (`kind = "gaussian"`): `observed` is the
posterior-weighted residual correlation, `expected` the correlation the
fitted covariance implies, and `statistic` the difference of their
Fisher transforms times `sqrt(effective_n - 3)`, referred to the normal
distribution. Two categorical or ordinal indicators (`"categorical"`):
`statistic` is the Pearson chi-square of the posterior-weighted two-way
table against the model's, Latent GOLD's bivariate residual times its
`df = (R - 1)(C - 1)`, with `observed` and `expected` the two tables'
Cramer's V; the reference ignores that the model's table is estimated
from the same data and can be conservative. A mixed pair or a pair with
a count indicator (`"mixed"`, `"count"`): `observed` is the
posterior-weighted correlation of the values (a categorical indicator
entering as its category number), `expected` zero, and no test.
`p_adjusted` applies `adjust` across the pairs with a p-value. A
weighted fit reports the descriptive columns and no tests.

## See also

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) for the
primary table, [`summary()`](https://rdrr.io/r/base/summary.html) for
every table,
[`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
for the classification tables with their checks, and the verbs
[`get_checks()`](https://pak.dynasite.org/latents/reference/get_checks.md),
[`get_entropy()`](https://pak.dynasite.org/latents/reference/get_entropy.md),
[`get_classification()`](https://pak.dynasite.org/latents/reference/get_classification.md),
[`get_avepp()`](https://pak.dynasite.org/latents/reference/get_avepp.md),
[`get_residuals()`](https://pak.dynasite.org/latents/reference/get_residuals.md),
[`get_replication()`](https://pak.dynasite.org/latents/reference/get_replication.md),
[`get_candidates()`](https://pak.dynasite.org/latents/reference/get_candidates.md),
[`get_best()`](https://pak.dynasite.org/latents/reference/get_best.md),
[`get_deltas()`](https://pak.dynasite.org/latents/reference/get_deltas.md),
[`get_knees()`](https://pak.dynasite.org/latents/reference/get_knees.md)
and
[`get_tree()`](https://pak.dynasite.org/latents/reference/get_tree.md).

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
