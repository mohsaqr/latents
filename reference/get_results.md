# Extract model estimates and analysis results

`get_results()` returns one table of a fitted model, an enumeration
grid, a comparison, a bootstrap test, a starting-value set or a
diagnostics result as a data frame, chosen by name with `what`. Without
`what` it returns the object's primary table (the same as
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)), and
with `what = "all"` a named list of the tables the object builds by
default. The catalogue below gives each table's rows, columns and
arguments. Many tables also have a verb of their own that takes any
object, such as
[`get_entropy()`](https://pak.dynasite.org/latents/reference/get_entropy.md),
[`get_avepp()`](https://pak.dynasite.org/latents/reference/get_avepp.md),
[`get_best()`](https://pak.dynasite.org/latents/reference/get_best.md)
or
[`get_tree()`](https://pak.dynasite.org/latents/reference/get_tree.md).

## Usage

``` r
get_results(x, what = NULL, ...)
```

## Arguments

- x:

  A fitted model of this package, an enumeration grid, a model
  comparison, a bootstrap comparison, a starting-value set, a
  diagnostics result, or the summary of any of them.

- what:

  Name of the table, or `"all"`. `NULL` (default) gives the primary
  table: the measurement model of a fit, the candidates of a grid, the
  test row of a bootstrap comparison, the entropy table of a diagnostics
  result. An unknown name is refused with the names the object offers.

- ...:

  Arguments of the requested table, named in its catalogue entry and
  described under *Table arguments*; a table refuses any other.

## Value

A data frame for one `what`, or a named list of them for `"all"`.
`"all"` builds each table with its defaults and leaves out the tables
this object cannot produce (a sequence table of a fit without `time`,
for example) and those that must be asked for (`"lmr"`, `"recovery"`),
so the names of the list say what was available. Every object also has
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html) (taking `level` and
`vcov_type` where the model has standard errors) and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).

## Tables of a profile, class or transition fit

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

  The ordinal and count measurement models of fits with such indicators.
  Take `data`.

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
  returns a recovery table (see *Recovery*). Rows with no observed
  indicator have `NA` assignments under `drop_empty = TRUE`, but keep a
  group class supported by other rows.

- `"classification"`:

  One row per level and class: `level`, `class`, `n_modal`,
  `proportion_modal`, `estimated_n`, `estimated_proportion`,
  `average_posterior`, `odds_correct_classification`. Takes `level`. See
  *Classification quality*.

- `"average_posteriors"`:

  For the units assigned to each class, their mean posterior for every
  class: `level`, `assigned_class`, `class`, `n_assigned`,
  `average_posterior`; rows sum to one within an assigned class. Takes
  `level`.

- `"classification_errors"`, `"bch_weights"`:

  The probability that a unit of one true class is assigned to another
  (`level`, `true_class`, `assigned_class`, `probability`), and the BCH
  weights that invert that matrix (`level`, `unit`, `assigned_class`,
  `class`, `weight`). Take `level`. On a covariate fit they are marginal
  over the units; not available for a noise fit.

- `"entropy"`:

  One row per level: `level`, `n_classes`, `n_units`, `entropy_sum`,
  `relative_entropy`.

- `"replication"`:

  One row: `n_starts`, `n_screened`, `n_converged`, `n_failed`,
  `n_best_replicated` (starts within `tolerance` of the best log
  likelihood), `share`, `tolerance`, `best_start`, `log_likelihood`.

- `"starts"`:

  One row per start: `start`, `log_likelihood`, `converged`,
  `iterations`, `error`, `boundary` where recorded, and under
  `search = "emEM"` `kind`, `stage` (`"final"` or `"screened"`) and
  `screened_log_likelihood`.

- `"residuals"`:

  Bivariate residuals of every indicator pair: `profile`, `indicator_1`,
  `indicator_2`, `kind`, `observed`, `expected`, `residual`,
  `effective_n`, `statistic`, `df`, `p_value`, `p_adjusted`. Takes
  `data`, `by` and `adjust`. See *Bivariate residuals*.

- `"gaussian_diagnostics"`:

  One row per Gaussian component: `context`, `profile`, `gaussian_mass`,
  `positive_rows`, `weight_effective_n`, `min_variance_ratio`,
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

  A `multilpa` fit: one row per group class and profile (`group_class`,
  `profile`, `probability`, `group_class_probability`), and one row per
  estimation stage (`stage`, `group_classes`, `fixed`, `log_likelihood`,
  `parameters`, `parameters_with_measurement`, `converged`; two rows for
  a
  [`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
  result).

- `"coefficients"`:

  A covariate fit: one row per membership coefficient (`level`
  `"profile"` or `"group"`, `outcome`, `term`, `parameter`, always
  `"logit"`, `estimate`); standard errors come from
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md).

- `"transitions"`, `"initial"`, `"sequence_lengths"`:

  A transition fit: one row per group class and ordered pair of profiles
  (`group_class`, `from`, `to`, `probability`, `expected_count`,
  `stable`, `estimated`, `origin_expected_count`, `supported`,
  `group_class_probability`; takes `estimated` and `stable`); one row
  per group class and profile (`group_class`, `profile`, `probability`,
  `prevalence`, `group_class_probability`); one row per group (`group`,
  `group_class`, `occasions`, `observations`, `complete`).

Under sampling weights every population quantity, here and in the grid
tables, is weighted: effective and estimated counts and proportions,
average posteriors, the odds of correct classification, entropy,
classification errors, BCH weights, the criteria and the tree's `mass`;
modal counts stay unweighted counts of the units of positive weight. A
unit of weight zero keeps its row and posteriors in the assignment and
posterior tables but enters no statistic.

## Tables of an enumeration grid

- `"candidates"`:

  One row per candidate: `candidate`, its class counts and `model`, log
  likelihood, parameter count, every criterion under both sample-size
  conventions, both entropies, the AvePP summaries of each level
  (`avepp_min`, `avepp_max`, `avepp_mean`, `avepp_below`, and the same
  with `group_`), convergence, boundary, replication, the Gaussian
  diagnostics, `warnings` and `error`. Failed candidates keep their row.
  Takes `cutoff` (default 0.8) for `avepp_below`; a class with no
  assigned unit counts as below it.
  [`get_candidates()`](https://pak.dynasite.org/latents/reference/get_candidates.md)
  drops the columns empty in every row.

- `"best"`:

  One row per criterion and convention: `criterion`, `convention`,
  `candidate` (the converged candidate with the lowest value),
  `n_profiles`, `n_group_classes`, `model`, `value`. `"criteria"` is its
  older name, still in `"all"`.

- `"classification"`:

  One row per fitted candidate, level and class: `candidate`, its class
  counts and `model`, the columns of a fit's `"classification"` table
  and the level's `relative_entropy`.

- `"deltas"`, `"knees"`:

  Model-selection evidence, taking `using` (default `"bic_groups"`);
  neither selects a model. `"deltas"` has one row per candidate:
  `candidate`, `n_profiles`, `n_group_classes`, `model`, `criterion`,
  `value`, `delta_profile` (the decrease from the converged candidate
  with one profile fewer; positive means this one is lower) with
  `from_profile`, and `delta_group_class` with `from_group_class`; `NA`
  where the neighbour is missing or unconverged. `"knees"` has one row
  per line of the grid (`direction` `"profiles"`, `"group_classes"` or
  `"parameters"`, `line`, `criterion`, `from`, `to`, `n_points`,
  `knee_candidate`, `distance`, `knee_without_largest`, `stable`,
  `note`): the candidate furthest below the chord of the line, a
  simplified Kneedle rule; a line needs three candidates for a knee and
  four for `stable`.

- `"tree"`:

  The flows of
  [`plot_tree()`](https://pak.dynasite.org/latents/reference/plot_tree.md),
  taking its arguments `n_group_classes`, `n_profiles` (`n_classes` on a
  regression grid), `model`, `candidates`, `from`, `to` and `level`
  (default `"individuals"`). One row per pair of classes in consecutive
  rows of a tree: `tree`, `level`, `from_candidate`, `from_class`,
  `from_share`, `to_candidate`, `to_class`, `to_share`, `mass` (the
  share of all units the two classes share, summed over the products of
  their posteriors), `share_of_from`, `share_of_to`.

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

A transition or family grid has `"candidates"` (with one entropy and no
AvePP summaries), `"best"` (taking `criterion`: `"bic"`, `"aic"` or
`"bic_individual"`) and `"tree"` (not in its `"all"`). The models of a
grid or a comparison are extracted with
[`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md).

## Tables of the other objects

A
[`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)
result: `"comparison"` (its table) and `"classification"` (one row per
model, level and class, for profile, class, covariate and transition
models). A
[`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
result: `"test"` (one row) and `"replicates"` (one row per simulated
data set). A starting-value set: `"profiles"`, `"responses"`,
`"profile_probabilities"` and `"covariances"`. A
[`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
result: `"entropy"`, `"classification"`, `"average_posteriors"`,
`"residuals"`, `"replication"` and `"checks"`.

The other model families return their own tables, the first being the
default. Tables with intervals take `level` and `vcov_type`
(`"observed"`, `"robust"` or `"opg"`; `NULL` reuses the type stored with
a regression fit).

- Mixture regressions:

  `"coefficients"` (one row per class and term with Wald statistics,
  `p_adjusted` by Benjamini-Hochberg over the non-intercepts, and
  `exp_estimate` where the link makes it an odds or rate ratio;
  intercept ratios are `NA`), `"classes"`, `"membership"` (logits and
  odds ratios), `"group_classes"` (two-level fits), `"fit"` (`bic` on
  the rows under observation nesting and on groups under group or
  two-level nesting, `bic_rows` on the rows), `"assignments"`,
  `"groups"`, `"fitted"` (ordinal observations, fitted values and
  residuals use category scores `1,...,C`), `"starts"` (`stage`
  `"screened"` or `"completed"`), `"classification"` and `"recovery"`
  (`by = "class"` or `"group_class"`). Identifier columns that collide
  with output names use an `id_` or `cluster_` prefix, with a unique
  suffix if needed.

- Growth mixtures:

  `"coefficients"`, `"classes"` (share with its interval,
  `expected_persons`, `assigned_persons`, `average_posterior`,
  `residual_sd`), `"membership"`, `"fit"`, `"assignments"` (one row per
  person: `class`, `probability`, one `posterior_*` column per class),
  `"starts"`, `"recovery"`, `"random"` (variances, standard deviations
  and correlations of the random effects, with intervals on the log and
  Fisher scales), `"random_effects"` (one row per person),
  `"trajectories"` (each class's mean curve with a band) and
  `"individual"` (each person's curve); the last two take `time`. A
  multilevel growth mixture adds `"group_classes"` and `"clusters"`, and
  its `"fit"` has `bic` on the clusters and `bic_persons` on the
  persons.

- Additive fits:

  `"parameters"` (class means, between- and within-group variances and
  class weights, with variance intervals on the log scale and weight
  intervals on the logit scale; takes `adjust`, `"none"` by default),
  `"group_classes"` (with `effective_groups`, the number of perfectly
  classified groups that would estimate the class weight as precisely;
  below 50 the class's intervals warn), `"groups"`, `"intercepts"`,
  `"conditional_intercepts"`, `"assignments"`, `"classification"`,
  `"fit"`, `"recovery"` and `"starts"`.

- Cross-level fits:

  `"profiles"` (the profiles, then the group classes' means of group
  means), `"group_classes"`, `"composition"` (profile shares within each
  class: model-implied for the full family, a posterior-weighted
  description for the restricted one), `"groups"`, `"assignments"`,
  `"fit"` and `"starts"`. The likelihood is a working likelihood that
  counts the ratings twice: its criteria favour group classes and the
  full form even with no group structure, and do not choose between
  them.

- General transition fits:

  An [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit
  with occasion-specific or covariate transitions, occasion-specific
  measurement, second order, a stayer class, weights, or ordinal or
  count indicators: `"transitions"`, `"initial"` and
  `"second_order_transitions"` (model-implied, averaged over the groups
  and applicable occasions), with `origin_expected_count` and
  `supported` describing transition-row information,
  `"transition_coefficients"`, `"initial_coefficients"` and
  `"second_order_coefficients"` (log odds with odds ratios),
  `"profiles"`, `"covariances"` (under `covariance_model = "full"`),
  `"responses"` (with `occasion`), `"ordinal"`, `"count_means"`,
  `"group_classes"`, `"assignments"`, `"posteriors"` (with `format`),
  `"entropy"`, `"information_criteria"` (with `format` and
  `definitions`; one entropy per level, as for the other transition
  fits), `"fit"` and `"starts"`. `profile` is `"profile_1"`, ... in
  `"profiles"`, `"responses"` and the transition tables, and the
  profile's number in the others.

- Pooled imputations:

  `"estimates"`, `"imputations"` (each imputation's aligned estimates)
  and `"fits"` (one row per imputation).

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
  population standard deviation. Weighted fits use a weighted mean and
  spread, excluding zero-weight rows; unweighted fits use the same
  population convention.

- `format`:

  For the posterior tables and `"sequences"`, `"long"` (default, one row
  per unit and class) or `"wide"` (one row per unit). For
  `"information_criteria"`, `"wide"` (default, one row) or `"long"` (one
  row per criterion and convention).

- `level`:

  `"individuals"`, `"groups"` or `"both"`; `NULL` (default) is `"both"`
  for a fit with group classes and `"individuals"` otherwise.

- `truth`:

  For `"assignments"`, column names of `data` holding known labels. See
  *Recovery*.

- `by`, `adjust`:

  For `"residuals"`: `"profile"` (default) or `"overall"`, and the
  multiplicity correction across the pairs (`"none"` by default, or a
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) method).

- `definitions`:

  For `"information_criteria"` in long format, `TRUE` adds each
  criterion's formula and reference.

- `estimated`, `stable`:

  For `"transitions"`, `TRUE` or `FALSE` keeps the rows estimated (or
  not) from observed occupancy, or the diagonal (or off-diagonal) cells.

## Recovery

`"assignments"` with `truth` cross-tabulates the model's labels against
known ones: one row per class and truth value, with `assignment` (the
label compared: `group_class` for a truth constant within every group of
a fit with group classes, `profile` otherwise), `class`, `truth`,
`value`, `n` and `proportion` (the share of units with that truth value
assigned to that class). Excluded units appear under `class = NA`.
Labels are arbitrary, so a truth value assigned mostly to a class with
another label is relabelling, not misclassification.

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
weight zero are in no sample size and no entropy. The long format with
`definitions = TRUE` reports each formula.

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
`p_adjusted` applies `adjust` (Benjamini-Hochberg by default, `"none"`
for the raw p-values) across the pairs with a p-value. A weighted fit
reports the descriptive columns and no tests.

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
