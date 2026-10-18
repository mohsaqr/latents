# Extract the result tables of a fitted model

Returns one tidy table of a fitted model or related object, named by
`what`, as a base `data.frame`; `what = "all"` returns every table the
object offers as a named list. The tables available depend on the class
of the object and are listed under *Tables*.

## Usage

``` r
get_results(x, what = NULL, ...)

# Default S3 method
get_results(x, what = NULL, ...)

# S3 method for class 'multilpa'
get_results(x, what = NULL, ...)

# S3 method for class 'multilpa_covariates'
get_results(x, what = NULL, ...)

# S3 method for class 'multilpa_transitions'
get_results(x, what = NULL, ...)

# S3 method for class 'multilpa_enumeration'
get_results(x, what = NULL, ...)

# S3 method for class 'latents_comparison'
get_results(x, what = NULL, ...)

# S3 method for class 'multilpa_bootstrap_lrt'
get_results(x, what = NULL, ...)

# S3 method for class 'multilpa_start'
get_results(x, what = NULL, ...)

# S3 method for class 'multilpa_diagnostics'
get_results(x, what = NULL, ...)

# S3 method for class 'summary_multilpa'
get_results(x, what = NULL, ...)

# S3 method for class 'summary_multilpa_covariates'
get_results(x, what = NULL, ...)

# S3 method for class 'summary_multilpa_transitions'
get_results(x, what = NULL, ...)

# S3 method for class 'summary_multilpa_enumeration'
get_results(x, what = NULL, ...)

# S3 method for class 'summary_multilpa_bootstrap_lrt'
get_results(x, what = NULL, ...)
```

## Arguments

- x:

  A fitted model of this package, an enumeration grid, a model
  comparison, a bootstrap comparison, a starting-value set, a
  classification-diagnostics object, or the summary of any of them.

- what:

  Name of the table to return, or `"all"` for every table the object
  offers. `NULL` (default) returns the primary table: the measurement
  model for a fit, the candidate grid for an enumeration, the test row
  for a bootstrap comparison and the entropy table for a diagnostics
  object. A name the object does not offer raises
  `latents_bad_argument`, listing the names it does offer.

- ...:

  Arguments for the requested table. Each table accepts only the
  arguments listed for it under *Tables*; any other argument raises
  `latents_bad_argument`.

## Value

A base `data.frame` for a single `what`, whose columns are given under
*Tables*, or a named list of such data frames for `what = "all"`.

## Details

Every object also has [`print()`](https://rdrr.io/r/base/print.html),
which shows a compact overview,
[`summary()`](https://rdrr.io/r/base/summary.html), which collects its
main tables (taking `level` and `vcov_type` where the model has standard
errors), and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), which
returns its primary table, the same as `get_results(x)`.

## Tables

Every fitted model offers these:

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
  `effective_proportion`. These differ from the modal counts in
  `"classification"` whenever entropy is below one.

- `"model"`:

  One row describing the fit: its dimensions, likelihood, parameter
  counts, `aic`, `bic_groups`, `bic_individual` and its convergence and
  boundary diagnostics. The columns differ by family (a covariate fit
  adds its predictor counts); `"information_criteria"` has the same
  columns for every family.

- `"posteriors"`:

  One row per individual and profile: `row`, `group`, `profile`,
  `posterior`, `modal`. Takes `format`.

- `"assignments"`:

  One row per observation: the columns the fit was built from (or those
  of `data`), then `profile`, `group_class` where the model has one,
  `uncertainty` (one minus the posterior of the assigned profile) and
  one `posterior_profile_*` column per profile, plus `posterior_noise`
  for a fit with a noise component (whose cases are profile 0). Takes
  `data` and `truth`; with `truth` it returns a recovery table instead
  (see *Recovery*).

- `"classification"`:

  Classification quality, one row per level and class: `level`, `class`,
  `n_modal`, `proportion_modal`, `estimated_n`, `estimated_proportion`,
  `average_posterior`, `odds_correct_classification`. Takes `level`.
  Under sampling weights `estimated_n` and `estimated_proportion` (and
  hence the odds of correct classification) are weighted posterior
  totals; the modal and average posterior columns describe the units
  supplied.

- `"average_posteriors"`:

  For the units modally assigned to each class, their mean posterior for
  every class: `level`, `assigned_class`, `class`, `n_assigned`,
  `average_posterior`. Rows sum to one within each assigned class. This
  conditions on the assigned class, whereas `"classification_errors"`
  conditions on the true class. Takes `level`.

- `"classification_errors"`:

  The probability that a unit in one true class is assigned to another:
  `level`, `true_class`, `assigned_class`, `probability`, summing to one
  within each `true_class`. Takes `level`. On a membership-covariate fit
  this is a marginal summary over the fitted units, not conditional on
  the covariates. A fit with a noise component raises
  `latents_unsupported_noise`, as
  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
  does, and `"all"` and
  [`summary()`](https://rdrr.io/r/base/summary.html) omit the table.

- `"bch_weights"`:

  The BCH inverse-error weights, one row per unit and class: `level`,
  `unit`, `assigned_class`, `class`, `weight`. Takes `level`; a weighted
  regression should use one level, not `"both"`. On a
  membership-covariate fit the weights invert the marginal error matrix
  and are not a conditional three-step correction. Not available for a
  fit with a noise component, as for `"classification_errors"`.

- `"entropy"`:

  One row per level: `level`, `n_classes`, `n_units`, `entropy_sum`,
  `relative_entropy`.

- `"residuals"`:

  Bivariate residuals for every indicator pair, testing within-profile
  independence: `profile`, `indicator_1`, `indicator_2`, `kind`,
  `observed`, `expected`, `residual`, `effective_n`, `statistic`, `df`,
  `p_value`, `p_adjusted`. Takes `data`, `by` and `adjust`.

- `"information_criteria"`:

  Every information criterion the package computes (see *Information
  criteria*). Takes `format` and `definitions`.

- `"sequences"`:

  The assignments in occasion order, one row per observation: `group`,
  `group_class`, `time`, `profile`. Takes `format`. Raises
  `latents_no_time` for a fit made without `time`.

- `"sequence_summary"`:

  One row per group class: `group_class`, `groups`, `observations`,
  `mean_length`, `median_length`, `shortest`, `longest`, `complete`,
  `gaps`. Raises `latents_no_time` for a fit made without `time`.

- `"starts"`:

  One row per EM start: `start`, `log_likelihood`, `converged`,
  `iterations`, `error` and, where the family records it, `boundary`.

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

Tables specific to one family or object:

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
  the convergence, boundary, replication, warning and error diagnostics.
  Failed candidates are kept with `NA` estimates and their error text.
  It takes `cutoff` (default `0.8`): `avepp_below` and
  `group_avepp_below` count the classes whose AvePP is below it. A class
  with no modally assigned unit has no AvePP; it is excluded from the
  minimum, maximum and mean and counted as below the cutoff. A level
  with one class has `NA` summaries. `"best"` has one row per
  information criterion and sample-size convention: `criterion`,
  `convention`, `candidate` (the converged candidate with the lowest
  value), `n_profiles`, `n_group_classes`, `model` and `value`;
  `what = "criteria"` is accepted as a deprecated name for it.
  `"classification"` has one row per fitted candidate, level and class:
  `candidate`, `n_profiles`, `n_group_classes`, `model`, then the
  columns of a fit's `"classification"` table and that level's
  `relative_entropy`. Failed candidates and single-class levels
  contribute no rows.

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
  simplified rule inspired by Kneedle, without its smoothing or
  sensitivity thresholds, so a very small positive distance can still
  name a knee. A line needs three candidates for a knee and four for
  `knee_without_largest` and `stable`; otherwise they are `NA` and
  `note` gives the reason.

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
  `"residuals"`:

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

  The tables above plus `"random"` (random-effect variances, standard
  deviations and correlations, with intervals on the log and Fisher-z
  scales), `"random_effects"` (one row per person), `"trajectories"`
  (each class's mean trajectory with a band) and `"individual"` (each
  person's predicted curve); the last two take `time`. A multilevel
  growth mixture adds `"group_classes"` and `"clusters"`, and its
  `"fit"` uses the number of clusters in `bic` (persons in
  `bic_persons`).

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

A named list of every table the object offers, in catalogue order, each
built with its default arguments; supplying any further argument raises
`latents_bad_argument`. A table this object cannot produce (for example
a sequence table for a fit made without `time`, or bivariate residuals
for a family without discrete group classes) is omitted, so the names of
the list show what was available. Only the conditions that mean the
table does not exist for this object are skipped; any other error is
raised.

## Arguments the tables take

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
to that class). Several truth columns are stacked.

A truth column that is constant within every group is compared against
`group_class`; one that varies within any group is compared against
`profile`. The `assignment` column records which comparison was made.
Group-level recovery counts each group once. Unused factor levels are
omitted. For a fit without discrete group classes every column is
compared against `profile`.

The model's labels are arbitrary, so mass off the diagonal indicates
relabelling, not failure.

## Information criteria

Let `q` be the number of free parameters, `n` the chosen sample size,
and `EN` the classification entropy of the level matching that
convention. The criteria are `deviance = -2L`, `aic = -2L + 2q`,
`bic = -2L + q log(n)`, `sabic = -2L + q log((n + 2) / 24)`,
`caic = -2L + q (log(n) + 1)`, `awe = -2(L - EN) + 2q (1.5 + log(n))`,
`icl = -2L + q log(n) + 2 EN`, `kic = -2L + 3(q + 1)`, and
`clc = -2L + 2 EN`. `clc`, like `icl` and `awe`, uses the entropy sum;
`tidyLPA` reports a `CLC` based on relative entropy, and the two values
are not comparable.

Every criterion that depends on a sample size is reported under both
multilevel conventions, `"groups"` and `"individuals"`, using
group-level posteriors under the first and individual-level posteriors
under the second. Use one convention consistently across compared
candidates. `deviance`, `aic` and `kic` carry
`convention = NA_character_`. `clc` depends on no sample size, but its
convention selects which level's classification entropy it penalises, so
it is reported once per convention. Individual sample sizes exclude rows
with no observed indicators, and group sample sizes exclude groups in
which every row has none (`n_informative_groups`); the classification
statistics of each level exclude the same rows and groups, as in Mplus.
Where the group classification entropy is undefined, group-level `awe`,
`icl` and `clc` are `NA`.

## Classification quality

The odds of correct classification for class `k` is
`(p / (1 - p)) / (r / (1 - r))`, where `p` is the average posterior in
the assigned class and `r` is the model-estimated class proportion.
Values near one indicate classification no better than the class
proportions alone; values of five or more are conventionally taken as
adequate separation. It is `NA_real_` when `average_posterior` is
missing, or when it or `estimated_proportion` is zero or one, so a
single-class solution always reports `NA_real_`. `average_posterior` is
`NA_real_` for a class with no modal members. `n_modal` and
`estimated_n` differ whenever entropy is below one.

`relative_entropy` is `1 - EN / (n log K)`, between zero and one, higher
meaning sharper classification; `entropy_sum` is the raw `EN` used by
the entropy-penalised criteria, lower meaning sharper.
`relative_entropy` is `NA_real_` for a level with a single class.

A fit with individual profiles but no discrete group classes reports
individuals only for `level = "both"` and raises
`latents_no_group_classes` for `level = "groups"`.

## Bivariate residuals

`"residuals"` tests, pair by pair, the assumption that indicators are
independent within a profile, weighting each observation by its
posterior. `data` must be in the order the model was fitted in; the
alignment is checked. A large residual indicates local dependence
between two indicators, which otherwise tends to be absorbed by an extra
profile. The p-values are approximate, and `p_adjusted` applies the
correction named by `adjust` across the pairs. A fit with sampling
weights uses them in the descriptive associations and `effective_n`;
design-based inference is not implemented, so `statistic`, `df`,
`p_value` and `p_adjusted` are `NA` for weighted fits.

## References

Schwarz, G. (1978). Estimating the dimension of a model. *Annals of
Statistics*, 6, 461–464. Sclove, S. L. (1987). Application of
model-selection criteria to some problems in multivariate analysis.
*Psychometrika*, 52, 333–343. Bozdogan, H. (1987). Model selection and
Akaike's information criterion. *Psychometrika*, 52, 345–370. Banfield,
J. D., & Raftery, A. E. (1993). Model-based Gaussian and non-Gaussian
clustering. *Biometrics*, 49, 803–821. Biernacki, C., & Govaert, G.
(1997). Using the classification likelihood to choose the number of
clusters. *Computing Science and Statistics*, 29, 451–457. Cavanaugh, J.
E. (1999). A large-sample model selection criterion based on Kullback's
symmetric divergence. *Statistics and Probability Letters*, 42, 333–343.
Biernacki, C., Celeux, G., & Govaert, G. (2000). Assessing a mixture
model for clustering with the integrated completed likelihood. *IEEE
Transactions on Pattern Analysis and Machine Intelligence*, 22, 719–725.
Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent
structure models with categorical variables. *Political Analysis*, 12,
3–27. Nagin, D. S. (2005). *Group-Based Modeling of Development*.
Harvard University Press. Vermunt, J. K. (2010). Latent class modeling
with covariates: two improved three-step approaches. *Political
Analysis*, 18, 450–469.

## See also

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), which
returns the primary table alone,
[`summary()`](https://rdrr.io/r/base/summary.html), which prints every
table, and
[`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
for the classification tables gathered with a condensed reading.

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
#> 2        1   lectures  0.4274032 0.8446922          0.9190714
#> 3        1 forum_read  0.6198249 0.4814742          0.6938834
#> 4        1 forum_post  0.5304889 0.6887894          0.8299334
#> 5        1 attendance  0.6490499 0.3841262          0.6197792
#> 6        2     browse -0.7655996 0.5284686          0.7269584
#> 7        2   lectures -0.6064852 0.5383733          0.7337393
#> 8        2 forum_read -0.8792135 0.3631131          0.6025887
#> 9        2 forum_post -0.7526882 0.4211624          0.6489703
#> 10       2 attendance -0.9206598 0.3736294          0.6112523
#>    mean_standard_error variance_standard_error
#> 1           0.02710160              0.02948401
#> 2           0.03248649              0.04210053
#> 3           0.02473030              0.02439626
#> 4           0.02953027              0.03481644
#> 5           0.02226482              0.01955857
#> 6           0.03121145              0.03221701
#> 7           0.03093059              0.03209025
#> 8           0.02621510              0.02279296
#> 9           0.02768122              0.02607601
#> 10          0.02642702              0.02278996
get_results(fit, what = "profile_probabilities")
#>   group_class profile probability group_class_probability
#> 1           1       1   0.7859317               0.6747028
#> 2           1       2   0.2140683               0.6747028
#> 3           2       1   0.1704509               0.3252972
#> 4           2       2   0.8295491               0.3252972
get_results(fit, what = "entropy")
#>         level n_classes n_units entropy_sum relative_entropy
#> 1 individuals         2    1422   60.939997        0.9381731
#> 2      groups         2     106    4.434664        0.9396428
names(get_results(fit, what = "all"))
#>  [1] "profiles"              "responses"             "covariances"          
#>  [4] "profile_probabilities" "counts"                "posteriors"           
#>  [7] "group_posteriors"      "assignments"           "classification"       
#> [10] "average_posteriors"    "classification_errors" "bch_weights"          
#> [13] "entropy"               "residuals"             "information_criteria" 
#> [16] "model"                 "stages"                "starts"               
#> [19] "data"                 

# `engagement` and `student_type` are the kinds each row was simulated from,
# which the model never saw. Both recovery checks are one call.
get_results(fit, what = "assignments", data = course_engagement,
         truth = c("engagement", "student_type"))
#>    assignment class        truth      value   n proportion
#> 1     profile     1   engagement disengaged  16 0.02735043
#> 2     profile     2   engagement disengaged 569 0.97264957
#> 3     profile     1   engagement    engaged 820 0.97968937
#> 4     profile     2   engagement    engaged  17 0.02031063
#> 5 group_class     1 student_type  committed  61 0.98387097
#> 6 group_class     2 student_type  committed   1 0.01612903
#> 7 group_class     1 student_type   wavering  10 0.22727273
#> 8 group_class     2 student_type   wavering  34 0.77272727
```
