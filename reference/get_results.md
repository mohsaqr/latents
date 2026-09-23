# Every table a fitted model holds, from one verb

One verb returns every tidy table this package computes, named by
`what`. The tables were once eleven separate verbs and seven
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) methods
with their own `what` catalogues; a reader had to know which verb owned
which table before they could ask for it. Here the object decides what
it can offer and `what` names it, so the same call shape reads a
measurement model, a classification diagnostic and a transition matrix.

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

  A fitted model of this package, an enumeration grid, a bootstrap
  comparison, a starting-value set, a classification-diagnostics object,
  or the summary of any of them.

- what:

  Which table to return, or `"all"` for every table the object can
  produce. The tables an object offers depend on its class and are
  listed under *Tables* below; asking for one it does not have raises
  `latents_bad_argument` naming the ones it does. `NULL`, the default,
  gives the object's primary table: the measurement model for a fit, the
  candidate grid for an enumeration, the test row for a bootstrap
  comparison, and the entropy summary for a diagnostics object.

- ...:

  Arguments for the requested table. Each table accepts only the
  arguments listed for it below, and any other name raises
  `latents_bad_argument` saying which that table takes, rather than
  being dropped on the way to a table that would then answer a different
  question.

## Value

A base `data.frame` for a single `what`, whose columns are given under
*Tables*, or a named list of such data frames for `what = "all"`.

## Tables

Every fitted model offers these:

- `"profiles"`:

  The Gaussian measurement model, one row per profile and continuous
  indicator: `profile`, `indicator`, `mean`, `variance`,
  `standard_deviation`. Zero rows when every indicator is categorical.
  Takes `scale`, and `data`, which adds `mean_standard_error` and
  `variance_standard_error`.

- `"responses"`:

  The categorical measurement model, one row per profile, categorical
  indicator and category: `profile`, `indicator`, `category`,
  `probability`, `threshold`. The threshold is
  `qlogis(P(y <= category))` and is `NA_real_` for each indicator's
  final category, where the cumulative probability is one. Takes `data`,
  which adds `probability_standard_error`.

- `"covariances"`:

  The within-profile residual covariance matrices, one row per profile
  and ordered pair of continuous indicators: `profile`, `indicator`,
  `indicator_2`, `covariance`. The diagonal parameterization carries no
  covariance array but does state a covariance matrix – the variances on
  the diagonal and exact zeros off it – and it is reported explicitly,
  so one shape answers both parameterizations.

- `"counts"`:

  Effective class memberships, one row per class at each level the fit
  has: `level`, `class`, `effective_count` (the summed posteriors) and
  `effective_proportion`. These differ from the modal counts in
  `"classification"` whenever entropy is below one.

- `"model"`:

  One row describing the whole fit: its dimensions, likelihood,
  parameter counts, `aic`, `bic_groups`, `bic_individual` and its
  convergence and boundary diagnostics. The columns are those of the
  family: a covariate fit reports its predictor counts, because padding
  those onto a two-level fit would leave the common case mostly `NA`.
  Comparing fits across families is what `"information_criteria"` is
  for, and that table does have one column set for every family.

- `"posteriors"`:

  One row per individual and profile: `row`, `group`, `profile`,
  `posterior`, `modal`. Takes `format`.

- `"assignments"`:

  One row per observation, the columns the fit was built from (or those
  of `data`) followed by `profile`, `group_class` where the model has
  one, `uncertainty`, and one `posterior_profile_*` column per profile.
  `uncertainty` is one minus the posterior of the assigned profile,
  which is what modal assignment discards. Takes `data` and `truth`;
  `truth` returns a recovery cross-tabulation instead, described below.

- `"classification"`:

  Classification quality, one row per class and level: `level`, `class`,
  `n_modal`, `proportion_modal`, `estimated_n`, `estimated_proportion`,
  `average_posterior`, `odds_correct_classification`. Takes `level`.

- `"average_posteriors"`:

  For the units assigned to each class, their mean posterior of
  belonging to every class: `level`, `assigned_class`, `class`,
  `n_assigned`, `average_posterior`. Each assigned class's rows sum to
  one. Conditions on the *assigned* class, where
  `"classification_errors"` conditions on the *true* one; the two are
  different numbers, not two spellings of one table. Takes `level`.

- `"classification_errors"`:

  The probability that a unit truly in one class is assigned to another:
  `level`, `true_class`, `assigned_class`, `probability`, summing to one
  within each `true_class`. Takes `level`. On a membership-covariate fit
  this is a marginal summary over the fitted units, not a
  classification-error adjustment conditional on covariates.

- `"bch_weights"`:

  The BCH inverse-error weights, one row per unit and class: `level`,
  `unit`, `assigned_class`, `class`, `weight`. Takes `level`; weight a
  regression with one level, not with `"both"`. On a
  membership-covariate fit these invert the marginal error matrix and
  are not a conditional three-step correction.

- `"entropy"`:

  One row per level: `level`, `n_classes`, `n_units`, `entropy_sum`,
  `relative_entropy`.

- `"residuals"`:

  Bivariate residuals for every indicator pair, testing the
  within-profile independence the measurement model assumes: `profile`,
  `indicator_1`, `indicator_2`, `kind`, `observed`, `expected`,
  `residual`, `effective_n`, `statistic`, `df`, `p_value`, `p_adjusted`.
  Takes `data`, `by` and `adjust`.

- `"information_criteria"`:

  Every likelihood-penalty criterion the package computes. Takes
  `format` and `definitions`.

- `"sequences"`:

  The assignments in occasion order, one row per observation: `group`,
  `group_class`, `time`, `profile`. Takes `format`. Raises
  `latents_no_time` for a fit made without `time`.

- `"sequence_summary"`:

  One row per group class, summarizing how much data it contributes:
  `group_class`, `groups`, `observations`, `mean_length`,
  `median_length`, `shortest`, `longest`, `complete`, `gaps`. Raises
  `latents_no_time` for a fit made without `time`.

- `"starts"`:

  One row per EM start: `start`, `log_likelihood`, `converged`,
  `iterations`, `error` and, where the family records it, `boundary`.

- `"data"`:

  The columns of the original data the model was fitted to, one row per
  observation in input order: the identifier, the occasion where there
  is one, and the indicators under their original names. Columns a model
  never saw are not here, which is why
  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
  and [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md)
  still take `data`.

Every family with discrete group classes adds:

- `"group_posteriors"`:

  One row per observed group and group class: `group`, `group_size`,
  `log_likelihood`, `group_class`, `posterior`, `modal`. Takes `format`.

And each family adds its own:

- `"profile_probabilities"`, `"stages"`:

  A `multilpa` fit. `"profile_probabilities"` has one row per group
  class and profile: `group_class`, `profile`, `probability`,
  `group_class_probability`. `"stages"` has one row per estimation
  stage: `stage`, `group_classes`, `fixed`, `log_likelihood`,
  `parameters`, `parameters_with_measurement`, `converged`. An ordinary
  fit has a single `"joint"` row and a
  [`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
  result has two; `fixed` is `NA_character_` for a stage that estimated
  every block, never a sentinel such as `"none"`.

- `"coefficients"`:

  A covariate fit. One row per estimated membership coefficient: `level`
  (`"profile"` or `"group"`), `outcome`, `term`, `parameter` (always
  `"logit"`), `estimate`. The standard errors are not here;
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  reports them.

- `"transitions"`, `"initial"`, `"sequence_lengths"`:

  A transition fit. `"transitions"` has one row per group class and
  ordered pair of profiles, with `group_class`, `from`, `to`,
  `probability`, `expected_count`, `stable`, `estimated` and
  `group_class_probability`; it takes `estimated` and `stable` to
  restrict it. `"initial"` has one row per group class and profile:
  `group_class`, `profile`, `probability`, `prevalence`,
  `group_class_probability`. `"sequence_lengths"` has one row per group:
  `group`, `group_class`, `occasions`, `observations`, `complete`.

- `"candidates"`, `"criteria"`:

  An enumeration grid. `"criteria"` has one row per information
  criterion, naming the candidate that minimises it. `"candidates"` has
  one row per candidate model, with its class counts, its covariance
  `model`, log likelihood, parameter count, every criterion under both
  sample-size conventions, both entropies, and the convergence,
  boundary, replication, warning and error diagnostics. Failed
  candidates are retained with `NA` estimates and their error text.

- `"test"`, `"replicates"`:

  A bootstrap comparison. `"test"` is the one-row result; `"replicates"`
  has one row per simulated dataset.

- `"covariances"`:

  A starting-value set. One row per profile and ordered pair of
  continuous indicators: `profile`, `indicator`, `indicator_2`,
  `covariance`.

## What `"all"` returns

A named list of every table the object can produce, in catalogue order,
built from the same definitions a single `what` uses, so the two can
never disagree. It takes no further arguments: every table is built with
its defaults, and supplying anything else raises `latents_bad_argument`.
A table this particular object cannot produce is left out rather than
erroring – a sequence table for a fit made without `time`, bivariate
residuals for a family with no discrete group classes – so the names of
the list say what was available. Only the classed refusals that mean
*this object has no such table* are skipped; anything else propagates.

## Arguments the tables take

These reach the requested table through `...`, and each table takes only
the ones listed for it under *Tables* above.

- `data`:

  A data frame with one row per observation of the fit, in the order it
  was fitted in. `NULL` uses the columns the fit carries. For
  `"profiles"` and `"responses"` it triggers
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  and adds a standard error beside every estimate.

- `scale`:

  For `"profiles"`, `"raw"` reports the estimates in input units and
  `"standardized"` divides each indicator's deviation from its grand
  mean by that indicator's observed standard deviation, which is what
  `plot(x, scale = "standardized")` draws. Every standard error beside
  an estimate is divided by the same constant, so the whole row is on
  one scale.

- `format`:

  For the posterior tables and `"sequences"`, `"long"` gives one row per
  unit and class and `"wide"` one row per unit with one column per
  class. For `"information_criteria"`, `"wide"` is the one-row reporting
  shape and `"long"` gives one row per criterion and sample-size
  convention.

- `level`:

  Which level of classification to report: `"individuals"`, `"groups"`
  or `"both"`. `NULL`, the default, gives `"both"` for a fit that has
  discrete group classes and `"individuals"` for one that does not,
  which is the same rule
  [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  follows.

- `truth`:

  For `"assignments"`, one or more column names of `data` holding a
  known label to check the assignments against. See *Recovery* below.

- `by`:

  For `"residuals"`, `"profile"` assesses each profile separately and
  `"overall"` pools them.

- `adjust`:

  For `"residuals"`, the multiplicity correction applied across the
  indicator pairs.

- `definitions`:

  For `"information_criteria"` in long format, `TRUE` adds the formula
  and the reference for each criterion.

- `estimated`, `stable`:

  For `"transitions"`, `TRUE` or `FALSE` restricts the table to the rows
  that were estimated from observed occupancy, or to the diagonal.
  `NULL` keeps every row.

## Recovery

`what = "assignments"` with `truth` returns a cross-tabulation of the
model's labels against known ones instead of the per-observation table:
one row per class and truth value, with columns `assignment` (which
model label the column was compared against), `class`, `truth` (the
column name), `value`, `n` and `proportion`, the share of the units
carrying that truth value that were assigned to that class. Several
truth columns may be named at once and are stacked.

Each truth column is compared against the level it describes, which is
read off the data rather than guessed: a column that takes one value
within every group is a property of the group and is compared against
`group_class`; a column that varies inside any group cannot be a
group-level label and is compared against `profile`. The `assignment`
column records which comparison was made, so the choice is never silent.
Group-level recovery counts each group once, even when groups have
different numbers of observations. Unused factor levels are omitted,
because their within-truth proportions are undefined. Every column is
compared against `profile` for a fit with no discrete group classes.

The labels are the model's own and carry no order, so a recovery table
whose mass sits off the diagonal is a relabelling, not a failure.

## Information criteria

Let `q` be the number of free parameters, `n` the chosen sample size,
and `EN` the classification entropy of the level matching that
convention. The criteria are `deviance = -2L`, `aic = -2L + 2q`,
`bic = -2L + q log(n)`, `sabic = -2L + q log((n + 2) / 24)`,
`caic = -2L + q (log(n) + 1)`, `awe = -2(L - EN) + 2q (1.5 + log(n))`,
`icl = -2L + q log(n) + 2 EN`, `kic = -2L + 3(q + 1)`, and
`clc = -2L + 2 EN`. Note that `clc` uses the entropy sum, as `icl` and
`awe` here do; `tidyLPA` reports a `CLC` built from relative entropy
instead, which is bounded by one and so penalizes almost nothing, and
the two numbers are not comparable.

Every criterion that depends on a sample size is reported under both
multilevel conventions, `"groups"` and `"individuals"`, and uses
group-level posteriors under the first and individual-level posteriors
under the second. That is a stated per-level choice, not a unique
multilevel definition: use one convention consistently across compared
candidates. `deviance`, `aic` and `kic` depend on no sample size and
carry `convention = NA_character_`. `clc` depends on none either, but
its convention selects which level's classification uncertainty it
penalizes, so it is reported once per convention. Individual sample
sizes exclude rows with no observed indicators. Where the group
classification entropy is undefined, group-level `awe`, `icl` and `clc`
are `NA`.

## Classification quality

The odds of correct classification for class `k` is
`(p / (1 - p)) / (r / (1 - r))`, where `p` is the average posterior in
the assigned class and `r` is the model-estimated class proportion.
Values near one indicate classification no better than the class
proportion alone; values of five or more are conventionally read as
adequate separation. It is `NA_real_` wherever it is undefined, that is
whenever `average_posterior` is missing, or it or `estimated_proportion`
is zero or one, so a single-class solution always reports `NA_real_`.
`average_posterior` is itself `NA_real_` for a class with no modal
members, where the mean is taken over nothing and is undefined rather
than zero. Modal assignment discards classification uncertainty, which
is why `n_modal` and `estimated_n` differ whenever entropy is below one.

`relative_entropy` is `1 - EN / (n log K)`, on the zero-to-one scale, so
higher is sharper; `entropy_sum` is the raw `EN` the entropy-penalized
criteria use, so lower is sharper. `relative_entropy` is `NA_real_` when
a level has a single class, where it is undefined rather than perfect.

A fit with individual profiles but no discrete group classes returns
individuals only for `level = "both"`, and raises
`latents_no_group_classes` for `level = "groups"`.

## Bivariate residuals

The measurement model assumes the indicators are independent within a
profile. `"residuals"` tests that pair by pair, weighting every
observation by the posterior the fit holds for it, so a frame in any
other row order would pair each observation with someone else's
posterior; the alignment is checked rather than assumed. A large
residual is evidence of local dependence between two indicators, which
the enumeration grid will otherwise absorb by asking for an extra
profile. The p-values are approximate, and `p_adjusted` applies the
multiplicity correction named by `adjust` across the pairs.

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
#> 1        1     browse  0.5395948 0.5898645          0.7680263
#> 2        1   lectures  0.4274286 0.8446803          0.9190649
#> 3        1 forum_read  0.6198475 0.4814627          0.6938751
#> 4        1 forum_post  0.5305152 0.6887698          0.8299216
#> 5        1 attendance  0.6490757 0.3841104          0.6197664
#> 6        2     browse -0.7655654 0.5284937          0.7269757
#> 7        2   lectures -0.6064745 0.5383778          0.7337424
#> 8        2 forum_read -0.8791779 0.3631402          0.6026112
#> 9        2 forum_post -0.7526676 0.4211809          0.6489846
#> 10       2 attendance -0.9206255 0.3736487          0.6112681
get_results(fit, what = "profile_probabilities")
#>   group_class profile probability group_class_probability
#> 1           1       1   0.1703746               0.3252304
#> 2           1       2   0.8296254               0.3252304
#> 3           2       1   0.7858826               0.6747696
#> 4           2       2   0.2141174               0.6747696
get_results(fit, what = "entropy")
#>         level n_classes n_units entropy_sum relative_entropy
#> 1 individuals         2    1422   60.943122        0.9381699
#> 2      groups         2     106    4.433469        0.9396590
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
#> 5 group_class     1 student_type  committed   1 0.01612903
#> 6 group_class     2 student_type  committed  61 0.98387097
#> 7 group_class     1 student_type   wavering  34 0.77272727
#> 8 group_class     2 student_type   wavering  10 0.22727273
```
