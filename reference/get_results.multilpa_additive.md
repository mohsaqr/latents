# Tables of an additive group-class fit

Tidy tables of a fit from `multilpa(family = "additive")`. Every table
is a base `data.frame`; the fit itself never needs to be reached into.

## Usage

``` r
# S3 method for class 'multilpa_additive'
get_results(
  x,
  what = "parameters",
  level = 0.95,
  vcov_type = c("observed", "robust", "opg"),
  adjust = .multilpa_p_adjust_methods,
  data = NULL,
  truth = NULL,
  ...
)

# S3 method for class 'multilpa_additive'
as.data.frame(x, row.names = NULL, optional = FALSE, what = "parameters", ...)

# S3 method for class 'multilpa_additive'
coef(object, scale = c("natural", "unconstrained"), ...)

# S3 method for class 'multilpa_additive'
vcov(
  object,
  type = c("observed", "robust", "opg"),
  scale = c("natural", "unconstrained"),
  ...
)

# S3 method for class 'multilpa_additive'
confint(object, parm, level = 0.95, type = c("observed", "robust", "opg"), ...)

# S3 method for class 'multilpa_additive'
logLik(object, ...)

# S3 method for class 'multilpa_additive'
nobs(object, ...)
```

## Arguments

- x:

  A `multilpa_additive` fit.

- what:

  Which table; `"all"` returns a named list of every table.

- level:

  Confidence level of the intervals in `"parameters"`.

- vcov_type:

  `"observed"` (inverse observed information), `"robust"` (CR0 sandwich
  over groups) or `"opg"` (inverse outer product of the group scores).

- adjust:

  p-value adjustment across the class means, one of
  [stats::p.adjust.methods](https://rdrr.io/r/stats/p.adjust.html);
  `"none"` by default, as elsewhere in latents.

- data, truth:

  For `what = "recovery"` only: the data the fit was made from and the
  name of a column, constant within groups, holding a known group
  classification.

- ...:

  Unused.

- row.names, optional:

  Unused; part of the generic.

- object:

  A `multilpa_additive` fit.

- scale:

  `"natural"` (means, variances, class weights) or `"unconstrained"`
  (means, log variances, class-1-baseline logits). Names of
  [`coef()`](https://rdrr.io/r/stats/coef.html) and the rows of
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) agree on either scale.

- type:

  Covariance type, as `vcov_type`.

- parm:

  Natural parameter names or positions; all by default.

## Value

A base `data.frame`, or a named list of them for `what = "all"`.

## Tables

- `parameters`:

  One row per natural parameter: `level` (`"between"` for class means
  and between-group variances, `"within"` for the shared within-group
  variances, `"group_class"` for class weights), `group_class`
  (`"shared"` for a parameter common to all classes), `indicator` (`NA`
  for weights), `parameter` (`"mean"`, `"variance"`, `"weight"`),
  `estimate`, `standard_error`, `statistic`, `p_value`, `p_adjusted`,
  `conf_low`, `conf_high`. Only means carry a Wald test. Intervals for
  variances are formed on the log scale and for weights on the logit
  scale. When inference is refused (boundary, unconverged, singular
  information) the standard-error columns are `NA` and a
  `latents_no_standard_errors` message says why.

- `group_classes`:

  One row per group class: `weight`, `count` (summed posterior),
  `n_assigned` (modal), `mean_posterior` among the groups assigned to
  it, and `effective_groups`: the number of perfectly classified groups
  that would estimate the class weight as precisely as this fit does,
  `w^2 (1 - w) / Var(w)`. It equals the class's expected group count
  when classes are perfectly separated and falls when a class is small
  or poorly separated. Below 50, intervals for class parameters are
  flagged with a `latents_weak_class` warning (threshold and its
  simulation evaluation: `validation/ADDITIVE_SIMULATION.md` in the
  source repository). `NA` for one class or when inference is refused.

- `groups`:

  One row per group (named by the `id` column): size `n`, modal
  `group_class`, its `posterior`, and `probability_group_class_h`.

- `intercepts`:

  One row per group and indicator: `n`, the observed group `average`,
  and the posterior mean and variance of the group's intercept, mixed
  over its class posterior (so the variance includes uncertainty about
  class membership). Conditional on the estimated parameters; not a
  statement about their sampling uncertainty.

- `conditional_intercepts`:

  As `intercepts`, one row per group, indicator and group class,
  conditional on that class, with the class `posterior`.

- `assignments`:

  One row per data row: `row`, the `id` column, the group's modal
  `group_class` and `posterior`, and `assignment_level = "group"`: every
  row inherits its group's class.

- `classification`:

  One row per assigned and group class: the `mean_posterior` of
  `group_class` among groups assigned to `assigned`.

- `fit`:

  One row: log likelihood, `n_parameters`, `aic`, `bic` (penalized by
  the number of groups), `bic_individual` (by rows), `n_groups`,
  `n_obs`, relative `entropy` over groups, the smallest class weight,
  `min_effective_groups`, `converged`, `boundary`, `kkt` and start
  replication.

- `recovery`:

  Needs `data` and `truth`: one row per modal group class and known
  value, with the number of groups `n` and their `share` of the assigned
  class. Not part of `what = "all"`.

- `starts`:

  One row per start: likelihood, convergence, iterations, number of
  between variances held at zero, error message of a failed start, and
  which start was selected.

## Examples

``` r
set.seed(1)
ratings <- data.frame(
  team = rep(seq_len(40), each = 6),
  climate = rep(rnorm(40, rep(c(-1, 1), each = 20), 0.5), each = 6) +
    rnorm(240),
  support = rep(rnorm(40, rep(c(-1, 1), each = 20), 0.5), each = 6) +
    rnorm(240))
fit <- multilpa(ratings, c("climate", "support"), "team",
                n_group_classes = 2, family = "additive", n_starts = 3,
                seed = 1)
get_results(fit)
#> Warning: Group classes supported by fewer than 50 effective groups (group_class_1: 19.7 effective of 20.4 expected groups, 96% of the information kept; group_class_2: 18.8 effective of 19.6 expected groups, 96% of the information kept). A small information share means poor separation; a small expected count means few groups. Wald intervals for class means and weights can be miscalibrated; see get_results(x, "group_classes").
#>          level   group_class indicator parameter    estimate standard_error
#> 1      between group_class_1   climate      mean -0.90841249     0.13365766
#> 2      between group_class_1   support      mean -0.91030607     0.13337034
#> 3      between group_class_2   climate      mean  1.07151874     0.12180149
#> 4      between group_class_2   support      mean  1.05245634     0.14155613
#> 5       within        shared   climate  variance  1.02556718     0.10255672
#> 6       within        shared   support  variance  1.14788037     0.11478804
#> 7      between group_class_1   climate  variance  0.17466828     0.11199090
#> 8      between group_class_1   support  variance  0.14745703     0.11201928
#> 9      between group_class_2   climate  variance  0.09965351     0.09154283
#> 10     between group_class_2   support  variance  0.18641044     0.12340834
#> 11 group_class group_class_1      <NA>    weight  0.51044770     0.08056876
#> 12 group_class group_class_2      <NA>    weight  0.48955230     0.08056876
#>    statistic      p_value   p_adjusted    conf_low  conf_high
#> 1  -6.796562 1.071455e-11 1.071455e-11 -1.17037668 -0.6464483
#> 2  -6.825401 8.768002e-12 8.768002e-12 -1.17170714 -0.6489050
#> 3   8.797255 1.402033e-18 1.402033e-18  0.83279221  1.3102453
#> 4   7.434905 1.046435e-13 1.046435e-13  0.77501142  1.3299013
#> 5         NA           NA           NA  0.84303180  1.2476256
#> 6         NA           NA           NA  0.94357511  1.3964223
#> 7         NA           NA           NA  0.04971126  0.6137242
#> 8         NA           NA           NA  0.03326823  0.6535837
#> 9         NA           NA           NA  0.01646529  0.6031368
#> 10        NA           NA           NA  0.05092766  0.6823179
#> 11        NA           NA           NA  0.35660584  0.6623357
#> 12        NA           NA           NA  0.33766432  0.6433942
get_results(fit, "group_classes")
#>     group_class    weight    count n_assigned mean_posterior effective_groups
#> 1 group_class_1 0.5104477 20.41791         20      0.9996703         19.65026
#> 2 group_class_2 0.4895523 19.58209         20      0.9787749         18.84587
get_results(fit, "parameters", vcov_type = "robust")
#> Warning: Group classes supported by fewer than 50 effective groups (group_class_1: 19.7 effective of 20.4 expected groups, 96% of the information kept; group_class_2: 18.9 effective of 19.6 expected groups, 96% of the information kept). A small information share means poor separation; a small expected count means few groups. Wald intervals for class means and weights can be miscalibrated; see get_results(x, "group_classes").
#>          level   group_class indicator parameter    estimate standard_error
#> 1      between group_class_1   climate      mean -0.90841249     0.13395472
#> 2      between group_class_1   support      mean -0.91030607     0.12998999
#> 3      between group_class_2   climate      mean  1.07151874     0.12302933
#> 4      between group_class_2   support      mean  1.05245634     0.14287131
#> 5       within        shared   climate  variance  1.02556718     0.09223567
#> 6       within        shared   support  variance  1.14788037     0.11855311
#> 7      between group_class_1   climate  variance  0.17466828     0.09511175
#> 8      between group_class_1   support  variance  0.14745703     0.08484683
#> 9      between group_class_2   climate  variance  0.09965351     0.07702866
#> 10     between group_class_2   support  variance  0.18641044     0.11534586
#> 11 group_class group_class_1      <NA>    weight  0.51044770     0.08049887
#> 12 group_class group_class_2      <NA>    weight  0.48955230     0.08049887
#>    statistic      p_value   p_adjusted    conf_low  conf_high
#> 1  -6.781489 1.189432e-11 1.189432e-11 -1.17095891 -0.6458661
#> 2  -7.002894 2.507293e-12 2.507293e-12 -1.16508177 -0.6555304
#> 3   8.709457 3.053329e-18 3.053329e-18  0.83038568  1.3126518
#> 4   7.366464 1.752136e-13 1.752136e-13  0.77243371  1.3324790
#> 5         NA           NA           NA  0.85982533  1.2232578
#> 6         NA           NA           NA  0.93752859  1.4054285
#> 7         NA           NA           NA  0.06007739  0.5078284
#> 8         NA           NA           NA  0.04774015  0.4554568
#> 9         NA           NA           NA  0.02190502  0.4533582
#> 10        NA           NA           NA  0.05543311  0.6268610
#> 11        NA           NA           NA  0.35673162  0.6622131
#> 12        NA           NA           NA  0.33778693  0.6432684
```
