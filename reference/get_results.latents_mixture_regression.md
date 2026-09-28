# Tables of a mixture-of-regressions fit

Tables of a mixture-of-regressions fit

## Usage

``` r
# S3 method for class 'latents_mixture_regression'
get_results(
  x,
  what = "coefficients",
  level = 0.95,
  vcov_type = NULL,
  data = NULL,
  truth = NULL,
  by = c("class", "group_class"),
  ...
)

# S3 method for class 'latents_mixture_regression'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  what = "coefficients",
  level = 0.95,
  vcov_type = NULL,
  ...
)

# S3 method for class 'latents_mixture_regression'
coef(object, ...)

# S3 method for class 'latents_mixture_regression'
vcov(object, type = NULL, ...)

# S3 method for class 'latents_mixture_regression'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'latents_mixture_regression'
logLik(object, ...)

# S3 method for class 'latents_mixture_regression'
nobs(object, ...)
```

## Arguments

- x:

  A fit from
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md).

- what:

  The table to return; see *Tables*. `"all"` returns every table in a
  named list.

- level:

  Confidence level for the intervals.

- vcov_type:

  `NULL` for the covariance stored with the fit, or `"observed"`,
  `"robust"` or `"opg"` to recompute it.

- data, truth:

  For `what = "recovery"`: the data frame the model was fitted to and
  the name of its column holding a known classification.

- by:

  For `what = "recovery"`: `"class"`, or `"group_class"` for the group
  classes of a two-level fit.

- ...:

  Unused.

- row.names, optional:

  Unused; part of the generic.

- object:

  A fit from
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md).

- type:

  Covariance type for [`vcov()`](https://rdrr.io/r/stats/vcov.html); as
  `vcov_type`.

- parm:

  Parameter names or positions; all by default.

## Value

A base `data.frame` as described under *Tables*, or a named list of them
for `what = "all"`.

## Tables

- `coefficients`:

  One row per class and regression term: `class` (`"common"` for a
  shared term), `term`, `estimate`, `std_error`, `statistic` (Wald z),
  `p_value`, `conf_low`, `conf_high`, `p_adjusted` (Benjamini-Hochberg
  across the non-intercept rows; `NA` for intercepts) and, for the
  binomial and Poisson families, `exp_estimate` (odds or rate ratio).

- `classes`:

  One row per class: `share` (mean posterior), `count` (summed
  posterior), `n_assigned` (modal assignment), `mean_posterior` (average
  posterior of the units assigned to it); `sigma` and `sigma_std_error`
  for the Gaussian family; `prior` and `prior_std_error` when membership
  has no covariates. Units are groups under `class_level = "group"`.

- `membership`:

  One row per multinomial-logit coefficient: `model` (`"class"` or
  `"group_class"`), `class`, `term`, `reference`, the Wald columns and
  `odds_ratio`. In the two-level model the class intercepts are per
  group class, `(Intercept):group_class_h`.

- `group_classes`:

  Two-level model: one row per group class and class, with the
  model-implied class `probability`, the `group_share`,
  `n_groups_assigned`, and `std_error` of the probability when
  membership has no covariates.

- `fit`:

  One row: log likelihood, parameter count, `aic`, `bic` (penalized by
  the number of independent units: rows, or groups when `id` defines
  them), `bic_rows` (always by rows, as flexmix reports it), `sabic`,
  `icl`, relative `entropy` (and `group_entropy`), the smallest class
  share, convergence and start replication.

- `assignments`:

  One row per data row used: `row` (its position in the supplied data),
  the `id` column, modal `class`, its `posterior`, and one
  `probability_class_k` column per class; two-level fits add
  `group_class` and `group_posterior`.

- `groups`:

  One row per group (`class_level = "group"` or the two-level model):
  its size, modal class or group class, and posterior probabilities.

- `fitted`:

  One row per data row: `observed`, the posterior-weighted `fitted`
  mean, `fitted_modal` (the modal class's mean), `residual` and `class`.

- `starts`:

  One row per start: its log likelihood, convergence, iterations,
  `stage` (`"screened"`: stopped after the 50-iteration screening
  because better starts existed; `"completed"`: run to convergence or
  `max_iter`), whether it degenerated, and which one was selected.
  `n_best_replicated` counts completed starts only.

- `classification`:

  One row per assigned (modal) class and class: the `mean_posterior`
  probability of `class` among the units assigned to `assigned`, and
  `n_assigned`. Its diagonal is each class's average posterior
  probability; values near one mean well-separated classes.

- `recovery`:

  Needs `data` and `truth`: the modal classes cross-tabulated against a
  known classification, one row per assigned class and true value, with
  the count `n` and its `share` of the assigned class. Group-level under
  `class_level = "group"` or `by = "group_class"`, where `truth` must be
  constant within groups.

## Examples

``` r
fit <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2,
                          n_starts = 3, seed = 1)
get_results(fit, "coefficients")
#>     class        term   estimate std_error statistic      p_value   conf_low
#> 1 class_1 (Intercept) 34.3545204 0.7192917 47.761595 0.000000e+00 32.9447345
#> 2 class_1       hours  4.5658711 0.1047310 43.596167 0.000000e+00  4.3606021
#> 3 class_2 (Intercept) 54.8951437 1.0927048 50.237855 0.000000e+00 52.7534818
#> 4 class_2       hours  0.8083711 0.1696013  4.766303 1.876367e-06  0.4759587
#>   conf_high   p_adjusted
#> 1 35.764306           NA
#> 2  4.771140 0.000000e+00
#> 3 57.036806           NA
#> 4  1.140784 1.876367e-06
get_results(fit, "fit")
#>     family     nesting n_classes n_group_classes n_observations n_groups
#> 1 gaussian observation         2               1            900       NA
#>   log_likelihood n_parameters      aic      bic bic_rows    sabic      icl
#> 1      -3144.007            7 6302.015 6335.632 6335.632 6313.401 7030.593
#>     entropy group_entropy smallest_share converged iterations n_starts
#> 1 0.4429906            NA      0.4360588      TRUE         42        4
#>   n_best_replicated vcov_type
#> 1                 3  observed
get_results(fit, "coefficients", vcov_type = "robust")
#>     class        term   estimate std_error statistic      p_value   conf_low
#> 1 class_1 (Intercept) 34.3545204 0.8102139 42.401790 0.000000e+00 32.7665303
#> 2 class_1       hours  4.5658711 0.1142614 39.959878 0.000000e+00  4.3419229
#> 3 class_2 (Intercept) 54.8951437 1.1172910 49.132361 0.000000e+00 52.7052937
#> 4 class_2       hours  0.8083711 0.1763408  4.584141 4.558558e-06  0.4627495
#>   conf_high   p_adjusted
#> 1 35.942510           NA
#> 2  4.789819 0.000000e+00
#> 3 57.084994           NA
#> 4  1.153993 4.558558e-06
get_results(fit, "recovery", data = study_hours, truth = "strategy")
#>   assigned strategy   n     share
#> 1  class_1     deep 448 0.7818499
#> 2  class_2     deep  49 0.1498471
#> 3  class_1  surface 125 0.2181501
#> 4  class_2  surface 278 0.8501529
```
