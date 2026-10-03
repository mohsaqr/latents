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
#> Regression coefficients (95% CI)
#> 
#> Class    Term       Estimate          95% CI      p
#> -------  ---------  --------  --------------  -----
#> Class 1  Intercept     34.35  [32.94, 35.76]  <.001
#> Class 1  hours          4.57  [ 4.36,  4.77]  <.001
#> Class 2  Intercept     54.90  [52.75, 57.04]  <.001
#> Class 2  hours          0.81  [ 0.48,  1.14]  <.001
get_results(fit, "fit")
#> Model fit
#> 
#> Family          gaussian
#> Nesting         observation
#> Classes         2
#> Group classes   1
#> Observations    900
#> Parameters      7
#> Log likelihood  -3144.01
#> AIC             6302.01
#> BIC             6335.63
#> SABIC           6313.40
#> ICL             7030.59
#> Entropy         0.443
#> Smallest class  43.6%
#> Converged       yes
get_results(fit, "coefficients", vcov_type = "robust")
#> Regression coefficients (95% CI)
#> 
#> Class    Term       Estimate          95% CI      p
#> -------  ---------  --------  --------------  -----
#> Class 1  Intercept     34.35  [32.77, 35.94]  <.001
#> Class 1  hours          4.57  [ 4.34,  4.79]  <.001
#> Class 2  Intercept     54.90  [52.71, 57.08]  <.001
#> Class 2  hours          0.81  [ 0.46,  1.15]  <.001
get_results(fit, "recovery", data = study_hours, truth = "strategy")
#> Recovery of a known classification
#> 
#> assigned  strategy    n  share
#> --------  --------  ---  -----
#> class_1   deep      448   0.78
#> class_2   deep       49   0.15
#> class_1   surface   125   0.22
#> class_2   surface   278   0.85
```
