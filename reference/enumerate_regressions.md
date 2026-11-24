# Fit mixture regression models with different numbers of classes, and compare them

`enumerate_regressions()` fits mixture regression models with different
numbers of latent classes and compares them. The same outcome model is
fitted across the requested class counts, so that researchers can assess
whether additional classes provide a useful representation of variation
in outcome-predictor relationships. The outcome distribution,
predictors, and other supplied model settings remain common across
candidates.

For supported two-level models, the function can also vary the number of
group classes and fit each requested combination of class counts. The
resulting table brings together information criteria, classification
measures, class sizes, convergence, and likelihood replication.
Individual fitted models are retained for examination of their
coefficients, predictions, and class assignments.

An optional parametric bootstrap likelihood-ratio test compares models
with consecutive numbers of regression classes while holding the number
of group classes fixed. This provides additional evidence about whether
the improvement in likelihood supports the more complex specification.

## Usage

``` r
enumerate_regressions(
  formula,
  data,
  n_classes = 1:4,
  n_group_classes = 1L,
  bootstrap = 0L,
  bootstrap_starts = 3L,
  seed = NULL,
  ...
)
```

## Arguments

- formula:

  A model formula `outcome ~ predictors`. Factors, interactions and
  [`offset()`](https://rdrr.io/r/stats/offset.html) are supported.

- data:

  A data frame.

- n_classes:

  Positive integer class counts to fit. Default `1:4`.

- n_group_classes:

  Positive integer group-class counts (two-level model; needs `id` in
  `...`), crossed with `n_classes`. Default 1, no group classes.

- bootstrap:

  Number of bootstrap replicates for the likelihood-ratio test of
  `k - 1` against `k` classes. `0` (the default) skips the test. A
  positive number adds the `blrt_` columns to the table; each replicate
  refits two models, so the run time grows with it. Not available with
  sampling weights.

- bootstrap_starts:

  Random starts for each bootstrap refit. Default 3.

- seed:

  Optional integer seed for the starts and the bootstrap; when given,
  the caller's random-number state is restored on exit.

- ...:

  Further arguments to
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md),
  used for every model:

  - `family`: `"gaussian"` (default), `"binomial"`, `"poisson"`,
    `"negative_binomial"` or `"ordinal"`;

  - `id` and `class_level`: `"observation"` (default, a class per row)
    or `"group"` (a class per unit of `id`);

  - `variance`: `"varying"` (default) or `"equal"`;

  - `common`, `membership`, `random` and the other model arguments;

  - `n_starts`, `max_iter`, `tol` and `missing` (`"error"` or `"omit"`).

  `n_classes`, `n_group_classes` and `vcov_type` are not accepted here:
  the class counts are this function's own arguments, and the models are
  fitted without a stored covariance matrix. A name that is not an
  argument of
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
  is refused.

## Value

An object of class `latents_regression_enumeration`. Its table (from
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
[`summary()`](https://rdrr.io/r/base/summary.html) or
`get_results(x, "fit")`) has one row per model, with these columns:

- `n_classes`, `n_group_classes`, `log_likelihood`, `n_parameters`;

- `aic`, `bic`, `bic_rows` (BIC with the number of rows as the sample
  size), `sabic` and `icl`: lower is better;

- `entropy`: between 0 and 1, higher meaning clearer classification;

- `smallest_share`: the proportion of the smallest class;

- `converged` and `n_best_replicated` (the number of starts that reached
  the best log likelihood);

- `note`: why a model could not be fitted; its other columns are `NA`;

- `best_bic`: `TRUE` for the rows with the smallest reported BIC,
  without filtering on convergence;

- with `bootstrap > 0`: `blrt_statistic`, `blrt_p_value`,
  `blrt_replicates` (the valid replicates), `blrt_n_invalid` (the
  failed, unconverged, separated or reversed-likelihood replicates) and
  `blrt_flagged` (a compatibility alias for `blrt_n_invalid`). The test
  of a row compares it with the row that has one class fewer. These
  columns are `NA` where the `k - 1` or `k` class model is missing or
  failed, and the p-value is `NA` when the observed fits fail validation
  or fewer than half of the requested replicates are valid.

`get_results(x, "model", n_classes = , n_group_classes = )` returns one
fitted model (a `latents_mixture_regression` or `latents_growth_mixture`
fit), `get_results(x, "tree")` the class tree and
`get_results(x, "all")` the named list of the `fit` and `tree` tables;
any other `what` raises `latents_bad_argument`.
[`print()`](https://rdrr.io/r/base/print.html) shows the main columns of
the table. [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
draws AIC, BIC, SABIC and ICL across the classes (see
[`plot_enumeration()`](https://pak.dynasite.org/latents/reference/plot_enumeration.md)
for its other views, and `plot_views(x)` for the views a regression grid
offers).

## Details

The optional parametric bootstrap likelihood-ratio test compares `k - 1`
with `k` classes at a fixed number of group classes. Outcomes are
simulated under the smaller model, and both models are refitted to each
simulated dataset. The p-value is `(1 + b) / (1 + V)`, where `b` is the
number of valid replicate statistics at least as large as the observed
statistic and `V` is the number of valid replicates.

A replicate is invalid if a refit fails, does not converge, produces a
separation or degenerate random-effect covariance warning, or yields a
negative likelihood-ratio statistic beyond numerical tolerance. The
p-value uses valid replicates only. Both valid and invalid counts are
reported, and `latents_failed_replicates` warns when any replicate is
invalid. The p-value is `NA` when fewer than half of the requested
replicates are valid. Excluding invalid replicates conditions the
reference distribution on successful estimation and may affect
calibration. Report the invalid proportion and examine the refit
diagnostics when interpreting the test.

## Interpretation

The `best_bic` column identifies the smallest reported BIC among the
fitted candidates; it does not exclude unconverged fits. Model
assessment should therefore consider convergence and likelihood
replication alongside information criteria, class sizes, and the meaning
of the estimated regression relationships. Bootstrap results also
require attention to unsuccessful refits, whose exclusion can affect the
test.

## Conditions

A model whose fit raises `latents_no_valid_start` or `latents_bad_data`
at some count has a row of `NA` values with the condition message in
`note` (and [`print()`](https://rdrr.io/r/base/print.html) says how
many); when every model fails with the same `latents_bad_data` error (a
column that does not exist), that error is raised instead. Other errors
stop the enumeration. `latents_bad_argument` for invalid `n_classes`,
`n_group_classes`, `bootstrap` or `bootstrap_starts`, for an invalid
`seed`, for an argument in `...` that
[`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
does not have, for `n_classes`, `n_group_classes` or `vcov_type` passed
through `...`, for a bootstrap test whose outcome is not a plain column
name, and from
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
for a model that is not in the grid or could not be fitted.
`latents_unsupported_weights` for a bootstrap test with sampling
weights. `latents_failed_replicates` is warned, naming the count, when
some replicates are invalid (the p-value then uses the valid ones), and
when the observed fits fail validation or fewer than half of the
replicates are valid (the p-value is then `NA`).

## Examples

``` r
classes <- enumerate_regressions(score ~ hours, data = study_hours,
                                 n_classes = 1:3, n_starts = 3, seed = 1)
classes
#> Mixture regression class enumeration
#> 
#>  n_classes n_group_classes log_likelihood n_parameters  bic  icl entropy
#>          1               1          -3245            3 6510 6510      NA
#>          2               1          -3144            7 6336 7031  0.4430
#>          3               1          -3133           11 6341 7271  0.5293
#>  smallest_share converged best_bic
#>         1.00000      TRUE    FALSE
#>         0.43587      TRUE     TRUE
#>         0.07512      TRUE    FALSE
#> 
#> One model: get_results(x, "model", n_classes = ).
get_results(classes, "model", n_classes = 2)
#> Mixture of gaussian regressions: 2 classes (one class per row)
#> score ~ hours
#> 900 rows | log likelihood -3144.0073 | BIC 6335.63 | entropy 0.443
#> Converged: TRUE | iterations: 49 | best likelihood reached by 2 of 2 completed starts (4 run)
#> 
#> Classes
#> 
#> Class    Share  Expected  Rows assigned  Avg. posterior  Residual SD
#> -------  -----  --------  -------------  --------------  -----------
#> Class 1  0.564    507.72            573           0.789         5.24
#> Class 2  0.436    392.28            327           0.829         6.89
#> 
#> Regression coefficients (95% CI)
#> 
#> Class    Term       Estimate          95% CI      p
#> -------  ---------  --------  --------------  -----
#> Class 1  Intercept     34.36  [32.95, 35.77]  <.001
#> Class 1  hours          4.57  [ 4.36,  4.77]  <.001
#> Class 2  Intercept     54.90  [52.76, 57.04]  <.001
#> Class 2  hours          0.81  [ 0.48,  1.14]  <.001
#> 
#> Every table: get_results(x, what = ), e.g. "classes", "membership", "fit", "assignments".
# \donttest{
# With the bootstrap likelihood-ratio test of one against two classes.
# `bootstrap` is small so the example runs quickly; use many more
# replicates for inference.
with_test <- enumerate_regressions(score ~ hours, data = study_hours,
                                   n_classes = 1:2, variance = "equal",
                                   n_starts = 3, seed = 1, bootstrap = 4)
as.data.frame(with_test)
#>   n_classes n_group_classes log_likelihood n_parameters      aic      bic
#> 1         1               1      -3244.958            3 6495.916 6510.323
#> 2         2               1      -3149.398            6 6310.797 6339.611
#>   bic_rows    sabic      icl   entropy smallest_share converged
#> 1 6510.323 6500.795 6510.323        NA        1.00000      TRUE
#> 2 6339.611 6320.556 7010.972 0.4619063        0.37308      TRUE
#>   n_best_replicated note best_bic blrt_statistic blrt_p_value blrt_replicates
#> 1                 4 <NA>    FALSE             NA           NA              NA
#> 2                 4 <NA>     TRUE        191.119          0.2               4
#>   blrt_n_invalid blrt_flagged
#> 1             NA           NA
#> 2              0            0
# }
```
