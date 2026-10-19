# Compare mixture regressions with different numbers of classes

Fits
[`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
for every requested number of classes (crossed with numbers of group
classes, for the two-level model) and returns their fit statistics in
one table. Optionally adds a parametric bootstrap likelihood-ratio test
of `k - 1` against `k` classes, at the same number of group classes
(McLachlan 1987; Nylund, Asparouhov and Muthen 2007): outcomes are
simulated from the `k - 1` class fit, both models are refitted to every
replicate, and the p-value is `(1 + b) / (1 + B)`, where `b` counts
replicate statistics at least as large as the observed one and `B` is
the number of replicates requested.

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

# S3 method for class 'latents_regression_enumeration'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'latents_regression_enumeration'
get_results(
  x,
  what = c("fit", "model"),
  n_classes = NULL,
  n_group_classes = 1L,
  ...
)

# S3 method for class 'latents_regression_enumeration'
summary(object, ...)

# S3 method for class 'latents_regression_enumeration'
print(x, digits = 4L, ...)
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

  Positive integer group-class counts (two-level model; needs `id`),
  crossed with `n_classes`. Default 1.

- bootstrap:

  Number of bootstrap replicates for the likelihood-ratio test; `0` (the
  default) skips it. Not available with sampling weights.

- bootstrap_starts:

  Random starts for each bootstrap refit. Default 3.

- seed:

  `NULL` (default) or an integer seed. The caller's random state is
  restored afterwards.

- ...:

  Further arguments to
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md),
  such as `family`, `id`, `class_level`, `common`, `membership` or
  `variance`. `n_classes`, `n_group_classes` and `vcov_type` are not
  accepted here.

- x:

  A `latents_regression_enumeration` object.

- row.names, optional:

  Not used; part of the generic.

- what:

  `"fit"` for the comparison table, `"model"` for one fitted model.

- object:

  A `latents_regression_enumeration` object.

- digits:

  Significant digits.

## Value

An object of class `latents_regression_enumeration`. Its table (from
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) or
`get_results(x, "fit")`) has one row per model, with columns
`n_classes`, `n_group_classes`, `log_likelihood`, `n_parameters`, `aic`,
`bic`, `bic_rows`, `sabic`, `icl`, `entropy`, `smallest_share`,
`converged`, `n_best_replicated`, `note` (why a model could not be
fitted), `best_bic` (`TRUE` for the minimum-BIC row), and, with
`bootstrap > 0`, `blrt_statistic`, `blrt_p_value`, `blrt_replicates`
(the successful replicates) and `blrt_flagged` (the failed, unconverged,
separated or reversed-likelihood replicates). The test columns are `NA`
where the `k - 1` or `k` class model is missing or failed, and the
p-value is `NA` if any requested replicate fails validation.
`get_results(x, "model", n_classes = , n_group_classes = )` returns one
fitted model. [`summary()`](https://rdrr.io/r/base/summary.html) returns
the same table as
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws AIC, BIC,
SABIC and ICL across the classes (see
[`plot_enumeration()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
for its other views, and `plot_views(x)` for the views a regression grid
offers).

## Conditions

A model whose fit raises `latents_no_valid_start` or `latents_bad_data`
at some count gets a row of `NA` with the condition message in `note`;
other errors stop the enumeration. Invalid `n_classes`,
`n_group_classes`, `bootstrap` or `bootstrap_starts` fail the argument
checks. `latents_bad_argument` for an invalid `seed`, for `n_classes`,
`n_group_classes` or `vcov_type` passed through `...`, for a bootstrap
test whose outcome is not a plain column name, and from
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
for a model that is not in the grid or could not be fitted.
`latents_unsupported_weights` for a bootstrap test with sampling
weights. `latents_failed_replicates` is warned, and the p-value is `NA`,
when the observed fits or any replicate fail validation.

## References

McLachlan, G. J. (1987). On bootstrapping the likelihood ratio test
statistic for the number of components in a normal mixture. *Applied
Statistics*, 36, 318–324.

Nylund, K. L., Asparouhov, T., & Muthen, B. O. (2007). Deciding on the
number of classes in latent class analysis and growth mixture modeling:
A Monte Carlo simulation study. *Structural Equation Modeling*, 14,
535–569.

## Examples

``` r
classes <- enumerate_regressions(score ~ hours, data = study_hours,
                                 n_classes = 1:3, n_starts = 3, seed = 1)
classes
#> Mixture regression class enumeration
#> 
#>  n_classes n_group_classes log_likelihood n_parameters  bic  icl entropy
#>          1               1          -3245            3 6510 6510      NA
#>          2               1          -3144            7 6336 7030  0.4432
#>          3               1          -3133           11 6341 7272  0.5292
#>  smallest_share converged best_bic
#>         1.00000      TRUE    FALSE
#>         0.43526      TRUE     TRUE
#>         0.07523      TRUE    FALSE
#> 
#> One model: get_results(x, "model", n_classes = ).
# \donttest{
with_test <- enumerate_regressions(score ~ hours, data = study_hours,
                                   n_classes = 1:2, n_starts = 3, seed = 1,
                                   bootstrap = 19)
#> Warning: Some bootstrap fits failed validation, so p_value is NA. summary() reports every replicate; improve fitting and rerun.
as.data.frame(with_test)
#>   n_classes n_group_classes log_likelihood n_parameters      aic      bic
#> 1         1               1      -3244.958            3 6495.916 6510.323
#> 2         2               1      -3144.007            7 6302.015 6335.632
#>   bic_rows    sabic      icl   entropy smallest_share converged
#> 1 6510.323 6500.795 6510.323        NA      1.0000000      TRUE
#> 2 6335.632 6313.401 7030.366 0.4431723      0.4352614      TRUE
#>   n_best_replicated note best_bic blrt_statistic blrt_p_value blrt_replicates
#> 1                 4 <NA>    FALSE             NA           NA              NA
#> 2                 4 <NA>     TRUE       201.9009           NA              14
#>   blrt_flagged
#> 1           NA
#> 2            5
# }
```
