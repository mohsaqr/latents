# Assess the number of classes in mixture regression models

`enumerate_regressions()` compares mixture regression models with
different numbers of latent classes. It fits the same outcome model
across the requested class counts, allowing researchers to assess
whether additional classes provide a useful representation of variation
in outcome–predictor relationships. The outcome distribution,
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

  Positive integer group-class counts (two-level model; needs `id`),
  crossed with `n_classes`. Default 1.

- bootstrap:

  Number of bootstrap replicates for the likelihood-ratio test; `0` (the
  default) skips it. Not available with sampling weights.

- bootstrap_starts:

  Random starts for each bootstrap refit. Default 3.

- seed:

  `NULL` (default) advances the caller's random stream. An integer seed
  restores the caller's `.Random.seed` afterwards. R does not store a
  cached Box-Muller normal draw in `.Random.seed`, so that hidden cached
  draw cannot be preserved.

- ...:

  Further arguments to
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md),
  such as `family`, `id`, `class_level`, `common`, `membership` or
  `variance`. `n_classes`, `n_group_classes` and `vcov_type` are not
  accepted here.

## Value

An object of class `latents_regression_enumeration`. Its table (from
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) or
`get_results(x, "fit")`) has one row per model, with columns
`n_classes`, `n_group_classes`, `log_likelihood`, `n_parameters`, `aic`,
`bic`, `bic_rows`, `sabic`, `icl`, `entropy`, `smallest_share`,
`converged`, `n_best_replicated`, `note` (why a model could not be
fitted), `best_bic` (`TRUE` for rows attaining the smallest reported
BIC, without filtering on convergence), and, with `bootstrap > 0`,
`blrt_statistic`, `blrt_p_value`, `blrt_replicates` (the valid
replicates), `blrt_n_invalid` (the failed, unconverged, separated or
reversed-likelihood replicates) and `blrt_flagged` (a compatibility
alias for `blrt_n_invalid`). The test columns are `NA` where the `k - 1`
or `k` class model is missing or failed, and the p-value is `NA` when
the observed fits fail validation or fewer than half of the requested
replicates are valid.
`get_results(x, "model", n_classes = , n_group_classes = )` returns one
fitted model, `get_results(x, "tree")` the class tree and
`get_results(x, "all")` the named list of the `fit` and `tree` tables;
any other `what` raises `latents_bad_argument`.
[`print()`](https://rdrr.io/r/base/print.html) shows the main columns of
the table. [`summary()`](https://rdrr.io/r/base/summary.html) returns
the same table as
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws AIC, BIC,
SABIC and ICL across the classes (see
[`plot_enumeration()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
for its other views, and `plot_views(x)` for the views a regression grid
offers).

## Details

The optional parametric bootstrap likelihood-ratio test compares `k - 1`
with `k` classes at a fixed number of group classes (McLachlan, 1987;
Nylund, Asparouhov and Muthen, 2007). Outcomes are simulated under the
smaller model, and both models are refitted to each simulated dataset.
The p-value is `(1 + b) / (1 + V)`, where `b` is the number of valid
replicate statistics at least as large as the observed statistic and `V`
is the number of valid replicates.

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
#>          2               1          -3144            7 6336 7031  0.4430
#>          3               1          -3133           11 6341 7271  0.5293
#>  smallest_share converged best_bic
#>         1.00000      TRUE    FALSE
#>         0.43587      TRUE     TRUE
#>         0.07512      TRUE    FALSE
#> 
#> One model: get_results(x, "model", n_classes = ).
# \donttest{
# A small bootstrap count is used only to limit example execution time.
with_test <- enumerate_regressions(score ~ hours, data = study_hours,
                                   n_classes = 1:2, n_starts = 3, seed = 1,
                                   bootstrap = 4)
#> Warning: 3 of 4 bootstrap replicates were invalid (a refit failed, did not converge or separated); fewer than half are valid, so blrt_p_value is NA.
as.data.frame(with_test)
#>   n_classes n_group_classes log_likelihood n_parameters      aic      bic
#> 1         1               1      -3244.958            3 6495.916 6510.323
#> 2         2               1      -3144.007            7 6302.015 6335.631
#>   bic_rows    sabic      icl   entropy smallest_share converged
#> 1 6510.323 6500.795 6510.323        NA      1.0000000      TRUE
#> 2 6335.631 6313.401 7030.539 0.4430333      0.4358709      TRUE
#>   n_best_replicated note best_bic blrt_statistic blrt_p_value blrt_replicates
#> 1                 4 <NA>    FALSE             NA           NA              NA
#> 2                 2 <NA>     TRUE        201.901           NA               1
#>   blrt_n_invalid blrt_flagged
#> 1             NA           NA
#> 2              3            3
# }
```
