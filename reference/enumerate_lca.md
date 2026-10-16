# Compare latent class models

Fits single-level latent class models with each number of classes in
`n_classes` and returns them in one table for comparison by information
criteria, entropy and the diagnostics of each fit. Every indicator is
categorical unless named in `ordinal` or `count`, as in
[`lca()`](https://pak.dynasite.org/latents/reference/lca.md). It is
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
with `id = NULL`, and it returns the same object; the number of classes
is its `n_profiles` column.

## Usage

``` r
enumerate_lca(data, vars, n_classes = 1:4, ...)
```

## Arguments

- data:

  A data frame with one row per observation.

- vars:

  Names of the categorical indicator columns.

- n_classes:

  Positive whole numbers of classes to compare.

- ...:

  Further arguments for
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md),
  such as `missing`, `n_starts`, `seed` and `min_probability`, except
  `id`, `n_group_classes` and `categorical`.

## Value

An object of class `multilpa_enumeration`, as
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
returns: [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
gives one row per candidate,
[`summary()`](https://rdrr.io/r/base/summary.html) the candidate each
criterion prefers,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) the criteria
against the number of classes, and
[`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md)
one fitted candidate, named by its `n_profiles`.

## Details

For observations nested in groups, use
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
with the grouping column as `id` and the indicators as `categorical`.

## Conditions

`latents_bad_argument` when `id`, `n_group_classes`, `categorical` or
`model` is passed.

## See also

[`lca()`](https://pak.dynasite.org/latents/reference/lca.md) to fit one
model;
[`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md)
for continuous indicators.

## Examples

``` r
models <- enumerate_lca(student_esm, c("time_with_friends",
                                       "on_social_media", "sports"),
                        n_classes = 1:3, n_starts = 2, seed = 1)
summary(models)
#> Class enumeration: 3 candidates, 3 converged, 0 failed to fit
#> 3 distinct candidate(s) are minimal under some criterion.
#> No candidate is selected automatically. Choose one convention and keep it.
#> The candidates table shows 13 of its 35 columns; get_results(x, what = "candidates") returns all of them.
#> 
#> -- candidates ------------------------------------------------------
#>  candidate n_profiles log_likelihood n_parameters  aic  bic icl_individual
#>   classes1          1          -3280            3 6566 6584           6584
#>   classes2          2          -3269            7 6552 6593           7454
#>   classes3          3          -3269           11 6560 6624           6651
#>  profile_entropy avepp_min avepp_mean avepp_below converged boundary
#>               NA        NA         NA          NA      TRUE    FALSE
#>           0.7594    0.8993     0.8993           1      TRUE    FALSE
#>           0.9953    0.9983     0.9992           0      TRUE    FALSE
#> 
#> -- criteria --------------------------------------------------------
#>  criterion  convention candidate n_profiles n_group_classes model value
#>        aic        <NA>  classes2          2               1  <NA>  6552
#>        kic        <NA>  classes2          2               1  <NA>  6562
#>        bic      groups  classes1          1               1  <NA>  6584
#>        bic individuals  classes1          1               1  <NA>  6584
#>      sabic      groups  classes2          2               1  <NA>  6571
#>      sabic individuals  classes2          2               1  <NA>  6571
#>       caic      groups  classes1          1               1  <NA>  6587
#>       caic individuals  classes1          1               1  <NA>  6587
#>        awe      groups  classes1          1               1  <NA>  6616
#>        awe individuals  classes1          1               1  <NA>  6616
#>    ... 4 more rows.  get_results(x, what = "criteria")
#> 
#> -- classification --------------------------------------------------
#>  candidate n_profiles n_group_classes model       level class n_modal
#>   classes2          2               1  <NA> individuals     1    2582
#>   classes2          2               1  <NA> individuals     2       0
#>   classes3          3               1  <NA> individuals     1    2054
#>   classes3          3               1  <NA> individuals     2     185
#>   classes3          3               1  <NA> individuals     3     343
#>  proportion_modal estimated_n estimated_proportion average_posterior
#>           1.00000      2322.1              0.89934            0.8993
#>           0.00000       259.9              0.10066                NA
#>           0.79551      2053.6              0.79535            0.9996
#>           0.07165       185.4              0.07182            0.9997
#>           0.13284       343.0              0.13283            0.9983
#>  odds_correct_classification relative_entropy
#>                          1.0           0.7594
#>                           NA           0.7594
#>                        573.9           0.9953
#>                      46383.9           0.9953
#>                       3920.3           0.9953
#> 
#> -- deltas ----------------------------------------------------------
#>  candidate n_profiles n_group_classes model  criterion value delta_profile
#>   classes1          1               1  <NA> bic_groups  6584            NA
#>   classes2          2               1  <NA> bic_groups  6593        -9.597
#>   classes3          3               1  <NA> bic_groups  6624       -31.210
#>  from_profile delta_group_class from_group_class
#>          <NA>                NA             <NA>
#>      classes1                NA             <NA>
#>      classes2                NA             <NA>
#> 
#> -- knees -----------------------------------------------------------
#>   direction           line  criterion     from       to n_points knee_candidate
#>    profiles  1 group class bic_groups classes1 classes3        3           <NA>
#>  parameters all candidates bic_groups classes1 classes3        3           <NA>
#>  distance knee_without_largest stable
#>        NA                 <NA>     NA
#>        NA                 <NA>     NA
#>                                         note
#>  the criterion does not fall along this line
#>  the criterion does not fall along this line
#> 
#> 5 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
```
