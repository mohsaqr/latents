# Compare latent class models

`enumerate_lca()` fits latent class models with different numbers of
classes to examine heterogeneity in categorical response patterns. It
brings together evidence about fit, classification, and estimation
stability to support selection of the number of classes. Each fitted
candidate remains available for examination of its class-specific
response probabilities and class proportions.

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

  Numbers of classes to compare, specified as positive integers.

- ...:

  Further arguments for
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md),
  such as `missing`, `n_starts`, `seed` and `min_probability`, except
  `id` and `n_group_classes`. Every indicator is categorical unless
  named in `ordinal` or `count`. `categorical` is accepted but has no
  further effect; its names must be in `vars`.

## Value

An object of class `multilpa_enumeration`, as
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
returns: [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
gives one row per candidate,
[`summary()`](https://rdrr.io/r/base/summary.html) the candidate
minimizing each criterion,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) the criteria
against the number of classes, and
[`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md)
one fitted candidate. The number of classes is reported in the
`n_profiles` column of the enumeration table and selected through the
`n_profiles` argument of
[`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md).

## Details

This is the single-level interface to
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
for categorical indicators. Every indicator is categorical unless named
in `ordinal` or `count`, as in
[`lca()`](https://pak.dynasite.org/latents/reference/lca.md).

For observations nested in groups, use
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
with the grouping column as `id` and the indicators as `categorical`.

## Conditions

`latents_bad_argument` when `id` or `n_group_classes` is passed, or when
`model` is passed (no indicator is continuous);
`latents_bad_categorical` when `categorical` names a column not in
`vars`.

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
#> Warning: The model is not identified: it has 11 free parameters, but its indicators' 8 response patterns leave only 7 free cell probabilities (-4 degrees of freedom). Its class sizes and response probabilities are not determined by the data; fit fewer classes or use more indicators.
summary(models)
#> Class enumeration: 3 candidates, 3 converged, 0 failed to fit
#> 2 distinct candidate(s) are minimal under some criterion.
#> No candidate is selected automatically. Choose one convention and keep it.
#> The candidates table shows 13 of its 37 columns; get_results(x, what = "candidates") returns all of them.
#> 
#> -- candidates ------------------------------------------------------
#>  candidate n_profiles log_likelihood n_parameters  aic  bic icl_individual
#>   classes1          1          -3280            3 6566 6584           6584
#>   classes2          2          -3270            7 6555 6596           6655
#>   classes3          3          -3269           11 6560 6624           6651
#>  profile_entropy avepp_min avepp_mean avepp_below converged boundary
#>               NA        NA         NA          NA      TRUE    FALSE
#>           0.9835    0.9942     0.9970           0      TRUE    FALSE
#>           0.9953    0.9983     0.9992           0      TRUE    FALSE
#> 
#> -- best ------------------------------------------------------------
#>  criterion  convention candidate n_profiles n_group_classes value
#>        aic        <NA>  classes2          2               1  6555
#>        kic        <NA>  classes2          2               1  6565
#>        bic      groups  classes1          1               1  6584
#>        bic individuals  classes1          1               1  6584
#>      sabic      groups  classes2          2               1  6574
#>      sabic individuals  classes2          2               1  6574
#>       caic      groups  classes1          1               1  6587
#>       caic individuals  classes1          1               1  6587
#>        awe      groups  classes1          1               1  6616
#>        awe individuals  classes1          1               1  6616
#>    ... 4 more rows.  get_results(x, what = "best")
#> 
#> -- criteria --------------------------------------------------------
#>  criterion  convention candidate n_profiles n_group_classes value
#>        aic        <NA>  classes2          2               1  6555
#>        kic        <NA>  classes2          2               1  6565
#>        bic      groups  classes1          1               1  6584
#>        bic individuals  classes1          1               1  6584
#>      sabic      groups  classes2          2               1  6574
#>      sabic individuals  classes2          2               1  6574
#>       caic      groups  classes1          1               1  6587
#>       caic individuals  classes1          1               1  6587
#>        awe      groups  classes1          1               1  6616
#>        awe individuals  classes1          1               1  6616
#>    ... 4 more rows.  get_results(x, what = "criteria")
#> 
#> -- classification --------------------------------------------------
#>  candidate n_profiles n_group_classes       level class n_modal
#>   classes2          2               1 individuals     1    1837
#>   classes2          2               1 individuals     2     745
#>   classes3          3               1 individuals     1    2054
#>   classes3          3               1 individuals     2     185
#>   classes3          3               1 individuals     3     343
#>  proportion_modal estimated_n estimated_proportion average_posterior
#>           0.71146      1841.1              0.71303            0.9998
#>           0.28854       740.9              0.28697            0.9942
#>           0.79551      2053.6              0.79535            0.9996
#>           0.07165       185.4              0.07182            0.9997
#>           0.13284       343.0              0.13283            0.9983
#>  odds_correct_classification relative_entropy
#>                       2505.7           0.9835
#>                        423.5           0.9835
#>                        573.9           0.9953
#>                      46383.9           0.9953
#>                       3920.3           0.9953
#> 
#> -- deltas ----------------------------------------------------------
#>  candidate n_profiles n_group_classes  criterion value delta_profile
#>   classes1          1               1 bic_groups  6584            NA
#>   classes2          2               1 bic_groups  6596        -12.42
#>   classes3          3               1 bic_groups  6624        -28.39
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
#> -- tree ------------------------------------------------------------
#>  tree       level from_candidate from_class from_share to_candidate to_class
#>  tree individuals       classes1          1      1.000     classes2        1
#>  tree individuals       classes1          1      1.000     classes2        2
#>  tree individuals       classes2          1      0.713     classes3        1
#>  tree individuals       classes2          2      0.287     classes3        1
#>  tree individuals       classes2          1      0.713     classes3        2
#>  tree individuals       classes2          2      0.287     classes3        2
#>  tree individuals       classes2          1      0.713     classes3        3
#>  tree individuals       classes2          2      0.287     classes3        3
#>  to_share    mass share_of_from share_of_to
#>   0.71303 0.71303       0.71303      1.0000
#>   0.28697 0.28697       0.28697      1.0000
#>   0.79535 0.55064       0.77225      0.6923
#>   0.79535 0.24471       0.85273      0.3077
#>   0.07182 0.05643       0.07914      0.7857
#>   0.07182 0.01539       0.05364      0.2143
#>   0.13283 0.10597       0.14861      0.7977
#>   0.13283 0.02687       0.09362      0.2023
#> 
#> 7 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
```
