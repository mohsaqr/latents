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
[`candidate_fit()`](https://pak.dynasite.org/latents/reference/candidate_fit.md)
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
#> The candidates table shows 9 of its 26 columns; get_results(x, what = "candidates") returns all of them.
#> 
#> -- candidates ------------------------------------------------------
#>  n_profiles log_likelihood n_parameters  aic  bic icl_individual
#>           1          -3280            3 6566 6584           6584
#>           2          -3269            7 6552 6593           7454
#>           3          -3269           11 6560 6624           6651
#>  profile_entropy converged boundary
#>               NA      TRUE    FALSE
#>           0.7594      TRUE    FALSE
#>           0.9953      TRUE    FALSE
#> 
#> -- criteria --------------------------------------------------------
#>  criterion  convention n_profiles n_group_classes model value
#>        aic        <NA>          2               1  <NA>  6552
#>        kic        <NA>          2               1  <NA>  6562
#>        bic      groups          1               1  <NA>  6584
#>        bic individuals          1               1  <NA>  6584
#>      sabic      groups          2               1  <NA>  6571
#>      sabic individuals          2               1  <NA>  6571
#>       caic      groups          1               1  <NA>  6587
#>       caic individuals          1               1  <NA>  6587
#>        awe      groups          1               1  <NA>  6616
#>        awe individuals          1               1  <NA>  6616
#>    ... 4 more rows.  get_results(x, what = "criteria")
#> 
#> 2 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
```
