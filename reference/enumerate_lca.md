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
#> 
#> ! 2 of 3 candidates reached their best log likelihood from one start only
#>     (classes2, classes3); refit those you consider with more starts.
#> 
#> -- lowest value of each criterion --------------------------------
#>  criterion  convention candidate value
#>        aic        <NA>  classes2  6555
#>        kic        <NA>  classes2  6565
#>        bic      groups  classes1  6584
#>        bic individuals  classes1  6584
#>      sabic      groups  classes2  6574
#>      sabic individuals  classes2  6574
#>       caic      groups  classes1  6587
#>       caic individuals  classes1  6587
#>        awe      groups  classes1  6616
#>        awe individuals  classes1  6616
#>        icl      groups  classes1  6584
#>        icl individuals  classes1  6584
#>        clc      groups  classes1  6560
#>        clc individuals  classes1  6560
#> 
#> -- knees of bic_groups ----------------------------------------
#>   direction           line knee_candidate knee_without_largest stable
#>    profiles  1 group class           <NA>                 <NA>     NA
#>  parameters all candidates           <NA>                 <NA>     NA
#> 
#> Tables: get_candidates(x), get_best(x), get_deltas(x), get_knees(x),
#>         get_tree(x), get_classification(x), get_entropy(x), get_checks(x).
#> Plots:  plot_criteria(x), plot_deltas(x), plot_avepp(x), plot_tree(x), plot_all(x).
```
