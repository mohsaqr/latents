# Fit latent class models with different numbers of classes, and compare them

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

  Numbers of classes to compare, specified as positive integers. Default
  `1:4`.

- ...:

  Further arguments for
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md),
  used for every candidate, such as `missing` (`"error"` or `"fiml"`),
  `n_starts`, `seed` and `min_probability`. `id` and `n_group_classes`
  are not accepted, and neither is `model`: a latent class model has no
  covariance structure. Every indicator is categorical unless named in
  `ordinal` or `count`. `categorical` is accepted but has no further
  effect; its names must be in `vars`.

## Value

An object of class `multilpa_enumeration`, as
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
returns.
[`get_candidates()`](https://pak.dynasite.org/latents/reference/get_candidates.md)
(also [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html))
gives one row per candidate with its log likelihood, information
criteria, entropy, convergence and replication;
[`get_best()`](https://pak.dynasite.org/latents/reference/get_best.md)
and [`summary()`](https://rdrr.io/r/base/summary.html) the candidate
with the lowest value of each criterion;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) the criteria
against the number of classes; and
[`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md)
one fitted candidate. The number of classes is in the `n_profiles`
column of the candidates table, and
[`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md)
selects it through its `n_profiles` argument. A candidate that could not
be fitted keeps its row, with `NA` values.

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
`vars`. A candidate with more free parameters than the table of response
patterns has free cells warns with class `latents_unidentified`. A
`latents_single_start` message names the candidates whose best log
likelihood was reached from one start only.

## See also

[`lca()`](https://pak.dynasite.org/latents/reference/lca.md) to fit one
model;
[`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md)
for continuous indicators.

## Examples

``` r
models <- enumerate_lca(student_esm, c("time_with_friends",
                                       "on_social_media", "tv_video_games",
                                       "sports"),
                        n_classes = 1:3, n_starts = 2, seed = 1)
#> 1 of 3 candidates reached their best log likelihood from one start only (classes3); refit those you consider with more starts.
summary(models)
#> Class enumeration: 3 candidates, 3 converged, 0 failed to fit
#> 3 distinct candidate(s) are minimal under some criterion.
#> No candidate is selected automatically. Choose one convention and keep it.
#> 
#> ! 1 of 3 candidates reached their best log likelihood from one start only
#>     (classes3); refit those you consider with more starts.
#> ! 1 candidates have a class with average posterior below 0.8 (classes2).
#> 
#> -- lowest value of each criterion --------------------------------
#>  criterion  convention candidate value
#>        aic        <NA>  classes2  9138
#>        kic        <NA>  classes2  9150
#>        bic      groups  classes2  9191
#>        bic individuals  classes2  9191
#>      sabic      groups  classes2  9163
#>      sabic individuals  classes2  9163
#>       caic      groups  classes2  9200
#>       caic individuals  classes2  9200
#>        awe      groups  classes1  9307
#>        awe individuals  classes1  9307
#>        icl      groups  classes3  9249
#>        icl individuals  classes3  9249
#>        clc      groups  classes3  9139
#>        clc individuals  classes3  9139
#> 
#> -- knees of bic_groups ----------------------------------------
#>   direction           line knee_candidate knee_without_largest stable
#>    profiles  1 group class       classes2                 <NA>     NA
#>  parameters all candidates       classes2                 <NA>     NA
#> 
#> Tables: get_candidates(x), get_best(x), get_deltas(x), get_knees(x),
#>         get_tree(x), get_classification(x), get_entropy(x), get_checks(x).
#> Plots:  plot_criteria(x), plot_deltas(x), plot_avepp(x), plot_tree(x), plot_all(x).
```
