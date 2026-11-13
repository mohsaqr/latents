# Compare latent profile models

`enumerate_lpa()` fits latent profile models with different numbers of
profiles and specified covariance structures. It brings together
evidence about fit, classification, and estimation stability to support
selection of the number of profiles and their covariance specification.
Each fitted candidate remains available for detailed examination.

## Usage

``` r
enumerate_lpa(data, vars, n_profiles = 1:4, model = "basic", ...)
```

## Arguments

- data:

  A data frame with one row per observation.

- vars:

  Names of the indicator columns.

- n_profiles:

  Numbers of profiles to compare, specified as positive integers.

- model:

  The covariance structures: `"basic"`, the default, for the four that
  combine equal or varying variances with covariances absent or present
  (`EEI`, `VVI`, `EEE`, `VVV`); `"all"` for all 14; or their
  three-letter codes. See
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md).

- ...:

  Further arguments for
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  such as `missing`, `n_starts`, `seed`, `boundary` and `categorical`,
  except `id` and `n_group_classes`.

## Value

An object of class `multilpa_enumeration`, as
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
returns: [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
gives one row per candidate,
[`summary()`](https://rdrr.io/r/base/summary.html) the candidate
minimizing each criterion,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) the criteria
against the number of profiles, and
[`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md)
one fitted candidate.

## Details

This is the single-level interface to
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md).
For observations nested in groups, use
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
with the grouping column as `id`, which also varies the number of group
classes.

## Conditions

`latents_bad_argument` when `id` or `n_group_classes` is passed, or when
`model` names an unknown structure.

## See also

[`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md) to fit one
model;
[`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md)
for categorical indicators;
[`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
to test one number of profiles against another.

## Examples

``` r
models <- enumerate_lpa(srl, c("cognitive_strategies", "intrinsic_value",
                               "self_efficacy", "self_regulation",
                               "test_anxiety"),
                        n_profiles = 1:3, n_starts = 2, seed = 1)
summary(models)
#> Class enumeration: 12 candidates, 12 converged, 0 failed to fit
#> 2 distinct candidate(s) are minimal under some criterion.
#> No candidate is selected automatically. Choose one convention and keep it.
#> 
#> ! 1 of 12 candidates reached their best log likelihood from one start only
#>     (profiles3_EEE); refit those you consider with more starts.
#> ! 1 candidates have a class with average posterior below 0.8 (profiles3_EEE).
#> 
#> -- lowest value of each criterion --------------------------------
#>  criterion  convention     candidate value
#>        aic        <NA> profiles3_VVV  1193
#>        kic        <NA> profiles3_VVV  1258
#>        bic      groups profiles1_EEE  1345
#>        bic individuals profiles1_EEE  1345
#>      sabic      groups profiles3_VVV  1226
#>      sabic individuals profiles3_VVV  1226
#>       caic      groups profiles1_EEE  1365
#>       caic individuals profiles1_EEE  1365
#>        awe      groups profiles1_EEE  1520
#>        awe individuals profiles1_EEE  1520
#>        icl      groups profiles1_EEE  1345
#>        icl individuals profiles1_EEE  1345
#>        clc      groups profiles1_EEE  1231
#>        clc individuals profiles1_EEE  1231
#> 
#> -- knees of bic_groups ----------------------------------------
#>   direction                line knee_candidate knee_without_largest stable
#>    profiles  1 group class, EEE           <NA>                 <NA>     NA
#>    profiles  1 group class, EEI  profiles2_EEI                 <NA>     NA
#>    profiles  1 group class, VVI  profiles2_VVI                 <NA>     NA
#>    profiles  1 group class, VVV           <NA>                 <NA>     NA
#>  parameters all candidates, EEE           <NA>                 <NA>     NA
#>  parameters all candidates, EEI  profiles2_EEI                 <NA>     NA
#>  parameters all candidates, VVI  profiles2_VVI                 <NA>     NA
#>  parameters all candidates, VVV           <NA>                 <NA>     NA
#> 
#> Tables: get_candidates(x), get_best(x), get_deltas(x), get_knees(x),
#>         get_tree(x), get_classification(x), get_entropy(x), get_checks(x).
#> Plots:  plot_criteria(x), plot_deltas(x), plot_avepp(x), plot_tree(x), plot_all(x).
plot(models)
```
