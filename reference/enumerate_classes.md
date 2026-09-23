# Enumerate numbers of individual profiles and group classes

Fits every requested combination and retains errors and warnings
alongside successful fits. Information criteria are descriptive: neither
the smallest BIC nor high entropy guarantees the correct number of
classes. Unconverged, boundary, and unreplicated fits are reported, not
silently selected.

## Usage

``` r
enumerate_classes(
  data,
  vars,
  id,
  n_profiles = 1:4,
  n_group_classes = 1:3,
  model = NULL,
  seed = NULL,
  ...
)
```

## Arguments

- data:

  Data frame.

- vars:

  Continuous indicator names.

- id:

  Group identifier column name.

- n_profiles:

  Positive integer profile counts to try.

- n_group_classes:

  Positive integer group-class counts to try.

- model:

  Covariance models to try, as the three-letter codes
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)'s
  `volume`, `shape` and `orientation` name: any of `"EII"`, `"VII"`,
  `"EEI"`, `"VEI"`, `"EVI"`, `"VVI"`, `"EEE"`, `"VEE"`, `"EVE"`,
  `"VVE"`, `"EEV"`, `"VEV"`, `"EVV"`, `"VVV"`. Naming them crosses the
  models with the class counts, which is the grid model-based clustering
  is usually selected over, and adds a `model` column to the candidate
  table. `NULL`, the default, fits whatever the other arguments already
  asked for, and the grid is the one this verb has always fitted.

- seed:

  Optional reproducible seed for each fit.

- ...:

  Further arguments to
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).

## Value

An object of class `multilpa_enumeration`. Read it with the verbs that
describe it rather than by reaching into it:
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives one
row per candidate model with every criterion and diagnostic,
[`summary()`](https://rdrr.io/r/base/summary.html) gives one row per
information criterion naming the candidate that minimises it,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws one
criterion across the grid, and
[`candidate_fit()`](https://pak.dynasite.org/latents/reference/candidate_fit.md)
returns the fitted model for one cell of the grid.

## See also

[`candidate_fit()`](https://pak.dynasite.org/latents/reference/candidate_fit.md)
to take one fitted model out of the grid,
[`summary.multilpa_enumeration()`](https://pak.dynasite.org/latents/reference/summary.multilpa_enumeration.md)
for the criterion-by-criterion comparison.

## Examples

``` r
set.seed(1)
d <- data.frame(g = rep(1:10, each = 10), y = rnorm(100))
candidates <- enumerate_classes(d, "y", "g", n_profiles = 1:2,
                              n_group_classes = 1, n_starts = 2, seed = 1)
as.data.frame(candidates)
#>   n_profiles n_group_classes model log_likelihood n_parameters      aic
#> 1          1               1   VVI      -130.6550            2 265.3100
#> 2          2               1   VVI      -130.5921            5 271.1842
#>        kic bic_groups bic_individual sabic_groups sabic_individual caic_groups
#> 1 270.3100   265.9152       270.5204     259.9237         264.2039    267.9152
#> 2 279.1842   272.6971       284.2100     257.7185         268.4188    277.6971
#>   caic_individual awe_groups awe_individual icl_groups icl_individual
#> 1        272.5204   276.5204       285.7307   265.9152       270.5204
#> 2        289.2100   299.2100       450.3309   272.6971       412.3050
#>   clc_groups clc_individual profile_entropy group_entropy converged boundary
#> 1   261.3100       261.3100              NA            NA      TRUE    FALSE
#> 2   261.1842       389.2792      0.07599015            NA      TRUE    FALSE
#>   n_best_replicated warnings error
#> 1                 2           <NA>
#> 2                 2           <NA>
summary(candidates)
#> Class enumeration: 2 candidates, 2 converged, 0 failed to fit
#> 2 distinct candidate(s) are minimal under some criterion.
#> No candidate is selected automatically. Choose one convention and keep it.
#> 
#> -- candidates ------------------------------------------------------
#>  n_profiles n_group_classes model log_likelihood n_parameters   aic   kic
#>           1               1   VVI         -130.7            2 265.3 270.3
#>           2               1   VVI         -130.6            5 271.2 279.2
#>  bic_groups bic_individual sabic_groups sabic_individual caic_groups
#>       265.9          270.5        259.9            264.2       267.9
#>       272.7          284.2        257.7            268.4       277.7
#>  caic_individual awe_groups awe_individual icl_groups icl_individual clc_groups
#>            272.5      276.5          285.7      265.9          270.5      261.3
#>            289.2      299.2          450.3      272.7          412.3      261.2
#>  clc_individual profile_entropy group_entropy converged boundary
#>           261.3              NA            NA      TRUE    FALSE
#>           389.3         0.07599            NA      TRUE    FALSE
#>  n_best_replicated warnings error
#>                  2           <NA>
#>                  2           <NA>
#> 
#> -- criteria --------------------------------------------------------
#>  criterion  convention n_profiles n_group_classes model value
#>        aic        <NA>          1               1   VVI 265.3
#>        kic        <NA>          1               1   VVI 270.3
#>        bic      groups          1               1   VVI 265.9
#>        bic individuals          1               1   VVI 270.5
#>      sabic      groups          2               1   VVI 257.7
#>      sabic individuals          1               1   VVI 264.2
#>       caic      groups          1               1   VVI 267.9
#>       caic individuals          1               1   VVI 272.5
#>        awe      groups          1               1   VVI 276.5
#>        awe individuals          1               1   VVI 285.7
#>    ... 4 more rows.  get_results(x, what = "criteria")
#> 
#> 2 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
```
