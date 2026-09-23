# Summarise a class-enumeration grid

Reports, for every information criterion the grid carries, which
candidate minimises it. Criteria disagree by construction: they differ
in how they penalise parameters and in whether they count individuals or
independent groups. Laying the minima side by side shows that
disagreement instead of hiding it behind one default. Nothing here
selects a model.

## Usage

``` r
# S3 method for class 'multilpa_enumeration'
summary(object, ...)
```

## Arguments

- object:

  An `multilpa_enumeration` result from
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md).

- ...:

  Reserved for compatibility with
  [`summary()`](https://rdrr.io/r/base/summary.html).

## Value

An object of class `summary_multilpa_enumeration`, with a `print` method
and an [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
accessor. [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
returns the candidate grid. `get_results(summary(object), "criteria")`
gives one row per information criterion, with its sample-size
convention, the class counts and covariance structure of the minimising
candidate, and its value. Only converged candidates are eligible; a
criterion with no converged candidate has `NA` in the candidate and
value columns.

## Examples

``` r
set.seed(1)
d <- data.frame(g = rep(1:10, each = 10), y = rnorm(100))
candidates <- enumerate_classes(d, "y", "g", n_profiles = 1:2,
                              n_group_classes = 1, n_starts = 2, seed = 1)
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
as.data.frame(summary(candidates))
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
```
