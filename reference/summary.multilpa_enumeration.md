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
#> Class enumeration: 4 candidates, 4 converged, 0 failed to fit
#> 3 distinct candidate(s) are minimal under some criterion.
#> No candidate is selected automatically. Choose one convention and keep it.
#> The candidates table shows 11 of its 26 columns; get_results(x, what = "candidates") returns all of them.
#> 
#> -- candidates ------------------------------------------------------
#>  n_profiles model log_likelihood n_parameters   aic bic_groups bic_individual
#>           1   EEI         -130.7            2 265.3      265.9          270.5
#>           2   EEI         -130.6            4 269.2      270.4          279.6
#>           1   VVI         -130.7            2 265.3      265.9          270.5
#>           2   VVI         -130.6            5 271.2      272.7          284.2
#>  icl_individual profile_entropy converged boundary
#>           270.5              NA      TRUE    FALSE
#>           355.6         0.45131      TRUE    FALSE
#>           270.5              NA      TRUE    FALSE
#>           411.6         0.08115      TRUE    FALSE
#> 
#> -- criteria --------------------------------------------------------
#>  criterion  convention n_profiles n_group_classes model value
#>        aic        <NA>          1               1   EEI 265.3
#>        kic        <NA>          1               1   EEI 270.3
#>        bic      groups          1               1   EEI 265.9
#>        bic individuals          1               1   EEI 270.5
#>      sabic      groups          2               1   VVI 257.7
#>      sabic individuals          1               1   EEI 264.2
#>       caic      groups          1               1   EEI 267.9
#>       caic individuals          1               1   EEI 272.5
#>        awe      groups          1               1   EEI 276.5
#>        awe individuals          1               1   EEI 285.7
#>    ... 4 more rows.  get_results(x, what = "criteria")
#> 
#> 2 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
as.data.frame(summary(candidates))
#>   n_profiles n_group_classes model log_likelihood n_parameters      aic
#> 1          1               1   EEI      -130.6550            2 265.3100
#> 2          2               1   EEI      -130.5791            4 269.1581
#> 3          1               1   VVI      -130.6550            2 265.3100
#> 4          2               1   VVI      -130.5920            5 271.1839
#>        kic bic_groups bic_individual sabic_groups sabic_individual caic_groups
#> 1 270.3100   265.9152       270.5204     259.9237         264.2039    267.9152
#> 2 276.1581   270.3684       279.5788     258.3855         266.9458    274.3684
#> 3 270.3100   265.9152       270.5204     259.9237         264.2039    267.9152
#> 4 279.1839   272.6968       284.2098     257.7182         268.4185    277.6968
#>   caic_individual awe_groups awe_individual icl_groups icl_individual
#> 1        272.5204   276.5204       285.7307   265.9152       270.5204
#> 2        283.5788   291.5788       386.0639   270.3684       355.6432
#> 3        272.5204   276.5204       285.7307   265.9152       270.5204
#> 4        289.2098   299.2098       449.6156   272.6968       411.5897
#>   clc_groups clc_individual profile_entropy group_entropy converged boundary
#> 1   261.3100       261.3100              NA            NA      TRUE    FALSE
#> 2   261.1581       337.2225      0.45131114            NA      TRUE    FALSE
#> 3   261.3100       261.3100              NA            NA      TRUE    FALSE
#> 4   261.1839       388.5639      0.08114782            NA      TRUE    FALSE
#>   n_best_replicated warnings error
#> 1                 2           <NA>
#> 2                 2           <NA>
#> 3                 2           <NA>
#> 4                 2           <NA>
```
