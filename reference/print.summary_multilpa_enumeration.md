# Print an enumeration summary

Print an enumeration summary

## Usage

``` r
# S3 method for class 'summary_multilpa_enumeration'
print(x, digits = 4L, rows = 10L, ...)
```

## Arguments

- x:

  A `summary_multilpa_enumeration` object.

- digits:

  Number of printed significant digits.

- rows:

  How many rows of each table to print. A longer table is shown to that
  depth, with its remaining row count and the
  [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
  call that returns it whole.

- ...:

  Passed to the underlying `data.frame` printing.

## Value

The summary, invisibly. Called for the side effect of printing the
candidate counts, the table of which candidate minimises each criterion,
and how many distinct candidates are minimal under some criterion.

## Examples

``` r
set.seed(1)
d <- data.frame(g = rep(seq_len(10), each = 10), y = rnorm(100))
candidates <- enumerate_classes(d, "y", "g", n_profiles = 1:2,
                                n_group_classes = 1, n_starts = 2, seed = 1)
print(summary(candidates))
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
