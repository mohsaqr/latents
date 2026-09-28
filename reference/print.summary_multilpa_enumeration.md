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
```
