# The best candidate under each criterion

The candidate with the lowest (or, by the criterion's convention,
highest) value of each information criterion. Several criteria may name
the same candidate; none selects a model on its own.

## Usage

``` r
get_best(x)
```

## Arguments

- x:

  An enumeration grid from
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md),
  [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md)
  or
  [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md).

## Value

A data frame. On a profile or class grid, one row per criterion and
convention:
`criterion, convention, candidate, n_profiles, n_group_classes, model, value`.
On a transition or family grid, the one row of the candidates table with
the lowest BIC. A regression grid raises `latents_bad_argument`, naming
the tables it has. Any other object raises `latents_not_a_fit`.

## See also

[`get_candidates()`](https://pak.dynasite.org/latents/reference/get_candidates.md),
[`get_deltas()`](https://pak.dynasite.org/latents/reference/get_deltas.md),
[`get_knees()`](https://pak.dynasite.org/latents/reference/get_knees.md).

## Examples

``` r
grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
                      n_profiles = 1:3, model = "EEI",
                      n_starts = 2, seed = 1)
get_best(grid)
#>    criterion  convention candidate n_profiles n_group_classes model    value
#> 1        aic        <NA> profiles3          3               1   EEI 2286.911
#> 2        kic        <NA> profiles3          3               1   EEI 2299.911
#> 3        bic      groups profiles3          3               1   EEI 2322.969
#> 4        bic individuals profiles3          3               1   EEI 2322.969
#> 5      sabic      groups profiles3          3               1   EEI 2291.262
#> 6      sabic individuals profiles3          3               1   EEI 2291.262
#> 7       caic      groups profiles3          3               1   EEI 2332.969
#> 8       caic individuals profiles3          3               1   EEI 2332.969
#> 9        awe      groups profiles2          2               1   EEI 2417.507
#> 10       awe individuals profiles2          2               1   EEI 2417.507
#> 11       icl      groups profiles2          2               1   EEI 2357.266
#> 12       icl individuals profiles2          2               1   EEI 2357.266
#> 13       clc      groups profiles2          2               1   EEI 2318.026
#> 14       clc individuals profiles2          2               1   EEI 2318.026
```
