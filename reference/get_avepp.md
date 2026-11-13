# Average posterior probability matrix

The average posterior probability of each class among the units assigned
to each class, at each level the fit classifies. For a fit, the table is
the one
[`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
computes with its defaults; for a grid, the tables of every candidate
are stacked with a first column `candidate` (a candidate whose
estimation failed has no row); for a diagnostics result, its own table.

## Usage

``` r
get_avepp(x)
```

## Arguments

- x:

  A fitted model of any family, an enumeration grid, or the result of
  [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md).

## Value

A data frame, one row per level, assigned class and class (and
candidate):
`level, assigned_class, class, n_assigned, average_posterior`. Any other
object raises `latents_not_a_fit`.

## See also

[`plot_avepp()`](https://pak.dynasite.org/latents/reference/plot_avepp.md)
draws it;
[`get_classification()`](https://pak.dynasite.org/latents/reference/get_classification.md),
[`get_entropy()`](https://pak.dynasite.org/latents/reference/get_entropy.md).

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
get_avepp(fit)
#>         level assigned_class class n_assigned average_posterior
#> 1 individuals              1     1        175      9.999434e-01
#> 2 individuals              1     2        175      5.663136e-05
#> 3 individuals              2     1         97      3.851364e-04
#> 4 individuals              2     2         97      9.996149e-01
```
