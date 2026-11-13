# Relative entropy of the classification

The entropy of the posterior class probabilities relative to a flat
posterior, at each level the fit classifies. For a fit, the table is the
one
[`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
computes with its defaults; for a grid, the tables of every candidate
are stacked with a first column `candidate` (a candidate whose
estimation failed has no row); for a diagnostics result, its own table.

## Usage

``` r
get_entropy(x)
```

## Arguments

- x:

  A fitted model of any family, an enumeration grid, or the result of
  [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md).

## Value

A data frame, one row per classification level (and candidate):
`level, n_classes, n_units, entropy_sum, relative_entropy`. Any other
object raises `latents_not_a_fit`.

## See also

[`plot_entropy()`](https://pak.dynasite.org/latents/reference/plot_entropy.md)
draws it;
[`get_classification()`](https://pak.dynasite.org/latents/reference/get_classification.md),
[`get_avepp()`](https://pak.dynasite.org/latents/reference/get_avepp.md).

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
get_entropy(fit)
#>         level n_classes n_units entropy_sum relative_entropy
#> 1 individuals         2     272   0.2262339           0.9988
#> 2      groups         1     272   0.0000000               NA
```
