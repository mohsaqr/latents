# Class sizes and classification certainty

Modal and estimated class sizes, average posterior probability and odds
of correct classification, at each level the fit classifies. For a fit,
the table is the one
[`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
computes with its defaults; for a grid, the tables of every candidate
are stacked with a first column `candidate` (a candidate whose
estimation failed has no row); for a diagnostics result, its own table.

## Usage

``` r
get_classification(x)
```

## Arguments

- x:

  A fitted model of any family, an enumeration grid, or the result of
  [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md).

## Value

A data frame, one row per level and class (and candidate):
`level, class, n_modal, proportion_modal, estimated_n, estimated_proportion, average_posterior, odds_correct_classification`.
Any other object raises `latents_not_a_fit`.

## See also

[`plot_sizes()`](https://pak.dynasite.org/latents/reference/plot_sizes.md)
and
[`plot_posteriors()`](https://pak.dynasite.org/latents/reference/plot_posteriors.md)
draw it;
[`get_avepp()`](https://pak.dynasite.org/latents/reference/get_avepp.md).

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
get_classification(fit)
#>         level class n_modal proportion_modal estimated_n estimated_proportion
#> 1 individuals     1     175        0.6433824   175.02745            0.6434833
#> 2 individuals     2      97        0.3566176    96.97255            0.3565167
#>   average_posterior odds_correct_classification
#> 1         0.9999434                    9782.752
#> 2         0.9996149                    4684.632
```
