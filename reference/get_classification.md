# Get how many units each class holds and how certain their assignment is

`get_classification()` returns, for every class, how many units it
holds, counted two ways, and how certain the assignment of its units is.

## Usage

``` r
get_classification(x)
```

## Arguments

- x:

  One of:

  - a fitted model of any family: from
    [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
    [`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
    [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md),
    [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) or
    [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md).
    The table is the one
    [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
    computes with its defaults;

  - an enumeration grid of any kind, from
    [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md),
    [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md),
    [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
    or
    [`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md).
    The tables of the fitted candidates are stacked, with a first column
    `candidate` that names the candidate of each row. A candidate whose
    estimation failed has no rows;

  - the result of
    [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md).
    The table is the one the result holds.

## Value

A data frame with one row per level and class (per candidate, level and
class for a grid) and the columns:

- `level`: `"individuals"` or `"groups"`;

- `class`: the number of the class; `0` is the noise component of a fit
  that has one;

- `n_modal`: the number of units whose most likely class this is;

- `proportion_modal`: `n_modal` as a share of the units;

- `estimated_n`: the sum of the units' posterior probabilities of the
  class;

- `estimated_proportion`: `estimated_n` as a share of the units;

- `average_posterior`: the mean posterior probability of the class among
  the units assigned to it;

- `odds_correct_classification`: the odds of `average_posterior` divided
  by the odds of `estimated_proportion`.

A fit without group classes has rows for the individuals only; a fit
that classifies groups only has rows for the groups only.

Errors: `latents_not_a_fit` when `x` is none of these, for example a
comparison, pooled imputations or a data frame.

## See also

[`plot_sizes()`](https://pak.dynasite.org/latents/reference/plot_sizes.md)
and
[`plot_posteriors()`](https://pak.dynasite.org/latents/reference/plot_posteriors.md)
draw it;
[`get_avepp()`](https://pak.dynasite.org/latents/reference/get_avepp.md)
for the full matrix of average posterior probabilities.

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
