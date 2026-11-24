# Get the relative entropy, from 0 to 1, of how sharply a model separates its classes

`get_entropy()` returns one number per classification level that says
how sharply the model separates its classes: the relative entropy, from
0 (every unit is equally likely to be in every class) to 1 (every unit
belongs to one class with certainty).

## Usage

``` r
get_entropy(x)
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

A data frame with one row per classification level (per candidate and
level for a grid) and the columns:

- `level`: `"individuals"` (the cases sorted into profiles or classes)
  or `"groups"` (the groups sorted into group classes);

- `n_classes`: the number of classes at that level;

- `n_units`: the number of classified units;

- `entropy_sum`: the sum of the units' entropies;

- `relative_entropy`: one minus the mean entropy of a unit divided by
  its largest possible value, the logarithm of `n_classes`. `NA` for a
  level with one class.

A profile, class or transition fit always has both rows; without group
classes its `"groups"` row has one class and `NA`. A fit of another
family has a row for each level it classifies.

Errors: `latents_not_a_fit` when `x` is none of these, for example a
comparison, pooled imputations or a data frame.

## See also

[`plot_entropy()`](https://pak.dynasite.org/latents/reference/plot_entropy.md)
draws the entropy of each unit;
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
