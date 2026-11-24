# Get the average posterior probability of every class among the units assigned to each class

`get_avepp()` returns the average posterior probability (AvePP) matrix
in long form: for the units assigned to each class, their mean posterior
probability of every class.

## Usage

``` r
get_avepp(x)
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

A data frame with one row per level, assigned class and class (and per
candidate for a grid) and the columns:

- `level`: `"individuals"` or `"groups"`;

- `assigned_class`: the class the units are assigned to (their most
  likely class);

- `class`: the class whose probability is averaged;

- `n_assigned`: the number of units assigned to `assigned_class`;

- `average_posterior`: the mean posterior probability of `class` among
  those units.

The rows in which `assigned_class` equals `class` are the diagonal of
the matrix: values near 1 mean that the class is well separated.

Errors: `latents_not_a_fit` when `x` is none of these, for example a
comparison, pooled imputations or a data frame.

## See also

[`plot_avepp()`](https://pak.dynasite.org/latents/reference/plot_avepp.md)
draws the matrix;
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
