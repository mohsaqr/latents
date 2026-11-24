# Get the converged candidate with the lowest value of each information criterion

`get_best()` returns, for each information criterion of an enumeration
grid, the converged candidate with the lowest value. Several criteria
may name the same candidate, and no criterion chooses a model on its
own.

## Usage

``` r
get_best(x)
```

## Arguments

- x:

  An enumeration grid from
  [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md),
  [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md)
  or
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md).
  A regression grid from
  [`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md)
  has no such table and is refused; its
  [`get_candidates()`](https://pak.dynasite.org/latents/reference/get_candidates.md)
  table marks the lowest BIC in `best_bic`.

## Value

A data frame.

On a profile or class grid it has one row per criterion and sample size,
with the columns:

- `criterion`: `aic`, `kic`, `bic`, `sabic`, `caic`, `awe`, `icl` or
  `clc`;

- `convention`: `"groups"` or `"individuals"`, the sample size the
  criterion uses; `NA` for AIC and KIC, which use none;

- `candidate`: the name of the converged candidate with the lowest
  value;

- `n_profiles`, `n_group_classes`, `model`: that candidate's numbers of
  classes and covariance structure;

- `value`: its value of the criterion.

On a transition grid or a grid of group-class families it is one row:
the row of `get_results(x, "candidates")` with the lowest BIC.

Errors: `latents_bad_argument` for a regression grid, naming the tables
it has; `latents_not_a_fit` when `x` is not an enumeration grid.

## See also

[`get_candidates()`](https://pak.dynasite.org/latents/reference/get_candidates.md)
for every candidate and criterion;
[`get_deltas()`](https://pak.dynasite.org/latents/reference/get_deltas.md),
[`get_knees()`](https://pak.dynasite.org/latents/reference/get_knees.md);
[`plot_criteria()`](https://pak.dynasite.org/latents/reference/plot_criteria.md).

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
