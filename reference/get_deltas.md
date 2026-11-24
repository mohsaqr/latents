# Get how much each information criterion improves as a class or a group class is added

`get_deltas()` returns, for every candidate of a profile or class grid,
how much lower its information criterion is than that of the candidate
with one profile fewer, and than that of the candidate with one group
class fewer. The criterion is `bic_groups`, BIC with the number of
groups as the sample size.

## Usage

``` r
get_deltas(x)
```

## Arguments

- x:

  A profile or class grid from
  [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md),
  [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md)
  or
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md).
  Transition grids, grids of group-class families and regression grids
  have no such table and are refused.

## Value

A data frame with one row per candidate and the columns:

- `candidate`, `n_profiles`, `n_group_classes`, `model`: the candidate,
  its numbers of classes and its covariance structure;

- `criterion`: the name of the criterion, `bic_groups`;

- `value`: the candidate's value of the criterion;

- `delta_profile`: the criterion of the candidate with one profile fewer
  (and the same group classes and structure) minus `value`. A positive
  number means that this candidate is lower;

- `from_profile`: the name of that smaller candidate;

- `delta_group_class`, `from_group_class`: the same for the candidate
  with one group class fewer.

A delta is `NA` where there is no such smaller candidate, or where
either of the two did not converge.

Errors: `latents_bad_argument` for a transition, family or regression
grid, naming the tables it has; `latents_not_a_fit` when `x` is not an
enumeration grid.

## See also

[`plot_deltas()`](https://pak.dynasite.org/latents/reference/plot_deltas.md)
draws the changes;
[`get_knees()`](https://pak.dynasite.org/latents/reference/get_knees.md),
[`get_best()`](https://pak.dynasite.org/latents/reference/get_best.md).

## Examples

``` r
grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
                      n_profiles = 1:3, model = "EEI",
                      n_starts = 2, seed = 1)
get_deltas(grid)
#>   candidate n_profiles n_group_classes model  criterion    value delta_profile
#> 1 profiles1          1               1   EEI bic_groups 3055.835            NA
#> 2 profiles2          2               1   EEI bic_groups 2354.601     701.23422
#> 3 profiles3          3               1   EEI bic_groups 2322.969      31.63182
#>   from_profile delta_group_class from_group_class
#> 1         <NA>                NA             <NA>
#> 2    profiles1                NA             <NA>
#> 3    profiles2                NA             <NA>
```
