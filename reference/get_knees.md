# Get the candidate after which adding classes improves each information criterion much less

`get_knees()` returns the knee of every line of a profile or class grid:
the candidate after which adding classes lowers the criterion much less
than before, and whether the knee stays the same when the largest
candidate is left out. The criterion is `bic_groups`, BIC with the
number of groups as the sample size.

## Usage

``` r
get_knees(x)
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

A data frame with one row per line and the columns:

- `direction`: `"profiles"` (along the number of profiles, at a fixed
  number of group classes), `"group_classes"` (along the number of group
  classes, at a fixed number of profiles) or `"parameters"` (every
  candidate, ordered by its number of parameters);

- `line`: what is held fixed along the line;

- `criterion`: the name of the criterion, `bic_groups`;

- `from`, `to`: the first and the last candidate of the line;

- `n_points`: the number of candidates on the line;

- `knee_candidate`: the candidate at the knee; `NA` when the line has
  none;

- `distance`: how far the knee lies below the straight line, in the
  rescaled units;

- `knee_without_largest`: the knee when the largest candidate is left
  out;

- `stable`: `TRUE` when the two knees are the same candidate; `NA` when
  the line has fewer than four candidates;

- `note`: why a knee or its stability is missing.

Errors: `latents_bad_argument` for a transition, family or regression
grid, naming the tables it has; `latents_not_a_fit` when `x` is not an
enumeration grid.

## Details

Both axes of a line are rescaled to run from 0 to 1. The knee is the
candidate that lies furthest below the straight line from the first
candidate to the last. A line needs three candidates to have a knee and
four for its stability to be assessed, and a line whose criterion does
not fall from its first to its last candidate has none.

## See also

[`get_deltas()`](https://pak.dynasite.org/latents/reference/get_deltas.md)
for the changes behind the lines;
[`get_best()`](https://pak.dynasite.org/latents/reference/get_best.md);
[`plot_deltas()`](https://pak.dynasite.org/latents/reference/plot_deltas.md).

## Examples

``` r
grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
                      n_profiles = 1:3, model = "EEI",
                      n_starts = 2, seed = 1)
get_knees(grid)
#>    direction           line  criterion      from        to n_points
#> 1   profiles  1 group class bic_groups profiles1 profiles3        3
#> 2 parameters all candidates bic_groups profiles1 profiles3        3
#>   knee_candidate  distance knee_without_largest stable
#> 1      profiles2 0.3230334                 <NA>     NA
#> 2      profiles2 0.3230334                 <NA>     NA
#>                                       note
#> 1 three candidates: stability not assessed
#> 2 three candidates: stability not assessed
```
