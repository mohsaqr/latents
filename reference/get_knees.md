# Where each criterion line bends

The knee of each line of the criterion over the number of classes, and
whether it holds when the largest candidate is left out.

## Usage

``` r
get_knees(x)
```

## Arguments

- x:

  An enumeration grid from
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md),
  [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md)
  or
  [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md).

## Value

A data frame, one row per line and direction:
`direction, line, criterion, from, to, n_points, knee_candidate, distance, knee_without_largest, stable, note`.
A regression, transition or family grid has no knees and raises
`latents_bad_argument`, naming the tables it has. Any other object
raises `latents_not_a_fit`.

## See also

[`get_deltas()`](https://pak.dynasite.org/latents/reference/get_deltas.md),
[`get_best()`](https://pak.dynasite.org/latents/reference/get_best.md).

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
