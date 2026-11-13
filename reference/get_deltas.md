# Change in the criterion with each added class

The decrease in the information criterion from the candidate with one
profile fewer, and from the one with one group class fewer (BIC over
groups).

## Usage

``` r
get_deltas(x)
```

## Arguments

- x:

  An enumeration grid from
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md),
  [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md)
  or
  [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md).

## Value

A data frame, one row per candidate:
`candidate, n_profiles, n_group_classes, model, criterion, value, delta_profile, from_profile, delta_group_class, from_group_class`.
A regression, transition or family grid has no deltas and raises
`latents_bad_argument`, naming the tables it has. Any other object
raises `latents_not_a_fit`.

## See also

[`plot_deltas()`](https://pak.dynasite.org/latents/reference/plot_deltas.md)
draws it;
[`get_knees()`](https://pak.dynasite.org/latents/reference/get_knees.md).

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
