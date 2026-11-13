# How the classes of each candidate split into the next

The flows
[`plot_tree()`](https://pak.dynasite.org/latents/reference/plot_tree.md)
draws: for each pair of consecutive candidates of a tree, how much of
each class goes to each class of the next.

## Usage

``` r
get_tree(x)
```

## Arguments

- x:

  An enumeration grid from
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md),
  [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md),
  [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md)
  or
  [`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md).

## Value

A data frame, one row per flow:
`tree, level, from_candidate, from_class, from_share, to_candidate, to_class, to_share, mass, share_of_from, share_of_to`.
A grid without a tree raises `latents_bad_argument`, naming the tables
it has. Any other object raises `latents_not_a_fit`.

## See also

[`plot_tree()`](https://pak.dynasite.org/latents/reference/plot_tree.md)
draws it.

## Examples

``` r
grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
                      n_profiles = 1:3, model = "EEI",
                      n_starts = 2, seed = 1)
get_tree(grid)
#>   tree       level from_candidate from_class from_share to_candidate to_class
#> 1 tree individuals      profiles1          1  1.0000000    profiles2        1
#> 2 tree individuals      profiles1          1  1.0000000    profiles2        2
#> 3 tree individuals      profiles2          1  0.6409952    profiles3        1
#> 4 tree individuals      profiles2          2  0.3590048    profiles3        1
#> 5 tree individuals      profiles2          1  0.6409952    profiles3        2
#> 6 tree individuals      profiles2          2  0.3590048    profiles3        2
#> 7 tree individuals      profiles2          1  0.6409952    profiles3        3
#> 8 tree individuals      profiles2          2  0.3590048    profiles3        3
#>    to_share         mass share_of_from  share_of_to
#> 1 0.6409952 6.409952e-01  6.409952e-01 1.000000e+00
#> 2 0.3590048 3.590048e-01  3.590048e-01 1.000000e+00
#> 3 0.1700851 1.674664e-01  2.612601e-01 9.846037e-01
#> 4 0.1700851 2.618684e-03  7.294286e-03 1.539631e-02
#> 5 0.3563994 1.362744e-05  2.125982e-05 3.823643e-05
#> 6 0.3563994 3.563858e-01  9.927047e-01 9.999618e-01
#> 7 0.4735155 4.735151e-01  7.387187e-01 9.999992e-01
#> 8 0.4735155 3.778992e-07  1.052630e-06 7.980715e-07
```
