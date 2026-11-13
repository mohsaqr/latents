# Every candidate of a grid

One row per candidate with its class counts, log likelihood, every
information criterion and quality measure, and its convergence, boundary
and replication status. Columns empty in every row are left out (for a
class grid, the covariance `model` and the Gaussian checks; `error` when
every candidate fitted); `get_results(x, "candidates")` keeps them.

## Usage

``` r
get_candidates(x)
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

A data frame, one row per candidate. On a profile or class grid:
candidate, n_profiles, n_group_classes, model, log_likelihood,
n_parameters, the criteria (`aic`, `bic_groups`, ...),
`profile_entropy`, the AvePP summaries, `converged`, `boundary`,
`n_best_replicated`, and the Gaussian checks where they apply. A
transition or family grid has one `entropy`, `bic` and `bic_individual`,
and no AvePP summaries; a regression grid has `n_classes`, `bic`,
`bic_rows`, `sabic`, `icl`, `entropy`, `smallest_share` and `best_bic`,
without `candidate` or `boundary`. Any other object raises
`latents_not_a_fit`.

## See also

[`get_best()`](https://pak.dynasite.org/latents/reference/get_best.md),
[`get_deltas()`](https://pak.dynasite.org/latents/reference/get_deltas.md),
[`plot_criteria()`](https://pak.dynasite.org/latents/reference/plot_criteria.md);
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
for the columns.

## Examples

``` r
grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
                      n_profiles = 1:3, model = "EEI",
                      n_starts = 2, seed = 1)
get_candidates(grid)
#>   candidate n_profiles n_group_classes model log_likelihood n_parameters
#> 1 profiles1          1               1   EEI      -1516.706            4
#> 2 profiles2          2               1   EEI      -1157.680            7
#> 3 profiles3          3               1   EEI      -1133.455           10
#>        aic      kic bic_groups bic_individual sabic_groups sabic_individual
#> 1 3041.412 3048.412   3055.835       3055.835     3043.152         3043.152
#> 2 2329.360 2339.360   2354.601       2354.601     2332.406         2332.406
#> 3 2286.911 2299.911   2322.969       2322.969     2291.262         2291.262
#>   caic_groups caic_individual awe_groups awe_individual icl_groups
#> 1    3059.835        3059.835   3090.258       3090.258   3055.835
#> 2    2361.601        2361.601   2417.507       2417.507   2357.266
#> 3    2332.969        2332.969   2483.786       2483.786   2397.728
#>   icl_individual clc_groups clc_individual profile_entropy avepp_min avepp_max
#> 1       3055.835   3033.412       3033.412              NA        NA        NA
#> 2       2357.266   2318.026       2318.026       0.9929302 0.9962846 0.9999908
#> 3       2397.728   2341.670       2341.670       0.8749101 0.8617788 0.9993315
#>   avepp_mean avepp_below converged boundary min_standardized_eigenvalue
#> 1         NA          NA      TRUE    FALSE                  1.00000000
#> 2  0.9981377           0      TRUE    FALSE                  0.10241010
#> 3  0.9306180           0      TRUE    FALSE                  0.05961025
#>   min_gaussian_mass n_best_replicated
#> 1         272.00000                 2
#> 2          97.64931                 2
#> 3          46.26315                 2
```
