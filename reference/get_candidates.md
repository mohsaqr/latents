# Get the table of every candidate model in an enumeration grid, with its criteria and diagnostics

`get_candidates()` returns the main table of an enumeration grid: one
row per candidate model, with its numbers of classes, its log
likelihood, its information criteria, its classification quality and
whether its estimation can be relied on.

## Usage

``` r
get_candidates(x)
```

## Arguments

- x:

  An enumeration grid of any kind, from
  [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md),
  [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md),
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  or
  [`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md).

## Value

A data frame with one row per candidate. The columns depend on the kind
of grid:

- a profile or class grid: `candidate` (its name), `n_profiles`,
  `n_group_classes`, `model` (the covariance structure),
  `log_likelihood`, `n_parameters`, the criteria (`aic`, `kic`, and
  `bic`, `sabic`, `caic`, `awe`, `icl` and `clc` with the suffix
  `_groups` or `_individual` for the sample size they use),
  `profile_entropy` and `group_entropy` (relative entropy at each
  level), `avepp_min`, `avepp_max`, `avepp_mean` and `avepp_below` (the
  average posterior probabilities of the classes and how many fall below
  the cutoff 0.8; with the prefix `group_` for the group classes),
  `converged`, `boundary`, `min_standardized_eigenvalue` and
  `min_gaussian_mass` (checks of Gaussian indicators) and
  `n_best_replicated` (the number of starts that reached the best log
  likelihood);

- a transition grid: `candidate`, `n_profiles`, `n_group_classes`,
  `log_likelihood`, `n_parameters`, `aic`, `bic`, `bic_individual`,
  `entropy`, `converged`, `boundary` and `n_best_replicated`;

- a grid of group-class families: `candidate`, `model`, `family`,
  `between_variance`, `n_group_classes`, `log_likelihood`,
  `n_parameters`, `aic`, `bic`, `bic_individual`, `entropy`,
  `smallest_class`, `min_effective_groups`, `converged`, `boundary` and
  `n_best_replicated`;

- a regression grid: `n_classes`, `n_group_classes`, `log_likelihood`,
  `n_parameters`, `aic`, `bic`, `bic_rows`, `sabic`, `icl`, `entropy`,
  `smallest_share`, `converged`, `n_best_replicated` and `best_bic`
  (`TRUE` for the candidate with the lowest BIC). It has no `candidate`
  or `boundary` column.

Errors: `latents_not_a_fit` when `x` is not an enumeration grid.

## Details

The table is `get_results(x, "candidates")` without the columns that are
empty in every row. On a class grid, for example, the covariance
structure `model` and the checks of Gaussian indicators are left out,
and the columns `warnings` and `error` appear only when a candidate has
a warning or failed.

## See also

[`get_best()`](https://pak.dynasite.org/latents/reference/get_best.md),
[`get_deltas()`](https://pak.dynasite.org/latents/reference/get_deltas.md)
and
[`get_knees()`](https://pak.dynasite.org/latents/reference/get_knees.md)
for what the criteria say;
[`plot_criteria()`](https://pak.dynasite.org/latents/reference/plot_criteria.md)
draws them;
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
for the full table.

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
