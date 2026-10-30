# Retrieve a fitted model from an enumeration or comparison

`get_candidate()` retrieves a specified fitted model from an enumeration
grid or model comparison for detailed inspection and subsequent
analysis. The returned object retains its original model class and
supports the summary, diagnostic, plotting, and inference methods
available for that model.

A candidate can be identified by its name or, for supported enumeration
grids, by its class counts and model specification. When several
candidates share the same class counts, the covariance structure or
candidate name distinguishes them. The function retrieves the stored fit
without repeating estimation.

## Usage

``` r
get_candidate(x, n_profiles, n_group_classes = 1L, model = NULL)

candidate_fit(x, n_profiles, n_group_classes = 1L, model = NULL)
```

## Arguments

- x:

  An enumeration grid from
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  (a profile, latent class, transition or group-class family grid), or a
  comparison from
  [`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md).
  On a family grid, give the name, or `n_group_classes` with `model`. A
  grid of mixture regressions from
  [`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md)
  is not supported (`latents_bad_argument`):
  `get_results(grid, "model", n_classes = )` returns one of its models.

- n_profiles:

  Number of individual profiles of the candidate, the candidate's name
  (`"profiles3_groups4"`), or its two class counts (`c(3, 4)`). For a
  comparison, the row label, as its `model` column shows it.

- n_group_classes:

  Number of group classes of the candidate. Default 1.

- model:

  The covariance structure of the candidate, needed when
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  crossed several structures and the class counts match several
  candidates; on a group-class family grid, the family model as the
  grid's `model` column shows it.

## Value

The fitted object for the selected candidate, with the class and methods
of its fitting function (for example, `multilpa` for a profile model).

## Details

A candidate is selected by its name, as the grid's `candidate` column
shows it (`get_candidate(grid, "profiles3_groups4")`; `classes` instead
of `profiles` when every indicator is categorical); by its class counts,
`get_candidate(grid, 3, 4)` or `get_candidate(grid, c(3, 4))`; or, on a
grid that crosses covariance structures, by its counts and `model`. From
a comparison, give the row label:
`get_candidate(comparison, "profiles3_groups4")`.

`candidate_fit()` is the former name of `get_candidate()`, kept for
compatibility. It takes the same arguments, forwards to
`get_candidate()` and gives no warning.

## Conditions

`latents_unknown_candidate` for a name or class counts not in the grid,
for counts that match several candidates (the message lists their
names), and for a label that is not a row of a comparison;
`latents_failed_candidate` when the candidate's fit raised an error, so
there is no model to return; `latents_bad_argument` for a name or pair
of counts given together with `n_group_classes` or `model`, for a value
that names no candidate, and for a comparison given no label;
`latents_incomplete_fit` for a comparison that no longer carries its
models. Counts that are not single positive integers fail the argument
checks.

## Examples

``` r
set.seed(1)
d <- data.frame(g = rep(1:10, each = 10), y = rnorm(100))
candidates <- enumerate_classes(d, "y", "g", n_profiles = 1:2,
                              n_group_classes = 1, n_starts = 2, seed = 1)
get_candidate(candidates, n_profiles = 2, n_group_classes = 1, model = "VVI")
#> Two-level latent profile analysis: 2 profiles, 1 group class
#> 100 individuals in 10 groups; varying diagonal residual covariance (VVI)
#> Log likelihood: -130.591962 | AIC: 271.184 | BIC (groups): 272.697
#> Converged: TRUE | iterations: 34 | best log likelihood reached by 2 of 2 starts (2 converged)
#> Smallest standardized Gaussian eigenvalue: 0.8201 (profiles, profile_2); Gaussian mass: 55.42; weight effective n: 92.77.
#> 
#> Class sizes   profile_1 44.6% (45)   profile_2 55.4% (55)
#> 
#> Means
#>  indicator profile_1 profile_2
#>          y    -0.193     0.352
#> 
#> Variances and standard errors: get_results(x, "profiles").
#> Every other table: get_results(x, what = ), or get_results(x, "all").
# This grid crosses two structures, so the names carry them:
candidates
#> Class enumeration: 4 candidate models (LPA)
#>              candidate n_profiles model log_likelihood n_parameters    aic
#>  profiles1_groups1_EEI          1   EEI        -130.66            2 265.31
#>  profiles2_groups1_EEI          2   EEI        -130.58            4 269.16
#>  profiles1_groups1_VVI          1   VVI        -130.66            2 265.31
#>  profiles2_groups1_VVI          2   VVI        -130.59            5 271.18
#>  bic_groups bic_individual icl_individual profile_entropy avepp_min avepp_max
#>      265.92         270.52         270.52              NA        NA        NA
#>      270.37         279.58         355.58           0.452     0.641     0.853
#>      265.92         270.52         270.52              NA        NA        NA
#>      272.70         284.21         411.59           0.081     0.638     0.641
#>  avepp_mean avepp_below converged boundary n_best_replicated
#>          NA          NA      TRUE    FALSE                 2
#>       0.747           1      TRUE    FALSE                 2
#>          NA          NA      TRUE    FALSE                 2
#>       0.639           2      TRUE    FALSE                 2
#> Every column, unrounded: as.data.frame(x).
#> No candidate is selected automatically. Compare one convention consistently.
get_candidate(candidates, "profiles2_groups1_VVI")
#> Two-level latent profile analysis: 2 profiles, 1 group class
#> 100 individuals in 10 groups; varying diagonal residual covariance (VVI)
#> Log likelihood: -130.591962 | AIC: 271.184 | BIC (groups): 272.697
#> Converged: TRUE | iterations: 34 | best log likelihood reached by 2 of 2 starts (2 converged)
#> Smallest standardized Gaussian eigenvalue: 0.8201 (profiles, profile_2); Gaussian mass: 55.42; weight effective n: 92.77.
#> 
#> Class sizes   profile_1 44.6% (45)   profile_2 55.4% (55)
#> 
#> Means
#>  indicator profile_1 profile_2
#>          y    -0.193     0.352
#> 
#> Variances and standard errors: get_results(x, "profiles").
#> Every other table: get_results(x, what = ), or get_results(x, "all").
compared <- compare_models(candidates, "profiles1_groups1_VVI",
                           "profiles2_groups1_VVI")
get_candidate(compared, "profiles2_groups1_VVI")
#> Two-level latent profile analysis: 2 profiles, 1 group class
#> 100 individuals in 10 groups; varying diagonal residual covariance (VVI)
#> Log likelihood: -130.591962 | AIC: 271.184 | BIC (groups): 272.697
#> Converged: TRUE | iterations: 34 | best log likelihood reached by 2 of 2 starts (2 converged)
#> Smallest standardized Gaussian eigenvalue: 0.8201 (profiles, profile_2); Gaussian mass: 55.42; weight effective n: 92.77.
#> 
#> Class sizes   profile_1 44.6% (45)   profile_2 55.4% (55)
#> 
#> Means
#>  indicator profile_1 profile_2
#>          y    -0.193     0.352
#> 
#> Variances and standard errors: get_results(x, "profiles").
#> Every other table: get_results(x, what = ), or get_results(x, "all").
```
