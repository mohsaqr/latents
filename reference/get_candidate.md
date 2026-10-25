# Take one fitted model out of an enumeration grid or a comparison

Returns one fitted model from an enumeration, selected by its name or by
its number of profiles and group classes (and its covariance structure
when the grid has several), or one model from a comparison, selected by
its row label.

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
  On a family grid, give the name, or `n_group_classes` with `model`.

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

The fitted model of that candidate, the same kind of object its fitting
function returns (for a profile grid, class `multilpa`), so the verbs
for a fitted model apply to it.

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
models. Counts that are not single positive whole numbers fail the
argument checks.

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
#> Converged: TRUE | iterations: 36 | best start: 1/2
#> Smallest standardized Gaussian eigenvalue: 0.8201 (profiles, profile_2); Gaussian mass: 55.42; weight effective n: 92.77.
#> 
#>  profile          y   count proportion
#>        1 -0.1928460 44.5831   0.445831
#>        2  0.3516482 55.4169   0.554169
#> 
#> Variances and standard errors: get_results(x, "profiles"). 
#> Every other table: get_results(x, what = ), or get_results(x, "all").
# This grid crosses two structures, so the names carry them:
candidates
#> Class enumeration: 4 candidate models (LPA)
#>              candidate n_profiles n_group_classes model log_likelihood
#>  profiles1_groups1_EEI          1               1   EEI      -130.6550
#>  profiles2_groups1_EEI          2               1   EEI      -130.5791
#>  profiles1_groups1_VVI          1               1   VVI      -130.6550
#>  profiles2_groups1_VVI          2               1   VVI      -130.5920
#>  n_parameters      aic bic_groups bic_individual icl_individual profile_entropy
#>             2 265.3100   265.9152       270.5204       270.5204              NA
#>             4 269.1581   270.3684       279.5788       355.4955      0.45237653
#>             2 265.3100   265.9152       270.5204       270.5204              NA
#>             5 271.1839   272.6968       284.2098       411.5897      0.08114782
#>  avepp_min avepp_max avepp_mean avepp_below converged boundary
#>         NA        NA         NA          NA      TRUE    FALSE
#>  0.6402627 0.8532616  0.7467621           1      TRUE    FALSE
#>         NA        NA         NA          NA      TRUE    FALSE
#>  0.6380900 0.6405463  0.6393181           2      TRUE    FALSE
#>  n_best_replicated
#>                  2
#>                  2
#>                  2
#>                  2
#> No candidate is selected automatically. Compare one convention consistently.
get_candidate(candidates, "profiles2_groups1_VVI")
#> Two-level latent profile analysis: 2 profiles, 1 group class
#> 100 individuals in 10 groups; varying diagonal residual covariance (VVI)
#> Log likelihood: -130.591962 | AIC: 271.184 | BIC (groups): 272.697
#> Converged: TRUE | iterations: 36 | best start: 1/2
#> Smallest standardized Gaussian eigenvalue: 0.8201 (profiles, profile_2); Gaussian mass: 55.42; weight effective n: 92.77.
#> 
#>  profile          y   count proportion
#>        1 -0.1928460 44.5831   0.445831
#>        2  0.3516482 55.4169   0.554169
#> 
#> Variances and standard errors: get_results(x, "profiles"). 
#> Every other table: get_results(x, what = ), or get_results(x, "all").
compared <- compare_models(candidates, "profiles1_groups1_VVI",
                           "profiles2_groups1_VVI")
get_candidate(compared, "profiles2_groups1_VVI")
#> Two-level latent profile analysis: 2 profiles, 1 group class
#> 100 individuals in 10 groups; varying diagonal residual covariance (VVI)
#> Log likelihood: -130.591962 | AIC: 271.184 | BIC (groups): 272.697
#> Converged: TRUE | iterations: 36 | best start: 1/2
#> Smallest standardized Gaussian eigenvalue: 0.8201 (profiles, profile_2); Gaussian mass: 55.42; weight effective n: 92.77.
#> 
#>  profile          y   count proportion
#>        1 -0.1928460 44.5831   0.445831
#>        2  0.3516482 55.4169   0.554169
#> 
#> Variances and standard errors: get_results(x, "profiles"). 
#> Every other table: get_results(x, what = ), or get_results(x, "all").
```
