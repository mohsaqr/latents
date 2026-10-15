# Take one fitted model out of an enumeration grid or a comparison

An enumeration keeps every candidate it fitted, and
[`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)
keeps every model it compared. This returns one of them, so that a
caller who wants to plot, summarise or test a particular candidate names
it by its class counts instead of indexing the grid by position.

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

  Number of individual profiles identifying the candidate, the
  candidate's name (`"profiles3_groups4"`), or its two class counts
  (`c(3, 4)`). For a comparison, the label of the row, as its `model`
  column shows it.

- n_group_classes:

  Number of group classes identifying the candidate.

- model:

  The covariance model identifying the candidate, needed when
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  was given several `model`s and the class counts alone name several
  candidates. Counts that match more than one candidate raise
  `latents_unknown_candidate`, listing their names.

## Value

The fitted model for that cell of the grid: an object of class
`multilpa`, the same kind of object
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
returns, so the verbs for a fitted model apply to it.

## Details

A candidate is named by its name,
`get_candidate(grid, "classes3_groups4")`, as the grid's `candidate`
column shows it; by its class counts, `get_candidate(grid, 3, 4)` or
`get_candidate(grid, c(3, 4))`; or, on a grid that crosses covariance
structures, by its counts and `model`. From a comparison, name the row:
`get_candidate(comparison, "classes3_groups4")`. The same model is
`grid$classes3_groups4`.

`candidate_fit()` is the name this function had in latents 0.9.21. It is
still accepted with the same arguments and without a warning, so code
written for that release keeps working; it forwards to
`get_candidate()`, which new code uses.

## Conditions

`latents_unknown_candidate` for a name or class counts not in the grid,
for counts that match several candidates (the message lists their
names), and for a label that is not a row of a comparison;
`latents_failed_candidate` when the candidate is in the grid but its fit
raised an error, so no model exists to return; `latents_bad_argument`
for a name or pair of counts given together with `n_group_classes` or
`model`, or for a value that names no candidate;
`latents_incomplete_fit` for a comparison that no longer carries its
models. Counts that are not whole positive numbers raise an error from
the argument checks. None of these returns `NULL`, because a `NULL`
would flow silently into whatever the caller did next.

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
#> 
#>  profile          y   count proportion
#>        1 -0.1928460 44.5831   0.445831
#>        2  0.3516482 55.4169   0.554169
#> 
#> Variances and standard errors: get_results(x, "profiles"). 
#> Every other table: get_results(x, what = ), or get_results(x, "all").
# This grid crosses two structures, so the names carry them:
candidates
#> Class enumeration: 4 candidate models
#>              candidate n_profiles n_group_classes model log_likelihood
#>  profiles1_groups1_EEI          1               1   EEI      -130.6550
#>  profiles2_groups1_EEI          2               1   EEI      -130.5791
#>  profiles1_groups1_VVI          1               1   VVI      -130.6550
#>  profiles2_groups1_VVI          2               1   VVI      -130.5920
#>  n_parameters      aic bic_groups bic_individual icl_individual profile_entropy
#>             2 265.3100   265.9152       270.5204       270.5204              NA
#>             4 269.1581   270.3684       279.5788       355.6432      0.45131114
#>             2 265.3100   265.9152       270.5204       270.5204              NA
#>             5 271.1839   272.6968       284.2098       411.5897      0.08114782
#>  group_entropy avepp_min avepp_max avepp_mean avepp_below group_avepp_min
#>             NA        NA        NA         NA          NA              NA
#>             NA 0.6413011 0.8527078  0.7470045           1              NA
#>             NA        NA        NA         NA          NA              NA
#>             NA 0.6380900 0.6405463  0.6393181           2              NA
#>  group_avepp_max group_avepp_mean group_avepp_below converged boundary
#>               NA               NA                NA      TRUE    FALSE
#>               NA               NA                NA      TRUE    FALSE
#>               NA               NA                NA      TRUE    FALSE
#>               NA               NA                NA      TRUE    FALSE
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
#> 
#>  profile          y   count proportion
#>        1 -0.1928460 44.5831   0.445831
#>        2  0.3516482 55.4169   0.554169
#> 
#> Variances and standard errors: get_results(x, "profiles"). 
#> Every other table: get_results(x, what = ), or get_results(x, "all").
```
