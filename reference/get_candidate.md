# Get one fitted model out of an enumeration grid or a model comparison

`get_candidate()` returns one fitted model from an enumeration grid or a
model comparison, for detailed inspection and further analysis. The
returned object keeps its original model class and supports the summary,
diagnostic, plotting and inference methods available for that model.

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
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md),
  [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md)
  or
  [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md)
  (a profile, latent class, transition or group-class family grid), or a
  comparison from
  [`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md).
  On a family grid, give the name, or `n_group_classes` with `model`. A
  grid of mixture regressions from
  [`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md)
  is not supported (`latents_bad_argument`):
  `get_results(grid, "model", n_classes = )` returns one of its models.

- n_profiles:

  Which candidate, in one of three forms:

  - a single number: the candidate's number of profiles (classes), used
    together with `n_group_classes` and, where needed, `model`;

  - a name: the candidate's name as the grid's `candidate` column shows
    it, such as `"profiles3_groups4"`;

  - two numbers: the number of profiles and the number of group classes,
    such as `c(3, 4)`.

  A name or a pair of numbers identifies the candidate on its own, so
  `n_group_classes` and `model` are then left out. For a comparison,
  this is the row label, as its `model` column shows it.

- n_group_classes:

  Number of group classes of the candidate. Default 1, which is also the
  value for every candidate of a single-level grid.

- model:

  `NULL` (default), or the covariance structure of the candidate as a
  single code such as `"VVI"`, needed when the grid crossed several
  structures and the class counts match several candidates. On a
  group-class family grid, the family model as the grid's `model` column
  shows it, such as `"additive_equal"`.

## Value

The fitted model of the selected candidate, with the class its fitting
function gives it: `multilpa` or `multilpa_covariates` from a profile or
class grid, `multilpa_transitions` or `multilpa_lta` from a transition
grid, `multilpa_additive` from a family grid. Every method and accessor
of that class applies. The stored fit is returned; nothing is estimated
again.

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
models. Counts that are not single positive integers, and a `model` that
is not a single string, also raise `latents_bad_argument`.

## Examples

``` r
# This grid crosses two structures, so the names carry them:
grid <- enumerate_lpa(faithful, c("eruptions", "waiting"),
                      n_profiles = 1:2, model = c("EEI", "VVI"),
                      n_starts = 2, seed = 1)
grid
#> Class enumeration: 4 candidate models (LPA)
#>      candidate n_profiles model log_likelihood n_parameters     aic bic_groups
#>  profiles1_EEI          1   EEI       -1516.71            4 3041.41    3055.83
#>  profiles2_EEI          2   EEI       -1157.68            7 2329.36    2354.60
#>  profiles1_VVI          1   VVI       -1516.71            4 3041.41    3055.83
#>  profiles2_VVI          2   VVI       -1147.81            9 2313.61    2346.06
#>  icl_individual profile_entropy avepp_min avepp_max avepp_mean avepp_below
#>         3055.83              NA        NA        NA         NA          NA
#>         2357.27           0.993     0.996         1      0.998           0
#>         3055.83              NA        NA        NA         NA          NA
#>         2346.52           0.999     1.000         1      1.000           0
#>  converged boundary n_best_replicated
#>       TRUE    FALSE                 2
#>       TRUE    FALSE                 2
#>       TRUE    FALSE                 2
#>       TRUE    FALSE                 2
#> Every column, unrounded: as.data.frame(x).
#> No candidate is selected automatically. Compare one convention consistently.
# By name:
get_candidate(grid, "profiles2_VVI")
#> Latent profile analysis: 2 profiles
#> 272 observations; varying diagonal residual covariance (VVI)
#> Log likelihood: -1147.806353 | AIC: 2313.613 | BIC: 2346.065
#> Converged: TRUE | iterations: 6 | best log likelihood reached by 2 of 2 starts (2 converged)
#> Smallest standardized Gaussian eigenvalue: 0.05419 (profiles, profile_2); Gaussian mass: 96.97; weight effective n: 97.02.
#> 
#> Class sizes   profile_1 64.3% (175)   profile_2 35.7% (97)
#> 
#> Means
#>  indicator profile_1 profile_2
#>  eruptions     4.291     2.038
#>    waiting    79.986    54.493
#> 
#> Variances and standard errors: get_results(x, "profiles").
#> Every other table: get_results(x, what = ), or get_results(x, "all").
# By number of profiles and structure:
get_candidate(grid, n_profiles = 2, model = "VVI")
#> Latent profile analysis: 2 profiles
#> 272 observations; varying diagonal residual covariance (VVI)
#> Log likelihood: -1147.806353 | AIC: 2313.613 | BIC: 2346.065
#> Converged: TRUE | iterations: 6 | best log likelihood reached by 2 of 2 starts (2 converged)
#> Smallest standardized Gaussian eigenvalue: 0.05419 (profiles, profile_2); Gaussian mass: 96.97; weight effective n: 97.02.
#> 
#> Class sizes   profile_1 64.3% (175)   profile_2 35.7% (97)
#> 
#> Means
#>  indicator profile_1 profile_2
#>  eruptions     4.291     2.038
#>    waiting    79.986    54.493
#> 
#> Variances and standard errors: get_results(x, "profiles").
#> Every other table: get_results(x, what = ), or get_results(x, "all").
# From a comparison, by row label:
compared <- compare_models(grid, "profiles1_VVI", "profiles2_VVI")
get_candidate(compared, "profiles2_VVI")
#> Latent profile analysis: 2 profiles
#> 272 observations; varying diagonal residual covariance (VVI)
#> Log likelihood: -1147.806353 | AIC: 2313.613 | BIC: 2346.065
#> Converged: TRUE | iterations: 6 | best log likelihood reached by 2 of 2 starts (2 converged)
#> Smallest standardized Gaussian eigenvalue: 0.05419 (profiles, profile_2); Gaussian mass: 96.97; weight effective n: 97.02.
#> 
#> Class sizes   profile_1 64.3% (175)   profile_2 35.7% (97)
#> 
#> Means
#>  indicator profile_1 profile_2
#>  eruptions     4.291     2.038
#>    waiting    79.986    54.493
#> 
#> Variances and standard errors: get_results(x, "profiles").
#> Every other table: get_results(x, what = ), or get_results(x, "all").
```
