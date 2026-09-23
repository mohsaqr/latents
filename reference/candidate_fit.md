# Take one fitted model out of an enumeration grid

An enumeration keeps every candidate it fitted. This returns one of
them, so that a caller who wants to plot, summarise or test a particular
candidate names it by its class counts instead of indexing the grid by
position.

## Usage

``` r
candidate_fit(x, n_profiles, n_group_classes = 1L, model = NULL)
```

## Arguments

- x:

  An `multilpa_enumeration` result from
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md).

- n_profiles:

  Number of individual profiles identifying the candidate.

- n_group_classes:

  Number of group classes identifying the candidate.

- model:

  The covariance model identifying the candidate, needed when
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  was given several `model`s and the class counts alone name several
  candidates. Naming counts that match more than one candidate raises
  `latents_unknown_candidate` listing the models it could have meant.

## Value

The fitted model for that cell of the grid: an object of class
`multilpa`, exactly as
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
returned it, with every verb of this package available on it.

## Conditions

`latents_unknown_candidate` when the requested class counts are not in
the grid, and `latents_failed_candidate` when they are in the grid but
that fit raised an error, so no model exists to return. Neither
situation returns `NULL`, because a `NULL` would flow silently into
whatever the caller did next.

## Examples

``` r
set.seed(1)
d <- data.frame(g = rep(1:10, each = 10), y = rnorm(100))
candidates <- enumerate_classes(d, "y", "g", n_profiles = 1:2,
                              n_group_classes = 1, n_starts = 2, seed = 1)
candidate_fit(candidates, n_profiles = 2, n_group_classes = 1)
#> Two-level latent profile analysis: 2 profiles, 1 group class
#> 100 individuals in 10 groups; varying diagonal residual covariance (VVI)
#> Log likelihood: -130.592099 | AIC: 271.184 | BIC (groups): 272.697
#> Converged: TRUE | iterations: 268 | best start: 1/2
#> 
#>  profile          y    count proportion
#>        1  0.3584555 53.72896  0.5372896
#>        2 -0.1808800 46.27104  0.4627104
#> 
#> Variances and standard errors: get_results(x, "profiles"). 
#> Every other table: get_results(x, what = ), or get_results(x, "all").
```
