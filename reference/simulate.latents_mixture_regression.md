# Simulate outcomes from a mixture-of-regressions fit

Draws class memberships from the fitted mixing model (group classes
first in the two-level model, one class per group under
`class_level = "group"`), then outcomes from each row's class
regression, at the fitted rows' predictors.

## Usage

``` r
# S3 method for class 'latents_mixture_regression'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- object:

  A fit from
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md).

- nsim:

  Number of simulated outcome vectors.

- seed:

  `NULL` or an integer; the caller's random state is restored.

- ...:

  Unused.

## Value

A base `data.frame` with one row per fitted row and one column per
simulation, `sim_1`, `sim_2`, ...; for a binomial fit with more than one
trial the columns hold success counts.

## Examples

``` r
fit <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2,
                          n_starts = 3, seed = 1)
head(simulate(fit, nsim = 2, seed = 1))
#>      sim_1    sim_2
#> 1 85.96382 84.66682
#> 2 43.05977 61.00305
#> 3 57.87184 62.04866
#> 4 51.20873 72.31212
#> 5 78.33330 84.48466
#> 6 59.89102 71.00138
```
