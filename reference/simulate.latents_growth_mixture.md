# Simulate outcomes from a growth mixture model

Draws new outcomes on the fitted design: each person's class from their
membership probabilities, their random effects from that class's
covariance, and their outcomes around their own trajectory with the
class's residual variance.

## Usage

``` r
# S3 method for class 'latents_growth_mixture'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- object:

  A `latents_growth_mixture` fit.

- nsim:

  Number of simulated outcome vectors.

- seed:

  `NULL` or an integer seed; the caller's random state is restored
  afterwards.

- ...:

  Unused.

## Value

A base `data.frame` with one row per fitted row and one column per
simulation, `sim_1`, `sim_2`, ...

## Examples

``` r
# \donttest{
fit <- mixture_regression(score ~ wave, growth_scores, n_classes = 3,
                          id = "student", class_level = "group",
                          random = "wave", random_covariance = "equal",
                          n_starts = 2, seed = 1)
head(simulate(fit, nsim = 2, seed = 1))
#>      sim_1    sim_2
#> 1 46.37865 49.53938
#> 2 54.89669 56.64919
#> 3 57.20172 56.30151
#> 4 65.90467 58.25365
#> 5 63.93110 61.80666
#> 6 65.30118 69.37129
# }
```
