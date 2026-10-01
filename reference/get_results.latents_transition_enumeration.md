# Tables of a transition-model enumeration

The result of `enumerate_classes(..., time = )`: one transition model
per number of profiles and group classes (and covariance structure, when
`model` names several).

## Usage

``` r
# S3 method for class 'latents_transition_enumeration'
get_results(
  x,
  what = c("candidates", "best"),
  criterion = c("bic", "aic", "bic_individual"),
  ...
)

# S3 method for class 'latents_transition_enumeration'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'latents_transition_enumeration'
print(x, ...)

# S3 method for class 'latents_transition_enumeration'
summary(object, ...)
```

## Arguments

- x:

  A `latents_transition_enumeration` result.

- what:

  `"candidates"` (one row per candidate) or `"best"` (the lowest value
  of `criterion`).

- criterion:

  `"bic"` (penalized by the number of groups), `"aic"` or
  `"bic_individual"`.

- ...:

  Unused.

- row.names, optional:

  Unused; part of the generic.

- object:

  A `latents_transition_enumeration` result.

## Value

A base `data.frame`, one row per candidate: `n_profiles`,
`n_group_classes`, `model`, `log_likelihood`, `n_parameters`, `aic`,
`bic`, `bic_individual`, relative profile `entropy`, `converged`,
`boundary`, `n_best_replicated`, `warnings` and `error`.

## Examples

``` r
grid <- enumerate_classes(course_engagement, c("browse", "lectures"),
                          "student", time = "sequence", n_profiles = 2:3,
                          n_group_classes = 1, n_starts = 2, seed = 1)
get_results(grid)
#>   n_profiles n_group_classes model log_likelihood n_parameters      aic
#> 1          2               1  <NA>      -3761.447           11 7544.894
#> 2          3               1  <NA>      -3747.521           20 7535.042
#>        bic bic_individual   entropy converged boundary n_best_replicated
#> 1 7574.192       7602.752 0.7290889      TRUE    FALSE                 2
#> 2 7588.311       7640.239 0.7627614      TRUE    FALSE                 1
#>   warnings error
#> 1           <NA>
#> 2           <NA>
get_results(grid, "best")
#>   n_profiles n_group_classes model log_likelihood n_parameters      aic
#> 1          2               1  <NA>      -3761.447           11 7544.894
#>        bic bic_individual   entropy converged boundary n_best_replicated
#> 1 7574.192       7602.752 0.7290889      TRUE    FALSE                 2
#>   warnings error
#> 1           <NA>
```
