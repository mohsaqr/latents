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
# The first 50 students keep the example quick
few_students <- subset(course_engagement, student <= 50)
grid <- enumerate_classes(few_students, c("browse", "lectures"),
                          "student", time = "sequence", n_profiles = 2:3,
                          n_group_classes = 1, n_starts = 2, seed = 1)
get_results(grid)
#>   n_profiles n_group_classes model log_likelihood n_parameters      aic
#> 1          2               1  <NA>      -1799.274           11 3620.548
#> 2          3               1  <NA>      -1790.352           20 3620.704
#>        bic bic_individual   entropy converged boundary n_best_replicated
#> 1 3641.580       3670.258 0.7752585      TRUE    FALSE                 2
#> 2 3658.944       3711.087 0.7016252      TRUE    FALSE                 2
#>   warnings error
#> 1           <NA>
#> 2           <NA>
get_results(grid, "best")
#>   n_profiles n_group_classes model log_likelihood n_parameters      aic     bic
#> 1          2               1  <NA>      -1799.274           11 3620.548 3641.58
#>   bic_individual   entropy converged boundary n_best_replicated warnings error
#> 1       3670.258 0.7752585      TRUE    FALSE                 2           <NA>
```
