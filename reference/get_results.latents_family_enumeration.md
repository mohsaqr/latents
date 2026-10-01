# Tables of a group-class family enumeration

The result of `enumerate_classes(family = )` with group-class families:
one fitted candidate per family, between-variance restriction and number
of group classes.

## Usage

``` r
# S3 method for class 'latents_family_enumeration'
get_results(
  x,
  what = c("candidates", "best"),
  criterion = c("bic", "aic", "bic_individual"),
  ...
)

# S3 method for class 'latents_family_enumeration'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'latents_family_enumeration'
print(x, ...)

# S3 method for class 'latents_family_enumeration'
summary(object, ...)
```

## Arguments

- x:

  A `latents_family_enumeration` result.

- what:

  `"candidates"` (one row per candidate) or `"best"` (the candidate with
  the lowest value of `criterion`).

- criterion:

  `"bic"` (penalized by the number of groups, the default), `"aic"` or
  `"bic_individual"`.

- ...:

  Unused.

- row.names, optional:

  Unused; part of the generic.

- object:

  A `latents_family_enumeration` result.

## Value

A base `data.frame`, one row per candidate: `model` (the code
[`candidate_fit()`](https://pak.dynasite.org/latents/reference/candidate_fit.md)
takes), `family`, `between_variance`, `n_group_classes`,
`log_likelihood`, `n_parameters`, `aic`, `bic`, `bic_individual`,
relative `entropy`, `smallest_class`, `min_effective_groups`,
`converged`, `boundary`, `n_best_replicated`, and the candidate's
`warnings` and `error`. With one group class every family is the same
model, so those rows share a likelihood.

## Examples

``` r
set.seed(1)
ratings <- data.frame(
  team = rep(seq_len(40), each = 6),
  climate = rep(rnorm(40, rep(c(-1, 1), each = 20), 0.5), each = 6) +
    rnorm(240))
grid <- enumerate_classes(ratings, "climate", "team",
                          family = c("additive", "dispersion"),
                          n_group_classes = 1:2, n_starts = 3, seed = 1)
get_results(grid)
#>            model     family between_variance n_group_classes log_likelihood
#> 1       additive   additive          varying               1      -383.9742
#> 2 additive_equal   additive            equal               1      -383.9742
#> 3     dispersion dispersion            equal               1      -383.9742
#> 4       additive   additive          varying               2      -378.4544
#> 5 additive_equal   additive            equal               2      -378.7540
#> 6     dispersion dispersion            equal               2      -383.9742
#>   n_parameters      aic      bic bic_individual   entropy smallest_class
#> 1            3 773.9485 779.0151       784.3904        NA      1.0000000
#> 2            3 773.9485 779.0151       784.3904        NA      1.0000000
#> 3            3 773.9485 779.0151       784.3904        NA      1.0000000
#> 4            6 768.9088 779.0421       789.7926 0.8381067      0.4155764
#> 5            5 767.5080 775.9524       784.9112 0.8429424      0.4997789
#> 6            5 777.9485 786.3929       795.3517 0.1877010      0.2506457
#>   min_effective_groups converged boundary n_best_replicated warnings error
#> 1                   NA      TRUE    FALSE                 3           <NA>
#> 2                   NA      TRUE    FALSE                 3           <NA>
#> 3                   NA      TRUE    FALSE                 3           <NA>
#> 4             5.197275      TRUE    FALSE                 1           <NA>
#> 5            15.574944      TRUE    FALSE                 3           <NA>
#> 6                   NA      TRUE    FALSE                 3           <NA>
get_results(grid, "best")
#>            model   family between_variance n_group_classes log_likelihood
#> 5 additive_equal additive            equal               2       -378.754
#>   n_parameters     aic      bic bic_individual   entropy smallest_class
#> 5            5 767.508 775.9524       784.9112 0.8429424      0.4997789
#>   min_effective_groups converged boundary n_best_replicated warnings error
#> 5             15.57494      TRUE    FALSE                 3           <NA>
```
