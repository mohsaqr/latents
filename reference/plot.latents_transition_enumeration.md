# Plot information criteria across transition models

The criterion against the number of profiles, one line per number of
group classes (and structure), told apart by colour, shape and line
type.

## Usage

``` r
# S3 method for class 'latents_transition_enumeration'
plot(
  x,
  criterion = c("bic", "aic", "bic_individual"),
  main = NULL,
  subtitle = NULL,
  ...
)
```

## Arguments

- x:

  A `latents_transition_enumeration` result.

- criterion:

  `"bic"`, `"aic"` or `"bic_individual"`.

- main, subtitle:

  Optional title and subtitle.

- ...:

  Unused.

## Value

A ggplot object. Raises `latents_missing_package` without ggplot2.
