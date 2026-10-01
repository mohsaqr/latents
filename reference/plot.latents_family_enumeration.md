# Plot information criteria across group-class families

One line per model code, criterion against the number of group classes;
models are told apart by colour, point shape and line type.

## Usage

``` r
# S3 method for class 'latents_family_enumeration'
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

  A `latents_family_enumeration` result.

- criterion:

  `"bic"`, `"aic"` or `"bic_individual"`.

- main, subtitle:

  Optional title and subtitle.

- ...:

  Unused.

## Value

A ggplot object. Raises `latents_missing_package` without ggplot2.
