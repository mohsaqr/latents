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
if (requireNamespace("ggplot2", quietly = TRUE)) plot(grid)
```
