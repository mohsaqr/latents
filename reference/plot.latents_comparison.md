# Plot a model comparison

Each model's AIC, BIC and sample-size adjusted BIC, lowest best; the
model with the smallest BIC is marked.

## Usage

``` r
# S3 method for class 'latents_comparison'
plot(x, main = NULL, subtitle = NULL, ...)
```

## Arguments

- x:

  A `latents_comparison` from
  [`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md).

- main, subtitle:

  Title and subtitle.

- ...:

  Unused.

## Value

A ggplot object. Raises `latents_missing_package` without ggplot2.
