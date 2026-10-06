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

## Examples

``` r
few_students <- subset(growth_scores, student <= 60)
trajectories <- mixture_regression(score ~ wave, few_students, n_classes = 2,
                                   id = "student", class_level = "group",
                                   seed = 1)
growth <- mixture_regression(score ~ wave, few_students, n_classes = 2,
                             id = "student", class_level = "group",
                             random = "intercept", seed = 1)
comparison <- compare_models(trajectories = trajectories, growth = growth)
if (requireNamespace("ggplot2", quietly = TRUE)) plot(comparison)
```
