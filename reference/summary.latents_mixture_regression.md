# Summarize a mixture-of-regressions fit

Summarize a mixture-of-regressions fit

## Usage

``` r
# S3 method for class 'latents_mixture_regression'
summary(object, level = 0.95, vcov_type = NULL, ...)

# S3 method for class 'summary_latents_mixture_regression'
print(x, digits = 4L, ...)
```

## Arguments

- object:

  A fit from
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md).

- level:

  Confidence level.

- vcov_type:

  As in
  [`get_results.latents_mixture_regression()`](https://pak.dynasite.org/latents/reference/get_results.latents_mixture_regression.md).

- ...:

  Unused.

- x:

  A `summary_latents_mixture_regression` object.

- digits:

  Significant digits.

## Value

An object of class `summary_latents_mixture_regression`: a named list of
the `fit`, `coefficients`, `classes`, `membership` and `group_classes`
tables, printed by its print method.
