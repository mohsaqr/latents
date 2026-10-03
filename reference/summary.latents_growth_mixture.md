# Summarize a growth mixture model

Summarize a growth mixture model

## Usage

``` r
# S3 method for class 'latents_growth_mixture'
summary(object, level = 0.95, vcov_type = NULL, ...)

# S3 method for class 'summary_latents_growth_mixture'
print(x, digits = 4L, ...)
```

## Arguments

- object:

  A `latents_growth_mixture` fit.

- level:

  Confidence level.

- vcov_type:

  As in
  [`get_results.latents_growth_mixture()`](https://pak.dynasite.org/latents/reference/get_results.latents_growth_mixture.md).

- ...:

  Unused.

- x:

  A `summary_latents_growth_mixture` object.

- digits:

  Significant digits.

## Value

An object of class `summary_latents_growth_mixture`: a list of the
`fit`, `classes`, `coefficients`, `random` and `membership` tables.
