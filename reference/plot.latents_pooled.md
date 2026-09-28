# Plot a pooled multiply imputed fit

One row per parameter: each imputation's aligned estimate as an open
circle and the pooled estimate as a filled circle with its interval, so
the between-imputation spread is visible next to the total uncertainty.

## Usage

``` r
# S3 method for class 'latents_pooled'
plot(x, level = c("profile", "group", "measurement"), main = NULL, ...)
```

## Arguments

- x:

  A `latents_pooled` object.

- level:

  Which parameter level to show: `"profile"` (the profile membership
  coefficients or probabilities), `"group"` or `"measurement"`.

- main:

  Plot title; `NULL` for the default.

- ...:

  Unused.

## Value

`x`, invisibly. Called for the side effect of drawing.
