# Pooled covariance of a multiply imputed fit

Pooled covariance of a multiply imputed fit

## Usage

``` r
# S3 method for class 'latents_pooled'
vcov(object, ...)
```

## Arguments

- object:

  A `latents_pooled` object.

- ...:

  Unused.

## Value

The total covariance matrix `U + (1 + 1/m) B` of the natural-scale
parameters, named as [`vcov()`](https://rdrr.io/r/stats/vcov.html) names
them for a single fit, or `NULL` when the imputations' fits carried no
covariance.
