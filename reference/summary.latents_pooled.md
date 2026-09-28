# Summarize a pooled multiply imputed fit

Summarize a pooled multiply imputed fit

## Usage

``` r
# S3 method for class 'latents_pooled'
summary(object, ...)
```

## Arguments

- object:

  A `latents_pooled` object.

- ...:

  Unused.

## Value

A one-row base `data.frame`: the number of imputations, the model,
whether every imputation converged, the range of their log likelihoods,
how many needed relabelling, and the largest and median fraction of
missing information.
