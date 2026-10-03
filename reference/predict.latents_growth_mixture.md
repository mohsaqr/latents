# Prediction for a growth mixture model (not implemented)

Predictions for new persons are not implemented yet. The fitted class
trajectories and each fitted person's own predicted curve are tables of
[`get_results.latents_growth_mixture()`](https://pak.dynasite.org/latents/reference/get_results.latents_growth_mixture.md):
`"trajectories"` and `"individual"`.

## Usage

``` r
# S3 method for class 'latents_growth_mixture'
predict(object, ...)
```

## Arguments

- object:

  A `latents_growth_mixture` fit.

- ...:

  Unused.

## Value

Never returns; raises `latents_unsupported_prediction`.
