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

## Examples

``` r
few_students <- subset(growth_scores, student <= 60)
fit <- mixture_regression(score ~ wave, few_students, n_classes = 2,
                          id = "student", class_level = "group",
                          random = "intercept", seed = 1)
# Prediction is refused; get_results(fit, "trajectories") has the curves
try(predict(fit))
#> Error : predict() is not implemented for growth mixture models yet. The class trajectories are get_results(x, "trajectories"), and each fitted person's class trajectory and own predicted curve are get_results(x, "individual").
```
