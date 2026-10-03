# Predict from a mixture-of-regressions fit

Predict from a mixture-of-regressions fit

## Usage

``` r
# S3 method for class 'latents_mixture_regression'
predict(
  object,
  newdata = NULL,
  type = c("response", "class_response", "posterior", "probabilities"),
  ...
)
```

## Arguments

- object:

  A fit from
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md).

- newdata:

  `NULL` for the fitted rows, or a data frame with the predictors (and
  the membership covariates and `id` the model uses).

- type:

  `"response"`: the mean outcome averaged over the classes with their
  prior probabilities (the covariates, but not the outcome, inform
  them). `"class_response"`: the mean outcome under every class, one row
  per row and class. `"posterior"`: the posterior class probabilities
  given the outcome, which `newdata` must then contain.
  `"probabilities"` (ordinal family only): the probability of every
  category under every class. For the ordinal family the "mean outcome"
  is the expected category score, the categories scored 1 to C in order.

- ...:

  Unused.

## Value

A base `data.frame`. For `"response"`: `row` and `fitted`. For
`"class_response"`: `row`, `class`, `prior` and `fitted`, one row per
data row and class. For `"posterior"`: the same columns as
`get_results(fit, "assignments")`. For `"probabilities"`: `row`,
`class`, `category`, `prior` and `probability`, one row per data row,
class and category. Raises `latents_bad_argument` for `"probabilities"`
with another family.

## Examples

``` r
fit <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2,
                          n_starts = 3, seed = 1)
new_students <- data.frame(hours = c(2, 10))
predict(fit, new_students, type = "class_response")
#>   row   class     prior   fitted
#> 1   1 class_1 0.5638921 43.48626
#> 2   2 class_1 0.5638921 80.01323
#> 3   1 class_2 0.4361079 56.51189
#> 4   2 class_2 0.4361079 62.97886
```
