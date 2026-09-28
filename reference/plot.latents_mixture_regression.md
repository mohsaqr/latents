# Plot a mixture-of-regressions fit

Plot a mixture-of-regressions fit

## Usage

``` r
# S3 method for class 'latents_mixture_regression'
plot(
  x,
  what = c("fitted", "coefficients", "posteriors"),
  predictor = NULL,
  level = 0.95,
  main = NULL,
  ...
)
```

## Arguments

- x:

  A fit from
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md).

- what:

  `"fitted"`: the outcome against one numeric predictor, rows marked by
  modal class (colour and symbol), each class's regression line drawn
  with the other predictors at their mean (or reference level).
  `"coefficients"`: every class's estimates with confidence intervals.
  `"posteriors"`: the distribution of each unit's largest posterior
  probability, by class.

- predictor:

  For `"fitted"`, the numeric predictor for the horizontal axis;
  defaults to the first one in the formula.

- level:

  Confidence level for `"coefficients"`.

- main:

  Plot title; `NULL` for the default.

- ...:

  Unused.

## Value

`x`, invisibly. Called for its plot.

## Examples

``` r
fit <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2,
                          n_starts = 3, seed = 1)
plot(fit)

plot(fit, what = "coefficients")
```
