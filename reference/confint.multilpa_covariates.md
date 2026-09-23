# Wald confidence intervals for a covariate fit

Wald confidence intervals for a covariate fit

## Usage

``` r
# S3 method for class 'multilpa_covariates'
confint(object, parm, level = 0.95, data = NULL, ...)
```

## Arguments

- object:

  A fitted `multilpa_covariates` model.

- parm:

  Optional parameter names or indices; defaults to all of them.

- level:

  Confidence level strictly between zero and one.

- data:

  Optional. The data frame the model was fitted to; when omitted it is
  rebuilt from the indicators, group index and designs the fit stores.
  Supplying it is the stronger check that the caller still holds that
  frame.

- ...:

  Passed to
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md),
  so `vcov_type = "robust"` and `step` reach it.

## Value

A two-column matrix of Wald intervals, one row per requested parameter,
named as [`coef()`](https://rdrr.io/r/stats/coef.html) names them.
Bounds are on the natural scale and are not constrained to respect a
variance's positivity or a probability's range.

## Examples

``` r
set.seed(5)
school <- rep(seq_len(16), each = 8)
high_class <- rep(rep(c(FALSE, TRUE), length.out = 16), each = 8)
x <- rnorm(128)
profile <- ifelse(
  runif(128) < plogis(-1 + 2 * high_class + 0.8 * x), 2L, 1L
)
example_data <- data.frame(
  school = school, x = x,
  y1 = rnorm(128, ifelse(profile == 2L, 2, -2), 0.7),
  y2 = rnorm(128, ifelse(profile == 2L, 1.5, -1.5), 0.7)
)
fit <- multilpa(example_data, c("y1", "y2"), "school", n_profiles = 2,
                n_group_classes = 2, profile_covariates = "x",
                n_starts = 2, seed = 1)
confint(fit, data = example_data)
#>                                                   2.5%       97.5%
#> measurement.mean.profile_1.y1                1.7699986  2.06818734
#> measurement.mean.profile_1.y2                1.2751549  1.62835603
#> measurement.mean.profile_2.y1               -2.0085544 -1.66222558
#> measurement.mean.profile_2.y2               -1.6679427 -1.29982964
#> measurement.variance.profile_1.y1            0.2660224  0.53250532
#> measurement.variance.profile_1.y2            0.3732428  0.74710769
#> measurement.variance.profile_2.y1            0.2942996  0.62673536
#> measurement.variance.profile_2.y2            0.3325405  0.70805306
#> profile.coefficient.profile_1.group_class_1 -1.3556252 -0.06080148
#> profile.coefficient.profile_1.group_class_2  0.4774802  2.10944992
#> profile.coefficient.profile_1.x              0.2950784  1.18530922
#> group.coefficient.group_class_1.(Intercept) -0.9400489  1.40374986
```
