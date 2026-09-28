# Covariance matrix of a covariate fit

Covariance matrix of a covariate fit

## Usage

``` r
# S3 method for class 'multilpa_covariates'
vcov(
  object,
  data = NULL,
  step = 1e-04,
  vcov_type = c("observed", "robust", "opg"),
  scale = c("natural", "unconstrained"),
  boundary = c("error", "fix"),
  ...
)
```

## Arguments

- object:

  A fitted `multilpa_covariates` model.

- data:

  Optional. The data frame the model was fitted to; when omitted it is
  rebuilt from the indicators, group index and designs the fit stores.
  Supplying it is the stronger check that the caller still holds that
  frame.

- step:

  Finite-difference step for the observed information.

- vcov_type:

  `"observed"`, `"robust"` or `"opg"`, as for
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md).

- scale:

  Which parameter scale the covariance is on, matching
  [`vcov.multilpa()`](https://pak.dynasite.org/latents/reference/vcov.multilpa.md).
  `"natural"`, the default, is the covariance of the estimates
  [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  report: variances and residual covariances in their own units, carried
  from the estimation scale by the delta method. `"unconstrained"` is
  the covariance on the scale the model is estimated on, with log
  variances and log-Cholesky coordinates. Means and membership
  coefficients are the same on both scales.

- boundary:

  `"error"` or `"fix"`, as for
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md):
  `"fix"` holds response probabilities on their bound, whose rows and
  columns are then zero.

- ...:

  Ignored, present for generic compatibility.

## Value

A square numeric matrix with one row and column per free parameter,
named as [`coef()`](https://rdrr.io/r/stats/coef.html) names them and
ordered measurement means, measurement variances or covariances, profile
logits, group logits.

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
vcov(fit, example_data, scale = "unconstrained")
#>                                             measurement.mean.profile_1.y1
#> measurement.mean.profile_1.y1                                5.786640e-03
#> measurement.mean.profile_1.y2                                2.179229e-07
#> measurement.mean.profile_2.y1                                2.254528e-07
#> measurement.mean.profile_2.y2                                1.077812e-07
#> measurement.log_variance.profile_1.y1                       -1.002598e-06
#> measurement.log_variance.profile_1.y2                       -7.547380e-07
#> measurement.log_variance.profile_2.y1                        8.950452e-07
#> measurement.log_variance.profile_2.y2                        1.299308e-07
#> profile.coefficient.profile_1.group_class_1                 -1.738891e-07
#> profile.coefficient.profile_1.group_class_2                 -8.213189e-07
#> profile.coefficient.profile_1.x                             -1.227484e-07
#> group.coefficient.group_class_1.(Intercept)                 -2.176109e-07
#>                                             measurement.mean.profile_1.y2
#> measurement.mean.profile_1.y1                                2.179229e-07
#> measurement.mean.profile_1.y2                                8.118727e-03
#> measurement.mean.profile_2.y1                                2.772009e-07
#> measurement.mean.profile_2.y2                                1.335890e-07
#> measurement.log_variance.profile_1.y1                       -1.028549e-06
#> measurement.log_variance.profile_1.y2                       -8.065083e-07
#> measurement.log_variance.profile_2.y1                        1.206184e-06
#> measurement.log_variance.profile_2.y2                        1.865066e-07
#> profile.coefficient.profile_1.group_class_1                 -2.310731e-07
#> profile.coefficient.profile_1.group_class_2                 -8.668341e-07
#> profile.coefficient.profile_1.x                             -1.542014e-07
#> group.coefficient.group_class_1.(Intercept)                 -2.331831e-07
#>                                             measurement.mean.profile_2.y1
#> measurement.mean.profile_1.y1                                2.254528e-07
#> measurement.mean.profile_1.y2                                2.772009e-07
#> measurement.mean.profile_2.y1                                7.805867e-03
#> measurement.mean.profile_2.y2                                2.650758e-07
#> measurement.log_variance.profile_1.y1                       -7.737157e-07
#> measurement.log_variance.profile_1.y2                       -7.151782e-07
#> measurement.log_variance.profile_2.y1                        2.602638e-06
#> measurement.log_variance.profile_2.y2                        5.823831e-07
#> profile.coefficient.profile_1.group_class_1                 -6.878763e-07
#> profile.coefficient.profile_1.group_class_2                 -7.373604e-07
#> profile.coefficient.profile_1.x                             -3.377541e-07
#> group.coefficient.group_class_1.(Intercept)                 -3.886793e-07
#>                                             measurement.mean.profile_2.y2
#> measurement.mean.profile_1.y1                                1.077812e-07
#> measurement.mean.profile_1.y2                                1.335890e-07
#> measurement.mean.profile_2.y1                                2.650758e-07
#> measurement.mean.profile_2.y2                                8.818734e-03
#> measurement.log_variance.profile_1.y1                       -3.188609e-07
#> measurement.log_variance.profile_1.y2                       -3.061991e-07
#> measurement.log_variance.profile_2.y1                        1.442898e-06
#> measurement.log_variance.profile_2.y2                        3.560253e-07
#> profile.coefficient.profile_1.group_class_1                 -4.266238e-07
#> profile.coefficient.profile_1.group_class_2                 -3.138273e-07
#> profile.coefficient.profile_1.x                             -1.959986e-07
#> group.coefficient.group_class_1.(Intercept)                 -2.301935e-07
#>                                             measurement.log_variance.profile_1.y1
#> measurement.mean.profile_1.y1                                       -1.002598e-06
#> measurement.mean.profile_1.y2                                       -1.028549e-06
#> measurement.mean.profile_2.y1                                       -7.737157e-07
#> measurement.mean.profile_2.y2                                       -3.188609e-07
#> measurement.log_variance.profile_1.y1                                2.899095e-02
#> measurement.log_variance.profile_1.y2                                3.920012e-06
#> measurement.log_variance.profile_2.y1                               -2.186616e-06
#> measurement.log_variance.profile_2.y2                                1.232158e-08
#> profile.coefficient.profile_1.group_class_1                          8.528079e-08
#> profile.coefficient.profile_1.group_class_2                          4.308017e-06
#> profile.coefficient.profile_1.x                                      2.874847e-07
#> group.coefficient.group_class_1.(Intercept)                          7.970725e-07
#>                                             measurement.log_variance.profile_1.y2
#> measurement.mean.profile_1.y1                                       -7.547380e-07
#> measurement.mean.profile_1.y2                                       -8.065083e-07
#> measurement.mean.profile_2.y1                                       -7.151782e-07
#> measurement.mean.profile_2.y2                                       -3.061991e-07
#> measurement.log_variance.profile_1.y1                                3.920012e-06
#> measurement.log_variance.profile_1.y2                                2.898848e-02
#> measurement.log_variance.profile_2.y1                               -2.552608e-06
#> measurement.log_variance.profile_2.y2                               -1.642739e-07
#> profile.coefficient.profile_1.group_class_1                          2.127237e-07
#> profile.coefficient.profile_1.group_class_2                          3.216488e-06
#> profile.coefficient.profile_1.x                                      3.121792e-07
#> group.coefficient.group_class_1.(Intercept)                          5.953360e-07
#>                                             measurement.log_variance.profile_2.y1
#> measurement.mean.profile_1.y1                                        8.950452e-07
#> measurement.mean.profile_1.y2                                        1.206184e-06
#> measurement.mean.profile_2.y1                                        2.602638e-06
#> measurement.mean.profile_2.y2                                        1.442898e-06
#> measurement.log_variance.profile_1.y1                               -2.186616e-06
#> measurement.log_variance.profile_1.y2                               -2.552608e-06
#> measurement.log_variance.profile_2.y1                                3.391306e-02
#> measurement.log_variance.profile_2.y2                                3.470680e-06
#> profile.coefficient.profile_1.group_class_1                         -4.005021e-06
#> profile.coefficient.profile_1.group_class_2                         -2.488502e-06
#> profile.coefficient.profile_1.x                                     -1.857712e-06
#> group.coefficient.group_class_1.(Intercept)                         -1.879375e-06
#>                                             measurement.log_variance.profile_2.y2
#> measurement.mean.profile_1.y1                                        1.299308e-07
#> measurement.mean.profile_1.y2                                        1.865066e-07
#> measurement.mean.profile_2.y1                                        5.823831e-07
#> measurement.mean.profile_2.y2                                        3.560253e-07
#> measurement.log_variance.profile_1.y1                                1.232158e-08
#> measurement.log_variance.profile_1.y2                               -1.642739e-07
#> measurement.log_variance.profile_2.y1                                3.470680e-06
#> measurement.log_variance.profile_2.y2                                3.389929e-02
#> profile.coefficient.profile_1.group_class_1                         -1.168408e-06
#> profile.coefficient.profile_1.group_class_2                         -1.354686e-07
#> profile.coefficient.profile_1.x                                     -4.912652e-07
#> group.coefficient.group_class_1.(Intercept)                         -5.228709e-07
#>                                             profile.coefficient.profile_1.group_class_1
#> measurement.mean.profile_1.y1                                             -1.738891e-07
#> measurement.mean.profile_1.y2                                             -2.310731e-07
#> measurement.mean.profile_2.y1                                             -6.878763e-07
#> measurement.mean.profile_2.y2                                             -4.266238e-07
#> measurement.log_variance.profile_1.y1                                      8.528079e-08
#> measurement.log_variance.profile_1.y2                                      2.127237e-07
#> measurement.log_variance.profile_2.y1                                     -4.005021e-06
#> measurement.log_variance.profile_2.y2                                     -1.168408e-06
#> profile.coefficient.profile_1.group_class_1                                1.091101e-01
#> profile.coefficient.profile_1.group_class_2                                2.142382e-02
#> profile.coefficient.profile_1.x                                           -1.246622e-02
#> group.coefficient.group_class_1.(Intercept)                                5.455471e-02
#>                                             profile.coefficient.profile_1.group_class_2
#> measurement.mean.profile_1.y1                                             -8.213189e-07
#> measurement.mean.profile_1.y2                                             -8.668341e-07
#> measurement.mean.profile_2.y1                                             -7.373604e-07
#> measurement.mean.profile_2.y2                                             -3.138273e-07
#> measurement.log_variance.profile_1.y1                                      4.308017e-06
#> measurement.log_variance.profile_1.y2                                      3.216488e-06
#> measurement.log_variance.profile_2.y1                                     -2.488502e-06
#> measurement.log_variance.profile_2.y2                                     -1.354686e-07
#> profile.coefficient.profile_1.group_class_1                                2.142382e-02
#> profile.coefficient.profile_1.group_class_2                                1.733277e-01
#> profile.coefficient.profile_1.x                                            5.582183e-03
#> group.coefficient.group_class_1.(Intercept)                                6.706351e-02
#>                                             profile.coefficient.profile_1.x
#> measurement.mean.profile_1.y1                                 -1.227484e-07
#> measurement.mean.profile_1.y2                                 -1.542014e-07
#> measurement.mean.profile_2.y1                                 -3.377541e-07
#> measurement.mean.profile_2.y2                                 -1.959986e-07
#> measurement.log_variance.profile_1.y1                          2.874847e-07
#> measurement.log_variance.profile_1.y2                          3.121792e-07
#> measurement.log_variance.profile_2.y1                         -1.857712e-06
#> measurement.log_variance.profile_2.y2                         -4.912652e-07
#> profile.coefficient.profile_1.group_class_1                   -1.246622e-02
#> profile.coefficient.profile_1.group_class_2                    5.582183e-03
#> profile.coefficient.profile_1.x                                5.157617e-02
#> group.coefficient.group_class_1.(Intercept)                   -4.059663e-03
#>                                             group.coefficient.group_class_1.(Intercept)
#> measurement.mean.profile_1.y1                                             -2.176109e-07
#> measurement.mean.profile_1.y2                                             -2.331831e-07
#> measurement.mean.profile_2.y1                                             -3.886793e-07
#> measurement.mean.profile_2.y2                                             -2.301935e-07
#> measurement.log_variance.profile_1.y1                                      7.970725e-07
#> measurement.log_variance.profile_1.y2                                      5.953360e-07
#> measurement.log_variance.profile_2.y1                                     -1.879375e-06
#> measurement.log_variance.profile_2.y2                                     -5.228709e-07
#> profile.coefficient.profile_1.group_class_1                                5.455471e-02
#> profile.coefficient.profile_1.group_class_2                                6.706351e-02
#> profile.coefficient.profile_1.x                                           -4.059663e-03
#> group.coefficient.group_class_1.(Intercept)                                3.575070e-01
```
