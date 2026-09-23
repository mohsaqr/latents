# Estimated parameters of a covariate fit

Estimated parameters of a covariate fit

## Usage

``` r
# S3 method for class 'multilpa_covariates'
coef(object, scale = c("natural", "unconstrained"), ...)
```

## Arguments

- object:

  A fitted `multilpa_covariates` model.

- scale:

  Natural coefficients, or the unconstrained coordinates the model is
  estimated on: log variances for a diagonal fit and log-Cholesky
  coordinates for a full-covariance fit. Matches
  [`coef.multilpa()`](https://pak.dynasite.org/latents/reference/coef.multilpa.md).

- ...:

  Ignored, present for generic compatibility.

## Value

A named numeric vector of every free parameter, in the order
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
reports them: measurement means, measurement variances or residual
covariances, profile logits, then group-class logits. Names are
`level.parameter.outcome.term`, the same four-part decomposition
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
reports as columns and the same spelling every fitted class in this
package uses; it is what keeps a mean and a variance on the same
indicator distinguishable, and the `parameter` part names the scale, so
a log variance is never served under a name that says `variance`. On the
natural scale variances and residual covariances are in their own units
and agree with
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)'s
`estimate` column exactly.

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
coef(fit)
#>               measurement.mean.profile_1.y1 
#>                                   1.9190929 
#>               measurement.mean.profile_1.y2 
#>                                   1.4517555 
#>               measurement.mean.profile_2.y1 
#>                                  -1.8353900 
#>               measurement.mean.profile_2.y2 
#>                                  -1.4838862 
#>           measurement.variance.profile_1.y1 
#>                                   0.3992639 
#>           measurement.variance.profile_1.y2 
#>                                   0.5601752 
#>           measurement.variance.profile_2.y1 
#>                                   0.4605175 
#>           measurement.variance.profile_2.y2 
#>                                   0.5202968 
#> profile.coefficient.profile_1.group_class_1 
#>                                  -0.7082134 
#> profile.coefficient.profile_1.group_class_2 
#>                                   1.2934651 
#>             profile.coefficient.profile_1.x 
#>                                   0.7401938 
#> group.coefficient.group_class_1.(Intercept) 
#>                                   0.2318505 
```
