# Extract multilevel LPA coefficients

Extract multilevel LPA coefficients

## Usage

``` r
# S3 method for class 'multilpa'
coef(object, scale = c("natural", "unconstrained"), ...)
```

## Arguments

- object:

  A fitted `multilpa` model.

- scale:

  Natural coefficients or unconstrained log variances (diagonal),
  log-Cholesky coordinates (full covariance), and baseline-category
  logits.

- ...:

  Reserved for generic compatibility.

## Value

A named numeric vector of every coefficient the model has, held ones
included: a block `fixed` held is part of the model and is reported
here, even though it has no standard error and no interval. Natural
coefficients include all mixing probabilities; unconstrained
coefficients exclude their reference categories. Names are
`level.parameter.outcome.term`, the same four-part decomposition
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
reports as columns and the same spelling every fitted class in this
package uses, so a name written for one fit means the same thing for
another. A parameter with no term – a group-class probability – carries
the first three parts only.
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
is the tidy form and the one to prefer.

## Examples

``` r
set.seed(3)
example_data <- data.frame(
  school = rep(seq_len(10), each = 8),
  score_a = stats::rnorm(80), score_b = stats::rnorm(80)
)
fit <- multilpa(example_data, c("score_a", "score_b"), "school",
                n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
coef(fit)
#>          measurement.mean.profile_1.score_a 
#>                                 -0.26522589 
#>          measurement.mean.profile_1.score_b 
#>                                  0.12143769 
#>          measurement.mean.profile_2.score_a 
#>                                  1.01708820 
#>          measurement.mean.profile_2.score_b 
#>                                 -0.51118996 
#>      measurement.variance.profile_1.score_a 
#>                                  0.60122256 
#>      measurement.variance.profile_1.score_b 
#>                                  1.03783944 
#>      measurement.variance.profile_2.score_a 
#>                                  0.04263245 
#>      measurement.variance.profile_2.score_b 
#>                                  0.65779843 
#> profile.probability.profile_1.group_class_1 
#>                                  0.80063012 
#> profile.probability.profile_2.group_class_1 
#>                                  0.19936988 
#>             group.probability.group_class_1 
#>                                  1.00000000 
```
