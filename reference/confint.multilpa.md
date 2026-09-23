# Wald confidence intervals for multilevel LPA coefficients

Wald confidence intervals for multilevel LPA coefficients

## Usage

``` r
# S3 method for class 'multilpa'
confint(object, parm, level = 0.95, data = NULL, ...)
```

## Arguments

- object:

  A fitted `multilpa` model.

- parm:

  Optional coefficient names or indices; defaults to every coefficient
  the fit estimated. Naming a coefficient that `fixed` held raises
  `latents_held_parameter`, because a held value has no interval.

- level:

  Confidence level strictly between zero and one.

- data:

  Optional, exactly as for
  [`vcov()`](https://rdrr.io/r/stats/vcov.html).

- ...:

  Additional arguments passed to
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md),
  including `method = "bootstrap"` for a structure the Wald path cannot
  chart.

## Value

A two-column matrix on the natural scale, one row per requested
coefficient and named as [`coef()`](https://rdrr.io/r/stats/coef.html)
names them. Wald intervals by default; with `method = "bootstrap"` the
percentile interval the replicates give, which is the same interval
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
reports rather than a normal approximation rebuilt from the bootstrap
standard error. Wald bounds are not clipped to the probability or
variance parameter space. For a fit made with `fixed`, the default rows
are the estimated coefficients only and the intervals are conditional on
the held measurement solution.

## Examples

``` r
set.seed(3)
example_data <- data.frame(
  school = rep(seq_len(10), each = 8),
  score_a = stats::rnorm(80), score_b = stats::rnorm(80)
)
fit <- multilpa(example_data, c("score_a", "score_b"), "school",
                n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
confint(fit, data = example_data)
#>                                                      2.5%       97.5%
#> measurement.mean.profile_1.score_a          -0.5083842790 -0.02206750
#> measurement.mean.profile_1.score_b          -0.1391045565  0.38197993
#> measurement.mean.profile_2.score_a           0.8825060016  1.15167041
#> measurement.mean.profile_2.score_b          -1.0313876060  0.00900769
#> measurement.variance.profile_1.score_a       0.3548112819  0.84763384
#> measurement.variance.profile_1.score_b       0.6676711274  1.40800776
#> measurement.variance.profile_2.score_a       0.0005662474  0.08469865
#> measurement.variance.profile_2.score_b       0.1139680930  1.20162877
#> profile.probability.profile_1.group_class_1  0.6693038676  0.93195638
#> profile.probability.profile_2.group_class_1  0.0680436177  0.33069613
#> group.probability.group_class_1              1.0000000000  1.00000000
```
