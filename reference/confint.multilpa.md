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
  including `method = "bootstrap"`.

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
#>                                                    2.5%        97.5%
#> measurement.mean.profile_1.score_a          -0.50822695 -0.021900709
#> measurement.mean.profile_1.score_b          -0.13911009  0.381943368
#> measurement.mean.profile_2.score_a           0.88251814  1.151708546
#> measurement.mean.profile_2.score_b          -1.03179187  0.008913755
#> measurement.variance.profile_1.score_a       0.39917942  0.905927249
#> measurement.variance.profile_1.score_b       0.72646619  1.482528823
#> measurement.variance.profile_2.score_a       0.01589038  0.114255437
#> measurement.variance.profile_2.score_b       0.28766115  1.503720763
#> profile.probability.profile_1.group_class_1  0.63826878  0.901493808
#> profile.probability.profile_2.group_class_1  0.09850619  0.361731225
#> group.probability.group_class_1              1.00000000  1.000000000
```
