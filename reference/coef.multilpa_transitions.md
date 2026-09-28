# Fitted parameters of a latent transition model

Fitted parameters of a latent transition model

## Usage

``` r
# S3 method for class 'multilpa_transitions'
coef(object, ...)
```

## Arguments

- object:

  A fitted `multilpa_transitions` model.

- ...:

  Reserved for compatibility with
  [`coef()`](https://rdrr.io/r/stats/coef.html).

## Value

A named numeric vector of every parameter on its natural scale, in the
order and with the names
[`vcov.multilpa_transitions()`](https://pak.dynasite.org/latents/reference/parameter_inference.multilpa_transitions.md)
and
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
use: profile means, then variances (one shared set under
`variance_model = "equal"`, the covariance matrices under
`covariance_model = "full"`), categorical response probabilities where
present, the group-class probabilities, each group class's initial
profile probabilities, and its transition probabilities row by row.

Names follow the package-wide `level.parameter.outcome.term` grammar
shared with
[`coef.multilpa()`](https://pak.dynasite.org/latents/reference/coef.multilpa.md),
so one pattern matches across fit classes: for example
`measurement.mean.profile_1.reading` and
`profile.transition_probability.profile_2.group_class_1:profile_1`,
whose term names the origin the move is from. `level`, `parameter` and
`outcome` never contain a dot, so everything after the third dot is the
term and the name parses back into the columns
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
reports even when an indicator name itself contains a dot.

Their standard errors come from
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
or
[`vcov.multilpa_transitions()`](https://pak.dynasite.org/latents/reference/parameter_inference.multilpa_transitions.md).
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
gives the same quantities as tidy tables, which is the form to prefer.

## Examples

``` r
fit <- lta(subset(course_engagement, student <= 40),
           c("browse", "lectures", "forum_read"), "student",
           n_profiles = 2, time = "sequence", n_starts = 2, seed = 1)
coef(fit)
#>                                measurement.mean.profile_1.browse 
#>                                                        0.5552352 
#>                              measurement.mean.profile_1.lectures 
#>                                                        0.4372193 
#>                            measurement.mean.profile_1.forum_read 
#>                                                        0.6089409 
#>                                measurement.mean.profile_2.browse 
#>                                                       -0.8706168 
#>                              measurement.mean.profile_2.lectures 
#>                                                       -0.6264837 
#>                            measurement.mean.profile_2.forum_read 
#>                                                       -0.9319730 
#>                            measurement.variance.profile_1.browse 
#>                                                        0.5490615 
#>                          measurement.variance.profile_1.lectures 
#>                                                        0.8253164 
#>                        measurement.variance.profile_1.forum_read 
#>                                                        0.4991390 
#>                            measurement.variance.profile_2.browse 
#>                                                        0.4601966 
#>                          measurement.variance.profile_2.lectures 
#>                                                        0.5736107 
#>                        measurement.variance.profile_2.forum_read 
#>                                                        0.3259582 
#>                                  group.probability.group_class_1 
#>                                                        1.0000000 
#>              profile.initial_probability.profile_1.group_class_1 
#>                                                        0.7421284 
#>              profile.initial_probability.profile_2.group_class_1 
#>                                                        0.2578716 
#> profile.transition_probability.profile_1.group_class_1:profile_1 
#>                                                        0.8906008 
#> profile.transition_probability.profile_2.group_class_1:profile_1 
#>                                                        0.1093992 
#> profile.transition_probability.profile_1.group_class_1:profile_2 
#>                                                        0.1896193 
#> profile.transition_probability.profile_2.group_class_1:profile_2 
#>                                                        0.8103807 
```
