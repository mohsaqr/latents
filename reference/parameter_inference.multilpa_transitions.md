# Wald inference for a latent transition model

Standard errors, tests and intervals for every parameter of an
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit: the
measurement model, the group-class shares, each group class's initial
profile distribution and its transition matrix. The scores follow the
Fisher identity (each is an expected count from the forward-backward
pass minus the count the fitted probabilities imply), and the observed
information is the numerical Jacobian of that analytic score.

## Usage

``` r
# S3 method for class 'multilpa_transitions'
parameter_inference(
  x,
  data = NULL,
  level = 0.95,
  step = 1e-04,
  vcov_type = c("observed", "robust", "opg"),
  adjust = .multilpa_p_adjust_methods,
  ...
)

# S3 method for class 'multilpa_transitions'
vcov(
  object,
  data = NULL,
  step = 1e-04,
  vcov_type = c("observed", "robust", "opg"),
  scale = c("natural", "unconstrained"),
  ...
)

# S3 method for class 'multilpa_transitions'
confint(object, parm, level = 0.95, data = NULL, ...)
```

## Arguments

- x:

  A fitted `multilpa_transitions` model from
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md).

- data:

  Optional. The data frame the model was fitted to. The fit stores its
  indicators, codes, groups and occasions, so it is not needed; when
  given it must reproduce them.

- level:

  Confidence level of the intervals.

- step:

  Finite-difference step for the observed information.

- vcov_type:

  `"observed"` (observed information), `"robust"` (the sandwich
  clustered on groups, which are the sequences) or `"opg"` (the outer
  product of the group scores).

- adjust:

  Multiplicity correction for `p_adjusted`, one of the methods
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) accepts.

- ...:

  Unused.

- object:

  A fitted `multilpa_transitions` model.

- scale:

  `"natural"` (probabilities) or `"unconstrained"` (logits).

- parm:

  Parameters to report, as names or positions of the natural
  coefficients; all by default.

## Value

A base `data.frame` with one row per natural parameter and the columns
of
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md):
`level`, `outcome`, `term`, `parameter`, `estimate`, `standard_error`,
`statistic`, `p_value`, `p_adjusted`, `conf_low` and `conf_high`, named
as [`coef()`](https://rdrr.io/r/stats/coef.html) names them. An initial
probability (`parameter = "initial_probability"`, `level = "profile"`)
has `outcome` the profile and `term` the group class; a transition
probability (`"transition_probability"`) has `outcome` the profile moved
to and `term` `group_class_h:profile_j`, the group class and the profile
moved from. Probabilities and variances carry no Wald test. Attributes
`covariance` (natural scale), `covariance_unconstrained`, `vcov_type`
and `level`.

## Conditions

`latents_no_converge` for an unconverged fit; `latents_boundary_fit`
when a variance, a response, initial or transition probability, or a
group-class share sits at its bound, or a transition row was never
informed (a profile no group occupies before its last occasion);
`latents_singular_information` when the information cannot be inverted;
`latents_bad_inference_data` when a supplied `data` does not reproduce
the fit; `latents_too_few_groups` for robust or OPG errors with no more
groups than parameters.

## Examples

``` r
fit <- lta(subset(course_engagement, student <= 40),
           c("browse", "lectures", "forum_read"), "student",
           n_profiles = 2, time = "sequence", n_starts = 2, seed = 1)
parameter_inference(fit)
#>          level       outcome                    term              parameter
#> 1  measurement     profile_1                  browse                   mean
#> 2  measurement     profile_1                lectures                   mean
#> 3  measurement     profile_1              forum_read                   mean
#> 4  measurement     profile_2                  browse                   mean
#> 5  measurement     profile_2                lectures                   mean
#> 6  measurement     profile_2              forum_read                   mean
#> 7  measurement     profile_1                  browse               variance
#> 8  measurement     profile_1                lectures               variance
#> 9  measurement     profile_1              forum_read               variance
#> 10 measurement     profile_2                  browse               variance
#> 11 measurement     profile_2                lectures               variance
#> 12 measurement     profile_2              forum_read               variance
#> 13       group group_class_1                    <NA>            probability
#> 14     profile     profile_1           group_class_1    initial_probability
#> 15     profile     profile_2           group_class_1    initial_probability
#> 16     profile     profile_1 group_class_1:profile_1 transition_probability
#> 17     profile     profile_2 group_class_1:profile_1 transition_probability
#> 18     profile     profile_1 group_class_1:profile_2 transition_probability
#> 19     profile     profile_2 group_class_1:profile_2 transition_probability
#>      estimate standard_error  statistic       p_value    p_adjusted    conf_low
#> 1   0.5552352     0.04015129  13.828579  1.713727e-43  1.713727e-43  0.47654014
#> 2   0.4372193     0.04891436   8.938464  3.946152e-19  3.946152e-19  0.34134888
#> 3   0.6089409     0.03909145  15.577343  1.037800e-54  1.037800e-54  0.53232308
#> 4  -0.8706168     0.05194498 -16.760366  4.756892e-63  4.756892e-63 -0.97242711
#> 5  -0.6264837     0.05713002 -10.965928  5.572617e-28  5.572617e-28 -0.73845648
#> 6  -0.9319730     0.04355606 -21.397095 1.421972e-101 1.421972e-101 -1.01734134
#> 7   0.5490615     0.04216140         NA            NA            NA  0.47234443
#> 8   0.8253164     0.06254016         NA            NA            NA  0.71140814
#> 9   0.4991390     0.03955163         NA            NA            NA  0.42733899
#> 10  0.4601966     0.04940331         NA            NA            NA  0.37287618
#> 11  0.5736107     0.06114026         NA            NA            NA  0.46546715
#> 12  0.3259582     0.03539700         NA            NA            NA  0.26346736
#> 13  1.0000000     0.00000000         NA            NA            NA  1.00000000
#> 14  0.7421284     0.07242704         NA            NA            NA  0.57817735
#> 15  0.2578716     0.07242704         NA            NA            NA  0.14199410
#> 16  0.8906008     0.01904335         NA            NA            NA  0.84732958
#> 17  0.1093992     0.01904335         NA            NA            NA  0.07727388
#> 18  0.1896193     0.03271224         NA            NA            NA  0.13357337
#> 19  0.8103807     0.03271224         NA            NA            NA  0.73793148
#>     conf_high
#> 1   0.6339303
#> 2   0.5330896
#> 3   0.6855587
#> 4  -0.7688065
#> 5  -0.5145109
#> 6  -0.8466047
#> 7   0.6382387
#> 8   0.9574634
#> 9   0.5830025
#> 10  0.5679657
#> 11  0.7068797
#> 12  0.4032710
#> 13  1.0000000
#> 14  0.8580059
#> 15  0.4218226
#> 16  0.9227261
#> 17  0.1526704
#> 18  0.2620685
#> 19  0.8664266
```
