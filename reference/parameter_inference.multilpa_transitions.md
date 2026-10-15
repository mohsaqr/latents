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
  method = c("wald", "bootstrap"),
  iter = 199L,
  n_starts = 10L,
  max_iter = 1000L,
  tol = 1e-10,
  seed = NULL,
  boundary = c("error", "fix")
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

- method:

  `"wald"` (the default) or `"bootstrap"`. The bootstrap resamples
  persons with replacement, refits the same specification, and matches
  each replicate's profiles and group classes to the original's;
  `standard_error` is the replicates' standard deviation, `conf_low` and
  `conf_high` their percentile interval, and `statistic`, `p_value` and
  `p_adjusted` are `NA`. It needs `data`.

- iter, n_starts, max_iter, tol, seed:

  Bootstrap resamples, the controls of each refit, and an optional seed
  (the caller's random state is restored).

- boundary:

  Accepted for the generic. Transition inference holds nothing at a
  bound: with no bound active `"fix"` gives the same table as `"error"`,
  and a fit with an active bound raises `latents_boundary_fit` under
  either.

- object:

  A fitted `multilpa_transitions` model.

- scale:

  `"natural"` (probabilities) or `"unconstrained"` (logits).

- ...:

  Unused.

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
groups than parameters; `latents_bad_argument` for
`method = "bootstrap"` without `data`.

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
#> 1   0.5552219     0.04015152  13.828167  1.723569e-43  1.723569e-43  0.47652636
#> 2   0.4372049     0.04891445   8.938153  3.957283e-19  3.957283e-19  0.34133430
#> 3   0.6089209     0.03909217  15.576546  1.050818e-54  1.050818e-54  0.53230168
#> 4  -0.8706338     0.05194466 -16.760794  4.722712e-63  4.722712e-63 -0.97244346
#> 5  -0.6264878     0.05713048 -10.965911  5.573650e-28  5.573650e-28 -0.73846148
#> 6  -0.9319806     0.04355619 -21.397204 1.418628e-101 1.418628e-101 -1.01734920
#> 7   0.5490699     0.04216212         NA            NA            NA  0.47235159
#> 8   0.8253232     0.06254069         NA            NA            NA  0.71141396
#> 9   0.4991604     0.03955361         NA            NA            NA  0.42735689
#> 10  0.4601900     0.04940240         NA            NA            NA  0.37287118
#> 11  0.5736153     0.06114120         NA            NA            NA  0.46547012
#> 12  0.3259593     0.03539715         NA            NA            NA  0.26346821
#> 13  1.0000000     0.00000000         NA            NA            NA  1.00000000
#> 14  0.7421403     0.07242601         NA            NA            NA  0.57818966
#> 15  0.2578597     0.07242601         NA            NA            NA  0.14198507
#> 16  0.8906157     0.01904225         NA            NA            NA  0.84734633
#> 17  0.1093843     0.01904225         NA            NA            NA  0.07726128
#> 18  0.1895985     0.03271126         NA            NA            NA  0.13355510
#> 19  0.8104015     0.03271126         NA            NA            NA  0.73795330
#>     conf_high
#> 1   0.6339174
#> 2   0.5330754
#> 3   0.6855402
#> 4  -0.7688241
#> 5  -0.5145141
#> 6  -0.8466121
#> 7   0.6382487
#> 8   0.9574712
#> 9   0.5830283
#> 10  0.5679571
#> 11  0.7068864
#> 12  0.4032725
#> 13  1.0000000
#> 14  0.8580149
#> 15  0.4218103
#> 16  0.9227387
#> 17  0.1526537
#> 18  0.2620467
#> 19  0.8664449
```
