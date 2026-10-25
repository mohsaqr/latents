# Standard errors and intervals for a latent transition model

Returns standard errors, tests and intervals for every parameter of an
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit: the
measurement model, the group-class shares, each group class's initial
profile distribution and its transition matrix. The scores follow the
Fisher identity (each is an expected count from the forward-backward
pass minus the count implied by the fitted probabilities), and the
observed information is the numerical Jacobian of that analytic score.

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
  max_iter = 5000L,
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

  Optional data frame the model was fitted to. Not needed for the Wald
  path, since the fit stores its indicators, codes, groups and
  occasions; when given it must reproduce them. Required for
  `method = "bootstrap"`.

- level:

  Confidence level of the intervals. Default `0.95`.

- step:

  Finite-difference step for the observed information.

- vcov_type:

  `"observed"` (default; observed information), `"robust"` (the sandwich
  clustered on groups, which are the sequences) or `"opg"` (the outer
  product of the group scores).

- adjust:

  Multiplicity correction for `p_adjusted`: `"none"` (default) or one of
  the [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html)
  methods.

- method:

  `"wald"` (default) or `"bootstrap"`. The bootstrap resamples persons
  with replacement, refits the same specification and matches each
  replicate's profiles and group classes to the original's;
  `standard_error` is the standard deviation of the replicates,
  `conf_low` and `conf_high` their percentile interval, and `statistic`,
  `p_value` and `p_adjusted` are `NA`.

- iter, n_starts, max_iter, tol, seed:

  Number of bootstrap resamples, the controls of each refit, and an
  optional seed (the caller's random state is restored).

- boundary:

  Accepted for consistency with the generic. Transition inference holds
  no parameter at a bound: without an active bound `"fix"` gives the
  same table as `"error"`, and a fit with an active bound raises
  `latents_boundary_fit` under either.

- object:

  A fitted `multilpa_transitions` model.

- scale:

  `"natural"` (probabilities) or `"unconstrained"` (logits).

- ...:

  For [`confint()`](https://rdrr.io/r/stats/confint.html), further
  arguments to
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md),
  such as `vcov_type`. [`vcov()`](https://rdrr.io/r/stats/vcov.html)
  takes none and warns (`latents_unused_argument`) about any it is
  given.

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
moved from. Probabilities and variances have no Wald test. The Wald
table has attributes `covariance` (natural scale),
`covariance_unconstrained`, `level`, `step`, `method` and `vcov_type`.

## Conditions

`latents_no_converge` for an unconverged fit; `latents_boundary_fit`
when a variance, a response, initial or transition probability, or a
group-class share is at its bound, or a transition row has no
information (a profile no group occupies before its last occasion);
`latents_singular_information` when the information cannot be inverted;
`latents_bad_inference_data` when a supplied `data` does not reproduce
the fit; `latents_too_few_groups` for robust or OPG errors with no more
groups than parameters; `latents_bad_argument` for
`method = "bootstrap"` without `data`; `latents_bootstrap_failed` when
fewer than two resamples give a usable fit, with a
`latents_bootstrap_dropped` warning for resamples that were dropped.

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
#>      estimate standard_error  statistic       p_value    p_adjusted   conf_low
#> 1   0.5552195     0.04015156  13.828093  1.725335e-43  1.725335e-43  0.4765239
#> 2   0.4372023     0.04891447   8.938097  3.959281e-19  3.959281e-19  0.3413317
#> 3   0.6089173     0.03909229  15.576403  1.053164e-54  1.053164e-54  0.5322979
#> 4  -0.8706368     0.05194461 -16.760870  4.716735e-63  4.716735e-63 -0.9724464
#> 5  -0.6264885     0.05713057 -10.965907  5.573909e-28  5.573909e-28 -0.7384623
#> 6  -0.9319820     0.04355621 -21.397222 1.418100e-101 1.418100e-101 -1.0173506
#> 7   0.5490715     0.04216225         NA            NA            NA  0.4723529
#> 8   0.8253244     0.06254078         NA            NA            NA  0.7114150
#> 9   0.4991643     0.03955397         NA            NA            NA  0.4273601
#> 10  0.4601888     0.04940225         NA            NA            NA  0.3728703
#> 11  0.5736161     0.06114138         NA            NA            NA  0.4654707
#> 12  0.3259595     0.03539718         NA            NA            NA  0.2634684
#> 13  1.0000000     0.00000000         NA            NA            NA  1.0000000
#> 14  0.7421425     0.07242582         NA            NA            NA  0.5781919
#> 15  0.2578575     0.07242582         NA            NA            NA  0.1419835
#> 16  0.8906184     0.01904205         NA            NA            NA  0.8473493
#> 17  0.1093816     0.01904205         NA            NA            NA  0.0772590
#> 18  0.1895947     0.03271108         NA            NA            NA  0.1335518
#> 19  0.8104053     0.03271108         NA            NA            NA  0.7379573
#>     conf_high
#> 1   0.6339151
#> 2   0.5330729
#> 3   0.6855368
#> 4  -0.7688272
#> 5  -0.5145146
#> 6  -0.8466133
#> 7   0.6382505
#> 8   0.9574727
#> 9   0.5830329
#> 10  0.5679556
#> 11  0.7068876
#> 12  0.4032728
#> 13  1.0000000
#> 14  0.8580165
#> 15  0.4218081
#> 16  0.9227410
#> 17  0.1526507
#> 18  0.2620427
#> 19  0.8664482
```
