# Relate latent class membership to external covariates

`r3step()` examines how external covariates are associated with
membership in latent classes identified by a previously fitted model. It
addresses whether the relative likelihood of belonging to particular
classes varies with characteristics that were not used to define them.
The measurement parameters remain fixed, preserving the class
definitions while estimating their associations with the covariates.

Because class membership is unobserved, a regression that treats
assigned classes as known outcomes can distort these associations. The
three-step approach incorporates estimated classification-error
probabilities into the membership regression, treating each assignment
as an imperfect indicator of latent membership. Covariates enter a
multinomial logistic model whose coefficients describe changes in the
log odds of membership in a class relative to a reference class, holding
the other covariates constant.

For supported multilevel models, the regression can concern individual
profiles or group classes. For supported transition models, a separate
specification relates covariates to initial profile membership and
movement between profiles over time. The function reports regression
coefficients, standard errors, confidence intervals, and hypothesis
tests, with options for clustered inference and adjustment for multiple
comparisons.

## Usage

``` r
r3step(
  x,
  data,
  covariates = character(),
  level = c("individuals", "groups"),
  ci_level = 0.95,
  vcov_type = c("robust", "observed"),
  adjust = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "none"),
  by_group_class = FALSE,
  initial_covariates = character(),
  transition_covariates = character()
)
```

## Arguments

- x:

  A covariate-free fit from
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  (`family = "profiles"`, the default) or from
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) with its
  default homogeneous first-order transitions and invariant measurement,
  or a trajectory model from
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
  that classifies persons (with or without `random`), whose person
  classes are then predicted. Fits with membership covariates and
  multilevel growth mixtures (`cluster`) raise
  `latents_unsupported_three_step`: the correction uses a fixed,
  unconditional classification-error matrix and independent units. The
  other families of
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  general transition fits of
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) and
  mixture regressions with observation classes raise
  `latents_unsupported_three_step` too, except that a covariate-free
  general transition fit takes `initial_covariates` and
  `transition_covariates`.

- data:

  The data frame carrying the covariates, in the fit's row order.
  Columns shared with the fit are checked row by row; a frame with too
  little identifying information raises a warning that alignment cannot
  be verified.

- covariates:

  Character vector of numeric columns, without missing or infinite
  values, that are not measurement indicators of `x`. For
  `level = "groups"` each must be constant within a group; for a
  trajectory model each must be constant within a person and not part of
  the trajectory regression (`latents_bad_covariate`).

- level:

  `"individuals"` (default) predicts profile membership, `"groups"`
  group-class membership. For a trajectory model the rows have
  `level = "persons"` and outcomes `class_1`, `class_2`, ...; `"groups"`
  is refused there with `latents_bad_argument`.

- ci_level:

  Confidence level for the intervals (default 0.95).

- vcov_type:

  `"robust"` (default) uses the sandwich clustered on the fit's groups
  with the small-sample (CR1) factor `G / (G - 1)` for `G` groups, as
  the default of
  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md).
  `"observed"` uses the observed information; for the pooled regression
  it treats every observation as independent, and with
  `by_group_class = TRUE` it is the model-based variance of the
  two-level likelihood. `"robust"` requires more independent groups than
  the regression has coefficients, because the group score contributions
  sum to zero at the estimate; otherwise it is refused with
  `latents_too_few_groups`, and `"observed"` remains available.

- adjust:

  Multiplicity correction across the covariate terms of every
  non-reference class (not the intercepts), passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html): `"BH"`
  (default), `"holm"`, `"hochberg"`, `"hommel"`, `"bonferroni"`, `"BY"`
  or `"none"`.

- by_group_class:

  `FALSE` (default) fits one pooled multinomial logit with a single
  intercept per profile. `TRUE` fits the two-level form of the one-step
  model: each group class has its own profile intercepts, the covariate
  slopes are shared across group classes, and a group's observations
  share one latent group class whose proportions are estimated jointly.
  Group classes represent between-group differences in profile
  membership that are not explained by the covariates. This changes the
  regression specification relative to the pooled analysis. Available
  only for `level = "individuals"` on a covariate-free
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  fit; refused with `latents_bad_argument` otherwise.

- initial_covariates, transition_covariates:

  For an [`lta()`](https://pak.dynasite.org/latents/reference/lta.md)
  fit with one group class: covariates of the profile at the first
  occasion (constant within a person) and of every transition (read at
  the destination occasion, so they may vary over occasions). Either one
  fits the three-step latent Markov model described below instead of the
  occasion-by-occasion profile regression; `covariates` must then be
  empty, and `level` and `by_group_class` do not apply. Refused with
  `latents_bad_argument` on a fit that is not a transition model.

## Value

A base `data.frame` with one row per non-reference class and term, and
the columns `level`, `outcome`, `term`, `estimate`, `standard_error`,
`statistic`, `p_value`, `p_adjusted`, `conf_low` and `conf_high`.
`p_value` is uncorrected; `p_adjusted` applies `adjust` (also recorded
in the `adjust` attribute) and is `NA` on the intercept rows. A
coefficient with zero standard error has `NA_real_` statistic and
p-value. Coefficients are log odds against the last class, as in
`multilpa(profile_covariates = )`; the reference class is named in the
`reference_class` attribute and the variance used in the `vcov_type`
attribute. With `by_group_class = TRUE` there is one intercept row per
group class (`group_class_1`, `group_class_2`, ...), followed by the
shared slopes and, at `level = "groups"`, the group-class logits against
the last group class; the result also carries the attributes
`by_group_class` and `log_likelihood` (the maximised log likelihood).

With `initial_covariates` or `transition_covariates`, one row per
coefficient of the step-three chain: `block` (`"initial"` or
`"transition"`), `from` (the origin profile of a transition, `NA` for an
initial logit), `to`, `term`, `estimate`, `standard_error`, `statistic`,
`p_value`, `p_adjusted`, `conf_low` and `conf_high`, with the attributes
`log_likelihood`, `classification_errors` (the step-two error matrix,
rows the true profile), `vcov_type` and `adjust`. Initial logits are
against the last profile and transition logits against staying, as in
`lta(initial_covariates =, transition_covariates =)`; `p_adjusted`
covers the covariate terms.

## Details

On an [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit,
`level = "individuals"` models the marginal profile at each occasion; it
does not regress initial-state or transition probabilities on the
covariates. `vcov_type = "robust"` clusters the occasion scores by
group; `"observed"` treats them as independent. On a trajectory model
the units are persons, each with one covariate value, and `"robust"` is
the sandwich over persons.

The three-step latent Markov model (`initial_covariates`,
`transition_covariates`) classifies each occasion from its own
indicators, using the transition fit's measurement model and its profile
shares pooled over occasions (the latent class model of the pooled
data), and assigns it to its modal profile. Step three is a hidden
Markov model for those assignments, with the classification errors fixed
at their step-two values and the initial and transition logits regressed
on the covariates (Di Mari, Oberski and Vermunt 2016). The transition
fit's own posteriors are not used for the assignments, because they draw
on a person's other occasions and their errors are therefore not
independent across occasions, as step three assumes. The step-one fit
must have no covariates and one group class
(`latents_unsupported_three_step` otherwise). `"robust"` uses a sandwich
covariance estimator with persons as independent units; `"observed"`
uses the observed information.

Because uncertainty in the estimated classification-error matrix is
omitted, interval coverage may be below the nominal level when class
separation is weak. Inspect the estimated error matrix with
`get_results(x, "classification_errors")`.

## Interpretation

The correction assumes that classification-error probabilities remain
constant across covariate values within each latent class. Reported
standard errors and confidence intervals treat those probabilities as
fixed and omit uncertainty from their estimation and from selecting the
class solution.

In the ordinary membership regression, the last class is the reference.
A positive coefficient indicates higher odds of membership in the
specified class relative to that reference; exponentiating the
coefficient gives the corresponding odds ratio. Transition regressions
instead compare movement to another profile with remaining in the origin
profile. These associations do not, by themselves, establish causal
relationships.

## References

Vermunt, J. K. (2010). Latent class modeling with covariates: two
improved three-step approaches. *Political Analysis*, 18, 450–469.
Asparouhov, T., & Muthen, B. (2014). Auxiliary variables in mixture
modeling: three-step approaches using Mplus. *Structural Equation
Modeling*, 21, 329–341. Di Mari, R., Oberski, D. L., & Vermunt, J. K.
(2016). Bias-adjusted three-step latent Markov modeling with covariates.
*Structural Equation Modeling*, 23, 649–660.

## See also

[`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
for a distal outcome, and `multilpa(profile_covariates = )` for the
one-step alternative that estimates everything jointly.

## Examples

``` r
set.seed(21)
g <- rep(seq_len(50), each = 12)
x <- rnorm(600)
truth <- 1L + as.integer(runif(600) < plogis(-0.3 + 1.2 * x))
example_data <- data.frame(
  g = g, x = x,
  a = rnorm(600, ifelse(truth == 2L, 1.2, -1.2)),
  b = rnorm(600, ifelse(truth == 2L, 1.2, -1.2))
)
fit <- multilpa(example_data, c("a", "b"), "g", n_profiles = 2,
                n_group_classes = 1, n_starts = 4, seed = 1)
r3step(fit, example_data, "x")
#>         level   outcome        term   estimate standard_error statistic
#> 1 individuals profile_1 (Intercept) -0.5703958      0.1191748 -4.786212
#> 2 individuals profile_1           x  1.2188658      0.1268352  9.609839
#>        p_value   p_adjusted   conf_low  conf_high
#> 1 1.699583e-06           NA -0.8039741 -0.3368175
#> 2 7.266244e-22 7.266244e-22  0.9702734  1.4674582

# Three-step regression of transition probabilities on covariates.
moves <- lta(course_engagement,
             vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
             id = "student", n_profiles = 2, time = "sequence", n_starts = 2, seed = 1)
r3step(moves, course_engagement, transition_covariates = "sequence")
#>        block      from        to        term     estimate standard_error
#> 1    initial      <NA> profile_1 (Intercept)  0.547155214     0.20869296
#> 2 transition profile_1 profile_2 (Intercept) -1.850705801     0.28098481
#> 3 transition profile_1 profile_2    sequence -0.003170648     0.03170193
#> 4 transition profile_2 profile_1 (Intercept) -1.630640960     0.30703090
#> 5 transition profile_2 profile_1    sequence  0.011857160     0.03415178
#>    statistic      p_value p_adjusted   conf_low   conf_high
#> 1  2.6218192 8.746182e-03         NA  0.1381245  0.95618591
#> 2 -6.5864976 4.503228e-11         NA -2.4014259 -1.29998569
#> 3 -0.1000143 9.203329e-01  0.9203329 -0.0653053  0.05896400
#> 4 -5.3109995 1.090256e-07         NA -2.2324105 -1.02887145
#> 5  0.3471901 7.284485e-01  0.9203329 -0.0550791  0.07879342
```
