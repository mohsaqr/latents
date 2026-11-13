# Relate latent class membership to external covariates

`r3step()` regresses membership in the classes of a fitted model on
covariates that were not used to define them, holding the measurement
model fixed. Treating the assigned classes as known outcomes distorts
the associations; the three-step approach instead treats each assignment
as an imperfect indicator of the true class, with the
classification-error probabilities of the fitted model built into a
multinomial logit. Coefficients are changes in the log odds of a class
against the reference class, the last one, holding the other covariates
constant; their exponentials are odds ratios. A transition fit can also
have its initial state and transitions regressed on covariates, a
separate specification. The correction assumes the classification-error
probabilities do not vary with the covariates within a class; the
standard errors treat them as fixed and ignore the choice of class
solution, so interval coverage can fall below the nominal level when the
classes are poorly separated; inspect the error matrix with
`get_results(x, "classification_errors")`. Associations are not causal
effects.

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

  A covariate-free fit without weights or a noise component: from
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  (the profile family), from
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md), or a
  trajectory model from
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
  that classifies persons. A default
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit takes
  `covariates` or the Markov covariates; one with occasion-specific or
  second-order transitions, occasion-specific measurement, or ordinal or
  count indicators takes the Markov covariates only. Refused: fits with
  membership covariates, multilevel growth mixtures (`cluster`), the
  other
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  families, observation-level mixture regressions and stayer-class fits.

- data:

  The data frame carrying the covariates, in the fit's row order;
  columns shared with the fit are checked row by row.

- covariates:

  Names of numeric columns, complete and finite, that are not indicators
  of `x`. For `level = "groups"` each must be constant within a group;
  for a trajectory model, constant within a person and not part of the
  trajectory regression.

- level:

  `"individuals"` (default) predicts profile membership, `"groups"`
  group-class membership. A trajectory model predicts its person
  classes, reported as `level = "persons"` with outcomes `class_1`,
  `class_2` and so on.

- ci_level:

  Confidence level of the intervals (default 0.95).

- vcov_type:

  `"robust"` (default): the sandwich clustered on the independent units
  (the fit's groups; the persons of a trajectory model or a transition
  chain), which needs more units than the regression has coefficients,
  `"observed"` remaining available otherwise. `"observed"`: the observed
  information, treating every observation as independent in the pooled
  regression; with `by_group_class = TRUE`, the model-based variance of
  the two-level likelihood.

- adjust:

  Multiplicity correction across the covariate terms of every
  non-reference class (not the intercepts), a
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) method;
  `"BH"` by default.

- by_group_class:

  `FALSE` (default): one pooled multinomial logit with one intercept per
  profile. `TRUE`: the two-level form of the one-step model, with
  profile intercepts per group class, slopes shared across group
  classes, and the group-class proportions estimated jointly; the group
  classes then capture between-group differences in membership the
  covariates do not explain. Only for `level = "individuals"` on a
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  fit.

- initial_covariates, transition_covariates:

  For an [`lta()`](https://pak.dynasite.org/latents/reference/lta.md)
  fit with one group class: covariates of the profile at the first
  occasion (constant within a person) and of every transition (read at
  the destination occasion, so they may vary over occasions). Either one
  fits the three-step latent Markov model (see *Transition fits*)
  instead of the occasion-by-occasion profile regression; `covariates`,
  `level` and `by_group_class` are refused then.

## Value

A data frame with one row per non-reference class and term: `level`,
`outcome`, `term`, `estimate`, `standard_error`, `statistic`, `p_value`,
`p_adjusted` (`NA` on the intercept rows), `conf_low` and `conf_high`.
Coefficients are log odds against the last class, named in the
`reference_class` attribute; `vcov_type` and `adjust` are attributes
too. With `by_group_class = TRUE` there is one intercept row per group
class and non-reference profile, then the shared slopes and, at
`level = "groups"`, the group-class logits, with the attributes
`by_group_class` and `log_likelihood`.

With `initial_covariates` or `transition_covariates`, one row per
coefficient of the step-three chain: `block` (`"initial"` or
`"transition"`), `from` (the origin profile, `NA` for an initial logit),
`to`, `term` and the estimate columns (`p_adjusted` on the covariate
terms), with the attributes `log_likelihood`, `classification_errors`
(the step-two error matrix, rows the true profile), `vcov_type` and
`adjust`. Initial logits are against the last profile, transition logits
against staying, as in
`lta(initial_covariates =, transition_covariates =)`.

## Transition fits

With `covariates` and `level = "individuals"`, an
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit is
treated occasion by occasion: the covariates predict the marginal
profile at each occasion, with the scores clustered by group under
`"robust"`. That is not a model of the initial state or the transitions.

With `initial_covariates` or `transition_covariates`, each occasion is
classified from its own indicators by the transition fit's measurement
model (with the profile shares pooled over occasions) and assigned to
its modal profile; step three is a hidden Markov model for those
assignments with the classification errors fixed and the initial and
transition logits regressed on the covariates. The transition fit's own
posteriors are not used, because they draw on a person's other occasions
and their errors would not be independent across occasions as step three
assumes. The step-one fit must have one group class and no covariates.
The step-three chain is homogeneous and first order, even when the
step-one transition structure varies over occasions or has higher order.
It uses one classification-error matrix pooled over occasions. Step-one
occasion-specific measurement is refused because its error probabilities
can change over occasions.

## Boundary estimates

A step-three logit with a coefficient exceeding 20 in centred, scaled
design coordinates, or an average class probability at or below 1e-6
over the contributing design rows, is treated as separated or on a
boundary. Wald inference is refused for that regression, under both
covariance types; changing the covariance estimator does not make the
boundary regular.

## See also

[`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
for a distal outcome; `multilpa(profile_covariates = )` for the one-step
alternative.

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
