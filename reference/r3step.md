# Covariates predicting class membership, corrected for misclassification

The R3STEP approach. A multinomial logit of class membership on
covariates, fitted after the measurement model rather than alongside it,
but treating the assigned class as an error-prone indicator of the true
one with the error rates held fixed at what step one found. Regressing
the modal class directly instead attenuates every coefficient, because
some units are in the wrong class and the covariate cannot explain why.

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

  A covariate-free fitted
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  or [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) model,
  or a trajectory model from
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
  (persons classified, with or without `random`), whose person classes
  are then predicted. A fit that already uses membership covariates
  cannot be corrected by the fixed, unconditional classification-error
  matrix used here, and a multilevel growth mixture (`cluster`) is
  refused because its persons are not independent; both raise
  `latents_unsupported_three_step`.

- data:

  The data frame carrying the covariates, in the fit's row order. Shared
  fit columns are checked row by row; a frame with too little
  identifying information warns that alignment cannot be verified.

- covariates:

  Character vector of numeric columns not used as measurement indicators
  in `x`. For `level = "groups"` each must be constant within a group;
  for a trajectory model each must be constant within a person and not
  part of the trajectory regression (`latents_bad_covariate`).

- level:

  `"individuals"` predicts profile membership, `"groups"` predicts
  group-class membership. A trajectory model classifies persons: its
  rows have `level = "persons"` and outcomes `class_1`, `class_2`, ...,
  and `"groups"` is refused with `latents_bad_argument`.

- ci_level:

  Confidence level for the intervals.

- vcov_type:

  `"robust"`, the default, uses the sandwich clustered on the fit's
  groups, with the small-sample (CR1) factor `G / (G - 1)` for `G`
  groups, which is the honest choice when the covariates are measured on
  observations nested inside them, and matches the default of
  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md).
  `"observed"` uses the observed information; for the pooled regression
  it treats every observation as independent, and with
  `by_group_class = TRUE` it is the model-based variance of the
  two-level likelihood. `"robust"` needs more independent groups than
  the regression has coefficients, because the group score contributions
  sum to zero at the estimate and so span at most one dimension fewer
  than there are groups; with too few it is refused with
  `latents_too_few_groups` rather than reporting a variance that is
  singular, or, with a single group, numerically zero. `"observed"`
  remains available there.

- adjust:

  Multiplicity correction applied across the covariate terms, passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). `"BH"`
  by default; `"none"` leaves the p-values uncorrected. The family is
  every covariate term of every non-reference class, which is the set of
  tests this call computes; the intercepts are not part of it.

- by_group_class:

  `FALSE`, the default, fits one pooled multinomial logit with a single
  intercept per profile. `TRUE` fits the two-level form of the one-step
  model instead: every group class has its own profile intercepts, the
  covariate slopes are shared across group classes, and a group's
  observations share one latent group class whose proportions are
  estimated alongside. Differences between groups that the covariates do
  not explain are then absorbed by the group classes rather than left in
  the residual, where the pooled regression attenuates the slopes and
  inflates their clustered standard errors. For `level = "individuals"`
  on a covariate-free
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  fit only; refused with `latents_bad_argument` otherwise.

- initial_covariates, transition_covariates:

  For an [`lta()`](https://pak.dynasite.org/latents/reference/lta.md)
  fit with one group class: covariates of the profile at the first
  occasion (constant within a person) and of every transition (read at
  the destination occasion, so they may change over occasions). Either
  one makes this the three-step latent Markov model described below, in
  place of the occasion-by-occasion profile regression of `covariates`,
  which must then be left empty; `level` and `by_group_class` do not
  apply. Refused with `latents_bad_argument` on a fit that is not a
  transition model.

## Value

A base `data.frame` with one row per non-reference class and term, and
the columns `level`, `outcome`, `term`, `estimate`, `standard_error`,
`statistic`, `p_value`, `p_adjusted`, `conf_low` and `conf_high`.
`p_value` is uncorrected and `p_adjusted` carries the correction named
by `adjust`, which is also recorded in the result's `adjust` attribute;
it is `NA` on the intercept rows, which are not part of the tested
family. A coefficient with zero standard error has an undefined Wald
statistic and p-value, reported as `NA_real_` rather than infinity or
zero. Coefficients are log odds against the final class, which is the
reference, matching `multilpa(profile_covariates = )`; that class is
named in the result's `reference_class` attribute, and the variance that
was used in its `vcov_type` attribute. With `by_group_class = TRUE` the
intercept rows are one per group class, named `group_class_1`,
`group_class_2`, ..., followed by the shared slopes and, at
`level = "groups"`, the group-class logits against the last group class;
the result also carries `by_group_class` and the maximised
`log_likelihood` as attributes.

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

The error matrix is held fixed rather than estimated jointly, which is
what makes this a three-step method and what keeps the covariates from
reshaping the classes. On a
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit,
`level = "individuals"` models the marginal profile at each occasion. It
does not regress initial-state or transition probabilities on the
covariates. `vcov_type = "robust"` clusters the occasion scores by
group; `"observed"` treats them as independent. On a trajectory model
the units are persons, each with one covariate value, and `"robust"` is
the sandwich over persons.

The three-step latent Markov model (`initial_covariates`,
`transition_covariates`) classifies each occasion from its own
indicators, with the transition fit's measurement and its profile shares
pooled over occasions (the latent class model of the pooled data), and
assigns it to its modal profile. Step three is a hidden Markov model for
those assignments, with the classification errors held at their step-two
values and the initial and transition logits regressed on the covariates
(Di Mari, Oberski and Vermunt 2016). The transition fit's own posteriors
are not used for the assignments: they borrow from a person's other
occasions, so their errors are not independent across occasions as step
three assumes. The step-one fit must have no covariates of its own and
one group class (`latents_unsupported_three_step` otherwise). It agrees
with Latent GOLD 6.1's `step3 modal ml` latent Markov model to its
printed precision on three designs (equivalence/latentgold-step3-markov
in the source repository); `"robust"` is the sandwich over persons and
`"observed"` the information Latent GOLD reports.

Over 60 replications of a two-profile design with a true log-odds slope
of 1.2 and intercept -0.3, this recovered the slope with a bias of
-0.002 and the intercept with a bias of +0.006, where a logistic
regression on the modal class was biased by -0.163 on the slope. The
correction costs some precision: the standard deviation of the slope
across replications was 0.127 against 0.100 for the naive fit. The
nominal 95% intervals covered the true slope in 95% of those
replications.

That coverage was measured where the profiles separate well. The error
matrix is still treated as known rather than estimated, so intervals
should be expected to run narrow where classification is poorer; check
`get_results(x, "classification_errors")` before relying on them.

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
#>         level   outcome        term  estimate standard_error statistic
#> 1 individuals profile_1 (Intercept) -0.570396      0.1191748 -4.786214
#> 2 individuals profile_1           x  1.218866      0.1268352  9.609839
#>        p_value   p_adjusted   conf_low  conf_high
#> 1 1.699572e-06           NA -0.8039742 -0.3368177
#> 2 7.266225e-22 7.266225e-22  0.9702734  1.4674582

# Covariates of the moves between profiles, three steps on from lta().
moves <- lta(course_engagement,
             vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
             id = "student", n_profiles = 2, time = "sequence", n_starts = 2, seed = 1)
r3step(moves, course_engagement, transition_covariates = "sequence")
#>        block      from        to        term     estimate standard_error
#> 1    initial      <NA> profile_1 (Intercept)  0.547153583     0.20869295
#> 2 transition profile_1 profile_2 (Intercept) -1.850702888     0.28098436
#> 3 transition profile_1 profile_2    sequence -0.003170733     0.03170190
#> 4 transition profile_2 profile_1 (Intercept) -1.630641206     0.30703086
#> 5 transition profile_2 profile_1    sequence  0.011857195     0.03415177
#>    statistic      p_value p_adjusted    conf_low   conf_high
#> 1  2.6218116 8.746377e-03         NA  0.13812293  0.95618424
#> 2 -6.5864980 4.503217e-11         NA -2.40142210 -1.29998367
#> 3 -0.1000171 9.203307e-01  0.9203307 -0.06530532  0.05896385
#> 4 -5.3110010 1.090247e-07         NA -2.23241063 -1.02887178
#> 5  0.3471912 7.284476e-01  0.9203307 -0.05507904  0.07879343
```
