# Compare latent classes on an external outcome

`three_step()` estimates differences in an external outcome between
latent classes defined by a previously fitted model. It addresses
whether classes distinguished by their indicator patterns also differ on
a variable that was not used to identify them. The measurement model
remains fixed, so the outcome comparison retains the original meaning of
the classes.

Latent class membership is inferred rather than observed. Assigning each
observation to its most probable class and comparing the resulting
groups treats those assignments as certain, which can distort estimated
outcome differences. Three-step methods address this problem by using
classification-error probabilities from the fitted model when estimating
class-specific outcome means. The default Bolck–Croon–Hagenaars
(`"bch"`) method applies correction weights; the alternative
maximum-likelihood (`"ml"`) method incorporates the classification-error
probabilities into a model for the outcome.

The analysis provides estimated outcome means within classes or
differences between pairs of classes, together with standard errors and
confidence intervals. Pairwise differences can also be tested against
zero, with adjustment for multiple comparisons. For supported multilevel
models, comparisons may concern individual profiles or group classes.
Uncorrected modal-assignment and posterior-weighted analyses are
available as additional methods.

## Usage

``` r
three_step(
  x,
  data,
  outcome,
  level = c("individuals", "groups"),
  method = c("bch", "proportional", "modal", "ml"),
  ci_level = 0.95,
  contrast = c("none", "pairs"),
  adjust = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "none"),
  vcov_type = c("cluster", "independent"),
  variance = c("varying", "equal")
)
```

## Arguments

- x:

  A covariate-free fit from
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  (`family = "profiles"`, the default) or from
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) with its
  default homogeneous first-order transitions and invariant measurement,
  or a trajectory model that classifies persons
  ([`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
  with `class_level = "group"`, with or without `random` effects) fitted
  without `membership` covariates. Fits with membership covariates and
  multilevel growth mixtures (`cluster`) raise
  `latents_unsupported_three_step`: the correction uses a fixed,
  unconditional classification-error matrix and independent units. The
  other families of
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fits from
  the general transition engine (occasion-varying, second-order or
  covariate transitions, occasion-specific measurement, a stayer class,
  ordinal or count indicators) and mixture regressions with observation
  classes raise `latents_unsupported_three_step` too.

- data:

  The data frame carrying the outcome, in the fit's row order. Columns
  shared with the fit are checked row by row; a frame with too little
  identifying information raises a warning that alignment cannot be
  verified. For a trajectory model it must be the data frame the model
  was fitted to.

- outcome:

  Name of a numeric outcome column, without missing or infinite values,
  that is not a measurement indicator of `x`. For `level = "groups"` it
  must be constant within each group, and for a trajectory model within
  each person.

- level:

  `"individuals"` (default) relates the outcome to profiles, `"groups"`
  to group classes. For a trajectory model the outcome is related to its
  person classes, reported with `level = "persons"`; `"groups"` is
  refused there with `latents_bad_argument`.

- method:

  `"bch"` (default) applies the Bolck-Croon-Hagenaars weights; `"ml"`
  fits the maximum-likelihood step-three model (see *Details*);
  `"proportional"` weights by the posterior probabilities without
  correction for classification error; `"modal"` uses the uncorrected
  modal assignment, for comparison.

- ci_level:

  Confidence level for the intervals (default 0.95).

- contrast:

  `"none"` (default) reports one outcome mean per class; `"pairs"`
  reports the difference between every pair of classes, with a test
  against zero.

- adjust:

  Multiplicity correction across the pairwise differences, passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html): `"BH"`
  (default), `"holm"`, `"hochberg"`, `"hommel"`, `"bonferroni"`, `"BY"`
  or `"none"`. Ignored for `contrast = "none"`.

- vcov_type:

  `"cluster"` (default) sums the influence contributions within the
  fit's groups and treats the groups as the independent units, as
  required when the outcome is measured on observations nested in
  groups. `"independent"` treats every observation as an independent
  unit, which is appropriate only when observations can be treated as
  independent. The choice is recorded in the result's `vcov_type`
  attribute.

- variance:

  For `method = "ml"`: `"varying"` (default) estimates a variance of the
  outcome in each class, `"equal"` one variance shared by the classes.
  Given with another method it raises `latents_bad_argument`.

## Value

A base `data.frame`. For `contrast = "none"`, one row per class with the
columns `level`, `method`, `class`, `estimate`, `standard_error`,
`conf_low`, `conf_high` and `effective_n`; class means are not tested.
For `method = "ml"`, a `variance` column (the outcome's estimated
variance in the class) comes before `effective_n`, which is then the
expected number of units in the class under the step-three model; the
`variance` attribute records the choice. For `contrast = "pairs"`, one
row per pair of classes with the columns `level`, `method`, `class`,
`reference_class`, `estimate` (the mean of `class` minus the mean of
`reference_class`), `standard_error`, `statistic`, `p_value`,
`p_adjusted`, `conf_low` and `conf_high`, as in
[`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md). A
contrast with zero or missing standard error has `NA_real_` statistic
and p-values.

The result's `vcov_type` attribute records the variance used; for
`contrast = "pairs"` the `adjust` attribute records the correction
applied to `p_adjusted`.

## Details

On an [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit,
`level = "individuals"` uses one row per occasion and relates the
outcome to the profile at that occasion, with the default
group-clustered variance accounting for the repeated occasions within a
group. This is a marginal profile-outcome analysis, not a model of
outcome effects on transitions or of a sequence-level outcome.

The cluster-robust variance is the sum of one outer product per
independent unit, scaled by the small-sample (CR1) factor `G / (G - 1)`
for `G` units (groups, or observations under
`vcov_type = "independent"`). These contributions sum to zero at the
estimate, so the variance has rank at most `G - 1`: it is zero with a
single group and singular with as many groups as classes. Both cases
raise `latents_too_few_groups`; `vcov_type = "independent"` gives the
unclustered variance. For the same reason a class carried by a single
unit (all its weight in one group, or under `vcov_type = "independent"`
in one observation; Kish's effective number of carrying units below 1.5)
has no standard error, whichever the method: its standard error, its
interval and those of the contrasts involving it are `NA`, with a
`latents_small_classes` warning.

With `method = "ml"` the assigned class is an indicator of the true
class whose error probabilities are those of step one, and the outcome
is normal within each true class. The log likelihood
`sum_i log sum_t pi_t D[t, W_i] phi(z_i; mu_t, sigma_t)` is maximised
over the class means and variances by EM followed by Newton steps, with
the class sizes `pi_t` fixed at the average posterior probabilities and
the classification errors `D` fixed, as in Latent GOLD's step-three
model for a dependent variable. Multiple initialisations are used
because the likelihood may have several local maxima. Equal-variance
fits additionally use bounded BFGS optimisation from quantile-based
initial values; accepted solutions must converge, have a sufficiently
small score and yield a positive-definite information matrix. These
searches do not guarantee the global maximum. Twenty random starts use a
fixed internal seed, and the caller's `.Random.seed` is restored. R's
hidden Box-Muller normal cache cannot be restored; see
[base::Random](https://rdrr.io/r/base/Random.html). Standard errors use
a sandwich covariance estimator, clustered as for the other methods. As
described above, a class supported by a single independent unit has no
standard error. Maximum-likelihood estimates depend on the assumed
within-class outcome distribution. Skewness or an incorrect
equal-variance constraint can bias estimated means; BCH is less
sensitive to these distributional assumptions (Bakk and Vermunt, 2016).

Because uncertainty in the estimated classification-error probabilities
is omitted, interval coverage can be below the nominal level,
particularly when classification is uncertain.

## Interpretation

Correction for classification error depends on the fitted measurement
model and assumes that the outcome and assigned class are independent
conditional on latent class membership. The maximum-likelihood method
additionally assumes a normal outcome distribution within each class.
Reported uncertainty is conditional on the estimated
classification-error probabilities and does not include uncertainty from
estimating those probabilities or selecting the class solution. Outcome
differences describe associations; they do not establish causal effects
of class membership.

## References

Vermunt, J. K. (2010). Latent class modeling with covariates: two
improved three-step approaches. *Political Analysis*, 18, 450–469. Bakk,
Z., Tekle, F. B., & Vermunt, J. K. (2013). Estimating the association
between latent class membership and external variables using
bias-adjusted three-step approaches. *Sociological Methodology*, 43,
272–311. Bakk, Z., & Vermunt, J. K. (2016). Robustness of stepwise
latent class modeling with continuous distal outcomes. *Structural
Equation Modeling*, 23, 20–31.

## Examples

``` r
set.seed(11)
group <- rep(seq_len(40), each = 10)
truth <- 1L + as.integer(runif(400) > 0.5)
example_data <- data.frame(
  g = group,
  a = rnorm(400, ifelse(truth == 2L, 1.1, -1.1)),
  b = rnorm(400, ifelse(truth == 2L, 1.1, -1.1)),
  y = rnorm(400, ifelse(truth == 2L, 10, 0))
)
fit <- multilpa(example_data, c("a", "b"), "g", n_profiles = 2,
                n_group_classes = 1, n_starts = 4, seed = 1)
three_step(fit, example_data, "y")
#> Distal outcome by class (bch, 95% CI)
#> 
#> Class    Mean            95% CI    SE  Effective n
#> -----  ------  ----------------  ----  -----------
#>     1  -0.069  [-0.463,  0.325]  0.20       182.74
#>     2   9.972  [ 9.531, 10.412]  0.22       170.37
three_step(fit, example_data, "y", contrast = "pairs")
#> Distal outcome: class differences (bch, 95% CI)
#> 
#> Class  vs  Difference         95% CI      p  p (BH)
#> -----  --  ----------  -------------  -----  ------
#>     2   1       10.04  [9.37, 10.72]  <.001   <.001
```
