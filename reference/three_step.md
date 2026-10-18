# Estimate a distal outcome's mean in each latent class

Estimates the mean of an outcome in each latent class of a fitted model,
or the difference between every pair of classes, by the three-step
approach, and returns a `data.frame`. The measurement model is fitted
first, units are assigned to classes, and the outcome is then related to
the assignments with a correction for classification error, so the
outcome does not influence the classes. `"bch"` weights give
approximately unbiased class means; regressing on the modal class
attenuates the differences (Vermunt 2010; Bakk and Vermunt 2016).

## Usage

``` r
three_step(
  x,
  data,
  outcome,
  level = c("individuals", "groups"),
  method = c("bch", "proportional", "modal"),
  ci_level = 0.95,
  contrast = c("none", "pairs"),
  adjust = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "none"),
  vcov_type = c("cluster", "independent")
)
```

## Arguments

- x:

  A covariate-free fitted
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  or [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) model,
  or a trajectory model that classifies persons
  ([`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
  with `class_level = "group"`, with or without `random` effects) fitted
  without `membership` covariates. Fits with membership covariates and
  multilevel growth mixtures (`cluster`) are refused with
  `latents_unsupported_three_step`: the correction uses a fixed,
  unconditional classification-error matrix and independent units.

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

  `"bch"` (default) applies the Bolck-Croon-Hagenaars weights;
  `"proportional"` weights by the posterior probabilities (attenuated);
  `"modal"` uses the uncorrected modal assignment, for comparison.

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
  unit, which is appropriate only without nesting (for example a fit
  with a single group). The choice is recorded in the result's
  `vcov_type` attribute.

## Value

A base `data.frame`. For `contrast = "none"`, one row per class with the
columns `level`, `method`, `class`, `estimate`, `standard_error`,
`conf_low`, `conf_high` and `effective_n`; class means are not tested.
For `contrast = "pairs"`, one row per pair of classes with the columns
`level`, `method`, `class`, `reference_class`, `estimate` (the mean of
`class` minus the mean of `reference_class`), `standard_error`,
`statistic`, `p_value`, `p_adjusted`, `conf_low` and `conf_high`, as in
[`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md). A
contrast with zero standard error has `NA_real_` statistic and p-values.

The result's `vcov_type` attribute records the variance used; for
`contrast = "pairs"` the `adjust` attribute records the correction
applied to `p_adjusted`.

## Details

The correction assumes the outcome is independent of the assigned class
given the true class. On an
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit,
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
single group and singular with as many groups as classes. Both cases are
refused with `latents_too_few_groups`; `vcov_type = "independent"` gives
the unclustered variance.

The classification-error matrix is treated as known, so the intervals
are too narrow. Over 60 replications of a two-profile design with 8%
misclassification and a true class difference of 10, `"bch"` had a bias
of -0.13 (against -1.28 for `"modal"` and -1.82 for `"proportional"`)
and its nominal 95% intervals covered the true difference in 85% of
replications; the intervals of the other two methods covered it in none.

## References

Vermunt, J. K. (2010). Latent class modeling with covariates: two
improved three-step approaches. *Political Analysis*, 18, 450–469. Bakk,
Z., & Vermunt, J. K. (2016). Robustness of stepwise latent class
modeling with continuous distal outcomes. *Structural Equation
Modeling*, 23, 20–31.

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
#>     1  -0.069  [-0.463,  0.326]  0.20       182.75
#>     2   9.972  [ 9.531, 10.412]  0.22       170.37
three_step(fit, example_data, "y", contrast = "pairs")
#> Distal outcome: class differences (bch, 95% CI)
#> 
#> Class  vs  Difference         95% CI      p  p (BH)
#> -----  --  ----------  -------------  -----  ------
#>     2   1       10.04  [9.36, 10.72]  <.001   <.001
```
