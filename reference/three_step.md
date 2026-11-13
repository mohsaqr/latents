# Compare latent classes on an external outcome

`three_step()` estimates how an outcome that was not used to define the
classes differs between them, holding the measurement model fixed so the
classes keep their meaning. Comparing classes by modal assignment treats
uncertain assignments as certain and biases the differences; the
three-step methods correct for the classification error of the fitted
model. The default is the Bolck-Croon-Hagenaars (`"bch"`) weighting; the
maximum-likelihood (`"ml"`) method models the outcome given the true
class instead. Both assume that the outcome and the assigned class are
independent given the true class; the differences are associations, not
effects of class membership.

## Usage

``` r
three_step(
  x,
  data,
  outcome,
  level = c("individuals", "groups"),
  method = c("bch", "proportional", "modal", "ml"),
  ci_level = 0.95,
  contrast = c("none", "pairs", "overall"),
  adjust = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "none"),
  vcov_type = c("cluster", "independent"),
  variance = c("varying", "equal"),
  step_one = c("fixed", "estimated")
)
```

## Arguments

- x:

  A covariate-free fit: from
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  (the profile family), from
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) with its
  default homogeneous first-order transitions and invariant measurement,
  or a trajectory model that classifies persons
  ([`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
  with `class_level = "group"`, with or without `random`) fitted without
  `membership` covariates. Refused: fits with membership covariates,
  multilevel growth mixtures (`cluster`), the other
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  families, observation-level mixture regressions, fits with weights or
  a noise component, and
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fits with
  occasion-specific or covariate transitions, second order, occasion-
  specific measurement, a stayer class or ordinal or count indicators.

- data:

  The data frame carrying the outcome, in the fit's row order; columns
  shared with the fit are checked row by row. For a trajectory model,
  the data the model was fitted to.

- outcome:

  Name of the outcome column, without missing or infinite values and not
  an indicator of `x`: numeric for class means, or a factor, character
  or logical column for the probability of each category in each class.
  For `level = "groups"` it must be constant within each group; for a
  trajectory model, within each person.

- level:

  `"individuals"` (default) relates the outcome to the profiles,
  `"groups"` to the group classes. A trajectory model relates it to its
  person classes, reported as `level = "persons"`.

- method:

  `"bch"` (default): Bolck-Croon-Hagenaars weights. `"ml"`: the
  maximum-likelihood step-three model. `"proportional"`: posterior
  weights without correction. `"modal"`: modal assignment without
  correction. See *Methods*.

- ci_level:

  Confidence level of the intervals (default 0.95).

- contrast:

  `"none"` (default): one outcome mean per class (one probability per
  class and category). `"pairs"`: the difference between every pair of
  classes, tested against zero. `"overall"`: one Wald test that every
  class has the same mean (the same category distribution).

- adjust:

  Multiplicity correction across the pairwise differences, a
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) method;
  `"BH"` by default. Ignored for `contrast = "none"` or `"overall"`.
  With a binary outcome, its two complementary category rows count as
  one test per class pair; with more categories, every reported category
  difference enters the family.

- vcov_type:

  `"cluster"` (default) treats the fit's groups as the independent units
  of the sandwich variance, as the nesting requires; `"independent"`
  treats every observation as independent.

- variance:

  For `method = "ml"` with a numeric outcome (refused otherwise):
  `"varying"` (default) estimates an outcome variance in each class,
  `"equal"` one shared by the classes.

- step_one:

  `"fixed"` (default) conditions on the step-one estimates.
  `"estimated"` adds the uncertainty of the classification-error
  probabilities (and, under `"ml"`, of the class sizes) to the
  step-three covariance; available for
  [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
  [`lca()`](https://pak.dynasite.org/latents/reference/lca.md) and
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  profile fits at either level when their Wald inference is available
  and their measurement was estimated jointly. Staged fits and fits
  holding a measurement block fixed are refused. See *Step-one
  uncertainty*.

## Value

A data frame. For `contrast = "none"`, one row per class (per class and
category for a categorical outcome, with `category` after `class`):
`level`, `method`, `class`, `estimate`, `standard_error`, `conf_low`,
`conf_high` and `effective_n`; means are not tested, and a probability's
interval is formed on the logit scale. Under `method = "ml"` with a
numeric outcome a `variance` column precedes `effective_n`, and under
`"ml"` `effective_n` is the expected number of units in the class. For
`contrast = "pairs"`, one row per pair (per pair and category): `level`,
`method`, `class`, `reference_class`, `estimate` (class minus
reference), `standard_error`, `statistic`, `p_value`, `p_adjusted`,
`conf_low` and `conf_high`, as in
[`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md). For
`contrast = "overall"`, one row: `level`, `method`, `statistic`, `df`
(`K - 1` for a mean, `(K - 1)(C - 1)` for `C` categories) and `p_value`.
The attributes `vcov_type` and `step_one` record the choices, `adjust`
on a pairs result and `variance` on a numeric ML result with
`contrast = "none"`.

## Methods

BCH gives each observation one weight per class, the row of the inverse
classification-error matrix for its assigned class, and takes weighted
means. With `method = "ml"` the assigned class is an indicator of the
true class with the step-one error probabilities, and the outcome is
normal within each true class (a multinomial distribution for a
categorical outcome), with the class sizes fixed at the average
posteriors, as in Latent GOLD's step-three model. The likelihood is
maximised from a fixed set of starts plus twenty seeded random ones, so
the result does not depend on the caller's seed; the global maximum is
not guaranteed. ML depends on the within-class outcome distribution:
skewness or a wrong equal-variance constraint can bias the means, and
BCH is less sensitive to that.

A categorical outcome is described by its probability in each class. BCH
weights can be negative, so with poorly separated classes a share can
fall outside 0 to 1; it is then reported without an interval or test,
with a warning. The ML probabilities stay in range; where the BCH shares
do not, the ML maximum usually sits on the bound and is reported without
standard errors, with a warning.

On an [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit,
`level = "individuals"` uses one row per occasion and relates the
outcome to the profile at that occasion, with the group-clustered
variance accounting for the repeated occasions. It is a marginal
profile-outcome analysis, not a model of outcome effects on transitions.

## Standard errors

The variance is a sandwich over the `G` independent units (the groups,
or the observations under `vcov_type = "independent"`) with the
`G / (G - 1)` small-sample factor. It needs more units than classes;
with as many or fewer, the analysis is refused. A class carried by one
unit (Kish's effective number of units carrying it below 1.5) has no
standard error: its interval and the contrasts involving it are `NA`,
with a warning.

## Step-one uncertainty

With `step_one = "fixed"` the intervals ignore the uncertainty of the
classification-error probabilities, which lowers their coverage when
classification is uncertain. `"estimated"` adds the first-order
contribution of the step-one covariance, with the data and modal
assignments held fixed; the estimates are unchanged, the standard errors
and intervals grow, and the tests use the corrected covariance. The
first-stage contribution uses its model-based observed information, even
when the step-three sandwich is clustered; it does not add a robust
first-stage covariance for residual dependence outside the fitted model.
The correction is a large-sample one: in a small, poorly separated
sample the intervals still cover less than the nominal level. The modal
method is unaffected, as it uses no posteriors. A staged fit or a fit
with a held measurement block does not carry the measurement covariance
needed by this correction; use `step_one = "fixed"` to report inference
conditional on that measurement. Neither option includes the uncertainty
of having selected the class solution.

## See also

[`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md) for
covariates of class membership;
[`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md).

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
# A categorical outcome: its probability in each class.
example_data <- transform(example_data, passed = y + rnorm(400, sd = 6) > 5)
three_step(fit, example_data, "passed")
#> Distal outcome by class (bch, 95% CI)
#> 
#> Class  Category  Probability        95% CI     SE  Effective n
#> -----  --------  -----------  ------------  -----  -----------
#>     1  FALSE            0.80  [0.73, 0.85]  0.030       182.74
#>     1  TRUE             0.20  [0.15, 0.27]  0.030       182.74
#>     2  FALSE            0.17  [0.12, 0.23]  0.030       170.37
#>     2  TRUE             0.83  [0.77, 0.88]  0.030       170.37
three_step(fit, example_data, "passed", contrast = "overall")
#> Distal outcome: equal across classes (bch, Wald)
#> 
#>  Wald  df      p
#> -----  --  -----
#> 99.11   1  <.001
```
