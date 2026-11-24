# Compare latent classes on an outcome not used to define them, correcting for classification error

`three_step()` estimates how an outcome that was not used to define the
classes differs between them: the mean of a numeric outcome in each
class, or the probability of each category of a categorical one. The
measurement model is held fixed, so the classes keep their meaning.

Comparing classes by modal assignment treats uncertain assignments as
certain and biases the differences. The three-step methods correct for
the classification error of the fitted model: step one is the fitted
model, step two assigns each unit to its most probable class and
measures how often each true class is assigned to another, and step
three relates the outcome to the classes with those error rates built
in. The default is the Bolck-Croon-Hagenaars (`"bch"`) weighting; the
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
  step_one = c("fixed", "estimated"),
  missing = c("error", "exclude")
)
```

## Arguments

- x:

  The fitted model whose classes are compared. It must have no
  membership covariates. Accepted:

  - a profile or class fit from
    [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
    [`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
    or
    [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
    (the default `family = "profiles"`);

  - an [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit
    with its default homogeneous first-order transitions, invariant
    measurement and continuous or categorical indicators;

  - a trajectory model that classifies persons:
    [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
    with `class_level = "group"`, with or without `random`, fitted
    without `membership` covariates and without `cluster`.

  Refused with `latents_unsupported_three_step`: a fit with membership
  covariates; the other
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  families (additive, dispersion and cross-level); an
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit with
  occasion-specific, second-order or covariate transitions, initial
  covariates, occasion-specific measurement, a stayer class, or ordinal
  or count indicators; a mixture regression that classifies
  observations; and a multilevel growth mixture (`cluster`). A fit with
  sampling weights raises `latents_unsupported_weights`, and one with a
  noise component `latents_unsupported_noise`.

- data:

  The data frame carrying the outcome, in the fit's row order; columns
  shared with the fit are checked row by row, and a frame that does not
  match raises `latents_bad_inference_data`. For a trajectory model, the
  data the model was fitted to.

- outcome:

  Name of the outcome column. It must not be an indicator of `x`, must
  have no infinite values, and must have no missing values unless
  `missing = "exclude"`. A numeric column gives a mean in each class; a
  factor, character or logical column gives the probability of each
  category in each class. For `level = "groups"` the outcome must be
  constant within each group; for a trajectory model, within each
  person.

- level:

  Which classes the outcome is related to:

  - `"individuals"` (default): the profiles (classes) of the
    observations;

  - `"groups"`: the group classes of a two-level fit.

  A trajectory model has person classes only: leave `level` at its
  default (the result reports `level = "persons"`); `"groups"` raises
  `latents_bad_argument`.

- method:

  How step three relates the outcome to the classes:

  - `"bch"` (default): Bolck-Croon-Hagenaars weights. Each unit gets one
    weight per class, taken from the inverse of the classification-error
    matrix, and the class means are weighted means;

  - `"ml"`: the maximum-likelihood step-three model. The assigned class
    is an indicator of the true class with the step-one error
    probabilities, and the outcome has a normal (numeric) or multinomial
    (categorical) distribution within each true class;

  - `"proportional"`: each unit is weighted by its posterior class
    probabilities, without correction for classification error;

  - `"modal"`: each unit counts fully in its most probable class,
    without correction.

  The last two are the uncorrected comparisons the correction is judged
  against. See *Methods*.

- ci_level:

  Confidence level of the intervals (default 0.95).

- contrast:

  What is reported:

  - `"none"` (default): the outcome mean of each class (for a
    categorical outcome, the probability of each category in each
    class), with a confidence interval and no test;

  - `"pairs"`: the difference between every pair of classes, each tested
    against zero;

  - `"overall"`: one Wald test that every class has the same mean (the
    same category distribution).

- adjust:

  Multiplicity correction of the p-values across the pairwise
  differences, giving `p_adjusted`. One of the
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) methods:

  - `"BH"` (default): Benjamini-Hochberg, controls the false discovery
    rate;

  - `"holm"`, `"hochberg"`, `"hommel"`, `"bonferroni"`: control the
    family-wise error rate;

  - `"BY"`: Benjamini-Yekutieli, controls the false discovery rate under
    any dependence;

  - `"none"`: no correction (`p_adjusted` equals `p_value`).

  Used only with `contrast = "pairs"`. With a binary outcome, its two
  complementary category rows count as one test per class pair; with
  more categories, every reported category difference enters the family.

- vcov_type:

  Which units the sandwich variance treats as independent:

  - `"cluster"` (default): the fit's groups (the persons of a transition
    fit or a trajectory model), as the nesting requires. For a
    single-level fit every observation is its own group;

  - `"independent"`: every observation, ignoring the nesting.

- variance:

  The outcome variances of the `"ml"` model for a numeric outcome:

  - `"varying"` (default): one variance in each class;

  - `"equal"`: one variance shared by the classes.

  Supplying it with another method or with a categorical outcome raises
  `latents_bad_argument`.

- step_one:

  Whether the standard errors include the uncertainty of the step-one
  estimates:

  - `"fixed"` (default): they condition on the step-one estimates;

  - `"estimated"`: the uncertainty of the classification-error
    probabilities (and, under `"ml"`, of the class sizes) is added to
    the step-three covariance. Available for
    [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
    [`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
    and
    [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
    fits at either level, when their Wald inference is available and
    their measurement was estimated jointly. A transition fit, a
    trajectory model, a staged fit and a fit holding a measurement block
    fixed are refused.

  See *Step-one uncertainty*.

- missing:

  What to do with units whose outcome is missing:

  - `"error"` (default): a missing outcome value raises
    `latents_bad_argument`;

  - `"exclude"`: the units without the outcome are left out of step
    three, as Latent GOLD's `missing excludeall` does. The step-one fit
    and its posteriors are unchanged, and the classification errors,
    weights and standard errors are those of the units analysed. This
    assumes the outcome is missing at random given the classes. At
    `level = "groups"` a group's outcome must be missing on all its rows
    or on none.

  An outcome missing for every unit is refused under either value.

## Value

A data frame of class `latents_table` (a `data.frame` that prints as a
formatted table). Its rows depend on `contrast`:

- `contrast = "none"`: one row per class, or per class and category for
  a categorical outcome. Columns: `level`, `method`, `class`, `category`
  (categorical outcomes only), `estimate` (the class mean, or the
  probability of the category in the class), `standard_error`,
  `conf_low`, `conf_high`, `variance` (the outcome variance in the
  class; `method = "ml"` with a numeric outcome only) and `effective_n`
  (Kish's effective number of units behind the class's weights; under
  `"ml"` the expected number of units in the class). Means are not
  tested, and a probability's interval is formed on the logit scale;

- `contrast = "pairs"`: one row per pair of classes, or per pair and
  category. Columns: `level`, `method`, `class`, `reference_class`,
  `category` (categorical outcomes only), `estimate` (class minus
  reference), `standard_error`, `statistic` (the estimate over its
  standard error), `p_value`, `p_adjusted`, `conf_low` and `conf_high`;

- `contrast = "overall"`: one row. Columns: `level`, `method`,
  `statistic` (the Wald chi-square), `df` (`K - 1` for a mean with `K`
  classes, `(K - 1)(C - 1)` for `C` categories) and `p_value`.

The attributes `vcov_type` and `step_one` record those choices; `adjust`
is recorded on a pairs result and `variance` on a numeric ML result with
`contrast = "none"`.

Errors: `latents_unsupported_three_step`, `latents_unsupported_weights`
and `latents_unsupported_noise` for a fit the correction does not cover
(see `x`); `latents_bad_outcome` for an outcome that is an indicator of
`x` (also of class `multilpa_indicator_reused`), that varies within a
group or person where one value is needed, or that takes a single
category (under `"ml"`, a single value); `latents_bad_data` when `data`
has no column named by `outcome`; `latents_bad_inference_data` when
`data` is not the fit's data; `latents_bad_argument` for any other
invalid argument; `latents_inseparable_classes` when the
classification-error matrix cannot be inverted or a single class is
asked for a contrast; and `latents_too_few_groups` when there are no
more independent units than classes. Warnings:
`latents_bch_out_of_range`, `latents_boundary` and
`latents_small_classes`, described under *Methods* and *Standard
errors*, and `latents_unconverged` when the `"ml"` model did not
converge.
[latents-conditions](https://pak.dynasite.org/latents/reference/latents-conditions.md)
describes each class.

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
with a `latents_bch_out_of_range` warning. The ML probabilities stay in
range; where the BCH shares do not, the ML maximum usually sits on the
bound and is reported without standard errors, with a `latents_boundary`
warning.

On an [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit,
`level = "individuals"` uses one row per occasion and relates the
outcome to the profile at that occasion, with the group-clustered
variance accounting for the repeated occasions. It is a marginal
profile-outcome analysis, not a model of outcome effects on transitions.

## Standard errors

The variance is a sandwich over the `G` independent units (the groups,
or the observations under `vcov_type = "independent"`) with the
`G / (G - 1)` small-sample factor. It needs more units than classes;
with as many or fewer, the analysis is refused
(`latents_too_few_groups`). A class carried by one unit (Kish's
effective number of units carrying it below 1.5) has no standard error:
its interval and the contrasts involving it are `NA`, with a
`latents_small_classes` warning.

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
[`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md);
`get_results(x, "classification_errors")` for the error matrix the
correction uses.

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
# The mean of `y` in each class, then the difference between the classes.
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
# An outcome missing for some units: they are left out of step three.
example_data <- transform(example_data, y_seen = replace(y, seq_len(25), NA))
three_step(fit, example_data, "y_seen", missing = "exclude")
#> Distal outcome by class (bch, 95% CI)
#> 
#> Class    Mean            95% CI    SE  Effective n
#> -----  ------  ----------------  ----  -----------
#>     1  -0.012  [-0.435,  0.412]  0.22       167.09
#>     2   9.949  [ 9.496, 10.402]  0.23       162.74
```
