# Tidy inference for a covariate fit

For a `multilpa(profile_covariates = )` fit, returns standard errors,
tests and intervals for every free parameter at all three levels: the
measurement model, the profile logits and the group-class logits.

Returns a table of the estimated parameters of a fitted model with their
standard errors, Wald tests and confidence intervals. By default the
standard errors are maximum-likelihood errors from the inverse observed
information, with delta-method Wald intervals; `method = "bootstrap"`
resamples groups instead. The estimates are local and asymptotic,
conditional on the selected numbers of classes, and do not account for
model selection. Boundary solutions and unidentified Hessians are
refused by the Wald path.

## Usage

``` r
# S3 method for class 'multilpa_additive'
parameter_inference(
  x,
  data = NULL,
  level = 0.95,
  step = 1e-04,
  vcov_type = c("observed", "robust", "opg"),
  adjust = .multilpa_p_adjust_methods,
  method = c("wald", "bootstrap"),
  ...
)

# S3 method for class 'multilpa_covariates'
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
  tol = 1e-08,
  seed = NULL,
  boundary = c("error", "fix")
)

# S3 method for class 'multilpa_cross_level'
parameter_inference(x, ...)

# S3 method for class 'latents_growth_mixture'
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
  tol = 1e-08,
  seed = NULL,
  boundary = c("error", "fix")
)

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
  tol = 1e-08,
  seed = NULL,
  boundary = c("error", "fix")
)

# S3 method for class 'multilpa'
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
  tol = 1e-08,
  seed = NULL,
  boundary = c("error", "fix")
)

# S3 method for class 'multilpa_lta'
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
```

## Arguments

- x:

  A converged `multilpa` fit with inactive variance bounds.

- data:

  Optional data frame the model was fitted to. When omitted it is
  rebuilt from the indicators, identifiers and occasions the fit stores.
  A supplied frame is checked exactly against the stored indicators and
  identifiers; for an older fit without stored indicators, the
  likelihood gives a weaker consistency check.

- level:

  Confidence level, strictly between zero and one. Default `0.95`.

- step:

  Positive finite-difference step for the analytic score derivative.
  Default `1e-4`.

- vcov_type:

  Covariance estimator for the Wald path. `"observed"` (default) inverts
  the observed information. `"robust"` returns the Huber-White sandwich
  `A^-1 B A^-1`, where `B` is the outer product of the per-group scores,
  multiplied by the CR1 small-sample factor `G / (G - 1)` for `G` groups
  (the default of
  [`sandwich::vcovCL()`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html),
  applied to every robust standard error in the package). Mplus's MLR
  applies no such factor, so these errors are `sqrt(G / (G - 1))` times
  Mplus's. Groups are the independent units: robust errors relax the
  assumption that the within-group model is correctly specified, not the
  assumption that groups are independent. `"opg"` inverts the outer
  product of the per-group scores (the BHHH estimator), as `glca`
  reports; it agrees with `"observed"` in large samples and can differ
  in small ones, where `"observed"` is usually more reliable. `"robust"`
  and `"opg"` need more groups than parameters. For a fit with sampling
  weights the default is `"robust"`, and any other value raises
  `latents_unsupported_weights`. `vcov_type` and `step` are not used by
  `method = "bootstrap"`.

- adjust:

  Multiplicity correction applied to `p_value` to give `p_adjusted`:
  `"none"` (default) or one of the
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) methods.
  The method applied is stored in the `adjust` attribute. The correction
  is taken over the tests the table reports; bounded parameters carry no
  test and are not counted.

- method:

  `"wald"` (default) inverts the observed information in the estimation
  coordinates; it is available for all fourteen covariance structures.
  EEI, VVI, EEE and VVV are estimated in log variances or log-Cholesky
  coordinates; the ten structures that constrain volume, shape or
  orientation across profiles are estimated in log volumes, log shapes
  with the determinant-one constraint built into the coordinates, and
  Cayley-transform orientations (or, for VEE and EVV, determinant-one
  log-Cholesky factors), so the information has one row per free
  parameter. The reported variances and covariances are obtained from
  these coordinates by the delta method, and their natural-scale
  covariance is singular, as for the class probabilities. `"bootstrap"`
  resamples groups with replacement, refits within the same covariance
  structure, realigns each refit's labels, and reports the percentile
  interval and the standard deviation of the replicates; it is also
  available for all fourteen structures and assumes independent groups,
  as `vcov_type = "robust"` does.

- ...:

  Further arguments; unused by most methods.

- iter:

  Number of resamples for `method = "bootstrap"`, at least two. Default
  `199`; a 95% percentile interval is more stable with 999 or more.

- n_starts, max_iter, tol:

  Passed to each bootstrap refit, as in
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).

- seed:

  Optional seed for the resampling: a whole number accepted by
  [`set.seed()`](https://rdrr.io/r/base/Random.html); any other value
  raises `latents_bad_argument` on either path. The Wald path uses no
  random numbers. The caller's random state is restored on exit.

- boundary:

  Handling of categorical response probabilities at their lower bound,
  `min_probability`, where the information is singular. `"error"`
  (default) raises `latents_boundary_fit`. `"fix"` holds those
  probabilities at the bound and reports every other parameter
  conditional on them, as Mplus does; the held probabilities are omitted
  from the table (or, for a covariate fit, reported with `NA` standard
  errors) and named in the `fixed_at_bound` attribute.

## Value

For a covariate fit, a base `data.frame` with one row per free parameter
and the same columns as for every other class: `level` (`"measurement"`,
`"profile"` or `"group"`), `outcome` (the profile or group class the
coefficient predicts, or the profile a measurement parameter belongs
to), `term`, `parameter` (`"mean"`, `"variance"`, `"covariance"`,
`"response"` or `"coefficient"`), `estimate`, `standard_error`,
`statistic`, `p_value`, `p_adjusted`, `conf_low` and `conf_high`.
Variances and covariances are in their natural units, with standard
errors obtained by the delta method from the log and log-Cholesky
coordinates on which they are estimated. Each categorical indicator
contributes one `"response"` row per profile, indicator and category,
named `indicator:category`: the category probability, with its standard
error from the multinomial logits and no Wald test.

A base `data.frame` with one row per reported parameter and the columns
`level` (`"measurement"`, `"profile"` or `"group"`), `outcome`, `term`,
`parameter`, `estimate`, `standard_error`, `statistic`, `p_value`,
`p_adjusted`, `conf_low` and `conf_high`. Estimates are on the natural
scale (variances in their own units, probabilities as probabilities). A
parameter whose null value lies on the boundary of its space (a
variance, a probability, a covariance diagonal) has `NA` `statistic`,
`p_value` and `p_adjusted`, but still has an interval. All class
probabilities are reported; their sum constraints make the natural
covariance singular. Wald intervals for probabilities are formed on the
logit scale and for variances (and covariance diagonals) on the log
scale, then back-transformed, so they stay within the parameter's range;
every other interval is `estimate +/- z * standard_error`.
[`confint()`](https://rdrr.io/r/stats/confint.html) returns the same
intervals. A probability at or below 1e-6 is treated as on its bound
(refused, or held by `boundary = "fix"`).

A fit made with `fixed` reports only the parameters it estimated: a held
measurement block contributes no row to the table or to the information
matrix. The table has one row per estimated natural coefficient:
`n_parameters` rows plus one for each set of probabilities whose
reference category the estimation scale drops. Standard errors,
intervals and p-values are conditional on the held measurement solution
and do not include its uncertainty (see
[`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)).
For a
[`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
result, `method = "bootstrap"` refits both stages in each resample, so
the table includes the measurement parameters and that uncertainty is
propagated (attribute `stages_resampled` is 2).

With `method = "bootstrap"` the table has the same columns and rows;
`estimate` is the fitted value, `standard_error` the standard deviation
of the replicates and `conf_low`/`conf_high` the percentile interval.
`statistic`, `p_value` and `p_adjusted` are `NA`. Its attributes are
`covariance` (the covariance of the replicates, returned by
[`vcov()`](https://rdrr.io/r/stats/vcov.html)), `level`, `method`,
`iter`, `n_valid`, `replicates` (the kept replicates, one row each),
`messages` (why a resample was dropped, `NA` where it was not),
`structure` and `fixed`; `covariance_unconstrained` is `NULL`.

For the Wald path, diagnostics of the whole fit are attributes:
`covariance` and `covariance_unconstrained` (natural and
estimation-scale covariance of the estimates), `hessian`, `gradient`,
`scaled_score`, `score_displacement` (the largest score times its
estimation-scale standard error, which the `latents_unconverged` warning
uses), `condition_ratio`, `level`, `step`, `method`, `vcov_type`,
`fixed` (the measurement blocks the fit held,
[`character()`](https://rdrr.io/r/base/character.html) when none),
`fixed_at_bound` (the probabilities held by `boundary = "fix"`) and
`adjust`. When `vcov_type` is `"robust"` or `"opg"`, `group_scores`
holds the groups-by-parameters score matrix and `scaling_correction` the
MLR scaling correction factor `tr(A^-1 B) / q` (with `B` without the CR1
factor), used for scaled likelihood-ratio difference tests; it is `NA`
otherwise.

## Details

Class-membership coefficients of a covariate fit are on the multinomial
logit scale, with the final profile and the final group class as
references, so each coefficient is a log odds against that reference.
The tests are Wald tests, which are invalid for a coefficient driven to
the boundary by separation or an empty class. A fit with a membership
logit above 20 in coordinates independent of the covariates' units and
origin (roughly, a slope per standard deviation of its covariate and an
intercept at the covariates' means) raises `latents_boundary_fit`; the
fit itself warns with `latents_extreme_coefficients`. Only
`method = "wald"` is implemented for a covariate fit; `"bootstrap"`
raises `latents_unsupported_inference`.

For a growth mixture
([`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
with `random`), the table has one row per natural parameter: `block`
(`"coefficient"`, `"membership"` and, for a multilevel growth mixture,
`"group_membership"`, `"residual_sd"`, `"random_variance"`,
`"random_sd"`, `"random_correlation"`, `"share"`), `class` (`"all"` for
a quantity shared by every class), `term`, `estimate`, `standard_error`,
`statistic`, `p_value`, `p_adjusted`, `conf_low` and `conf_high`. Wald
rows use the fit's covariance (`vcov_type`) with the delta method;
intervals of standard deviations and variances are formed on the log
scale, of correlations on the Fisher scale and of shares on the logit
scale, and only coefficients and logits are tested against zero.
`method = "bootstrap"` resamples persons, refits each resample from
`n_starts` starts and matches its classes to the original by their
posterior memberships of the original persons; it is refused for a
multilevel growth mixture (`latents_unsupported_inference`). In a
simulation of 1000 datasets per design, Wald intervals covered
0.936-0.965 with well-separated classes and 0.885-0.943 with weakly
separated ones, whose standard errors were 6-30% too small
(validation/GROWTH_SIMULATION.md in the source repository); the
bootstrap is preferable for weakly separated classes. A class whose
random-effect covariance is degenerate raises `latents_boundary_fit`
under `boundary = "error"`; under `"fix"` its random-effect rows have no
standard error or interval. `data`, when given, must be the fitted data;
`step` is not used (the information is the derivative of the analytic
score).

For a general transition fit (an
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit with
covariate, occasion-varying or second-order transitions), the table has
one row per transition, initial and second-order logit: `block`,
`group_class`, `from`, `to`, `term`, `estimate`, `standard_error`,
`statistic`, `p_value`, `p_adjusted`, `conf_low`, `conf_high` (and
`previous` for a second-order fit). `method = "bootstrap"` resamples
persons, refits, and matches every replicate's profiles (on the
measurement) and group classes to the original's, rebasing the initial
logits when the reference profile moves; it needs `data`, and works
where Wald inference is refused (for example full covariance
structures). Transition inference fixes no parameter at a bound: without
an active bound `boundary = "fix"` gives the same table as `"error"`,
and a fit with an active bound raises `latents_boundary_fit` under
either.

## Computation

The observed information takes two evaluations of the score over all the
data for each free parameter, which makes it the slowest part of Wald
inference. For profile, class, covariate and transition models it is
computed once per R session for a fit, its data and `step`, and reused
by later calls, including those made by the plots that draw intervals
and by
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md).
Any change to the fit or its data computes it anew, and the result is
the same either way.

## Conditions

Errors: `latents_no_converge` for an unconverged fit;
`latents_boundary_fit` when a variance, response probability or mixing
probability is at its bound; `latents_singular_information` when the
observed information cannot be inverted; `latents_no_free_parameters`
when `fixed` held every parameter; `latents_bad_inference_data` when a
supplied `data` does not reproduce the fit; `latents_too_few_groups` for
`vcov_type = "robust"` or `"opg"` when the fit has no more independent
groups than free parameters. The Wald path raises
`latents_unsupported_inference` for a fit made with `multilpa(prior = )`
(a posterior mode, where the score is not zero), and either path raises
`latents_unsupported_noise` for a fit with a noise component. A fit
whose score is still far from zero gives a `latents_unconverged`
warning.

`method = "bootstrap"` also raises `latents_unsupported_inference` for
an additive or covariate fit, and for a fit holding a measurement block
through `multilpa(fixed = )`, since resampling the data does not
resample the held values (a
[`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
result is bootstrapped by refitting both stages instead), and
`latents_bootstrap_failed` when fewer than two resamples gave a usable
fit; its message gives the first reason. Resamples that fail or do not
converge are dropped with a `latents_bootstrap_dropped` warning stating
how many.

## Examples

``` r
# Two separated profiles, membership driven by `x` and by the group's class.
set.seed(5)
school <- rep(seq_len(16), each = 8)
high_class <- rep(rep(c(FALSE, TRUE), length.out = 16), each = 8)
x <- rnorm(128)
profile <- ifelse(
  runif(128) < plogis(-1 + 2 * high_class + 0.8 * x), 2L, 1L
)
example_data <- data.frame(
  school = school, x = x,
  y1 = rnorm(128, ifelse(profile == 2L, 2, -2), 0.7),
  y2 = rnorm(128, ifelse(profile == 2L, 1.5, -1.5), 0.7)
)
fit <- multilpa(example_data, c("y1", "y2"), "school", n_profiles = 2,
                n_group_classes = 2, profile_covariates = "x",
                n_starts = 2, seed = 1)
parameter_inference(fit, example_data)
#>          level       outcome          term   parameter   estimate
#> 1  measurement     profile_1            y1        mean  1.9190929
#> 2  measurement     profile_1            y2        mean  1.4517555
#> 3  measurement     profile_2            y1        mean -1.8353900
#> 4  measurement     profile_2            y2        mean -1.4838862
#> 5  measurement     profile_1            y1    variance  0.3992639
#> 6  measurement     profile_1            y2    variance  0.5601752
#> 7  measurement     profile_2            y1    variance  0.4605175
#> 8  measurement     profile_2            y2    variance  0.5202968
#> 9      profile     profile_1 group_class_1 coefficient -0.7082123
#> 10     profile     profile_1 group_class_2 coefficient  1.2934666
#> 11     profile     profile_1             x coefficient  0.7401937
#> 12       group group_class_1   (Intercept) coefficient  0.2318526
#>    standard_error  statistic       p_value    p_adjusted   conf_low   conf_high
#> 1      0.07606996  25.227998 1.975161e-140 1.975161e-140  1.7699986  2.06818734
#> 2      0.09010398  16.112002  2.101105e-58  2.101105e-58  1.2751549  1.62835603
#> 3      0.08835082 -20.773888  7.457039e-96  7.457039e-96 -2.0085544 -1.66222558
#> 4      0.09390811 -15.801469  3.039500e-56  3.039500e-56 -1.6679427 -1.29982964
#> 5      0.06798158         NA            NA            NA  0.2859751  0.55743188
#> 6      0.09537546         NA            NA            NA  0.4012345  0.78207700
#> 7      0.08480659         NA            NA            NA  0.3209912  0.66069214
#> 8      0.09579579         NA            NA            NA  0.3626853  0.74640116
#> 9      0.33031825  -2.144030  3.203048e-02  3.203048e-02 -1.3556241 -0.06080037
#> 10     0.41632648   3.106856  1.890883e-03  1.890883e-03  0.4774817  2.10945146
#> 11     0.22710388   3.259274  1.116978e-03  1.116978e-03  0.2950783  1.18530914
#> 12     0.59791891   0.387766  6.981892e-01  6.981892e-01 -0.9400469  1.40375215
set.seed(42)
example_data <- data.frame(
  school = rep(seq_len(10), each = 10),
  score_a = rnorm(100), score_b = rnorm(100)
)
fit <- multilpa(example_data, c("score_a", "score_b"), "school",
                n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
parameter_inference(fit)
#>          level       outcome          term   parameter    estimate
#> 1  measurement     profile_1       score_a        mean  0.16569033
#> 2  measurement     profile_1       score_b        mean -0.07267264
#> 3  measurement     profile_2       score_a        mean -2.35114132
#> 4  measurement     profile_2       score_b        mean -0.35258109
#> 5  measurement     profile_1       score_a    variance  0.78587057
#> 6  measurement     profile_1       score_b    variance  0.71248993
#> 7  measurement     profile_2       score_a    variance  0.22391803
#> 8  measurement     profile_2       score_b    variance  2.46889204
#> 9      profile     profile_1 group_class_1 probability  0.94708605
#> 10     profile     profile_2 group_class_1 probability  0.05291395
#> 11       group group_class_1          <NA> probability  1.00000000
#>    standard_error  statistic      p_value   p_adjusted    conf_low  conf_high
#> 1      0.10691280  1.5497708 1.211965e-01 1.211965e-01 -0.04385490  0.3752356
#> 2      0.09184042 -0.7912926 4.287733e-01 4.287733e-01 -0.25267656  0.1073313
#> 3      0.38614364 -6.0887740 1.137786e-09 1.137786e-09 -3.10796894 -1.5943137
#> 4      0.79370695 -0.4442207 6.568830e-01 6.568830e-01 -1.90821812  1.2030559
#> 5      0.14710614         NA           NA           NA  0.54452229  1.1341915
#> 6      0.10799758         NA           NA           NA  0.52936541  0.9589631
#> 7      0.22273970         NA           NA           NA  0.03186888  1.5732992
#> 8      1.68065393         NA           NA           NA  0.65022241  9.3743738
#> 9      0.03538773         NA           NA           NA  0.81767987  0.9861939
#> 10     0.03538773         NA           NA           NA  0.01380615  0.1823201
#> 11     0.00000000         NA           NA           NA  1.00000000  1.0000000

# Many tests in one table: name the correction, do not apply one by stealth.
parameter_inference(fit, adjust = "BH")
#>          level       outcome          term   parameter    estimate
#> 1  measurement     profile_1       score_a        mean  0.16569033
#> 2  measurement     profile_1       score_b        mean -0.07267264
#> 3  measurement     profile_2       score_a        mean -2.35114132
#> 4  measurement     profile_2       score_b        mean -0.35258109
#> 5  measurement     profile_1       score_a    variance  0.78587057
#> 6  measurement     profile_1       score_b    variance  0.71248993
#> 7  measurement     profile_2       score_a    variance  0.22391803
#> 8  measurement     profile_2       score_b    variance  2.46889204
#> 9      profile     profile_1 group_class_1 probability  0.94708605
#> 10     profile     profile_2 group_class_1 probability  0.05291395
#> 11       group group_class_1          <NA> probability  1.00000000
#>    standard_error  statistic      p_value   p_adjusted    conf_low  conf_high
#> 1      0.10691280  1.5497708 1.211965e-01 2.423931e-01 -0.04385490  0.3752356
#> 2      0.09184042 -0.7912926 4.287733e-01 5.716977e-01 -0.25267656  0.1073313
#> 3      0.38614364 -6.0887740 1.137786e-09 4.551145e-09 -3.10796894 -1.5943137
#> 4      0.79370695 -0.4442207 6.568830e-01 6.568830e-01 -1.90821812  1.2030559
#> 5      0.14710614         NA           NA           NA  0.54452229  1.1341915
#> 6      0.10799758         NA           NA           NA  0.52936541  0.9589631
#> 7      0.22273970         NA           NA           NA  0.03186888  1.5732992
#> 8      1.68065393         NA           NA           NA  0.65022241  9.3743738
#> 9      0.03538773         NA           NA           NA  0.81767987  0.9861939
#> 10     0.03538773         NA           NA           NA  0.01380615  0.1823201
#> 11     0.00000000         NA           NA           NA  1.00000000  1.0000000

# A structure that constrains the shape across profiles (VEI) is charted in
# its own free coordinates, so it has Wald standard errors too.
shaped <- multilpa(example_data, c("score_a", "score_b"), "school",
                   n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1,
                   volume = "varying", shape = "equal", orientation = "axis")
parameter_inference(shaped)
#>          level       outcome          term   parameter    estimate
#> 1  measurement     profile_1       score_a        mean  0.25022879
#> 2  measurement     profile_1       score_b        mean -0.06203038
#> 3  measurement     profile_2       score_a        mean -1.35411736
#> 4  measurement     profile_2       score_b        mean -0.24959732
#> 5  measurement     profile_1       score_a    variance  0.70077548
#> 6  measurement     profile_1       score_b    variance  0.70904425
#> 7  measurement     profile_2       score_a    variance  1.31232456
#> 8  measurement     profile_2       score_b    variance  1.32780927
#> 9      profile     profile_1 group_class_1 probability  0.86429738
#> 10     profile     profile_2 group_class_1 probability  0.13570262
#> 11       group group_class_1          <NA> probability  1.00000000
#>    standard_error  statistic    p_value p_adjusted    conf_low  conf_high
#> 1      0.13357063  1.8733818 0.06101568 0.06101568 -0.01156483 0.51202242
#> 2      0.09872161 -0.6283364 0.52978357 0.52978357 -0.25552118 0.13146041
#> 3      0.72483385 -1.8681762 0.06173751 0.06173751 -2.77476561 0.06653089
#> 4      0.41275474 -0.6047110 0.54537106 0.54537106 -1.05858176 0.55938711
#> 5      0.13997195         NA         NA         NA  0.47376290 1.03656549
#> 6      0.11795100         NA         NA         NA  0.51176924 0.98236413
#> 7      0.61223539         NA         NA         NA  0.52593539 3.27453859
#> 8      0.48454076         NA         NA         NA  0.64940766 2.71490090
#> 9      0.10218086         NA         NA         NA  0.53592034 0.97231994
#> 10     0.10218086         NA         NA         NA  0.02768006 0.46407966
#> 11     0.00000000         NA         NA         NA  1.00000000 1.00000000

# The bootstrap resamples schools and refits inside the same family.
parameter_inference(shaped, method = "bootstrap", iter = 10, n_starts = 1,
                    seed = 1)
#> Warning: Profile 2: 3.0 Gaussian membership mass for 3 free means and variances. This is a conservative diagnostic, not an identification test; inspect the Gaussian covariance diagnostics.
#>          level       outcome          term   parameter    estimate
#> 1  measurement     profile_1       score_a        mean  0.25022879
#> 2  measurement     profile_1       score_b        mean -0.06203038
#> 3  measurement     profile_2       score_a        mean -1.35411736
#> 4  measurement     profile_2       score_b        mean -0.24959732
#> 5  measurement     profile_1       score_a    variance  0.70077548
#> 6  measurement     profile_1       score_b    variance  0.70904425
#> 7  measurement     profile_2       score_a    variance  1.31232456
#> 8  measurement     profile_2       score_b    variance  1.32780927
#> 9      profile     profile_1 group_class_1 probability  0.86429738
#> 10     profile     profile_2 group_class_1 probability  0.13570262
#> 11       group group_class_1          <NA> probability  1.00000000
#>    standard_error statistic p_value p_adjusted    conf_low conf_high
#> 1       0.3673629        NA      NA         NA  0.03448968 1.0490512
#> 2       0.4769439        NA      NA         NA -0.30194200 1.1224914
#> 3       1.0307447        NA      NA         NA -2.68649048 0.1362221
#> 4       0.5430616        NA      NA         NA -0.69622362 0.9697939
#> 5       0.2765594        NA      NA         NA  0.20977054 1.0503136
#> 6       0.2198624        NA      NA         NA  0.15771978 0.7749002
#> 7       0.5999388        NA      NA         NA  0.03880418 1.7155546
#> 8       0.5841325        NA      NA         NA  0.04243513 1.6601809
#> 9       0.2848561        NA      NA         NA  0.16211188 0.9682255
#> 10      0.2848561        NA      NA         NA  0.03177455 0.8378881
#> 11      0.0000000        NA      NA         NA  1.00000000 1.0000000
```
