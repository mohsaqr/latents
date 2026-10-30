# Standard errors, confidence intervals and tests for model parameters

`parameter_inference()` quantifies sampling uncertainty in the parameter
estimates of a fitted mixture model. It provides standard errors and
confidence intervals for quantities such as indicator means, response
probabilities, membership coefficients, and transition probabilities,
according to the model family. These results describe the precision of
the estimates and support interpretation of the fitted model.

Two approaches are available for supported models. Wald inference uses
an estimated parameter covariance matrix to obtain standard errors,
confidence intervals, and tests for eligible parameters. Bootstrap
inference repeatedly resamples independent observational units and
refits the model, using variation across the resulting estimates to
calculate standard errors and percentile confidence intervals.

The function returns a table of parameter estimates and their
uncertainty measures. The available parameters, covariance estimators,
and inference methods depend on the fitted model. Wald tests are omitted
where an ordinary test against zero is inappropriate, including
probabilities and variance parameters. Bootstrap results contain
standard errors and confidence intervals but no hypothesis-test
p-values.

## Usage

``` r
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
```

## Arguments

- x:

  A converged fit from
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
  [`lca()`](https://pak.dynasite.org/latents/reference/lca.md) or
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
  (with or without covariates),
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md), the
  additive family, or a growth mixture. A fit with an active bound
  raises `latents_boundary_fit` and an unconverged fit
  `latents_no_converge`. Mixture regressions raise
  `latents_bad_argument` (their inference is in
  [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md),
  [`confint()`](https://rdrr.io/r/stats/confint.html) and
  [`vcov()`](https://rdrr.io/r/stats/vcov.html)), and the cross-level
  families raise `latents_unsupported_inference`.

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

  Covariance estimator for the Wald method. `"observed"` (default)
  inverts the observed information. `"robust"` returns the Huber-White
  sandwich `A^-1 B A^-1`, where `B` is the outer product of the
  per-group scores, multiplied by the CR1 small-sample factor
  `G / (G - 1)` for `G` groups (the default of `sandwich::vcovCL()`,
  applied to every robust standard error in the package). Mplus's MLR
  applies no such factor, so these errors are `sqrt(G / (G - 1))` times
  Mplus's. Groups are the independent units: robust errors relax the
  assumption that the within-group model is correctly specified, not the
  assumption that groups are independent. `"opg"` inverts the outer
  product of the per-group scores (the BHHH estimator). The OPG and
  observed-information estimators are asymptotically equivalent under a
  correctly specified regular model, but can differ in finite samples.
  `"robust"` and `"opg"` need more groups than parameters; `"opg"` also
  needs group scores of full rank, since it inverts their outer product,
  while the sandwich does not. For a fit with sampling weights the
  default is `"robust"`, and any other value raises
  `latents_unsupported_weights`. For a fit with `centering = "person"`
  or `scaling = "person"` the default is also `"robust"`, and
  `"observed"` or `"opg"` raise `latents_unsupported_inference`.
  Transforming observations within a unit induces dependence; the
  sandwich estimator aggregates score contributions at the
  independent-unit level. Under `centering = "grand"` or
  `scaling = "grand"` the sample means and standard deviations the
  indicators were transformed with are estimated jointly with the model,
  as estimating equations, and their variability is added to every
  covariance type: a centred mean `mu_k - xbar` has the standard error
  of that difference, not that of `mu_k`. `vcov_type` and `step` are not
  used by `method = "bootstrap"`.

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
  these coordinates by the delta method. Constraints can make the
  natural-scale covariance singular; the sum constraints on class
  probabilities also imply a singular covariance. `"bootstrap"`
  resamples groups with replacement, refits within the same covariance
  structure, realigns each refit's labels, and reports the percentile
  interval and the standard deviation of the replicates; it is also
  available for all fourteen structures and assumes independent groups,
  as `vcov_type = "robust"` does.

- iter:

  Number of resamples for `method = "bootstrap"`, at least two. Default
  `199`; a 95% percentile interval is more stable with 999 or more.

- n_starts, max_iter, tol:

  Passed to each bootstrap refit, as in
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).

- seed:

  Optional seed for the resampling: a whole number accepted by
  [`set.seed()`](https://rdrr.io/r/base/Random.html); any other value
  raises `latents_bad_argument` with either method. The Wald method uses
  no random numbers. A supplied integer seed restores the caller's
  `.Random.seed` on exit; `NULL` uses its current stream. R does not
  expose the Box-Muller cached normal variate, so that cache cannot be
  restored; see [base::Random](https://rdrr.io/r/base/Random.html).

- boundary:

  Handling of categorical response probabilities at their lower bound,
  `min_probability`, where the information is singular. `"error"`
  (default) raises `latents_boundary_fit`. `"fix"` holds those
  probabilities at the bound and reports every other parameter
  conditional on them, as Mplus does; the held probabilities are omitted
  from the table (or, for a covariate fit, reported with `NA` standard
  errors) and named in the `fixed_at_bound` attribute.

## Value

A base `data.frame` with one row per reported parameter and the columns
`level` (`"measurement"`, `"profile"` or `"group"`), `outcome`, `term`,
`parameter`, `estimate`, `standard_error`, `statistic`, `p_value`,
`p_adjusted`, `conf_low` and `conf_high`. Estimates are on the natural
scale (variances in squared indicator units, probabilities on `[0, 1]`).
A parameter whose null value lies on the boundary of its space (a
variance, a probability, a covariance diagonal, a count mean) has `NA`
`statistic`, `p_value` and `p_adjusted`, but still has an interval. All
class probabilities are reported; their sum constraints make the natural
covariance singular. Wald intervals for probabilities are formed on the
logit scale and for variances (and covariance diagonals) and count means
on the log scale, then back-transformed, so they stay within the
parameter's range; every other interval is
`estimate +/- z * standard_error`.
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

For the Wald method, diagnostics of the whole fit are attributes:
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

For a covariate fit (`multilpa(profile_covariates = )`), the table
covers every free parameter at all three levels: the measurement model,
the profile logits and the group-class logits. `outcome` is the profile
or group class a coefficient predicts, or the profile a measurement
parameter belongs to, and `parameter` is `"mean"`, `"variance"`,
`"covariance"`, `"response"` or `"coefficient"`. Each categorical
indicator contributes one `"response"` row per profile, indicator and
category, named `indicator:category`: the category probability, with its
standard error from the multinomial logits and no Wald test.
Class-membership coefficients are on the multinomial logit scale, with
the final profile and the final group class as references, so each
coefficient is a log odds against that reference. The tests are Wald
tests, which are invalid for a coefficient driven to the boundary by
separation or an empty class. A fit with a membership logit above 20 in
coordinates independent of the covariates' units and origin (roughly, a
slope per standard deviation of its covariate and an intercept at the
covariates' means) raises `latents_boundary_fit`; the fit itself warns
with `latents_extreme_coefficients`. Only `method = "wald"` is
implemented for a covariate fit; `"bootstrap"` raises
`latents_unsupported_inference`.

For a latent transition model from
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md), the table
covers the measurement model, the group-class shares, each group class's
initial profile distribution and its transition matrix. An initial
probability (`parameter = "initial_probability"`, `level = "profile"`)
has `outcome` the profile and `term` the group class; a transition
probability (`"transition_probability"`) has `outcome` the profile moved
to and `term` `group_class_h:profile_j`, the group class and the profile
moved from. The scores follow the Fisher identity and the observed
information is the numerical Jacobian of that analytic score;
`vcov_type = "robust"` clusters on the groups, which are the sequences.
`method = "bootstrap"` resamples persons, refits and matches each
replicate's profiles and group classes to the original's; it needs
`data` (`latents_bad_argument` otherwise). A transition row with no
information (a profile no group occupies before its last occasion)
raises `latents_boundary_fit`.

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
posterior memberships of the original persons; it is unavailable for a
multilevel growth mixture (`latents_unsupported_inference`). Weak class
separation can reduce the accuracy of Wald intervals; bootstrap
inference provides an alternative assessment of sampling uncertainty. A
class whose random-effect covariance is degenerate, or (with correlated
residuals) a residual standard deviation or residual process held on its
bound, raises `latents_boundary_fit` under `boundary = "error"`, naming
which; under `"fix"` those rows have no standard error or interval.
`data`, when given, must be the fitted data; `step` is not used (the
information is the derivative of the analytic score).

For a general transition fit (an
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit with
covariate, occasion-varying or second-order transitions), the table has
one row per transition, initial and second-order logit: `block`,
`group_class`, `from`, `to`, `term`, `estimate`, `standard_error`,
`statistic`, `p_value`, `p_adjusted`, `conf_low`, `conf_high`,
`estimated` (and `previous` for a second-order fit). `estimated` is
`FALSE` for a logit the likelihood does not identify: a move into an
occasion at which no row observes an indicator, under
`transitions = "occasion"` (or every second-order logit when nothing is
observed from occasion 3 on). Such a logit is not counted in
`n_parameters`, AIC or BIC, and its standard error and test are `NA`.
`method = "bootstrap"` resamples persons, refits, and matches every
replicate's profiles (on the measurement) and group classes to the
original's, rebasing the initial logits when the reference profile
moves; it needs `data`, and works where Wald inference is refused (for
example full covariance structures).

Transition inference (both kinds of
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit) fixes
no parameter at a bound: without an active bound `boundary = "fix"`
gives the same table as `"error"`, and a fit with an active bound raises
`latents_boundary_fit` under either.

## Interpretation

Inference is conditional on the specified model, including its number of
classes and measurement assumptions. The reported uncertainty does not
account for selecting that specification from several candidates. Wald
inference requires an adequately converged and locally identified
solution; estimates near parameter boundaries may require special
treatment. Bootstrap inference also depends on successful refitting and
consistent matching of class labels across resamples.

## Computation

Numerical differentiation of the score requires two score evaluations
per free parameter. For profile, class, covariate and transition models,
the information matrix is cached within an R session for the same fitted
model, data and `step`. Subsequent inference, plotting and
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
calls reuse this calculation; changes to its inputs trigger
recalculation.

## Conditions

Errors: `latents_no_converge` for an unconverged fit;
`latents_boundary_fit` when a variance, response probability or mixing
probability is at its bound; `latents_singular_information` when the
observed information cannot be inverted; `latents_no_free_parameters`
when `fixed` held every parameter; `latents_bad_inference_data` when a
supplied `data` does not reproduce the fit; `latents_too_few_groups` for
`vcov_type = "robust"` or `"opg"` when the fit has no more independent
groups than free parameters. The Wald method raises
`latents_unsupported_inference` for a fit made with `multilpa(prior = )`
(a posterior mode, where the score is not zero), and either method
raises `latents_unsupported_noise` for a fit with a noise component. A
fit whose score is still far from zero gives a `latents_unconverged`
warning. An invalid `level`, `step`, `iter`, `n_starts`, `max_iter` or
`tol` raises `latents_bad_argument`, and a `data` that is not a data
frame `latents_bad_data`. A fitted mixture regression raises
`latents_bad_argument` naming `get_results(x, "coefficients")`, which
reports its standard errors and tests; an object that is not a fitted
model raises `latents_not_a_fit`.

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
set.seed(42)
example_data <- data.frame(
  school = rep(seq_len(10), each = 10),
  score_a = rnorm(100), score_b = rnorm(100)
)
fit <- multilpa(example_data, c("score_a", "score_b"), "school",
                n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
parameter_inference(fit)
#>          level       outcome          term   parameter    estimate
#> 1  measurement     profile_1       score_a        mean  0.16571533
#> 2  measurement     profile_1       score_b        mean -0.07265950
#> 3  measurement     profile_2       score_a        mean -2.35097507
#> 4  measurement     profile_2       score_b        mean -0.35274803
#> 5  measurement     profile_1       score_a    variance  0.78583377
#> 6  measurement     profile_1       score_b    variance  0.71247734
#> 7  measurement     profile_2       score_a    variance  0.22400619
#> 8  measurement     profile_2       score_b    variance  2.46857665
#> 9      profile     profile_1 group_class_1 probability  0.94707314
#> 10     profile     profile_2 group_class_1 probability  0.05292686
#> 11       group group_class_1          <NA> probability  1.00000000
#>    standard_error  statistic      p_value   p_adjusted    conf_low  conf_high
#> 1      0.10689493  1.5502636 1.210782e-01 1.210782e-01 -0.04379489  0.3752255
#> 2      0.09183129 -0.7912281 4.288109e-01 4.288109e-01 -0.25264553  0.1073265
#> 3      0.38600600 -6.0905143 1.125485e-09 1.125485e-09 -3.10753293 -1.5944172
#> 4      0.79331210 -0.4446523 6.565711e-01 6.565711e-01 -1.90761117  1.2021151
#> 5      0.14706811         NA           NA           NA  0.54453908  1.1340503
#> 6      0.10798809         NA           NA           NA  0.52936710  0.9589262
#> 7      0.22275459         NA           NA           NA  0.03190175  1.5729162
#> 8      1.68018293         NA           NA           NA  0.65027166  9.3712690
#> 9      0.03537887         NA           NA           NA  0.81774059  0.9861813
#> 10     0.03537887         NA           NA           NA  0.01381871  0.1822594
#> 11     0.00000000         NA           NA           NA  1.00000000  1.0000000

# Apply Benjamini-Hochberg adjustment to the reported tests.
parameter_inference(fit, adjust = "BH")
#>          level       outcome          term   parameter    estimate
#> 1  measurement     profile_1       score_a        mean  0.16571533
#> 2  measurement     profile_1       score_b        mean -0.07265950
#> 3  measurement     profile_2       score_a        mean -2.35097507
#> 4  measurement     profile_2       score_b        mean -0.35274803
#> 5  measurement     profile_1       score_a    variance  0.78583377
#> 6  measurement     profile_1       score_b    variance  0.71247734
#> 7  measurement     profile_2       score_a    variance  0.22400619
#> 8  measurement     profile_2       score_b    variance  2.46857665
#> 9      profile     profile_1 group_class_1 probability  0.94707314
#> 10     profile     profile_2 group_class_1 probability  0.05292686
#> 11       group group_class_1          <NA> probability  1.00000000
#>    standard_error  statistic      p_value   p_adjusted    conf_low  conf_high
#> 1      0.10689493  1.5502636 1.210782e-01 2.421565e-01 -0.04379489  0.3752255
#> 2      0.09183129 -0.7912281 4.288109e-01 5.717478e-01 -0.25264553  0.1073265
#> 3      0.38600600 -6.0905143 1.125485e-09 4.501941e-09 -3.10753293 -1.5944172
#> 4      0.79331210 -0.4446523 6.565711e-01 6.565711e-01 -1.90761117  1.2021151
#> 5      0.14706811         NA           NA           NA  0.54453908  1.1340503
#> 6      0.10798809         NA           NA           NA  0.52936710  0.9589262
#> 7      0.22275459         NA           NA           NA  0.03190175  1.5729162
#> 8      1.68018293         NA           NA           NA  0.65027166  9.3712690
#> 9      0.03537887         NA           NA           NA  0.81774059  0.9861813
#> 10     0.03537887         NA           NA           NA  0.01381871  0.1822594
#> 11     0.00000000         NA           NA           NA  1.00000000  1.0000000

# Wald inference for a VEI model with equal shape across profiles.
shaped <- multilpa(example_data, c("score_a", "score_b"), "school",
                   n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1,
                   volume = "varying", shape = "equal", orientation = "axis")
parameter_inference(shaped)
#>          level       outcome          term   parameter    estimate
#> 1  measurement     profile_1       score_a        mean -1.35412149
#> 2  measurement     profile_1       score_b        mean -0.24961049
#> 3  measurement     profile_2       score_a        mean  0.25022652
#> 4  measurement     profile_2       score_b        mean -0.06202866
#> 5  measurement     profile_1       score_a    variance  1.31232687
#> 6  measurement     profile_1       score_b    variance  1.32780761
#> 7  measurement     profile_2       score_a    variance  0.70077856
#> 8  measurement     profile_2       score_b    variance  0.70904523
#> 9      profile     profile_1 group_class_1 probability  0.13570105
#> 10     profile     profile_2 group_class_1 probability  0.86429895
#> 11       group group_class_1          <NA> probability  1.00000000
#>    standard_error  statistic    p_value p_adjusted    conf_low conf_high
#> 1       0.7248330 -1.8681841 0.06173642 0.06173642 -2.77476808 0.0665251
#> 2       0.4127571 -0.6047394 0.54535214 0.54535214 -1.05859949 0.5593785
#> 3       0.1335693  1.8733842 0.06101535 0.06101535 -0.01156440 0.5120174
#> 4       0.0987215 -0.6283197 0.52979456 0.52979456 -0.25551923 0.1314619
#> 5       0.6122361         NA         NA         NA  0.52593663 3.2745424
#> 6       0.4845408         NA         NA         NA  0.64940624 2.7149001
#> 7       0.1399720         NA         NA         NA  0.47376579 1.0365683
#> 8       0.1179509         NA         NA         NA  0.51177031 0.9823648
#> 9       0.1021789         NA         NA         NA  0.02768012 0.4640724
#> 10      0.1021789         NA         NA         NA  0.53592760 0.9723199
#> 11      0.0000000         NA         NA         NA  1.00000000 1.0000000

# Resample schools and refit the same covariance structure.
parameter_inference(shaped, method = "bootstrap", iter = 10, n_starts = 1,
                    seed = 1)
#> Warning: Profile 2: 3.0 Gaussian membership mass for 3 free means and variances. This is a conservative diagnostic, not an identification test; inspect the Gaussian covariance diagnostics.
#>          level       outcome          term   parameter    estimate
#> 1  measurement     profile_1       score_a        mean -1.35412149
#> 2  measurement     profile_1       score_b        mean -0.24961049
#> 3  measurement     profile_2       score_a        mean  0.25022652
#> 4  measurement     profile_2       score_b        mean -0.06202866
#> 5  measurement     profile_1       score_a    variance  1.31232687
#> 6  measurement     profile_1       score_b    variance  1.32780761
#> 7  measurement     profile_2       score_a    variance  0.70077856
#> 8  measurement     profile_2       score_b    variance  0.70904523
#> 9      profile     profile_1 group_class_1 probability  0.13570105
#> 10     profile     profile_2 group_class_1 probability  0.86429895
#> 11       group group_class_1          <NA> probability  1.00000000
#>    standard_error statistic p_value p_adjusted    conf_low conf_high
#> 1       1.0307711        NA      NA         NA -2.68649048 0.1362212
#> 2       0.5431321        NA      NA         NA -0.69622746 0.9699979
#> 3       0.3674272        NA      NA         NA  0.03449099 1.0492880
#> 4       0.4769342        NA      NA         NA -0.30187374 1.1224894
#> 5       0.5999731        NA      NA         NA  0.03880418 1.7156408
#> 6       0.5841342        NA      NA         NA  0.04243513 1.6602573
#> 7       0.2765998        NA      NA         NA  0.20977145 1.0503225
#> 8       0.2198904        NA      NA         NA  0.15764115 0.7748999
#> 9       0.2848955        NA      NA         NA  0.03177455 0.8379474
#> 10      0.2848955        NA      NA         NA  0.16205256 0.9682255
#> 11      0.0000000        NA      NA         NA  1.00000000 1.0000000

# A covariate fit: membership driven by `x` and by the group's class.
set.seed(5)
school <- rep(seq_len(16), each = 8)
high_class <- rep(rep(c(FALSE, TRUE), length.out = 16), each = 8)
x <- rnorm(128)
profile <- ifelse(
  runif(128) < plogis(-1 + 2 * high_class + 0.8 * x), 2L, 1L
)
covariate_data <- data.frame(
  school = school, x = x,
  y1 = rnorm(128, ifelse(profile == 2L, 2, -2), 0.7),
  y2 = rnorm(128, ifelse(profile == 2L, 1.5, -1.5), 0.7)
)
covariate_fit <- multilpa(covariate_data, c("y1", "y2"), "school",
                          n_profiles = 2, n_group_classes = 2,
                          profile_covariates = "x", n_starts = 2, seed = 1)
parameter_inference(covariate_fit, covariate_data)
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

# Inference for a latent transition model, including transition probabilities.
moves <- lta(subset(course_engagement, student <= 40),
             c("browse", "lectures", "forum_read"), "student",
             n_profiles = 2, time = "sequence", n_starts = 2, seed = 1)
parameter_inference(moves)
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
#> 4  -0.8706368     0.05194461 -16.760870  4.716736e-63  4.716736e-63 -0.9724464
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
