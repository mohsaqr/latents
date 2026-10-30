# Parametric bootstrap likelihood-ratio test for nested mixture models

The parametric bootstrap likelihood-ratio test evaluates whether a more
flexible mixture model provides an improvement in fit that is
inconsistent with a simpler, nested model. `bootstrap_lrt()` compares
two fitted models using twice the difference between their maximized log
likelihoods. In class enumeration, this comparison assesses the evidence
for adding a latent profile or group class.

The function generates datasets from the fitted null model and refits
both candidate models to each dataset. These simulated comparisons
approximate the distribution of the likelihood-ratio statistic under the
null hypothesis. The observed statistic is then evaluated against this
distribution to obtain a bootstrap p-value. This approach avoids relying
on the usual chi-squared approximation, which is generally inappropriate
when comparing different numbers of mixture components.

The result includes the observed test statistic, bootstrap p-value,
Monte Carlo standard error, and diagnostics for the simulated refits.
The function supports profile and class models, latent transition
models, and group-class model families, subject to the nesting and
compatibility requirements below.

## Usage

``` r
bootstrap_lrt(
  null_model,
  alternative_model,
  data = NULL,
  iter = 199L,
  n_starts = NULL,
  max_iter = NULL,
  tol = NULL,
  seed = NULL
)
```

## Arguments

- null_model:

  A converged fitted model representing the null hypothesis. Supported
  inputs are profile or class models, with or without membership
  covariates; transition models from
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md); and
  group-class models fitted through
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).
  See *Supported comparisons* for nesting requirements.

- alternative_model:

  A compatible, more flexible model fitted to the same data, within
  which `null_model` is nested.

- data:

  The data frame used to fit both models. Required for transition
  models, whose simulations retain the observed groups, occasions and
  covariates. For ordinary profile and class models, omission
  reconstructs the data from the null fit, or from the alternative for
  membership-covariate models. For group-class models, the stored
  sufficient statistics are used when `data` is omitted. Supplied data
  are checked against the fitted models.

- iter:

  Number of simulated datasets; default 199. Profile, class and
  group-class comparisons require at least two replicates; transition
  comparisons accept any positive integer. Larger values reduce Monte
  Carlo uncertainty in the estimated p-value.

- n_starts:

  Number of starts for each simulated fit. With `NULL` (the default),
  each model retains its original number of starts and search strategy.
  A supplied number applies to both models.

- max_iter:

  Maximum iterations for each simulated fit. With `NULL` (the default),
  each model retains its original iteration limit.

- tol:

  Convergence tolerance for each simulated fit. With `NULL` (the
  default), each model retains its original tolerance.

- seed:

  Optional integer seed accepted by
  [`set.seed()`](https://rdrr.io/r/base/Random.html). The caller's
  random-number state is restored on exit.

## Value

An object of class `multilpa_bootstrap_lrt` containing the observed
statistic, bootstrap p-value, Monte Carlo standard error, numbers of
valid and invalid replicates, and replicate-level diagnostics. The
`fixed` component identifies any measurement blocks held constant.

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns a
one-row test summary. `get_results(result, what = "replicates")` returns
one row per simulated replicate.
[`summary()`](https://rdrr.io/r/base/summary.html) combines the test
summary and replicate diagnostics;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) displays the
simulated reference distribution and observed statistic.

## Interpretation

A small p-value provides evidence against the null model in favour of
the specified alternative, subject to the model assumptions and the
quality of estimation. It does not establish that the alternative has
the correct number of classes or fits the data adequately in an absolute
sense. Interpretation should consider convergence, replication across
starting values, and the frequency of invalid bootstrap replicates.
Increasing the number of replicates reduces simulation uncertainty but
does not resolve systematic estimation problems.

## Supported comparisons

For ordinary profile and class models, with or without membership
covariates, the models must differ by exactly one observation-level
profile or class, or by exactly one group class. They must use the same
observations, group layout, covariance structure, centering and scaling,
missing-data handling, and `drop_empty` setting. Grand-mean centering is
repeated in each simulated refit. Person-centred or person-scaled fits
are not supported, because the model does not specify a distribution for
the group baselines removed by these transformations.

Transition models fitted with
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) support
comparisons in which the null specification is nested within the
alternative. Supported differences include class counts and restrictions
on measurement, transition probabilities, Markov order, and covariate
effects. Both models must use compatible indicators, grouping, occasion
layout, and measurement distributions. The original `data` must be
supplied for transition comparisons.

Group-class models fitted through the `family` argument of
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
support comparisons across class counts and parameter restrictions. The
null must be nested within, and distinct from, the alternative: the
alternative must have at least as many group classes and cannot impose
additional restrictions on the class-specific means or variance
components. Models from different families cannot be combined
arbitrarily.

For membership-covariate models, covariates are held at their observed
values when simulating data. Group classes and profiles are drawn using
the fitted membership probabilities. Both models must use the same
profile covariates, group covariates and slope specification, except
that a one-group-class null may be compared with an alternative that
adds group covariates or slopes that vary by group class.

For models fitted with `fixed`, including those from
[`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md),
each replicate retains the same fixed blocks at the same values. The
resulting p-value is conditional on that measurement solution;
uncertainty in its estimation is not propagated. Pairs for which nesting
cannot be established are rejected.

## Bootstrap calculation

The observed statistic is
`2 * (logLik(alternative_model) - logLik(null_model))`. Datasets are
simulated under the fitted null model, retaining the group structure and
other design features required by the model family. Both models are
refitted to each simulated dataset and the statistic is recalculated. By
default, each model retains its original search controls; `n_starts`,
`max_iter` and `tol` can override the corresponding settings for both
models.

The bootstrap p-value is `(1 + b) / (V + 1)`, where `b` is the number of
valid simulated statistics at least as large as the observed statistic
and `V` is the number of valid replicates. When all replicates are
valid, `V = iter`. The reported Monte Carlo standard error summarizes
uncertainty from the finite number of simulated replicates.

## Estimation failures

Both original models must have converged. A simulated replicate is
invalid if a refit fails, does not converge, or produces an
alternative-model likelihood below the null-model likelihood beyond
numerical tolerance. The result reports `n_valid` and `n_invalid`;
`latents_failed_replicates` warns when any replicate is invalid. The
p-value uses only valid replicates and is `NA` when fewer than half of
the requested replicates are valid.

Excluding invalid replicates conditions the reference distribution on
successful estimation and may affect calibration. Report the invalid
proportion and examine the replicate diagnostics before interpreting the
test.

For membership-covariate models, a simulated refit is considered
converged for this test when its likelihood has converged, even if a
membership logit is still changing. This can occur with an empty or
separated class in an alternative model. The `logits_settled` column
records whether the logits also converged. The original covariate models
must be fully converged.

## Missing-data assumptions

Fits estimated with `missing = "fiml"` are supported when both models
use that setting. Each simulated dataset retains the original
missing-value pattern before refitting. Calibration requires missingness
independent of the indicators and latent classes, conditional on any
fixed covariates. A fixed missingness pattern does not reproduce a
missing-at-random (MAR) mechanism that depends on observed indicators or
a missing-not-at-random (MNAR) mechanism. Valid FIML estimation under
MAR does not by itself establish validity of this bootstrap under MAR.

## Conditions

For ordinary profile and class comparisons, `latents_bad_nesting`
indicates incompatible class counts, fixed measurement blocks, covariate
specifications, or parameter counts. `latents_incomparable_models`
indicates incompatible data or fitting specifications. Transition and
group-class comparisons use their family-specific nesting checks and
raise `latents_bad_argument` for unsupported nesting requests.
Transition models also raise `latents_bad_argument` when `data` is
omitted.

`latents_unsupported_bootstrap` identifies unsupported person centering
or scaling. Models with a prior, sampling weights or a noise component
raise `latents_unsupported_prior`, `latents_unsupported_weights` or
`latents_unsupported_noise`, respectively. `latents_bad_data` identifies
unsupported missing indicators, and `latents_bad_inference_data`
identifies supplied data that do not reproduce the original fitting data
or likelihood.

`latents_no_converge` indicates that an original model has not
converged; for ordinary profile and class models it also identifies
active variance bounds. `latents_reversed_likelihood` indicates that the
alternative has a lower likelihood than the null beyond numerical
tolerance. Invalid seeds or incompatible model families raise
`latents_bad_argument`. `latents_failed_replicates` warns about invalid
simulated replicates; see *Estimation failures* for their treatment.

## References

Nylund, K. L., Asparouhov, T., & Muthen, B. O. (2007). Deciding on the
number of classes in latent class analysis and growth mixture modeling:
A Monte Carlo simulation study. *Structural Equation Modeling*, 14(4),
535–569.

## See also

[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
for fitting candidate models;
[`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)
for information-criterion comparisons.

## Examples

``` r
set.seed(1)
example_data <- data.frame(
  school = rep(seq_len(10), each = 10),
  score = rnorm(100, rep(c(-2, 2), each = 50))
)
smaller <- multilpa(example_data, "score", "school", n_profiles = 1,
                    n_group_classes = 1, n_starts = 2, seed = 1)
larger <- multilpa(example_data, "score", "school", n_profiles = 2,
                   n_group_classes = 1, n_starts = 2, seed = 1)
# `iter` is small so the example runs quickly; use many more for inference.
comparison <- bootstrap_lrt(smaller, larger, iter = 9, n_starts = 1,
                            max_iter = 2000, tol = 1e-6, seed = 1)
comparison
#> Parametric bootstrap likelihood-ratio comparison
#> Null: 1 profile; alternative: 2 profiles
#> Observed statistic: 44.198110
#> p-value: 0.1 (Monte Carlo SE 0.0949) from 9 of 9 valid replicates (0 invalid)
```
