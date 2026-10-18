# Parametric bootstrap likelihood-ratio test of two nested models

Tests a model against one with one more individual profile or one more
group class. Indicators are simulated from the null model with the
observed group sizes, both models are refitted to every simulated
dataset, and the observed likelihood-ratio statistic is compared with
the simulated ones. The p-value is `(1 + b) / (iter + 1)`, where `b`
counts simulated statistics at least as large as the observed one; no
chi-square reference distribution is used. This is a native parametric
bootstrap, not an implementation of Mplus TECH14.

## Usage

``` r
bootstrap_lrt(
  null_model,
  alternative_model,
  data = NULL,
  iter = 199L,
  n_starts = 10L,
  max_iter = 1000L,
  tol = 1e-08,
  seed = NULL
)
```

## Arguments

- null_model:

  The smaller model: a converged
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  fit, with or without membership covariates.

- alternative_model:

  The larger model, fitted to the same data.

- data:

  Optional. The data both models were fitted to, used to verify both
  fitted likelihoods. When omitted it is rebuilt from what the null
  model stores (from the alternative for covariate fits).

- iter:

  Number of simulated datasets, at least 2. Default 199; use many for
  inference.

- n_starts:

  Number of random starts for each simulated fit. Default 10.

- max_iter:

  Maximum EM iterations for each simulated fit. Default 1000.

- tol:

  Relative likelihood convergence tolerance. Default `1e-8`.

- seed:

  Optional seed: a whole number
  [`set.seed()`](https://rdrr.io/r/base/Random.html) accepts. The
  caller's random-number state is restored on exit.

## Value

An object of class `multilpa_bootstrap_lrt`, holding the observed
statistic, the p-value, its Monte Carlo standard error, the measurement
blocks held fixed in every fit, and one record per replicate.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives the
single-row test result, `get_results(what = "replicates")` one row per
simulated replicate, [`summary()`](https://rdrr.io/r/base/summary.html)
the test beside the replicate diagnostics, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) the simulated
null distribution with the observed statistic marked. Failed starts,
boundary flags and likelihood replication should be inspected before the
result is interpreted.

## Details

The two models must differ by exactly one individual profile or one
group class, with the observations, group layout, covariance structure,
centering and scaling, missing-data handling and `drop_empty` the same.
Grand-mean centering is repeated in every simulated refit.
Person-centred or person-scaled fits are refused, because the model
specifies no distribution for the group baselines that transformation
removes. If any replicate fails validation (fails, does not converge, or
has a reversed likelihood), the p-value is `NA`; no replicate is
dropped.

Fits made with `missing = "fiml"` are supported when both models were:
every replicate is given the observed data's missing cells before it is
refitted. Calibration requires missingness independent of the indicators
and latent classes, conditional on any fixed covariates. A fixed mask
does not reproduce a MAR mechanism that depends on observed indicators,
or an MNAR mechanism, and FIML estimation under MAR does not by itself
validate this bootstrap under MAR; the size study covers MCAR only.

Membership-covariate fits
(`multilpa(profile_covariates = , group_covariates = )`) are supported:
the covariates are held at their observed values, group classes are
drawn from the fitted group logits and profiles from each row's profile
logits. Both models must use the same profile covariates, group
covariates and slope specification, except that a one-group-class null,
which cannot carry group covariates or slopes by group class, may be
compared with an alternative that adds them. When `data` is omitted it
is rebuilt from the alternative, whose columns include the null's. A
covariate replicate counts as valid when its likelihood has converged
even if a membership logit is still changing, as happens for an empty or
separated class in an over-fitted alternative, because the statistic
uses only the maximised likelihood. The replicate table's
`logits_settled` column records whether both refits' logits had also
converged. The original models must be fully converged.

Models fitted with `fixed`, including those from
[`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md),
are supported: every replicate is refitted with the same blocks held at
the same values. The p-value is then conditional on that measurement
solution, whose own uncertainty is not propagated. A constrained pair
whose nesting cannot be established is refused.

Two transition fits
([`lta()`](https://pak.dynasite.org/latents/reference/lta.md)) or two
group-class family fits are also accepted; each is handled by its own
family's bootstrap.

## Conditions

`latents_bad_nesting` when the models do not differ by exactly one class
at one level; when only one of them holds measurement blocks fixed, when
they hold different blocks or the same blocks at different values, or
when they hold a measurement fixed while differing in the number of
profiles; when the alternative does not estimate more free parameters
than the null; and when two covariate models use different covariates or
slope specifications (other than a one-group-class null).
`latents_incomparable_models` when the models differ in their
observations, grouping, covariance structure, centering, scaling,
missing-data handling or `drop_empty`. `latents_unsupported_bootstrap`
for person-centred or person-scaled fits. `latents_unsupported_prior`
for fits estimated with a prior, whose likelihoods are evaluated at
posterior modes rather than ML estimates. `latents_unsupported_weights`
for fits with sampling weights; `latents_unsupported_noise` for fits
with a noise component. `latents_bad_data` when the data have missing
indicators and the models were not fitted with `missing = "fiml"`;
`latents_bad_inference_data` when `data` does not reproduce a fitted
likelihood or the original observations and their order.
`latents_no_converge` when either original model has not converged or
has an active variance bound; `latents_reversed_likelihood` when the
alternative's log likelihood is below the null's by more than `tol` can
explain; `latents_bad_argument` for an invalid `seed`, or for a
transition fit paired with a fit of another family.
`latents_failed_replicates` is a warning, given when some replicate
fails validation; the p-value is then `NA`.

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
#> p-value: 0.1 (Monte Carlo SE 0.0949) from 9 of 9 valid replicates
```
