# Parametric bootstrap likelihood-ratio comparison

Simulates complete indicators under the null model while preserving
observed group sizes, refits both models, and compares their likelihood
differences. Models must differ by exactly one individual profile or one
group class, with the other count, covariance structure and centering
mode fixed. Grand-mean centering is repeated in every simulated refit.
Person-centred fits are refused because this model does not specify a
generative distribution for the group baselines removed by that
transformation. This is a native parametric bootstrap, not an
implementation of Mplus TECH14. It does not use a chi-square reference
distribution. Any failed/nonconverged or reversed replicate makes the
p-value NA, avoiding silent deletion of difficult fits.

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

  Smaller, converged
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  model on complete data.

- alternative_model:

  Larger model fitted to exactly the same data.

- data:

  Optional. The data both models were fitted to, used to verify both
  fitted likelihoods; when omitted it is rebuilt from what the null
  model stores.

- iter:

  Number of simulated datasets (at least two; use many for inference).

- n_starts:

  Number of starts for each simulated fit.

- max_iter:

  Maximum EM iterations for each simulated fit.

- tol:

  Relative likelihood convergence tolerance.

- seed:

  Optional seed: any whole number
  [`set.seed()`](https://rdrr.io/r/base/Random.html) accepts. The
  caller's random-number state is restored on exit.

## Value

An object of class `multilpa_bootstrap_lrt`, carrying the observed
statistic, the finite-simulation corrected p-value, its Monte Carlo
standard error, the measurement blocks held fixed in every fit, and one
record per replicate. Read it with the verbs that describe it:
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives the
single-row test result, `get_results(what = "replicates")` one row per
simulated replicate, [`summary()`](https://rdrr.io/r/base/summary.html)
the test beside the replicate diagnostics, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) the simulated
null distribution with the observed statistic marked. Inspect failed
starts, boundary flags and likelihood replication before interpreting
results.

## Details

Models fitted with `fixed`, including those from
[`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md),
are supported: every replicate is refitted with the same blocks held at
the same values, so the simulated statistics compare the same two
constrained models as the observed statistic does. The p-value is then
conditional on that measurement solution, which is treated as known and
whose own uncertainty is not propagated. A constrained pair whose
nesting cannot be established is refused rather than given a p-value.

## Conditions

`latents_bad_nesting` when only one of the two models holds measurement
blocks fixed, when they hold different blocks, when they hold the same
blocks at different values, when they hold a measurement fixed while
differing in the number of profiles (the held blocks then have different
shapes, so the null is not a restriction of the alternative), or when
the alternative does not estimate more free parameters than the null.
`latents_failed_replicates` is warned when some replicate fails
validation, and the p-value is `NA`. `latents_unsupported_bootstrap`
refuses person-centred fits, for which this simulator has no
group-baseline distribution.

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
comparison <- bootstrap_lrt(smaller, larger, iter = 9, n_starts = 2,
                            max_iter = 2000, tol = 1e-6, seed = 1)
comparison
#> Parametric bootstrap likelihood-ratio comparison
#> Null: 1 profiles, 1 group classes; alternative: 2 profiles, 1 group classes
#> Observed statistic: 44.198110
#> p-value: 0.1 (Monte Carlo SE 0.0949) from 9 of 9 valid replicates
```
