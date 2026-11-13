# Parametric bootstrap likelihood-ratio test for nested mixture models

`bootstrap_lrt()` tests whether a more flexible model improves on a
simpler one nested in it, typically one more profile or group class. The
statistic is twice the difference in maximised log likelihoods; its null
distribution is simulated by generating data from the fitted null model
and refitting both models to each data set, since the chi-squared
approximation does not hold when the number of mixture components
changes. A small p-value is evidence against the null in favour of this
alternative; it does not say the alternative has the right number of
classes.

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

  A converged fit under the null: a profile or class model with or
  without membership covariates, an
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit, or a
  group-class family fit from
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).
  See *Supported comparisons*.

- alternative_model:

  A converged, more flexible model fitted to the same data, in which
  `null_model` is nested.

- data:

  The data frame both models were fitted to. Required for transition
  models, whose simulations keep the observed groups, occasions and
  covariates; otherwise taken from the fits when omitted, and checked
  against them when supplied.

- iter:

  Number of simulated data sets (default 199; at least 2). More
  replicates reduce the Monte Carlo error of the p-value but do not cure
  estimation problems.

- n_starts, max_iter, tol:

  Controls for each simulated refit; `NULL` (default) keeps each model's
  own controls and search strategy. A supplied value applies to both
  models.

- seed:

  Optional integer seed. The caller's random-number state is restored on
  exit.

## Value

A `multilpa_bootstrap_lrt` object: the observed statistic, the bootstrap
p-value and its Monte Carlo standard error, `n_valid` and `n_invalid`
replicates, `fixed` (the measurement blocks held in every refit) and
per-replicate diagnostics, including how many starts of each refit
replicated its best likelihood.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the one-row test summary, `get_results(result, "replicates")` one row
per replicate, [`summary()`](https://rdrr.io/r/base/summary.html) both,
and [`plot()`](https://rdrr.io/r/graphics/plot.default.html) the
simulated distribution with the observed statistic.

## Supported comparisons

Both models must be of one kind: profile or class models, transition
models, or group-class family fits. Profile and class models, with or
without membership covariates, must differ by exactly one profile
(class) or exactly one group class, and share their observations, group
layout, covariance structure, centring and scaling, missing-data
handling and `drop_empty`. `centering = "grand"` is repeated in each
refit; person-centred or person-scaled fits are refused, because the
model gives the group means removed by those transformations no
distribution. Covariate models must share their covariates and slope
specification (a one-group-class null may be compared with an
alternative that adds group covariates or group-class-specific slopes);
covariates are held at their observed values when simulating, and both
originals must be fully converged, logits included.

Transition models may differ in class counts or in restrictions on
measurement, transitions, Markov order or covariate effects, with
compatible indicators, grouping, occasion layout and measurement
distributions. Group-class family fits may differ in the number of group
classes or in restrictions on the means and variance components, the
alternative being the less restricted one. A profile or class fit with
`fixed` (including a
[`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
result) keeps the same held blocks in every replicate, so the p-value is
conditional on that measurement solution; a covariate fit's held blocks
are re-estimated in each replicate. Fits with a prior, sampling weights
or a noise component are refused.

## Calculation

The p-value is `(1 + b) / (V + 1)`, with `b` the number of valid
simulated statistics at least as large as the observed one and `V` the
number of valid replicates. A replicate is invalid when a refit fails or
does not converge, or when the alternative's likelihood falls below the
null's beyond numerical tolerance; such replicates are left out with a
warning, and the p-value is `NA` when fewer than half of the requested
replicates are valid. Leaving them out conditions the reference
distribution on successful estimation. Report the invalid share, and
look at whether the refits replicated their likelihoods across starts.
For covariate models a refit counts as converged when its likelihood has
converged even if a membership logit is still moving, as happens with an
empty or separated class in the alternative; `logits_settled` records
whether the logits converged too.

Fits with `missing = "fiml"` are supported when both models use it; each
simulated data set keeps the original missing-value pattern. That is
calibrated only when missingness is independent of the indicators and
the classes given any fixed covariates: a fixed pattern does not
reproduce an MAR mechanism that depends on observed indicators, nor an
MNAR one, so valid FIML estimation does not by itself make this test
valid.

## See also

[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
for fitting the candidates;
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
