# Test whether a model with one more class improves on a nested model, with a parametric bootstrap likelihood ratio

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
  covariates. For the other models `NULL` (default) uses the data stored
  in the fits; a supplied data frame is checked against them.

- iter:

  Number of simulated data sets (default 199; at least 2). More
  replicates reduce the Monte Carlo error of the p-value but do not cure
  estimation problems. Both models are refitted to every data set, so
  the run time grows with `iter` times the cost of the two fits.

- n_starts, max_iter, tol:

  Controls for each simulated refit; `NULL` (default) keeps each model's
  own controls and search strategy. A supplied value applies to both
  models.

- seed:

  Optional integer seed. The caller's random-number state is restored on
  exit.

## Value

An object of class `multilpa_bootstrap_lrt`, for every kind of model
compared. Printing it shows the two models, the observed statistic and
the p-value. Its tables:

- `as.data.frame(result)`: the one-row test summary, with the class
  counts of both models, `statistic` (twice the difference in log
  likelihoods), `p_value`, `monte_carlo_se`, `iter`, `n_valid` and
  `n_invalid` (the replicates used and left out) and `fixed` (the
  measurement blocks held in every refit);

- `get_results(result, "replicates")`: one row per simulated data set,
  with its `statistic`, whether it is `valid`, whether a refit was at a
  bound, how many starts of each refit replicated its best likelihood,
  and any warning or error;

- `summary(result)`: both tables.

`plot(result)` draws the simulated distribution with the observed
statistic.

## Conditions

Errors, raised before anything is simulated: `latents_bad_argument` for
two models of different kinds, an invalid `iter`, `n_starts`, `max_iter`
or `tol`, or transition models without `data`; `latents_bad_nesting` for
profile or class models that do not differ by exactly one profile or one
group class; `latents_incomparable_models` for models that differ in
their observations, grouping, covariance structure, centring, scaling,
missing-data handling or `drop_empty`; `latents_no_converge` for an
original model that did not converge or has a variance at its bound;
`latents_reversed_likelihood` when the alternative has the lower log
likelihood; `latents_bad_inference_data` for `data` that are not the
rows the models were fitted to; `latents_unsupported_weights`,
`latents_unsupported_prior`, `latents_unsupported_noise` and
`latents_unsupported_bootstrap` for fits with sampling weights, a prior,
a noise component, or person centring or scaling.

A `latents_failed_replicates` warning says how many replicates were
invalid and left out; see *Calculation*.

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
indicators <- c("eruptions", "waiting")
smaller <- lpa(faithful, indicators, n_profiles = 1, n_starts = 2, seed = 1)
larger <- lpa(faithful, indicators, n_profiles = 2, n_starts = 2, seed = 1)
# `iter` is small so the example runs quickly; use many more for inference.
comparison <- bootstrap_lrt(smaller, larger, iter = 9, seed = 1)
comparison
#> Parametric bootstrap likelihood-ratio comparison
#> Null: 1 profile; alternative: 2 profiles
#> Observed statistic: 737.798948
#> p-value: 0.1 (Monte Carlo SE 0.0949) from 9 of 9 valid replicates (0 invalid)
as.data.frame(comparison)
#>   null_profiles null_group_classes alternative_profiles
#> 1             1                  1                    2
#>   alternative_group_classes statistic p_value monte_carlo_se iter n_valid
#> 1                         1  737.7989     0.1     0.09486833    9       9
#>   n_invalid fixed null_family alternative_family
#> 1         0  <NA>    profiles           profiles
```
