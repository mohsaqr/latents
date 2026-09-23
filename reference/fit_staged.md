# Fit a two-level model in stages, holding the measurement solution

Estimates the measurement model once, on its own, and then estimates the
group-class structure with that measurement held fixed. This is the
staged workflow used when the measurement solution is meant to be
decided before the grouping question is asked, so that adding group
classes cannot move the profiles they are meant to describe.

## Usage

``` r
fit_staged(
  data,
  vars,
  id,
  n_profiles,
  n_group_classes = 2L,
  measurement = NULL,
  variance_model = c("varying", "equal"),
  n_starts = 10L,
  max_iter = 1000L,
  tol = 1e-08,
  min_variance = 1e-06,
  seed = NULL,
  missing = c("error", "fiml"),
  covariance_model = c("diagonal", "full"),
  categorical = character(),
  min_probability = 1e-10,
  time = NULL
)
```

## Arguments

- data:

  A data frame containing indicators and a group identifier.

- vars:

  Unique character vector of continuous indicator column names.

- id:

  Name of the observed group identifier column. Character, factor, or
  numeric identifiers are supported; missing identifiers are not.
  Required here: this model has a second level by construction, so it
  has no single-level form and does not take
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)'s
  `id = NULL`.

- n_profiles:

  Positive integer number of individual profiles.

- n_group_classes:

  Number of latent group classes to estimate in the second stage. Must
  be at least two; one would leave nothing to estimate.

- measurement:

  Optional fitted `multilpa` model to use as the first stage, so that a
  measurement solution already chosen and inspected is carried forward
  rather than refitted. Its profiles, indicators and measurement options
  must match the ones requested here, and it must have one group class:
  a multi-class fit already estimates group structure and is not a
  measurement-only first stage. When `NULL`, the first stage is fitted
  here with one group class.

- variance_model:

  Either `"varying"` (profile-specific variances or covariance matrices)
  or `"equal"` (shared across profiles).

- n_starts:

  Number of EM starts, used for both stages.

- max_iter:

  Nonnegative integer maximum number of EM updates per start.
  `max_iter = 0` performs no update. With `start`, it returns the model
  evaluated at exactly those values, whatever `n_starts` is, so that
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) scores a parameter
  set supplied from elsewhere rather than the best of some random
  initializations that were never asked for. Without `start` there is
  nothing to evaluate at, so the random initializations are scored and
  the highest is returned.

- tol:

  Positive relative log-likelihood tolerance. Convergence requires
  absolute change no greater than
  `tol * (1 + abs(previous log likelihood))`.

- min_variance:

  Positive lower bound on each variance, or each covariance eigenvalue
  for full covariance, in squared input units. This defines a
  constrained maximum-likelihood problem. Bound-active estimates are
  explicitly reported and generate a warning.

- seed:

  Optional random seed: any whole number
  [`set.seed()`](https://rdrr.io/r/base/Random.html) accepts, negative
  ones included. With a supplied seed, the caller's random-number state
  is restored on exit.

- missing:

  `"error"` rejects missing indicators; `"fiml"` maximizes the
  observed-data likelihood under an ignorable missingness mechanism
  (MAR). Missing indicators are integrated out, not filled in for
  likelihood fitting.

- covariance_model:

  `"diagonal"` assumes conditional independence; `"full"` estimates
  within-profile residual covariances.

- categorical:

  Character vector naming indicators to treat as categorical. Each is
  modelled by unrestricted, profile-specific response probabilities over
  its observed categories, which is the latent class measurement model.
  Binary, ordinal and unordered indicators are all handled by the same
  unrestricted parameterization; numeric, integer, logical, character
  and factor columns are accepted. Indicators not named here stay
  Gaussian, so naming a subset fits a mixed-mode model.

- min_probability:

  Positive lower bound on every categorical response probability,
  defining a constrained maximum-likelihood problem in the same way
  `min_variance` does for Gaussian indicators.

- time:

  Optional name of a column giving each observation's position within
  its group, such as a wave, occasion or course number. The model does
  not use it; it is stored so that `get_results(x, "sequences")`,
  `get_results(x, "sequence_summary")` and `plot(what = "sequences")`
  can read the assignments back in order. Values must be complete and
  unique within each group.

## Value

A `multilpa` object for the second stage, so that every accessor,
diagnostic and method works on it unchanged. Its measurement parameters
are exactly the first stage's. It additionally carries `fixed`, naming
the held blocks; `n_parameters`, counting only the parameters this stage
estimated; and `n_parameters_with_measurement`, which adds the held
measurement back. Use
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
with `what = "stages"` for a tidy two-row summary of both stages.

## Details

First-stage uncertainty is **not** propagated. The second stage treats
the measurement solution as known, so its standard errors, information
criteria and likelihood-ratio comparisons are conditional on that
solution and are narrower than they would be if the measurement had been
estimated jointly. This is a property of staging itself, not of this
implementation, and it is the reason both parameter counts are reported:
compare staged fits with one another using `n_parameters`, and compare a
staged fit with a jointly estimated one using
`n_parameters_with_measurement`, remembering that the staged likelihood
is not the joint maximum and the comparison is descriptive.

## See also

[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
with `fixed` for finer control over which blocks are held, and
[`starting_values()`](https://pak.dynasite.org/latents/reference/starting_values.md)
with `what = "measurement"` for the values a stage hands on.

## Examples

``` r
set.seed(7)
example_data <- data.frame(
  school = rep(seq_len(30), each = 8),
  score_a = stats::rnorm(240), score_b = stats::rnorm(240)
)
staged <- fit_staged(example_data, c("score_a", "score_b"), "school",
                     n_profiles = 2, n_group_classes = 2, n_starts = 2,
                     seed = 1)
get_results(staged, what = "stages")
#>         stage group_classes            fixed log_likelihood parameters
#> 1 measurement             1             <NA>      -671.9067          9
#> 2  membership             2 means, variances      -670.9827          3
#>   parameters_with_measurement converged
#> 1                           9      TRUE
#> 2                          11      TRUE
get_results(staged, what = "profile_probabilities")
#>   group_class profile  probability group_class_probability
#> 1           1       1 9.999951e-01                0.844672
#> 2           1       2 4.946711e-06                0.844672
#> 3           2       1 8.991468e-01                0.155328
#> 4           2       2 1.008532e-01                0.155328
```
