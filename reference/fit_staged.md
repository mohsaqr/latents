# Estimate group classes with fixed individual profile definitions

`fit_staged()` estimates a multilevel mixture model in which observed
groups are distinguished by their composition of individual profiles. It
separates the estimation of profile definitions from the estimation of
group classes, allowing between-group differences to be studied while
retaining a common measurement model.

The first stage estimates the individual profiles using a model with one
group class. These profiles are defined by their indicator
distributions, including means and covariance parameters for continuous
indicators and response probabilities for categorical indicators. The
second stage holds these measurement parameters fixed and estimates the
profile probabilities within each group class, together with the
proportions of the group classes. Group classes therefore differ in
their profile composition while the profiles retain the same measurement
definitions.

A compatible fitted measurement model may also be supplied, allowing
established profile definitions to be reused. The returned model
contains the second-stage estimates and retains the first-stage fit, so
both stages can be inspected.

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
  min_variance = NULL,
  seed = NULL,
  missing = c("error", "fiml"),
  covariance_model = c("diagonal", "full"),
  categorical = character(),
  min_probability = 1e-10,
  time = NULL,
  search = NULL,
  n_screen = NULL
)
```

## Arguments

- data:

  A data frame containing the indicators and the group identifier.

- vars:

  Unique character vector of indicator column names. Indicators named in
  `categorical` use categorical response distributions; the remaining
  indicators are continuous and use Gaussian distributions.

- id:

  Name of the observed group identifier column. Character, factor and
  numeric identifiers are supported; missing identifiers are not.
  Required; `id = NULL` is not accepted.

- n_profiles:

  Positive integer number of individual profiles.

- n_group_classes:

  Number of latent group classes in the second stage. Must be at least
  two; the default is 2.

- measurement:

  An optional fitted `multilpa` model supplying the first-stage
  measurement estimates. With `NULL` (the default), the first stage is
  estimated from `data`. A supplied model must have one group class and
  matching profiles, indicators and measurement options. Models with
  `weights`, `centering`, `scaling`, `ordinal` or `count` indicators,
  `direct_effects`, or `local_dependence` between continuous indicators
  are not supported. Supported covariance structures are EEI, VVI, EEE
  and VVV. An incompatible measurement model raises `latents_bad_stage`.

- variance_model:

  `"varying"` (the default; profile-specific variances or covariance
  matrices) or `"equal"` (shared across profiles).

- n_starts:

  Number of EM starts used in each estimated stage. When `measurement`
  is supplied, only the second stage is estimated.

- max_iter:

  Nonnegative integer maximum number of EM updates per start in each
  estimated stage (default 1000). Zero evaluates initial values without
  performing EM updates.

- tol:

  Positive convergence tolerance (default `1e-8`). EM has converged when
  one step (or SQUAREM cycle) changes the log likelihood by no more than
  `tol` times the number of observed indicator values (each non-missing
  cell in a row of positive weight, with a minimum count of one). The
  tolerance is therefore expressed per observed indicator value.

- min_variance:

  Lower bound on each continuous-indicator variance, or on each
  covariance eigenvalue under full covariance. `NULL` (the default) uses
  `1e-6` times the smallest observed sample variance of the continuous
  indicators. A numeric value specifies an absolute bound in squared
  indicator units. Estimates at the bound raise a warning. In the second
  stage, measurement parameters are held at their first-stage values.

- seed:

  Optional random seed: any integer
  [`set.seed()`](https://rdrr.io/r/base/Random.html) accepts, including
  negative integers. When supplied, the caller's random-number state is
  restored on exit.

- missing:

  `"error"` (the default) rejects missing indicators; `"fiml"` maximises
  the observed-data likelihood under an ignorable missingness mechanism
  (MAR), integrating over missing indicators.

- covariance_model:

  `"diagonal"` (the default) assumes conditional independence; `"full"`
  estimates within-profile residual covariances.

- categorical:

  Character vector naming indicators to treat as categorical. Each is
  modelled by unrestricted, profile-specific response probabilities over
  its observed categories (the latent class measurement model); binary,
  ordinal and unordered indicators are handled by the same
  parameterization. Numeric, integer, logical, character and factor
  columns are accepted. Indicators not named here stay Gaussian, so
  naming a subset fits a mixed-mode model.

- min_probability:

  Positive lower bound on every categorical response probability
  (default `1e-10`), constraining the maximum-likelihood problem as
  `min_variance` does for Gaussian indicators.

- time:

  Optional name of a column giving each observation's position within
  its group, such as a wave, occasion or course number. It is not used
  in estimation; it is stored for `get_results(x, "sequences")`,
  `get_results(x, "sequence_summary")` and `plot(what = "sequences")`.
  Values must be complete and unique within each group, and are ordered
  as [`sort()`](https://rdrr.io/r/base/sort.html) orders them: numbers
  numerically, a factor by its levels, and character strings
  alphabetically. Numbers stored as text whose alphabetical order is not
  their numeric order (`"2"`, `"10"`), and numbered labels sharing their
  text around the number whose alphabetical order is not the numbers'
  (`"T2"`, `"T10"`; `"week 2"`, `"week 10"`), raise `latents_bad_time`
  (also `latents_bad_data`); give the order with `factor(levels = )` or
  zero-pad the numbers. Zero-padded numbers and ISO dates are accepted.

- search:

  Initialisation strategy for each estimated stage: `NULL` (the
  default), `"kmeans"` or `"emEM"`. The default is `"emEM"` for models
  containing continuous indicators with a shared covariance matrix (EEI
  or EEE), and `"kmeans"` otherwise. `"kmeans"` uses clustering-based
  starts; `"emEM"` screens a pool of starts with short EM runs before
  completing `n_starts` selected searches. The second stage uses the
  first-stage measurement estimates as supplied starting values and
  holds them fixed. See
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  for the screening and start-selection procedure.

- n_screen:

  Number of starts screened under `search = "emEM"`, at least
  `n_starts`. The default is five times `n_starts`. Screening is skipped
  in the second stage when `n_starts = 1`, because a starting
  measurement solution is supplied. A value with `search = "kmeans"`, or
  one below `n_starts`, raises `latents_bad_argument`.

## Value

A second-stage `multilpa` object with measurement parameters equal to
the first-stage estimates. The `fixed` component identifies the
parameter blocks held constant. `n_parameters` counts the parameters
estimated in the second stage; `n_parameters_with_measurement` also
includes the first-stage measurement parameters.
`get_results(fit, what = "stages")` returns a two-row summary of the
estimation stages. Standard result-extraction and diagnostic methods for
`multilpa` objects are available.

## Details

Information criteria and likelihood-ratio comparisons condition on the
fixed measurement model.

With `parameter_inference(fit, method = "bootstrap")`, groups are
resampled and both stages are re-estimated. This procedure propagates
first-stage estimation uncertainty and reports standard errors and
intervals for the measurement parameters as well as the second-stage
parameters. A supplied `measurement` model is re-estimated in each
resample using its fitting options.

For staged models with the same fixed measurement specification, use
`n_parameters` when comparing conditional information criteria.
Comparisons with jointly estimated models can include the measurement
parameters through `n_parameters_with_measurement`, but remain
descriptive because the staged estimates do not generally maximize the
joint likelihood.

## Interpretation

The second-stage estimates are conditional on the fixed measurement
model. Default Wald standard errors omit uncertainty from estimating
that model; bootstrap inference can propagate this uncertainty by
resampling groups and re-estimating both stages. Staged estimates do not
generally maximise the joint likelihood, so comparisons with jointly
estimated models require care.

## See also

[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
for selective constraints through `fixed`;
[`starting_values()`](https://pak.dynasite.org/latents/reference/starting_values.md)
with `what = "measurement"` to extract measurement parameters for a
subsequent fit.

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
#> 1 measurement             1             <NA>      -671.8928          9
#> 2  membership             2 means, variances      -671.8928          3
#>   parameters_with_measurement converged
#> 1                           9      TRUE
#> 2                          11      TRUE
get_results(staged, what = "profile_probabilities")
#>   group_class profile probability group_class_probability
#> 1           1       1   0.9782787                     0.5
#> 2           1       2   0.0217213                     0.5
#> 3           2       1   0.9782787                     0.5
#> 4           2       2   0.0217213                     0.5
```
