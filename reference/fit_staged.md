# Estimate group classes with fixed individual profile definitions

`fit_staged()` fits a two-level profile model in two stages. The first
stage estimates the profiles with one group class; the second holds
their measurement parameters (means and covariances of continuous
indicators, response probabilities of categorical ones) fixed and
estimates the profile probabilities within each group class and the
group-class proportions. The profiles therefore keep one definition
while the group classes differ in composition. A fitted measurement
model can be supplied instead of the first stage, to reuse established
profiles. The result holds the second-stage estimates and keeps the
first-stage fit.

Information criteria and likelihood-ratio comparisons condition on the
held measurement. Staged fits that hold the same measurement are
compared on `n_parameters`. A comparison with a jointly estimated fit
uses `n_parameters_with_measurement` and stays descriptive, since staged
estimates do not maximise the joint likelihood. Wald standard errors
omit the uncertainty of the first stage.
`parameter_inference(fit, method = "bootstrap")` resamples groups and
re-estimates both stages (a supplied same-sample `measurement` model
with its own options), so that uncertainty is propagated and the
measurement parameters get intervals too. A measurement fitted to
different observations, or with unknown provenance, is refused by this
bootstrap; resampling the current data cannot represent uncertainty from
an external measurement sample.

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

  A data frame with the indicators and the group identifier.

- vars:

  Names of the indicator columns. Those not named in `categorical` are
  continuous.

- id:

  Name of the group identifier column (character, factor or numeric,
  without missing values). Required.

- n_profiles:

  Number of individual profiles.

- n_group_classes:

  Number of group classes in the second stage, at least two (default 2).

- measurement:

  `NULL` (default) estimates the first stage from `data`; otherwise a
  fitted `multilpa` model supplying the measurement estimates, with one
  group class and the same profiles, indicators and measurement options.
  Its covariance structure must be EEI, VVI, EEE or VVV, and it must not
  have `weights`, `centering`, `scaling`, `ordinal` or `count`
  indicators, `direct_effects`, `local_dependence` or `noise`.

- variance_model:

  `"varying"` (default): each profile has its own variances or
  covariance matrix; `"equal"`: one shared by all profiles.

- n_starts:

  Number of EM starts in each estimated stage (default 10). With
  `measurement`, only the second stage is estimated.

- max_iter:

  Maximum EM updates per start in each stage (default 1000); zero
  evaluates the starting values without updating them.

- tol:

  Convergence tolerance (default `1e-8`), per observed indicator value:
  EM stops when a step changes the log likelihood by no more than `tol`
  times the number of observed indicator values.

- min_variance:

  Lower bound on every variance (on every covariance eigenvalue under
  full covariance), as in
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  applied in both stages; a supplied `measurement` must satisfy it.

- seed:

  Optional integer seed. The caller's random-number state is restored on
  exit.

- missing:

  `"error"` (default) rejects missing indicator values; `"fiml"`
  maximises the observed-data likelihood under MAR, integrating the
  missing values out.

- covariance_model:

  `"diagonal"` (default): the indicators are independent within a
  profile; `"full"`: within-profile covariances are estimated.

- categorical:

  Indicators modelled by profile-specific response probabilities over
  their observed categories, the latent class measurement model. Binary,
  ordered and unordered indicators share this parameterization. Naming a
  subset of `vars` fits a mixed-mode model.

- min_probability:

  Lower bound on every categorical response probability (default
  `1e-10`).

- time:

  Optional column giving each row's occasion within its group, used only
  to order the rows for `get_results(x, "sequences")`,
  `"sequence_summary"` and
  [`plot_sequences()`](https://pak.dynasite.org/latents/reference/plot_sequences.md).
  Values must be unique within a group and sort in their intended order
  (numbers numerically, factors by level, text alphabetically, so
  `"T2"`, `"T10"` is refused: use a factor or zero-pad).

- search:

  `NULL` (default), `"kmeans"` or `"emEM"` for each estimated stage, as
  under
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md):
  the default is `"emEM"` for continuous indicators with a shared
  covariance matrix (EEI or EEE) and `"kmeans"` otherwise. The second
  stage starts from the first-stage measurement and holds it.

- n_screen:

  For `search = "emEM"` only, the size of the screened pool (default
  five times `n_starts`, at least `n_starts`). The second stage skips
  screening when `n_starts = 1`.

## Value

A second-stage `multilpa` object whose measurement parameters equal the
first-stage estimates; `fixed` names the held blocks, `n_parameters`
counts the second-stage parameters and `n_parameters_with_measurement`
adds the measurement ones. `get_results(fit, what = "stages")`
summarises the two stages in two rows. Every table, plot and diagnostic
of a `multilpa` fit applies.

## See also

[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
for holding chosen blocks through `fixed`;
`starting_values(fit, what = "measurement")` to extract measurement
parameters for a later fit.

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
