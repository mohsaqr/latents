# Fit a two-level profile model in two stages, holding the profiles fixed while the group classes are estimated

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
  without missing values). Required: a staged fit is a two-level model.

- n_profiles:

  Number of individual profiles.

- n_group_classes:

  Number of group classes in the second stage, at least two (default 2).

- measurement:

  Where the profiles come from:

  - `NULL` (default): the first stage is estimated from `data`, with one
    group class;

  - a fitted `multilpa` model: its measurement estimates are used and
    only the second stage is estimated. It must have one group class and
    the same profiles, indicators and measurement options as this call.
    Its covariance structure must be EEI, VVI, EEE or VVV, and it must
    not have `weights`, `centering`, `scaling`, `ordinal` or `count`
    indicators, `direct_effects`, `local_dependence` or `noise`.

- variance_model, covariance_model:

  The covariance structure of the continuous indicators, as under
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md):
  `variance_model` is `"varying"` (default, variances per profile) or
  `"equal"` (shared by all profiles), and `covariance_model` is
  `"diagonal"` (default, no within-profile covariances) or `"full"`.
  These give the structures VVI, EEI, VVV and EEE; a staged fit takes no
  other structure.

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

  What to do with missing indicator values, in both stages:

  - `"error"` (default): stop, naming the indicators that have them;

  - `"fiml"`: keep every row and fit each one on the indicators it
    observes (full-information maximum likelihood), valid when the
    missingness depends only on observed values.

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

  How the starting values of each estimated stage are searched, as under
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md):

  - `"kmeans"`: every one of the `n_starts` starts is run to
    convergence;

  - `"emEM"`: a larger pool of `n_screen` starts is screened with a
    short EM run each, and only the first start and the best of the
    rest, `n_starts` in all, are run to convergence;

  - `NULL` (default): `"emEM"` for continuous indicators with a shared
    covariance matrix (EEI or EEE), `"kmeans"` otherwise.

  The second stage starts from the first-stage measurement and holds it.

- n_screen:

  For `search = "emEM"` only, the size of the screened pool (default
  five times `n_starts`, at least `n_starts`). The second stage skips
  screening when `n_starts = 1`.

## Value

A fitted model of class `multilpa`: the second-stage fit, whose
measurement parameters equal the first-stage estimates. Every table,
plot and diagnostic of a
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
fit applies. What a staged fit adds:

- `get_results(fit, what = "stages")`: two rows, `measurement` (the
  first stage) and `membership` (the second), each with its number of
  group classes, the held blocks (`fixed`), `log_likelihood`,
  `parameters`, `parameters_with_measurement` and `converged`;

- `get_results(fit, what = "profile_probabilities")`: the profile
  proportions within each group class, the result of the second stage;

- `n_parameters` counts the second-stage parameters only, which the
  information criteria use; `n_parameters_with_measurement` adds the
  held measurement parameters.

## Conditions

`latents_bad_argument` for `n_group_classes` below two, `id = NULL`, or
a `measurement` that is not a `multilpa` fit; `latents_bad_stage` for a
`measurement` with group classes, with a different number of profiles,
other indicators, categories or covariance structure than this call, or
with an option the second stage does not use (listed under
`measurement`); `latents_unsupported_noise` for a `measurement` with a
noise component. Both stages raise the conditions of
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).

## See also

[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
for holding chosen blocks through `fixed`;
`starting_values(fit, what = "measurement")` to extract measurement
parameters for a later fit.

## Examples

``` r
# Stage 1: profiles of course sessions. Stage 2: classes of students that
# differ in how often their sessions are in each profile.
staged <- fit_staged(course_engagement, c("browse", "lectures"),
                     id = "student", n_profiles = 2, n_group_classes = 2,
                     n_starts = 2, seed = 1)
get_results(staged, what = "stages")
#>         stage group_classes            fixed log_likelihood parameters
#> 1 measurement             1             <NA>      -3894.214          9
#> 2  membership             2 means, variances      -3800.754          3
#>   parameters_with_measurement converged
#> 1                           9      TRUE
#> 2                          11      TRUE
get_results(staged, what = "profile_probabilities")
#>   group_class profile probability group_class_probability
#> 1           1       1   0.7840790               0.6747345
#> 2           1       2   0.2159210               0.6747345
#> 3           2       1   0.1792617               0.3252655
#> 4           2       2   0.8207383               0.3252655
```
