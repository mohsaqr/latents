# Fit a two-level latent profile model

Fits individual Gaussian profiles nested within observed groups. A
discrete latent group class determines the profile proportions. Profile
means and variances are shared across group classes (measurement
invariance), and indicators are independent conditional on individual
profile membership by default. Full residual covariance and
observed-data maximum likelihood for missing indicators are available.

## Usage

``` r
multilpa(
  data,
  vars,
  id,
  n_profiles,
  n_group_classes = 2L,
  profile_covariates = character(),
  group_covariates = character(),
  variance_model = c("varying", "equal"),
  n_starts = 10L,
  max_iter = 1000L,
  tol = 1e-08,
  min_variance = 1e-06,
  seed = NULL,
  start = NULL,
  missing = c("error", "fiml"),
  covariance_model = c("diagonal", "full"),
  categorical = character(),
  min_probability = 1e-10,
  time = NULL,
  fixed = character(),
  centering = c("none", "person", "grand"),
  volume = NULL,
  shape = NULL,
  orientation = NULL
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

  `id` has no default: omitting it raises `latents_bad_argument`,
  because a forgotten grouping would otherwise be fitted as a different
  model without saying so. Passing `id = NULL` explicitly fits a
  **single-level** model: the observations are treated as independent,
  each row is its own unit, and `n_group_classes` becomes one. That fit
  raises a `latents_single_level` warning, since this package exists for
  the two-level model. This is the ordinary Gaussian or latent-class
  mixture that the two-level model reduces to, and every verb of this
  package works on it. The unit column is fabricated internally as
  `.observation`; it is not returned by `get_results(x, "data")`, and a
  `data` that already has a column of that name raises
  `latents_bad_data`. Asking for more than one group class without an
  `id` raises `latents_bad_argument`, because one observation per unit
  leaves no composition for a second-level class to differ in.

- n_profiles:

  Positive integer number of individual profiles.

- n_group_classes:

  Positive integer number of latent group classes. With `id = NULL` it
  is one, and naming anything else is an error.

- profile_covariates, group_covariates:

  Names of numeric columns of `data` predicting individual-profile and
  group-class membership through multinomial logits, with the final
  class as reference. Naming either one fits the one-step covariate
  model and returns a `multilpa_covariates` object: profile slopes are
  shared across group classes, profile intercepts differ by group class,
  and a `group_covariate` must be constant within each group. Covariates
  enter in their supplied units, so centre or scale them beforehand if
  that is what you want. This is one-step maximum likelihood, not a
  regression on assigned classes;
  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
  and [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md)
  are the staged alternatives. The covariate path supports neither
  `start`, nor `missing = "fiml"`, nor `fixed`, and refuses them by name
  rather than ignoring them. Use one for an ordinary, pooled LPA.

- variance_model:

  Either `"varying"` (profile-specific variances or covariance matrices)
  or `"equal"` (shared across profiles).

- n_starts:

  Positive integer number of EM starts. When `start` is supplied, it
  supplies the first start; remaining starts are random initializations.
  It is ignored when `max_iter = 0` and `start` is supplied: that call
  evaluates the supplied parameters and nothing else, so exactly one
  start is run whatever `n_starts` says.

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

- start:

  Optional list of `means` and `variances` (profiles by indicators),
  `profile_probabilities` (group classes by profiles), and
  `group_probabilities` (vector). Starting probabilities must be
  positive. For full covariance, supply `covariances` (indicators by
  indicators by profiles); `variances` may be omitted or must match
  their diagonals. A categorical model also takes
  `response_probabilities`, a list of one profiles-by-categories matrix
  per indicator named in `categorical`. That list is read by position
  unless it is labelled: name its elements after the indicators, or its
  columns after the categories, and each block is matched to the
  indicator and category its labels name, so a differently ordered
  `categorical` or a differently ordered set of factor levels cannot
  attach a distribution to the wrong item. A label naming an indicator
  or a category this fit does not have raises `latents_bad_start` rather
  than being aligned by position.

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

- fixed:

  Character vector naming measurement blocks to hold at the values
  `start` supplies, instead of estimating them: any of `"means"`,
  `"variances"` and `"response_probabilities"`, or `"measurement"` for
  every block the model has. Under `covariance_model = "full"`,
  `"variances"` holds the residual covariance matrices. A held block
  stays exactly as supplied, in every restart, and stops counting
  towards `n_parameters`, so this is a different model rather than a
  different starting point for the same one. `start` must carry the
  named blocks; `starting_values(fit, what = "measurement")` produces
  them, and
  [`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
  wraps the whole two-stage workflow in one call. The mixing parameters
  cannot be held: they are what a fixed-measurement fit is for.

- centering:

  How to centre the continuous indicators before fitting. `"none"`, the
  default, fits them as supplied. `"person"` subtracts each group's own
  mean from its rows, so a value reads as a deviation from that unit's
  average and the profiles become profiles of *change* rather than of
  level: this is the within-person, person-mean-centred design
  (Quintana, 2021; Voelkle, Brose, Schmiedek, & Lindenberger, 2014).
  `"grand"` subtracts one mean per indicator, which moves the origin
  without touching the within-group structure. The offsets are kept on
  the fit, so `get_results(x, "data")` still returns the columns you
  supplied and every verb that checks row alignment still checks it.
  Centring removes exactly the between-unit variation, so `"person"`
  refuses with `latents_bad_data` when it leaves an indicator constant —
  which is what happens when a unit has one observation of it. With
  `"person"` the group classes become types of *change pattern*, not
  types of unit.

- volume, shape, orientation:

  The covariance structure, in the three pieces it is made of. Each
  profile's covariance decomposes as
  `Sigma_k = lambda_k * D_k * A_k * D_k'`: a *volume*
  `lambda_k = |Sigma_k|^(1/d)`, an *orientation* `D_k` of eigenvectors,
  and a *shape* `A_k`, diagonal with determinant one. Constraining the
  three across profiles gives the fourteen models `mclust` names with
  three letters, and this package fits all of them.

  `volume` is `"equal"` or `"varying"`. `shape` is `"equal"`,
  `"varying"` or `"spherical"`, the last making every indicator's spread
  equal within a profile, which leaves no orientation to constrain.
  `orientation` is `"axis"` (axis-parallel, a diagonal covariance),
  `"equal"` (one orientation shared by every profile) or `"varying"`.
  Each is `NULL` by default, which follows `variance_model` and
  `covariance_model`, so a call that names none of them fits exactly
  what it always did.

  |               |               |               |                    |
  |---------------|---------------|---------------|--------------------|
  | `volume`      | `shape`       | `orientation` | model              |
  | equal         | spherical     | —             | EII                |
  | varying       | spherical     | —             | VII                |
  | equal/varying | equal/varying | axis          | EEI, VEI, EVI, VVI |
  | equal/varying | equal/varying | equal         | EEE, VEE, EVE, VVE |
  | equal/varying | equal/varying | varying       | EEV, VEV, EVV, VVV |

  The estimates are those of Celeux and Govaert (1995); the two models
  with a shared orientation and a free shape, EVE and VVE, have no
  closed form and use the minorize-maximize step of Browne and
  McNicholas (2014). Parameter counts match `mclust`'s own for all
  fourteen.

  Anything other than EEI, VVI, EEE or VVV is maximized across every
  profile at once, so it cannot be combined with a held `variances`
  block, and `parameter_inference(method = "wald")` refuses it with
  `latents_unsupported_inference`: the free coordinates are log
  variances, which is the wrong chart for a constrained volume, shape or
  orientation. `parameter_inference(method = "bootstrap")` reports all
  fourteen: it resamples groups and refits inside the same family, so it
  needs no chart.

## Value

An `multilpa` object containing `means`, `variances`, optional
`covariances` (indicators by indicators by profiles),
`profile_probabilities`, `group_probabilities`, posterior matrices,
classifications, log likelihood, information criteria, restart
diagnostics, and convergence history. Individual rows retain their input
order; groups retain first-occurrence order. `bic` and `bic_groups` use
the observed group count; `bic_individual` uses `n_informative`, the
number of rows carrying at least one observed indicator. These are
alternative conventions, not interchangeable criteria. `boundary`
identifies variance bounds; `small_classes` flags effective memberships
below one. The fit itself carries no standard errors:
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
computes them from the fit and the data it was fitted to, and
[`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
tests nested models. No guarantee of global optimality is given,
whatever `n_starts` is used. Read the tidy form with
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html); `what`
selects which table.

## Details

Infinite and constant observed indicators are rejected. Each indicator
must have at least two distinct observed values. No rows are silently
dropped. In FIML mode, fully missing individuals contribute no direct
measurement likelihood but receive posterior probabilities from their
group's information. Observation counts retain these individuals; the
individual-level BIC excludes them through `n_informative`, while the
group-level BIC uses all observed groups. Missing patterns can prevent
parameter identification; no general identification guarantee is made.
Initialization alone uses indicator-mean filling. EM uses conditional
Gaussian sufficient statistics and observed marginal densities. More
than one group class requires more than one profile and at least one
group with multiple individuals; these checks are necessary but do not
establish identification. The highest finite likelihood across starts is
returned, even if that start did not converge; inspect `converged` and
`starts`. Profile and group-class labels are arbitrary.

## References

Vermunt, J. K. (2003). Multilevel latent class models. Sociological
Methodology, 33, 213–239. doi:10.1111/j.0081-1750.2003.t01-1-00131.x.

## See also

[`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
for a model in which every indicator is categorical.

## Examples

``` r
set.seed(7)
# Two kinds of school, differing only in how often a pupil scores highly.
example_data <- data.frame(
  school = rep(seq_len(24), each = 10),
  school_type = rep(c("mixed", "high"), each = 120)
)
example_data$high <- rbinom(240, 1L,
  ifelse(example_data$school_type == "high", 0.8, 0.2))
example_data$score_a <- rnorm(240, mean = 2 * example_data$high)
example_data$score_b <- rnorm(240, mean = 2 * example_data$high)

fit <- multilpa(example_data, c("score_a", "score_b"), "school",
                n_profiles = 2, n_group_classes = 2, n_starts = 4,
                seed = 42)
summary(fit)
#> Multilevel LPA: 2 profiles and 2 group classes
#> Individuals: 240; groups: 24; parameters: 11; converged: TRUE
#> Log likelihood: -784.609205; AIC: 1591.218
#> BIC (groups): 1604.177; BIC (individuals): 1629.505
#> Best likelihood replicated in 4/4 starts (absolute tolerance 0.000786).
#> 
#> -- profiles --------------------------------------------------------
#>  profile indicator     mean variance standard_deviation
#>        1   score_a  1.82743   1.4076             1.1864
#>        1   score_b  1.96280   0.9206             0.9595
#>        2   score_a  0.07696   0.9562             0.9779
#>        2   score_b -0.04159   1.0240             1.0119
#> 
#> -- responses -------------------------------------------------------
#>    (no rows)
#> 
#> -- covariances -----------------------------------------------------
#>  profile indicator indicator_2 covariance
#>        1   score_a     score_a     1.4076
#>        1   score_b     score_a     0.0000
#>        1   score_a     score_b     0.0000
#>        1   score_b     score_b     0.9206
#>        2   score_a     score_a     0.9562
#>        2   score_b     score_a     0.0000
#>        2   score_a     score_b     0.0000
#>        2   score_b     score_b     1.0240
#> 
#> -- profile_probabilities -------------------------------------------
#>  group_class profile probability group_class_probability
#>            1       1      0.8101                   0.516
#>            1       2      0.1899                   0.516
#>            2       1      0.1717                   0.484
#>            2       2      0.8283                   0.484
#> 
#> -- counts ----------------------------------------------------------
#>        level class effective_count effective_proportion
#>  individuals     1          120.27               0.5011
#>  individuals     2          119.73               0.4989
#>       groups     1           12.38               0.5160
#>       groups     2           11.62               0.4840
#> 
#> -- posteriors ------------------------------------------------------
#>  row group profile posterior modal
#>    1     1       1 0.5306745  TRUE
#>    2     1       1 0.0004077 FALSE
#>    3     1       1 0.0169366 FALSE
#>    4     1       1 0.0073025 FALSE
#>    5     1       1 0.0813076 FALSE
#>    6     1       1 0.0197391 FALSE
#>    7     1       1 0.0274373 FALSE
#>    8     1       1 0.9720359  TRUE
#>    9     1       1 0.1178553 FALSE
#>   10     1       1 0.1704138 FALSE
#>    ... 470 more rows.  get_results(x, what = "posteriors")
#> 
#> -- group_posteriors ------------------------------------------------
#>  group group_size log_likelihood group_class posterior modal
#>      1         10         -34.41           1 9.332e-03 FALSE
#>      2         10         -29.65           1 1.057e-03 FALSE
#>      3         10         -37.14           1 4.329e-01 FALSE
#>      4         10         -30.63           1 1.013e-03 FALSE
#>      5         10         -31.67           1 1.539e-01 FALSE
#>      6         10         -30.04           1 9.911e-06 FALSE
#>      7         10         -33.47           1 1.742e-04 FALSE
#>      8         10         -32.82           1 7.818e-04 FALSE
#>      9         10         -33.93           1 1.409e-03 FALSE
#>     10         10         -32.16           1 3.257e-04 FALSE
#>    ... 38 more rows.  get_results(x, what = "group_posteriors")
#> 
#> -- assignments -----------------------------------------------------
#>  school score_a  score_b profile group_class uncertainty posterior_profile_1
#>       1  0.4453  2.20269       1           2   0.4693255           0.5306745
#>       1  1.5699 -2.30855       2           2   0.0004077           0.0004077
#>       1  0.6884 -0.05664       2           2   0.0169366           0.0169366
#>       1 -0.1776  0.06284       2           2   0.0073025           0.0073025
#>       1  0.7292  0.71023       2           2   0.0813076           0.0813076
#>       1  1.5333 -0.59231       2           2   0.0197391           0.0197391
#>       1  0.5066  0.29852       2           2   0.0274373           0.0274373
#>       1  2.0333  2.64254       1           2   0.0279641           0.9720359
#>       1 -1.4676  2.11280       2           2   0.1178553           0.1178553
#>       1  1.0192  0.91778       2           2   0.1704138           0.1704138
#>  posterior_profile_2
#>              0.46933
#>              0.99959
#>              0.98306
#>              0.99270
#>              0.91869
#>              0.98026
#>              0.97256
#>              0.02796
#>              0.88214
#>              0.82959
#>    ... 230 more rows.  get_results(x, what = "assignments")
#> 
#> -- classification --------------------------------------------------
#>        level class n_modal proportion_modal estimated_n estimated_proportion
#>  individuals     1     123           0.5125      120.27               0.5011
#>  individuals     2     117           0.4875      119.73               0.4989
#>       groups     1      12           0.5000       12.38               0.5160
#>       groups     2      12           0.5000       11.62               0.4840
#>  average_posterior odds_correct_classification
#>             0.9275                       12.73
#>             0.9470                       17.97
#>             0.9818                       50.49
#>             0.9498                       20.18
#> 
#> -- average_posteriors ----------------------------------------------
#>        level assigned_class class n_assigned average_posterior
#>  individuals              1     1        123           0.92747
#>  individuals              1     2        123           0.07253
#>  individuals              2     1        117           0.05296
#>  individuals              2     2        117           0.94704
#>       groups              1     1         12           0.98176
#>       groups              1     2         12           0.01824
#>       groups              2     1         12           0.05017
#>       groups              2     2         12           0.94983
#> 
#> -- classification_errors -------------------------------------------
#>        level true_class assigned_class probability
#>  individuals          1              1     0.94849
#>  individuals          1              2     0.05151
#>  individuals          2              1     0.07451
#>  individuals          2              2     0.92549
#>       groups          1              1     0.95138
#>       groups          1              2     0.04862
#>       groups          2              1     0.01884
#>       groups          2              2     0.98116
#> 
#> -- bch_weights -----------------------------------------------------
#>        level unit assigned_class class   weight
#>  individuals    1              1     1  1.05894
#>  individuals    1              1     2 -0.05894
#>  individuals    2              2     1 -0.08526
#>  individuals    2              2     2  1.08526
#>  individuals    3              2     1 -0.08526
#>  individuals    3              2     2  1.08526
#>  individuals    4              2     1 -0.08526
#>  individuals    4              2     2  1.08526
#>  individuals    5              2     1 -0.08526
#>  individuals    5              2     2  1.08526
#>    ... 518 more rows.  get_results(x, what = "bch_weights")
#> 
#> -- entropy ---------------------------------------------------------
#>        level n_classes n_units entropy_sum relative_entropy
#>  individuals         2     240      40.873           0.7543
#>       groups         2      24       1.901           0.8857
#> 
#> -- residuals -------------------------------------------------------
#>    profile indicator_1 indicator_2     kind observed expected residual
#>  profile_1     score_a     score_b gaussian  0.11599        0  0.11599
#>  profile_2     score_a     score_b gaussian -0.08285        0 -0.08285
#>  effective_n statistic df p_value p_adjusted
#>        120.3    1.2618 NA  0.2070     0.2070
#>        119.7   -0.8971 NA  0.3697     0.3697
#> 
#> -- information_criteria --------------------------------------------
#>  log_likelihood n_parameters  aic  kic bic_groups bic_individual sabic_groups
#>          -784.6           11 1591 1605       1604           1630         1570
#>  sabic_individual caic_groups caic_individual awe_groups awe_individual
#>              1595        1615            1641       1676           1805
#>  icl_groups icl_individual clc_groups clc_individual
#>        1608           1711       1573           1651
#> 
#> -- model -----------------------------------------------------------
#>  n_observations n_informative n_groups n_profiles n_group_classes centering
#>             240           240       24          2               2      none
#>  covariance_structure n_parameters n_parameters_with_measurement log_likelihood
#>                   VVI           11                            11         -784.6
#>   aic bic_groups bic_individual converged iterations boundary small_classes
#>  1591       1604           1630      TRUE         21    FALSE         FALSE
#>  best_start n_best_replicated
#>           1                 4
#> 
#> -- stages ----------------------------------------------------------
#>  stage group_classes fixed log_likelihood parameters
#>  joint             2  <NA>         -784.6         11
#>  parameters_with_measurement converged
#>                           11      TRUE
#> 
#> -- starts ----------------------------------------------------------
#>  start log_likelihood converged iterations error boundary
#>      1         -784.6      TRUE         21  <NA>    FALSE
#>      2         -784.6      TRUE         24  <NA>    FALSE
#>      3         -784.6      TRUE         21  <NA>    FALSE
#>      4         -784.6      TRUE         24  <NA>    FALSE
#> 
#> -- data ------------------------------------------------------------
#>  school score_a  score_b
#>       1  0.4453  2.20269
#>       1  1.5699 -2.30855
#>       1  0.6884 -0.05664
#>       1 -0.1776  0.06284
#>       1  0.7292  0.71023
#>       1  1.5333 -0.59231
#>       1  0.5066  0.29852
#>       1  2.0333  2.64254
#>       1 -1.4676  2.11280
#>       1  1.0192  0.91778
#>    ... 230 more rows.  get_results(x, what = "data")
#> 
#> 19 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
as.data.frame(fit)
#>   profile indicator        mean  variance standard_deviation
#> 1       1   score_a  1.82742675 1.4076126          1.1864285
#> 2       1   score_b  1.96280211 0.9206189          0.9594889
#> 3       2   score_a  0.07696489 0.9562409          0.9778757
#> 4       2   score_b -0.04158577 1.0239945          1.0119261
get_results(fit, what = "profile_probabilities")
#>   group_class profile probability group_class_probability
#> 1           1       1   0.8101391               0.5159839
#> 2           1       2   0.1898609               0.5159839
#> 3           2       1   0.1716925               0.4840161
#> 4           2       2   0.8283075               0.4840161
```
