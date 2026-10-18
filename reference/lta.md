# Fit a latent transition model

Fits latent profiles to repeated observations of the same group,
estimates the probabilities of moving between profiles from one occasion
to the next, and returns a `multilpa_transitions` (or `multilpa_lta`)
object. The measurement model is that of
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
(Gaussian, categorical or mixed indicators) and is invariant across
occasions by default, so a profile has the same meaning at every
occasion. With more than one group class, each class has its own initial
distribution and its own transition matrix, so the group classes are a
mixture over transition patterns.
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
is the cross-sectional counterpart.

## Usage

``` r
lta(
  data,
  vars,
  id,
  n_profiles,
  time,
  n_group_classes = 1L,
  variance_model = c("varying", "equal"),
  n_starts = 10L,
  max_iter = 5000L,
  tol = 1e-10,
  min_variance = 1e-06,
  seed = NULL,
  missing = c("error", "fiml"),
  covariance_model = c("diagonal", "full"),
  categorical = character(),
  min_probability = 1e-10,
  occasions = c("observed", "grid"),
  select_start = c("likelihood", "converged"),
  transitions = c("homogeneous", "occasion"),
  transition_covariates = character(),
  initial_covariates = character(),
  measurement = c("invariant", "occasion"),
  order = 1L,
  model = NULL,
  mover_stayer = FALSE,
  weights = NULL,
  ordinal = character(),
  count = character(),
  count_model = c("poisson", "negative_binomial"),
  count_dispersion = c("varying", "equal"),
  drop_empty = TRUE
)
```

## Arguments

- data:

  A data frame containing the indicators and the group identifier.

- vars:

  Unique character vector of indicator column names. Columns not named
  in `categorical`, `ordinal` or `count` are continuous indicators.

- id:

  Name of the observed group identifier column. Character, factor and
  numeric identifiers are supported; missing identifiers are not.
  Required; `id = NULL` is not accepted.

- n_profiles:

  Integer number of individual profiles, at least two
  (`latents_bad_transition` otherwise).

- time:

  Name of the column giving each observation's occasion within its
  group. Required. Values must be complete and unique within each group.
  They are ordered as [`sort()`](https://rdrr.io/r/base/sort.html)
  orders them: numbers numerically, a factor by its levels, and
  character strings alphabetically. Numbers stored as text whose
  alphabetical order is not their numeric order (`"2"`, `"10"`) raise
  `latents_bad_time` (also `latents_bad_data`); zero-padded numbers and
  ISO dates are accepted.

- n_group_classes:

  Positive integer number of latent group classes (default 1). One gives
  an ordinary latent transition model; more than one fits a mixture of
  transition models over groups.

- variance_model:

  `"varying"` (the default; profile-specific variances or covariance
  matrices) or `"equal"` (shared across profiles).

- n_starts:

  Positive integer number of EM starts (default 10). The first start is
  deterministic: Ward's hierarchical clustering of the standardized
  indicators (on at most 2000 evenly spaced rows), each cluster
  described by its own means and spread, as in mclust's initialization.
  The remaining starts are k-means from random centres. A supplied
  `start` replaces the first start, and the remaining starts are random.
  With `start` and `max_iter = 0`, exactly one start is run whatever
  `n_starts` is.

- max_iter:

  Nonnegative integer maximum number of EM updates per start (default
  5000, higher than the 1000 of
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  to match the tighter `tol`). A fit that stops short warns with
  `latents_unconverged`. `max_iter = 0` updates nothing and evaluates
  the starting values.

- tol:

  Positive relative log-likelihood tolerance (default `1e-10`, tighter
  than the `1e-8` of
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)).
  Convergence requires an absolute change no greater than
  `tol * (1 + abs(previous log likelihood))`. The group-class shares are
  re-estimated from posteriors one EM step behind them; at `1e-8` the
  two differed by up to about `2e-4`, at `1e-10` by about `2e-5`, for
  roughly 1.5 times the iterations.

- min_variance:

  Positive lower bound on each variance, or on each covariance
  eigenvalue under full covariance, in squared input units (default
  `1e-6`); the fit is then constrained maximum likelihood. Estimates at
  the bound are reported and raise a warning.

- seed:

  Optional random seed: any whole number
  [`set.seed()`](https://rdrr.io/r/base/Random.html) accepts, including
  negative ones. When supplied, the caller's random-number state is
  restored on exit.

- missing:

  `"error"` (the default) rejects missing indicators; `"fiml"` maximizes
  the observed-data likelihood under an ignorable missingness mechanism
  (MAR), integrating missing indicators out rather than filling them in.
  This also holds with membership covariates; the covariates themselves
  must be complete (`latents_bad_data`).

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

- occasions:

  How a group's occasions are placed on the transition grid.
  `"observed"` (the default) numbers each group's own observations
  consecutively, so a transition always links two adjacent observations.
  `"grid"` places them on the grid of every position seen in the data,
  so a group that skips a position still makes a transition across the
  gap and contributes no measurement information at it. The two agree
  when every group is observed at every position. Under either, a
  group's first observation is its first occasion: occasion-specific
  transitions and measurement count occasions from each group's own
  start, not from the first position in the data.

- select_start:

  Which start the fit reports. `"likelihood"` (the default) takes the
  highest log likelihood across starts, preferring a converged start
  only among those within a relative `1e-10` of it. `"converged"` takes
  the highest log likelihood among the converged starts when at least
  one converged, and among all starts otherwise. The two differ when an
  unconverged start exceeds a converged one by a negligible amount, as
  happens when a membership logit drifts along a flat ridge of the
  likelihood; only the converged start is a maximum on which standard
  errors can be based.
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  `get_results(fit, "starts")` show every start.

- transitions:

  `"homogeneous"` (the default) uses one transition matrix for every
  occasion; `"occasion"` estimates a separate one for each move
  (occasion 1 to 2, 2 to 3, ...).

- transition_covariates:

  Names of columns whose values shift the transition probabilities
  through a multinomial logit per origin profile: the log odds of moving
  to each other profile rather than staying. A column may change over
  occasions; the value at the destination occasion is used. Factors are
  expanded to indicator columns. With `occasions = "grid"`, a move into
  a position the group skipped has no row to read the covariates from,
  and the fit raises `latents_bad_data`; use `occasions = "observed"`,
  or add the skipped rows with their covariate values (their indicators
  may be missing under `missing = "fiml"`).

- initial_covariates:

  Names of columns that shift the initial profile distribution (a
  multinomial logit, reference: the last profile), read at each group's
  first occasion.

- measurement:

  `"invariant"` (the default) keeps the profiles' means, variances and
  response probabilities equal across occasions; `"occasion"` estimates
  them separately at each occasion.

- order:

  `1` (the default) or `2`. With `2`, the profile at each occasion from
  the third on depends on the profiles at the two previous occasions (a
  second-order Markov chain); the first move stays first order.

- model:

  A covariance structure as
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  takes it (`"EEI"`, `"VVI"`, `"VEI"`, ... ; mclust's codes). Standard
  errors are given for every structure with the extensions (covariate,
  occasion-varying or second-order transitions, occasion-specific
  measurement), and for EEI, VVI, EEE and VVV in the homogeneous model,
  which refuses the other ten with `latents_unsupported_inference`.

- mover_stayer:

  `FALSE` (the default), or `TRUE` to add a class of stayers: groups
  that never change profile (identity transitions) and have their own
  initial profile distribution (Goodman's mover-stayer model).
  `n_group_classes` then counts the mover classes; the stayer class is
  reported as `"stayers"`.

- weights:

  `NULL`, or the name of a column of `data` holding one sampling weight
  per `id` (constant over its occasions). As in
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md):
  pseudo maximum likelihood with the weights scaled to sum to the number
  of `id` units, sandwich standard errors by default, and integer
  weights equal to repeating each sequence. A weighted fit uses the
  general transition engine and returns a `multilpa_lta`, whose
  model-implied `"transitions"` and `"initial"` tables average the
  groups with their weights.

- ordinal, count, count_model, count_dispersion:

  Ordinal and count indicators, as in
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md);
  with `measurement = "occasion"` their parameters are estimated per
  occasion. Read them with `get_results(fit, "ordinal")` and
  `get_results(fit, "count_means")`. A fit with them uses the general
  transition engine.

- drop_empty:

  `TRUE` (the default) or `FALSE`: how the reported statistics treat
  rows whose every indicator is missing, and groups consisting only of
  such rows. These contribute nothing to the likelihood, so the
  estimates are the same either way, and their posterior is the prior.
  With `TRUE` they are excluded, as in Mplus, from the individual and
  group sample sizes of the information criteria, from relative entropy,
  AvePP and the entropy term of ICL, AWE and CLC, from the
  classification tables and diagnostic plots, and from the effective
  class counts. With `FALSE` every row and group is counted in those
  statistics (the results of latents 0.9.37 and earlier), except in the
  individual sample size of a
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  fit, which excludes such rows either way.
  [`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)
  and
  [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
  do not compare fits made with different `drop_empty`. The setting is
  stored on the fit and passed on by enumerations and refits; it has no
  effect on data without such rows.

## Value

A `multilpa_transitions` object containing `means`, `variances`,
optional `covariances` and `response_probabilities`,
`initial_probabilities` (group classes by profiles),
`transition_probabilities` (profiles by profiles by group classes, rows
indexing the profile moved from), `group_probabilities`, posterior
matrices, classifications, log likelihood, information criteria, restart
diagnostics and convergence history.
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
returns the tidy transition table and the other tables,
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
Wald or bootstrap standard errors, and
[`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
compares nested fits. Global optimality is not guaranteed.

With `transitions = "occasion"`, `transition_covariates`,
`initial_covariates`, `measurement = "occasion"`, `order = 2`,
`mover_stayer = TRUE`, `weights`, `ordinal` or `count`, the result is a
`multilpa_lta` object, read with
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md):
transition and initial logit coefficients with Wald standard errors
(observed, robust or OPG, from analytic scores), model-implied
transition probabilities per occasion, and profiles per occasion.
Missing indicators (`missing = "fiml"`) and covariance structures are
supported; diagonal fits are finished by a quasi-Newton search on the
exact likelihood, others by EM. Agreement with Mplus (User's Guide
examples 8.13, 8.14; occasion-specific thresholds), depmixS4 and LMest
is recorded in `equivalence/lta-extensions/` of the source repository.

## Details

By default transitions are first order and homogeneous over occasions:
the probability of moving from one profile to another does not depend on
the occasion or on earlier profiles. Measurement parameters are shared
across occasions and across group classes. A profile that no group
occupies before its final occasion leaves its transition row without
information; the row is then uniform by construction, not estimated.
Such a row raises a warning when the model is fitted, is flagged by the
`estimated` column of the transition table and is listed by
`get_results(fit, "transitions", estimated = FALSE)`. Profile labels are
arbitrary and are not comparable across fits without alignment.

`max_iter = 0` updates nothing: the starting values are evaluated and
returned with the expectation they imply, and the fit reports
`converged = FALSE` after zero iterations. Every reported quantity,
including `expected_count` and the implied profile prevalence, is then
read off that single expectation.

## References

Collins, L. M., & Lanza, S. T. (2010). Latent class and latent
transition analysis. Wiley.

Vermunt, J. K. (2003). Multilevel latent class models. Sociological
Methodology, 33, 213–239. doi:10.1111/j.0081-1750.2003.t01-1-00131.x.

## See also

[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
for the fitted transition probabilities and for the assignments in
order, and
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
for the cross-sectional model this one shares its measurement parameters
with.

## Examples

``` r
# Students take their courses in their own order, which `sequence` records,
# so the same enrolments that fit a cross-sectional model fit a transition
# one. Engagement mostly persists from one course to the next.
fit <- lta(
  course_engagement,
  vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
  id = "student", n_profiles = 2, time = "sequence", n_starts = 2, seed = 1
)
get_results(fit, what = "transitions")
#>   group_class from to probability expected_count stable estimated
#> 1           1    1  1   0.8613981      670.56011   TRUE      TRUE
#> 2           1    1  2   0.1386019      107.89512  FALSE      TRUE
#> 3           1    2  1   0.1850265       99.46008  FALSE      TRUE
#> 4           1    2  2   0.8149735      438.08469   TRUE      TRUE
#>   group_class_probability
#> 1                       1
#> 2                       1
#> 3                       1
#> 4                       1
summary(fit)
#> Latent transition model: 2 profiles and 1 group class
#> Observations: 1422; groups: 106; up to 15 occasions (unbalanced, observed grid)
#> Parameters: 23; converged: TRUE
#> Log likelihood: -8335.287419; AIC: 16716.575
#> BIC (groups): 16777.834; BIC (individuals): 16837.551
#> Best likelihood replicated in 2/2 starts (absolute tolerance 0.00834).
#> 
#> -- profiles --------------------------------------------------------
#>  profile  indicator    mean variance standard_deviation mean_standard_error
#>        1     browse  0.5405   0.5929             0.7700             0.02688
#>        1   lectures  0.4266   0.8376             0.9152             0.03197
#>        1 forum_read  0.6173   0.4860             0.6972             0.02454
#>        1 forum_post  0.5297   0.6877             0.8293             0.02902
#>        1 attendance  0.6436   0.3916             0.6258             0.02207
#>        2     browse -0.7720   0.5139             0.7168             0.03027
#>        2   lectures -0.6093   0.5456             0.7387             0.03095
#>        2 forum_read -0.8814   0.3582             0.5985             0.02531
#>        2 forum_post -0.7565   0.4181             0.6466             0.02731
#>        2 attendance -0.9189   0.3780             0.6148             0.02598
#>  variance_standard_error
#>                  0.02930
#>                  0.04135
#>                  0.02453
#>                  0.03411
#>                  0.01971
#>                  0.03088
#>                  0.03266
#>                  0.02152
#>                  0.02552
#>                  0.02257
#> 
#> -- responses -------------------------------------------------------
#>    (no rows)
#> 
#> -- covariances -----------------------------------------------------
#>  profile  indicator indicator_2 covariance
#>        1     browse      browse     0.5929
#>        1   lectures      browse     0.0000
#>        1 forum_read      browse     0.0000
#>        1 forum_post      browse     0.0000
#>        1 attendance      browse     0.0000
#>        1     browse    lectures     0.0000
#>        1   lectures    lectures     0.8376
#>        1 forum_read    lectures     0.0000
#>        1 forum_post    lectures     0.0000
#>        1 attendance    lectures     0.0000
#>    ... 40 more rows.  get_results(x, what = "covariances")
#> 
#> -- transitions -----------------------------------------------------
#>  group_class from to probability expected_count stable estimated
#>            1    1  1      0.8614         670.56   TRUE      TRUE
#>            1    1  2      0.1386         107.90  FALSE      TRUE
#>            1    2  1      0.1850          99.46  FALSE      TRUE
#>            1    2  2      0.8150         438.08   TRUE      TRUE
#>  group_class_probability
#>                        1
#>                        1
#>                        1
#>                        1
#> 
#> -- initial ---------------------------------------------------------
#>  group_class profile probability prevalence group_class_probability
#>            1       1      0.6258     0.5882                       1
#>            1       2      0.3742     0.4118                       1
#> 
#> -- counts ----------------------------------------------------------
#>        level class effective_count effective_proportion
#>  individuals     1           836.4               0.5882
#>  individuals     2           585.6               0.4118
#>       groups     1           106.0               1.0000
#> 
#> -- posteriors ------------------------------------------------------
#>  row group profile posterior modal
#>    1     1       1 1.000e+00  TRUE
#>    2     1       1 9.971e-01  TRUE
#>    3     1       1 1.387e-05 FALSE
#>    4     1       1 8.981e-07 FALSE
#>    5     1       1 1.263e-06 FALSE
#>    6     1       1 4.927e-06 FALSE
#>    7     1       1 4.797e-06 FALSE
#>    8     1       1 2.279e-06 FALSE
#>    9     1       1 4.709e-07 FALSE
#>   10     1       1 2.102e-06 FALSE
#>    ... 2834 more rows.  get_results(x, what = "posteriors")
#> 
#> -- group_posteriors ------------------------------------------------
#>  group group_size log_likelihood group_class posterior modal
#>      1         13         -59.47           1         1  TRUE
#>      2         13         -71.93           1         1  TRUE
#>      3         14         -75.38           1         1  TRUE
#>      4         13         -66.30           1         1  TRUE
#>      5         14         -79.40           1         1  TRUE
#>      6         12         -80.17           1         1  TRUE
#>      7         15         -88.32           1         1  TRUE
#>      8         14         -81.68           1         1  TRUE
#>      9         14         -82.52           1         1  TRUE
#>     10         15         -82.01           1         1  TRUE
#>    ... 96 more rows.  get_results(x, what = "group_posteriors")
#> 
#> -- assignments -----------------------------------------------------
#>  student sequence browse lectures forum_read forum_post attendance profile
#>        1        1   0.73    -0.08       0.30       0.67       0.72       1
#>        1        2   0.68     0.53      -0.02      -0.55       0.70       1
#>        1        3  -1.51    -0.51      -0.91      -0.80      -0.97       2
#>        1        4   0.64    -1.02      -1.62      -1.20      -1.26       2
#>        1        5  -0.94    -0.37      -0.71      -0.28      -1.52       2
#>        1        6  -0.36    -0.46      -0.77      -0.43      -1.34       2
#>        1        7   0.38    -1.85      -0.84      -1.09      -1.12       2
#>        1        8  -1.44     0.66      -0.67      -1.67      -1.00       2
#>        1        9  -0.76    -0.43      -1.04      -0.69      -1.37       2
#>        1       10   0.26    -0.95      -0.51      -1.89      -1.45       2
#>  group_class uncertainty posterior_profile_1 posterior_profile_2
#>            1   1.174e-05           1.000e+00           1.174e-05
#>            1   2.865e-03           9.971e-01           2.865e-03
#>            1   1.387e-05           1.387e-05           1.000e+00
#>            1   8.981e-07           8.981e-07           1.000e+00
#>            1   1.263e-06           1.263e-06           1.000e+00
#>            1   4.927e-06           4.927e-06           1.000e+00
#>            1   4.797e-06           4.797e-06           1.000e+00
#>            1   2.279e-06           2.279e-06           1.000e+00
#>            1   4.709e-07           4.709e-07           1.000e+00
#>            1   2.102e-06           2.102e-06           1.000e+00
#>    ... 1412 more rows.  get_results(x, what = "assignments")
#> 
#> -- classification --------------------------------------------------
#>        level class n_modal proportion_modal estimated_n estimated_proportion
#>  individuals     1     837           0.5886       836.4               0.5882
#>  individuals     2     585           0.4114       585.6               0.4118
#>       groups     1     106           1.0000       106.0               1.0000
#>  average_posterior odds_correct_classification
#>             0.9915                       81.98
#>             0.9890                      128.23
#>             1.0000                          NA
#> 
#> -- average_posteriors ----------------------------------------------
#>        level assigned_class class n_assigned average_posterior
#>  individuals              1     1        837           0.99153
#>  individuals              1     2        837           0.00847
#>  individuals              2     1        585           0.01101
#>  individuals              2     2        585           0.98899
#>       groups              1     1        106           1.00000
#> 
#> -- classification_errors -------------------------------------------
#>        level true_class assigned_class probability
#>  individuals          1              1    0.992296
#>  individuals          1              2    0.007704
#>  individuals          2              1    0.012105
#>  individuals          2              2    0.987895
#>       groups          1              1    1.000000
#> 
#> -- bch_weights -----------------------------------------------------
#>        level unit assigned_class class   weight
#>  individuals    1              1     1  1.00786
#>  individuals    1              1     2 -0.00786
#>  individuals    2              1     1  1.00786
#>  individuals    2              1     2 -0.00786
#>  individuals    3              2     1 -0.01235
#>  individuals    3              2     2  1.01235
#>  individuals    4              2     1 -0.01235
#>  individuals    4              2     2  1.01235
#>  individuals    5              2     1 -0.01235
#>  individuals    5              2     2  1.01235
#>    ... 2940 more rows.  get_results(x, what = "bch_weights")
#> 
#> -- entropy ---------------------------------------------------------
#>        level n_classes n_units entropy_sum relative_entropy
#>  individuals         2    1422        38.3           0.9611
#>       groups         1     106         0.0               NA
#> 
#> -- residuals -------------------------------------------------------
#>    profile indicator_1 indicator_2     kind observed expected residual
#>  profile_1    lectures  attendance gaussian  0.23622        0  0.23622
#>  profile_2  forum_read  attendance gaussian  0.23302        0  0.23302
#>  profile_1  forum_read  attendance gaussian  0.21662        0  0.21662
#>  profile_1  forum_post  attendance gaussian  0.20818        0  0.20818
#>  profile_2      browse  attendance gaussian  0.19680        0  0.19680
#>  profile_1      browse  attendance gaussian  0.16968        0  0.16968
#>  profile_2    lectures  attendance gaussian  0.11963        0  0.11963
#>  profile_2  forum_post  attendance gaussian  0.11548        0  0.11548
#>  profile_2      browse  forum_read gaussian  0.10229        0  0.10229
#>  profile_2  forum_read  forum_post gaussian  0.06942        0  0.06942
#>  effective_n statistic df   p_value p_adjusted
#>        836.4     6.951 NA 3.638e-12  3.638e-12
#>        585.6     5.730 NA 1.005e-08  1.005e-08
#>        836.4     6.354 NA 2.097e-10  2.097e-10
#>        836.4     6.099 NA 1.068e-09  1.068e-09
#>        585.6     4.813 NA 1.487e-06  1.487e-06
#>        836.4     4.946 NA 7.574e-07  7.574e-07
#>        585.6     2.901 NA 3.714e-03  3.714e-03
#>        585.6     2.800 NA 5.112e-03  5.112e-03
#>        585.6     2.478 NA 1.322e-02  1.322e-02
#>        585.6     1.678 NA 9.328e-02  9.328e-02
#>    ... 10 more rows.  get_results(x, what = "residuals")
#> 
#> -- information_criteria --------------------------------------------
#>  log_likelihood n_parameters   aic   kic bic_groups bic_individual sabic_groups
#>           -8335           23 16717 16743      16778          16838        16705
#>  sabic_individual caic_groups caic_individual awe_groups awe_individual
#>             16764       16801           16861      16954          17150
#>  icl_groups icl_individual clc_groups clc_individual
#>       16778          16914      16671          16747
#> 
#> -- model -----------------------------------------------------------
#>  n_observations n_informative n_groups n_profiles n_group_classes centering
#>            1422          1422      106          2               1      none
#>  scaling covariance_structure n_parameters n_parameters_with_measurement
#>     none                 <NA>           23                            23
#>  log_likelihood   aic bic_groups bic_individual converged iterations boundary
#>           -8335 16717      16778          16838      TRUE          6    FALSE
#>  small_classes best_start n_best_replicated weights
#>          FALSE          1                 2    <NA>
#> 
#> -- sequences -------------------------------------------------------
#>  group group_class time profile
#>      1           1    1       1
#>      1           1    2       1
#>      1           1    3       2
#>      1           1    4       2
#>      1           1    5       2
#>      1           1    6       2
#>      1           1    7       2
#>      1           1    8       2
#>      1           1    9       2
#>      1           1   10       2
#>    ... 1412 more rows.  get_results(x, what = "sequences")
#> 
#> -- sequence_summary ------------------------------------------------
#>  group_class groups observations mean_length median_length shortest longest
#>            1    106         1422       13.42            14       10      15
#>  complete gaps
#>        24    0
#> 
#> -- sequence_lengths ------------------------------------------------
#>  group group_class occasions observations complete
#>      1           1        13           13    FALSE
#>      2           1        13           13    FALSE
#>      3           1        14           14    FALSE
#>      4           1        13           13    FALSE
#>      5           1        14           14    FALSE
#>      6           1        12           12    FALSE
#>      7           1        15           15     TRUE
#>      8           1        14           14    FALSE
#>      9           1        14           14    FALSE
#>     10           1        15           15     TRUE
#>    ... 96 more rows.  get_results(x, what = "sequence_lengths")
#> 
#> -- starts ----------------------------------------------------------
#>  start log_likelihood converged iterations error
#>      1          -8335      TRUE          6  <NA>
#>      2          -8335      TRUE          7  <NA>
#> 
#> -- data ------------------------------------------------------------
#>  student sequence browse lectures forum_read forum_post attendance
#>        1        1   0.73    -0.08       0.30       0.67       0.72
#>        1        2   0.68     0.53      -0.02      -0.55       0.70
#>        1        3  -1.51    -0.51      -0.91      -0.80      -0.97
#>        1        4   0.64    -1.02      -1.62      -1.20      -1.26
#>        1        5  -0.94    -0.37      -0.71      -0.28      -1.52
#>        1        6  -0.36    -0.46      -0.77      -0.43      -1.34
#>        1        7   0.38    -1.85      -0.84      -1.09      -1.12
#>        1        8  -1.44     0.66      -0.67      -1.67      -1.00
#>        1        9  -0.76    -0.43      -1.04      -0.69      -1.37
#>        1       10   0.26    -0.95      -0.51      -1.89      -1.45
#>    ... 1412 more rows.  get_results(x, what = "data")
#> 
#> 22 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
```
