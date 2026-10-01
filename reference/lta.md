# Latent transition analysis

Fits latent profiles to repeated observations of the same group and
estimates the probabilities of moving between them from one occasion to
the next. The measurement model is that of
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md) –
Gaussian, categorical or mixed indicators, invariant across occasions –
so profiles keep the same meaning at every occasion and a change of
profile is a change of state rather than a change of definition. With
more than one group class, each class carries its own initial
distribution and its own transition matrix, which is how groups that
differ in their dynamics are separated from groups that differ only in
where they start.

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
  max_iter = 1000L,
  tol = 1e-08,
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
  model = NULL
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

  Integer number of individual profiles, at least two.

- time:

  Name of the column giving each observation's occasion within its
  group. Required: this model is defined by the ordering. Values must be
  complete and unique within each group.

- n_group_classes:

  Positive integer number of latent group classes. One gives an ordinary
  single-level latent transition model; more than one fits a mixture of
  transition models over groups.

- variance_model:

  Either `"varying"` (profile-specific variances or covariance matrices)
  or `"equal"` (shared across profiles).

- n_starts:

  Positive integer number of EM starts. The first is deterministic:
  Ward's hierarchical clustering of the standardized indicators (on at
  most 2000 evenly spaced rows), each cluster described by its own means
  and spread, which is the start mclust's name refers to and finds
  solutions whose profiles differ in shape. The rest are k-means from
  random centres. When `start` is supplied, it supplies the first start
  instead; remaining starts are random initializations. It is ignored
  when `max_iter = 0` and `start` is supplied: that call evaluates the
  supplied parameters and nothing else, so exactly one start is run
  whatever `n_starts` says.

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
  likelihood fitting. This holds with membership covariates too; the
  covariates themselves must be complete (`latents_bad_data`).

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

- occasions:

  How a group's occasions are placed on the transition grid.
  `"observed"` numbers each group's own observations consecutively, so a
  transition always links two adjacent observations. `"grid"` places
  them on the grid of every position seen in the data, so a group that
  skips a position still consumes a transition across the gap and
  contributes no measurement information at it. The two agree whenever
  every group is observed at every position.

- select_start:

  Which start the fit reports. `"likelihood"`, the default, takes the
  highest log likelihood across starts, preferring a converged start
  only among those tied with it to within a relative `1e-10`.
  `"converged"` takes the highest log likelihood among the converged
  starts whenever at least one converged, and falls back to every start
  when none did. Use it when an unconverged start edges out a converged
  one by a negligible amount, which happens when a membership logit
  creeps along a flat ridge of the likelihood: the converged start is
  then a maximum that the standard errors can work with, and the other
  is not. [`summary()`](https://rdrr.io/r/base/summary.html) and
  `get_results(fit, "starts")` show every start either way.

- transitions:

  `"homogeneous"` (the default) uses one transition matrix for every
  occasion; `"occasion"` estimates a separate one for each move
  (occasion 1 to 2, 2 to 3, ...).

- transition_covariates:

  Names of columns whose values shift the transition probabilities
  through a multinomial logit per origin profile: the log odds of moving
  to each other profile rather than staying. A column may change over
  occasions; the value at the destination occasion is used. Factors are
  expanded to indicator columns.

- initial_covariates:

  Names of columns that shift the initial profile distribution (a
  multinomial logit, reference: the last profile), read at each group's
  first occasion.

- measurement:

  `"invariant"` (the default) keeps the profiles' means, variances and
  response probabilities equal across occasions; `"occasion"` estimates
  them separately at each occasion, so a profile is defined by its
  position in the sequence of transitions rather than by one fixed
  measurement model.

- order:

  `1` (the default) or `2`. With `2`, the profile at each occasion from
  the third on depends on the profiles at the two previous occasions (a
  second-order Markov chain); the first move stays first order.

- model:

  A covariance structure as
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  takes it (`"EEI"`, `"VVI"`, `"VEI"`, ... ; mclust's codes). Available
  with homogeneous first-order transitions and invariant measurement;
  standard errors for structures other than EEI, VVI, EEE and VVV are
  refused (`latents_unsupported_inference`).

## Value

A `multilpa_transitions` object containing `means`, `variances`,
optional `covariances` and `response_probabilities`,
`initial_probabilities` (group classes by profiles),
`transition_probabilities` (profiles by profiles by group classes, rows
indexing the profile moved from), `group_probabilities`, posterior
matrices, classifications, log likelihood, information criteria, restart
diagnostics and convergence history. Use
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
for the tidy transition table and for the other tables. No standard
errors, likelihood-ratio tests or guarantees of global optimality are
given for this model family.

With `transitions = "occasion"`, `transition_covariates`,
`initial_covariates`, `measurement = "occasion"` or `order = 2` the
result is a `multilpa_lta` object instead, read with
[`get_results.multilpa_lta()`](https://pak.dynasite.org/latents/reference/get_results.multilpa_lta.md):
transition and initial logit coefficients with Wald standard errors
(observed, robust or OPG, from analytic scores), model-implied
transition probabilities per occasion, profiles per occasion. These need
complete indicators and diagonal covariances. EM is finished by a
quasi-Newton search on the exact likelihood. Agreement with Mplus
(User's Guide examples 8.13, 8.14; occasion-specific thresholds),
depmixS4 and LMest is recorded in `equivalence/lta-extensions/` of the
source repository.

## Details

This is latent transition analysis (LTA); with more than one group class
it is a mixture over transition patterns, so the classes are
trajectories rather than states.
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
is the cross-sectional counterpart.

By default transitions are first order and homogeneous over occasions:
the probability of moving from one profile to another does not depend on
the occasion or on earlier profiles. Measurement parameters are shared
across occasions and across group classes. A profile that no group
occupies before its final occasion leaves its transition row without
information; the row is then uniform by construction rather than
estimated. Such a row is warned about when the model is fitted, is
flagged by the `estimated` column of the transition table and is listed
by `get_results(fit, "transitions", estimated = FALSE)`. Profile labels
are arbitrary and are not comparable across fits without alignment.

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
#> 1           1    1  1   0.8614039      670.56162   TRUE      TRUE
#> 2           1    1  2   0.1385961      107.89439  FALSE      TRUE
#> 3           1    2  1   0.1850203       99.45941  FALSE      TRUE
#> 4           1    2  2   0.8149797      438.08458   TRUE      TRUE
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
#>            1    1  2      0.1386         107.89  FALSE      TRUE
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
#>    4     1       1 8.982e-07 FALSE
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
#>            1   2.864e-03           9.971e-01           2.864e-03
#>            1   1.387e-05           1.387e-05           1.000e+00
#>            1   8.982e-07           8.982e-07           1.000e+00
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
#>             0.9890                      128.22
#>             1.0000                          NA
#> 
#> -- average_posteriors ----------------------------------------------
#>        level assigned_class class n_assigned average_posterior
#>  individuals              1     1        837          0.991531
#>  individuals              1     2        837          0.008469
#>  individuals              2     1        585          0.011015
#>  individuals              2     2        585          0.988985
#>       groups              1     1        106          1.000000
#> 
#> -- classification_errors -------------------------------------------
#>        level true_class assigned_class probability
#>  individuals          1              1    0.992295
#>  individuals          1              2    0.007705
#>  individuals          2              1    0.012104
#>  individuals          2              2    0.987896
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
#>  profile_1    lectures  attendance gaussian  0.23623        0  0.23623
#>  profile_2  forum_read  attendance gaussian  0.23302        0  0.23302
#>  profile_1  forum_read  attendance gaussian  0.21662        0  0.21662
#>  profile_1  forum_post  attendance gaussian  0.20818        0  0.20818
#>  profile_2      browse  attendance gaussian  0.19679        0  0.19679
#>  profile_1      browse  attendance gaussian  0.16968        0  0.16968
#>  profile_2    lectures  attendance gaussian  0.11963        0  0.11963
#>  profile_2  forum_post  attendance gaussian  0.11548        0  0.11548
#>  profile_2      browse  forum_read gaussian  0.10229        0  0.10229
#>  profile_2  forum_read  forum_post gaussian  0.06942        0  0.06942
#>  effective_n statistic df   p_value p_adjusted
#>        836.4     6.951 NA 3.637e-12  3.637e-12
#>        585.6     5.730 NA 1.005e-08  1.005e-08
#>        836.4     6.354 NA 2.096e-10  2.096e-10
#>        836.4     6.099 NA 1.068e-09  1.068e-09
#>        585.6     4.813 NA 1.487e-06  1.487e-06
#>        836.4     4.946 NA 7.573e-07  7.573e-07
#>        585.6     2.901 NA 3.714e-03  3.714e-03
#>        585.6     2.800 NA 5.113e-03  5.113e-03
#>        585.6     2.478 NA 1.322e-02  1.322e-02
#>        585.6     1.678 NA 9.329e-02  9.329e-02
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
#>  covariance_structure n_parameters n_parameters_with_measurement log_likelihood
#>                  <NA>           23                            23          -8335
#>    aic bic_groups bic_individual converged iterations boundary small_classes
#>  16717      16778          16838      TRUE          6    FALSE         FALSE
#>  best_start n_best_replicated
#>           2                 2
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
#>      1          -8335      TRUE          4  <NA>
#>      2          -8335      TRUE          6  <NA>
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
