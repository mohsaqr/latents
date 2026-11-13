# Fit a latent transition model

Latent transition analysis describes stability and change in latent
states across repeated occasions. The states (profiles) are patterns in
continuous or discrete indicators; transition probabilities give the
chance of staying in a state or moving to another between successive
occasions. `lta()` estimates the measurement and transition parameters
jointly.

By default the chain is first order and time-homogeneous: the state at
an occasion depends on the state before it, with the same transition
matrix at every occasion. The measurement parameters are shared across
occasions and across group classes, so a state means the same thing
throughout and group classes differ only in their initial distribution
and transitions. Options add covariate-dependent transitions,
occasion-specific transitions or measurement, second-order dependence, a
mover-stayer class, sampling weights and ordinal or count indicators;
any of these changes the result class (see *Value*).

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
  min_variance = NULL,
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
  drop_empty = TRUE,
  search = c("kmeans", "emEM"),
  n_screen = NULL,
  acceleration = c("squarem", "none"),
  cores = 1L
)
```

## Arguments

- data:

  A data frame with the indicators and the group identifier.

- vars:

  Names of the indicator columns. Those not named in `categorical`,
  `ordinal` or `count` are continuous.

- id:

  Name of the group identifier column (character, factor or numeric,
  without missing values). Required.

- n_profiles:

  Number of profiles, at least two.

- time:

  Column giving each row's occasion within its group. Required; values
  must be unique within a group and sort in their intended order
  (numbers numerically, factors by level, text alphabetically, so
  `"T2"`, `"T10"` is refused; use a factor or zero-pad).

- n_group_classes:

  Number of group classes (default 1). More than one fits a mixture of
  transition models over the units.

- variance_model:

  `"varying"` (default): each profile has its own variances or
  covariance matrix; `"equal"`: one shared by all profiles.

- n_starts:

  Number of EM starts (default 10). The first is deterministic (Ward's
  clustering), the rest random; the reported start need not have
  converged, so check `converged`. See *Starts and convergence* under
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).

- max_iter:

  Maximum EM updates per start (default 5000). The fit warns when its
  reported start stopped short. `max_iter = 0` evaluates the starting
  values without updating them and reports `converged = FALSE`.

- tol:

  Convergence tolerance (default `1e-10`): EM stops when a step changes
  the log likelihood by no more than `tol` times the number of observed
  indicator values.

- min_variance:

  Lower bound on every variance (on every covariance eigenvalue under
  full covariance). `NULL` (default) uses `1e-6` times the smallest
  sample variance of the continuous indicators, on the scale fitted; a
  number is an absolute bound in squared indicator units. An estimate at
  the bound is a boundary solution and warns.

- seed:

  Optional integer seed. The caller's random-number state is restored on
  exit.

- missing:

  `"error"` (default) rejects missing indicator values; `"fiml"`
  maximizes the observed-data likelihood under MAR. Covariates must be
  complete either way. See *Occasions without data*.

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

- occasions:

  How a group's occasions are placed on the transition grid.
  `"observed"` (default) numbers each group's own rows consecutively, so
  a transition always links two adjacent rows. `"grid"` places them on
  the grid of every position seen in the data, so a group that skips a
  position still makes a transition across the gap and contributes no
  measurement at it. The two agree when every group is observed at every
  position. Under either, a group's first row is its first occasion:
  occasion-specific transitions and measurement count occasions from
  each group's own start, not from the first position in the data.

- select_start:

  Which start is reported. `"likelihood"` (default): the highest log
  likelihood, preferring a converged start among those tied within
  tolerance. `"converged"`: the best converged start when there is one.
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  `get_results(fit, "starts")` report every start.

- transitions:

  `"homogeneous"` (default): one transition matrix for every occasion.
  `"occasion"`: a separate matrix for each move (occasion 1 to 2, 2 to
  3, and so on). See *Occasions without data*.

- transition_covariates:

  Columns whose values shift the transition probabilities through a
  multinomial logit per origin profile (the log odds of moving to each
  other profile rather than staying). A column may change over
  occasions; the value at the destination occasion is used. Factors are
  expanded to indicator columns after dropping unused levels. With
  `occasions = "grid"`, a skipped position has no row to read the
  covariates from and is refused: use `"observed"`, or add the skipped
  rows with their covariate values and missing indicators under FIML.

- initial_covariates:

  Columns that shift the initial profile distribution (a multinomial
  logit, last profile as reference), read at each group's first
  occasion. Factors are expanded as for `transition_covariates`.

- measurement:

  `"invariant"` (default) keeps the profiles' means, variances and
  response probabilities equal across occasions; `"occasion"` estimates
  them separately at each occasion.

- order:

  `1` (default) or `2`. With `2`, the profile at each occasion from the
  third on depends on the profiles at the two previous occasions; the
  first transition stays first order.

- model:

  A covariance structure as
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  takes it (mclust's codes, such as `"VVI"`). A `multilpa_transitions`
  fit gives standard errors for EEI, VVI, EEE and VVV only; a
  `multilpa_lta` fit (see *Value*) gives them for every structure.

- mover_stayer:

  `TRUE` adds a class of stayers: units that never change profile
  (identity transitions) with their own initial distribution, the
  mover-stayer model. `n_group_classes` then counts the mover classes;
  the stayer class is reported as `"stayers"`.

- weights:

  Name of a column holding one sampling weight per `id` (constant over
  its occasions): pseudo maximum likelihood with sandwich standard
  errors, as under
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).
  The model-implied `"transitions"` and `"initial"` tables then average
  the units with their weights.

- ordinal, count, count_model, count_dispersion:

  Ordinal and count indicators, as under
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md);
  with `measurement = "occasion"` their parameters are estimated per
  occasion. Read with `get_results(fit, "ordinal")` and
  `get_results(fit, "count_means")`.

- drop_empty:

  `TRUE` (default) leaves rows with no observed indicator, and groups
  made only of such rows, out of the sample sizes, entropies and
  classification tables, as Mplus does, keeping their rows with `NA` in
  the assignment tables. `FALSE` counts them in the entropies, the
  tables and the number of groups, so the criteria over groups and those
  using the entropy change; the number of rows in the individual-level
  criteria stays the informative rows for profile and class models. For
  transition models, `FALSE` also counts all positive-weight rows in the
  individual-level BIC. Estimation is the same either way; fits with
  different settings are not compared.

- search:

  `"kmeans"` (default) or `"emEM"`, as under
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).
  Unlike
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  `lta()` does not pick the default by covariance structure.

- n_screen:

  For `search = "emEM"` only, the size of the screened pool (default
  five times `n_starts`, at least `n_starts`).

- acceleration:

  `"squarem"` (default) accelerates EM with SQUAREM, as under
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md);
  `"none"` runs plain EM.

- cores:

  Number of processor cores for the starts (default 1); with `seed` the
  result is identical to one core. To fit several models, give `cores`
  to
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  instead.

## Value

A `multilpa_transitions` object with `means`, `variances` (and
`covariances` or `response_probabilities`), `initial_probabilities`
(group classes by profiles), `transition_probabilities` (profiles by
profiles by group classes, rows indexing the profile moved from),
`group_probabilities`, the posterior matrices and classifications,
`log_likelihood`, the information criteria, the starts table and the
convergence history.
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
returns its tables,
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
Wald or bootstrap standard errors and
[`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
compares nested fits.

With `transitions = "occasion"`, `transition_covariates`,
`initial_covariates`, `measurement = "occasion"`, `order = 2`,
`mover_stayer = TRUE`, `weights`, `ordinal` or `count`, the result is a
`multilpa_lta` object: its tables are transition and initial logit
coefficients with Wald standard errors, model-implied transition
probabilities per occasion (with an `estimated` column), profiles per
occasion and the assignments.

Profile labels are arbitrary and not comparable across fits without
alignment.

## Occasions without data

A transition the data never observe is not estimated. A profile that no
unit of a group class occupies before the unit's last occasion has no
move out of it: its row of the transition matrix is uniform by
construction, flagged `estimated = FALSE`, not counted in
`n_parameters`, AIC or BIC, and reported with a warning. On a
`multilpa_transitions` fit,
`get_results(fit, "transitions", estimated = FALSE)` lists such rows.
Nonzero origins with fewer than one expected outgoing move trigger a
weak- support warning. Transition tables show `origin_expected_count`
and `supported`; Wald inference is withheld. These nonzero rows remain
formally estimated and counted in the criteria: this is a support
diagnostic, not a declaration of structural non-identification.
Homogeneous rows pool support across their applicable moves; occasion-
specific rows use that occasion's count. The averaged second-order table
reports the total count, but `supported` is false if any separately
fitted occasion has fewer than one expected outgoing move.

Under FIML, three more things are unidentified and flagged the same way:
under `transitions = "occasion"`, a move into an occasion at which no
row observes an indicator (after the last observed occasion it does not
enter the likelihood; before it, only through its product with the move
out); under `order = 2`, the second-order logits when nothing is
observed from occasion 3 on; under `measurement = "occasion"`, the
measurement block of an empty occasion.

For second-order occasion-specific chains, an empty interior occasion
does not use the first-order product rule: remaining logits stay counted
and the parameter count is an upper bound. A second-order fit's
first-order table describes only the move into occasion 2. The
second-order table averages over groups and applicable occasions, with
move-frequency weights.

## Short panels and bootstrap support

Without initial or transition covariates, two occasions cannot identify
multiple group classes or a separate stayer class. Three occasions
cannot identify multiple freely moving group classes with
occasion-specific transitions, or multiple classes in a second-order
chain. These specifications are refused; convergence cannot resolve an
identification ridge. Specifications whose chain coordinates outnumber
the latent path probabilities are refused as well.

Bootstrap resamples whose occasions or factor levels change the
coefficient set are recorded as failed resamples. Intervals use retained
resamples only, with a warning and recorded reasons, and are conditional
on those resamples supporting the requested coefficients. With
`drop_empty = FALSE`, both engines count every positive-weight row in
the individual BIC; the default counts rows that observe an indicator.
Initial or transition probabilities at the effective numerical bound
also set `boundary = TRUE` in both engines; bootstrap likelihood-ratio
tests refuse such inputs, as they already do for general transition
fits.

## See also

[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
for the transition table and the assignments;
[`plot_transitions()`](https://pak.dynasite.org/latents/reference/plot_transitions.md)
and
[`plot_sequences()`](https://pak.dynasite.org/latents/reference/plot_sequences.md);
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
for the cross-sectional model whose measurement this one shares;
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
to compare numbers of profiles.

## Examples

``` r
# Model engagement transitions across each student's course sequence.
fit <- lta(
  course_engagement,
  vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
  id = "student", n_profiles = 2, time = "sequence", n_starts = 2, seed = 1
)
get_results(fit, what = "transitions")
#>   group_class from to probability expected_count origin_expected_count
#> 1           1    1  1   0.8613986      670.56049              778.4556
#> 2           1    1  2   0.1386014      107.89506              778.4556
#> 3           1    2  1   0.1850266       99.46003              537.5444
#> 4           1    2  2   0.8149734      438.08441              537.5444
#>   supported stable estimated group_class_probability
#> 1      TRUE   TRUE      TRUE                       1
#> 2      TRUE  FALSE      TRUE                       1
#> 3      TRUE  FALSE      TRUE                       1
#> 4      TRUE   TRUE      TRUE                       1
summary(fit)
#> Latent transition model: 2 profiles and 1 group class
#> Observations: 1422; groups: 106; up to 15 occasions (unbalanced, observed grid)
#> Parameters: 23; converged: TRUE
#> Log likelihood: -8335.287419; AIC: 16716.575
#> BIC (groups): 16777.834; BIC (individuals): 16837.551
#> Best likelihood replicated in 2/2 starts (absolute tolerance 0.00711).
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
#> -- gaussian_diagnostics --------------------------------------------
#>   context   profile gaussian_mass positive_rows weight_effective_n
#>  profiles profile_1         836.4          1422              847.0
#>  profiles profile_2         585.6          1422              596.3
#>  min_variance_ratio min_standardized_eigenvalue boundary
#>              0.4007                      0.4007    FALSE
#>              0.3664                      0.3664    FALSE
#> 
#> -- transitions -----------------------------------------------------
#>  group_class from to probability expected_count origin_expected_count supported
#>            1    1  1      0.8614         670.56                 778.5      TRUE
#>            1    1  2      0.1386         107.90                 778.5      TRUE
#>            1    2  1      0.1850          99.46                 537.5      TRUE
#>            1    2  2      0.8150         438.08                 537.5      TRUE
#>  stable estimated group_class_probability
#>    TRUE      TRUE                       1
#>   FALSE      TRUE                       1
#>   FALSE      TRUE                       1
#>    TRUE      TRUE                       1
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
#>  average_posterior odds_correct_classification
#>             0.9915                       81.98
#>             0.9890                      128.23
#> 
#> -- average_posteriors ----------------------------------------------
#>        level assigned_class class n_assigned average_posterior
#>  individuals              1     1        837          0.991531
#>  individuals              1     2        837          0.008469
#>  individuals              2     1        585          0.011015
#>  individuals              2     2        585          0.988985
#> 
#> -- classification_errors -------------------------------------------
#>        level true_class assigned_class probability
#>  individuals          1              1    0.992296
#>  individuals          1              2    0.007704
#>  individuals          2              1    0.012104
#>  individuals          2              2    0.987896
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
#>    ... 2834 more rows.  get_results(x, what = "bch_weights")
#> 
#> -- entropy ---------------------------------------------------------
#>        level n_classes n_units entropy_sum relative_entropy
#>  individuals         2    1422        38.3           0.9611
#>       groups         1     106         0.0               NA
#> 
#> -- replication -----------------------------------------------------
#>  n_starts n_screened n_converged n_failed n_best_replicated share tolerance
#>         2          0           2        0                 2     1   0.00711
#>  best_start log_likelihood
#>           1          -8335
#> 
#> -- residuals -------------------------------------------------------
#>    profile indicator_1 indicator_2     kind observed expected residual
#>  profile_1    lectures  attendance gaussian  0.23622        0  0.23622
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
#>        836.4     6.951 NA 3.638e-12  7.276e-11
#>        585.6     5.730 NA 1.005e-08  5.023e-08
#>        836.4     6.354 NA 2.096e-10  2.096e-09
#>        836.4     6.099 NA 1.068e-09  7.122e-09
#>        585.6     4.813 NA 1.487e-06  4.955e-06
#>        836.4     4.946 NA 7.573e-07  3.029e-06
#>        585.6     2.901 NA 3.714e-03  1.061e-02
#>        585.6     2.800 NA 5.112e-03  1.278e-02
#>        585.6     2.478 NA 1.322e-02  2.938e-02
#>        585.6     1.678 NA 9.329e-02  1.866e-01
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
#>           -8335 16717      16778          16838      TRUE          9    FALSE
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
#>      1          -8335      TRUE          9  <NA>
#>      2          -8335      TRUE          9  <NA>
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
#> 23 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
```
