# Fit a two-level latent profile model

Multilevel latent profile analysis describes heterogeneity among
observations and among the groups they are nested in: measurements
within persons, students within schools. Observation-level profiles are
patterns in the indicators; group-level classes differ in how prevalent
each profile is among their members. A profile has the same measurement
parameters in every group class, so group classes differ only in
composition.

`multilpa()` estimates the measurement parameters, the profile
probabilities within each group class, the group-class proportions and
the posterior memberships at both levels. The default is Gaussian
indicators with diagonal within-profile covariances. The reported start
need not have converged: read *Starts and convergence* before
interpreting a fit.

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
  min_variance = NULL,
  seed = NULL,
  start = NULL,
  missing = c("error", "fiml"),
  profile_slopes = c("shared", "group_class"),
  covariance_model = c("diagonal", "full"),
  categorical = character(),
  min_probability = 1e-10,
  time = NULL,
  fixed = character(),
  centering = c("none", "person", "grand"),
  volume = NULL,
  shape = NULL,
  orientation = NULL,
  model = NULL,
  select_start = c("likelihood", "converged"),
  prior = NULL,
  noise = FALSE,
  acceleration = c("squarem", "none"),
  family = c("profiles", "additive", "dispersion", "additive_dispersion",
    "restricted_cross_level", "full_cross_level"),
  between_variance = c("varying", "equal"),
  weights = NULL,
  ordinal = character(),
  count = character(),
  count_model = c("poisson", "negative_binomial"),
  count_dispersion = c("varying", "equal"),
  scaling = c("none", "person", "grand"),
  drop_empty = TRUE,
  cores = 1L,
  local_dependence = character(),
  direct_effects = character(),
  search = NULL,
  n_screen = NULL,
  boundary = c("allow", "exclude")
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
  without missing values). `id = NULL` fits a single-level model in
  which each row is its own unit; `n_group_classes` must then be one.

- n_profiles:

  Number of individual profiles.

- n_group_classes:

  Number of group classes (default 2).

- profile_covariates, group_covariates:

  Complete numeric columns of `data` that predict profile and
  group-class membership through multinomial logits (last class as
  reference), estimated in one step and entering in their own units. A
  group covariate must be constant within its group. See *Covariates*.

- variance_model:

  `"varying"` (default): each profile has its own variances or
  covariance matrix; `"equal"`: one shared by all profiles.

- n_starts:

  Number of EM starts (default 10). The first is deterministic (Ward's
  clustering), the rest random. See *Starts and convergence*.

- max_iter:

  Maximum EM updates per start (default 1000). `max_iter = 0` evaluates
  the starts without updating them; with `start`, that is the model at
  exactly those values.

- tol:

  Convergence tolerance (default `1e-8`), per observed indicator value:
  EM stops when a step changes the log likelihood by no more than `tol`
  times the number of observed indicator values.

- min_variance:

  Lower bound on every variance (on every covariance eigenvalue under
  full covariance). `NULL` (default) uses `1e-6` times the smallest
  sample variance of the continuous indicators, on the scale fitted; a
  number is an absolute bound in squared indicator units. An estimate at
  the bound is a boundary solution and warns.

- seed:

  Optional integer seed. The caller's random-number state is restored on
  exit.

- start:

  Optional starting values, replacing the first start: a list of
  `means`, `variances` (profiles by indicators), `profile_probabilities`
  (group classes by profiles), `group_probabilities`, and where the
  model has them `covariances` and `response_probabilities` (one
  profiles-by-categories matrix per categorical indicator).
  [`starting_values()`](https://pak.dynasite.org/latents/reference/starting_values.md)
  produces one from a fit.

- missing:

  `"error"` (default) rejects missing indicator values; `"fiml"`
  maximizes the observed-data likelihood under MAR, integrating the
  missing values out. Covariates must be complete either way. See
  *Missing data*.

- profile_slopes:

  `"shared"` (default): one slope per profile logit, common to all group
  classes. `"group_class"`: separate slopes in each group class, a
  cross-level interaction with `n_group_classes` times as many slope
  parameters.

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

- fixed:

  Measurement blocks held at the values in `start` instead of estimated:
  `"means"`, `"variances"` (EEI, VVI, EEE and VVV only),
  `"response_probabilities"`, or `"measurement"` for all; the mixing
  parameters cannot be held. A held block is not counted in
  `n_parameters`.
  [`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
  runs the two-stage workflow in one call.

- centering:

  How to centre the continuous indicators before fitting. `"none"`
  (default) fits them as supplied. `"person"` subtracts each group's own
  mean, so the profiles describe within-unit change and the group
  classes are types of change pattern. `"grand"` subtracts one mean per
  indicator. See *Scale*.

- volume, shape, orientation:

  The covariance structure by its three components: `volume` `"equal"`
  or `"varying"`; `shape` `"equal"`, `"varying"` or `"spherical"`;
  `orientation` `"axis"` (no covariances), `"equal"` or `"varying"`.
  `NULL` (default) follows `variance_model` and `covariance_model`. The
  combinations are the fourteen mclust models:

  |               |               |               |                    |
  |---------------|---------------|---------------|--------------------|
  | `volume`      | `shape`       | `orientation` | model              |
  | equal         | spherical     | \-            | EII                |
  | varying       | spherical     | \-            | VII                |
  | equal/varying | equal/varying | axis          | EEI, VEI, EVI, VVI |
  | equal/varying | equal/varying | equal         | EEE, VEE, EVE, VVE |
  | equal/varying | equal/varying | varying       | EEV, VEV, EVV, VVV |

- model:

  The covariance structure as one of those fourteen codes, such as
  `"VVI"`, instead of the other structure arguments; the codes
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  takes.

- select_start:

  Which start is reported. `"likelihood"` (default): the highest log
  likelihood, preferring a converged start among those within the
  convergence tolerance of it. `"converged"`: the best converged start
  when there is one. See *Starts and convergence*.

- prior:

  `NULL` (default) for maximum likelihood, or a
  [`prior_control()`](https://pak.dynasite.org/latents/reference/prior_control.md)
  object for maximum a posteriori estimation of the Gaussian means and
  covariances under mclust's conjugate prior, with its data-dependent
  defaults. See *Prior and noise*.

- noise:

  `TRUE` adds a noise component of constant density over the data's
  hypervolume, which absorbs outliers and scatter no Gaussian profile
  explains. See *Prior and noise*.

- acceleration:

  `"squarem"` (default) accelerates EM with SQUAREM, which never lowers
  the likelihood; `"none"` runs plain EM. Fits with `prior` or
  membership covariates always use plain EM.

- family:

  `"profiles"` (default) is the model described above. The others have
  group classes without individual profiles (`"additive"`,
  `"dispersion"`, `"additive_dispersion"`) or estimate profiles and
  group classes from the same indicators (`"restricted_cross_level"`,
  `"full_cross_level"`). See *Model families*.

- between_variance:

  For the additive, additive-dispersion and cross-level families only
  (an error with the others): between-group variances per group class
  (`"varying"`, default) or one set shared by all classes (`"equal"`).
  The dispersion family takes `"equal"` only.

- weights:

  Name of a numeric column holding a sampling weight per independent
  unit: per group in a two-level fit (constant within the group), per
  row in a single-level one. See *Weights*.

- ordinal:

  Indicators modelled as ordinal (an ordered factor or integers): an
  adjacent-category logit with category intercepts shared by all
  profiles and one location per profile, Latent GOLD's default ordinal
  model, with `(K - 1) + (C - 1)` parameters per indicator against
  `K * (C - 1)` for `categorical` (`K` profiles, `C` categories). An
  unused factor level warns; drop it with
  [`droplevels()`](https://rdrr.io/r/base/droplevels.html).

- count:

  Indicators holding non-negative whole numbers, modelled as Poisson
  with one mean per profile, or negative binomial under `count_model`.

- count_model:

  `"poisson"` (default) or `"negative_binomial"` (NB2, variance
  `mu + alpha * mu^2`) for counts more dispersed than a Poisson within a
  profile. A dispersion estimated at zero is a boundary solution.

- count_dispersion:

  One negative-binomial dispersion per profile (`"varying"`, default) or
  one shared by all (`"equal"`).

- scaling:

  How to scale the continuous indicators before fitting, at the same
  levels as `centering`: `"none"` (default), `"grand"` (each indicator
  divided by its standard deviation over the sample) or `"person"` (by
  each group's own standard deviation). With the matching centring these
  are z-scores. See *Scale*.

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

- cores:

  Number of processor cores for the starts (default 1); the result is
  identical to one core, and a covariate fit runs on one. To fit several
  models, give `cores` to
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  instead.

- local_dependence:

  Pairs of indicators that stay associated within each class, written
  `"item1:item2"`: a within-class covariance for two continuous
  indicators, a joint distribution for two categorical ones (Latent
  GOLD's `item1 <-> item2 | Cluster`). See *Indicator types*.

- direct_effects:

  Direct effects of covariates on continuous indicators, written
  `"indicator ~ covariate"` (several joined with `+`): the covariate
  shifts the indicator's mean by the same amount in every profile, so
  what it explains is not read as a profile difference; `means` are then
  at covariate zero. See *Indicator types*.

- search:

  `"kmeans"` runs `n_starts` starts to convergence; `"emEM"` screens a
  pool of `n_screen` starts with a short EM run each and runs the first
  start and the best of the rest, `n_starts` in all, to convergence.
  `NULL` (default) is `"emEM"` when continuous indicators share one
  covariance matrix across profiles (EII, EEI or EEE) and `"kmeans"`
  otherwise.

- n_screen:

  For `search = "emEM"` only, the size of the screened pool (default
  five times `n_starts`, at least `n_starts`).

- boundary:

  Whether a start with a variance at `min_variance` may be reported:
  `"allow"` (default) does, as Mplus does, with a warning; `"exclude"`
  chooses among the other starts, and reports the best with a warning
  when every start is at the bound.

## Value

A `multilpa` object: the measurement parameters, the profile
probabilities within each group class, the group-class proportions, the
posteriors and classifications at both levels, `log_likelihood`, the
information criteria (`bic` and `bic_groups` on the number of groups,
`bic_individual` on the number of rows with an observed indicator), the
starts table, `converged` and `boundary`.
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
returns its tables,
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
its standard errors and
[`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
tests it against a nested fit.

## Starts and convergence

The first start is Ward's clustering of the standardized indicators; the
others are k-means from random centres, or under `search = "emEM"` a
screened pool that also holds random partitions.
`get_results(fit, "starts")` lists every start. The reported start is
the one `select_start` picks; `converged` says whether it converged, and
[`summary()`](https://rdrr.io/r/base/summary.html) says how many starts
replicated its log likelihood. A replicated likelihood does not
establish a global optimum. With profile-specific variances the highest
likelihood can belong to a nearly degenerate profile on a few
observations: a variance at `min_variance` warns, `boundary = "exclude"`
sets such starts aside, and `get_results(fit, "gaussian_diagnostics")`
describes each profile's support. Wald inference is refused for a fit
that did not converge or has a parameter at its bound.

More than one group class needs more than one profile and a group with
more than one row; a model of categorical and ordinal indicators with
more free parameters than its table has free cells warns. Neither check
establishes identification. Labels are arbitrary.

## Missing data

Under `missing = "fiml"` EM uses the observed-data likelihood; only the
initialization fills missing values with means. A row with no observed
indicator contributes no measurement likelihood but receives posteriors
from its group, and `drop_empty` governs how it enters the reported
statistics. Missing patterns can prevent identification.

## Covariates

Naming `profile_covariates` or `group_covariates` fits the one-step
covariate model (class `multilpa_covariates`), with profile intercepts
per group class, slopes shared unless `profile_slopes = "group_class"`,
and Wald standard errors from
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md).
With `start` and `fixed = "measurement"` from a fit without covariates,
only the logits are estimated: the two-step estimator, whose criteria
and standard errors condition on the held measurement.
[`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
and [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md)
are the staged alternatives. A covariate model does not take `prior`,
`noise`, `direct_effects` or a continuous `local_dependence` pair.

## Indicator types

Continuous indicators are Gaussian within profiles; `categorical`,
`ordinal` and `count` indicators are independent of each other and of
the continuous ones given the profile. Ordinal and count indicators work
with FIML, weights, two-level fits, covariates, inference,
[`predict()`](https://rdrr.io/r/stats/predict.html) and
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md), but not
with `start`, `fixed`, `prior`, `noise` or the other families. A
categorical `local_dependence` pair becomes one indicator over the
combinations of its categories; an unobserved combination sits on the
`min_probability` floor and still counts as a parameter.
`local_dependence` is not combined with the structure arguments
(`covariance_model`, `model`, `volume`, `shape`, `orientation`), `time`,
`prior`, `noise` or held variances, and a continuous pair not with
covariates either; `direct_effects` needs a diagonal covariance and
complete covariates, and is not combined with membership covariates,
continuous pairs, `time`, `prior`, `noise`, `fixed`, `start`,
`centering` or `scaling`. The tables `"ordinal"`, `"count_means"`,
`"local_dependence"` and `"direct_effects"` of
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
hold the estimates.

## Scale

Centring and scaling are applied before fitting and stored on the fit,
so
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
returns the supplied data and
[`predict()`](https://rdrr.io/r/stats/predict.html) transforms new units
the same way. Estimates, likelihood and criteria are on the fitted
scale, so only fits with the same centring and scaling are comparable.
Grand scaling changes only the units under a structure unchanged by
rescaling each indicator (EEI, VEI, EVI, VVI, EEE, VEE, EVV, VVV) and
changes the solution under the others. Person centring needs each group
observed more than once on each indicator, and refuses singly observed
indicators even when other groups have repeats; person scaling also
needs variation within the group, and a group observed exactly twice is
flagged as uninformative. Grand offsets and scales use sampling weights
when supplied. Grand scaling divides by the population standard
deviation (observed weight in the denominator); person scaling uses the
within-group sample standard deviation.

## Weights

The fit maximizes the pseudo log likelihood `sum(w_j * log(L_j))`, the
weights scaled to sum to the number of units of positive weight as Mplus
does, so integer weights act as repeated units. A unit of weight zero
changes nothing, enters no statistic, and keeps its row and posterior in
the assignment tables. Standard errors are sandwich estimates
(`vcov_type = "robust"`); the bootstrap resamples units with their
weights. Weights work with every family, covariates, categorical
indicators and FIML;
[`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md),
[`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md),
[`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md),
`prior` and `noise` refuse a weighted fit.

## Model families

Three group-class-only families: each group has a Gaussian intercept per
indicator drawn from its class's distribution, and its members vary
around it. `"additive"`: classes differ in means (and in between-group
variances under `between_variance = "varying"`), with one within-group
variance per indicator; `"dispersion"`: in within-group variances;
`"additive_dispersion"`: in both. They take only `data`, `vars`, `id`,
`n_group_classes`, `between_variance`, `n_starts`, `max_iter`, `tol`,
`min_variance`, `seed`, `weights`, `search` and `n_screen`, need
complete continuous indicators, and return a `multilpa_additive` fit
with [`summary()`](https://rdrr.io/r/base/summary.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html), Wald
inference,
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
by `family` and
[`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
across numbers of group classes. They are experimental: no bootstrap
inference, FIML, covariates or discrete indicators.

Two cross-level families estimate profiles and group classes from the
same indicators, each group's indicator means being between-level
indicators of its class: `"restricted_cross_level"` holds the profile
prevalences equal across classes and estimates the levels separately;
`"full_cross_level"` lets them differ. They take `data`, `vars`, `id`,
`n_profiles`, `n_group_classes`, `variance_model`, `between_variance`,
`n_starts`, `max_iter`, `tol`, `min_variance`, `seed` and `weights`,
have no standard errors, and return a `multilpa_cross_level` fit whose
likelihood is computed on the same scale only between cross-level fits
of the same data. It counts the ratings twice, once directly and once
through their group means, so its criteria favour group classes and the
full form even on data with no group structure: they do not choose the
number of group classes or the form, and a fit with group classes warns
with class `latents_cross_level_criteria`.

## Prior and noise

`prior` is available for EII, VII, EEI, VEI, EVI, VVI, EEE, EEV, VEV and
VVV with complete continuous indicators and without `fixed` or
covariates; the mixing proportions have no prior. As in mclust, the
likelihood and criteria are unpenalized values at the posterior mode;
Wald inference is not defined there, and the bootstrap refits each
resample with the same prior.

`noise` is available for complete continuous indicators with one group
class, under any structure, with or without `prior`, and not with
`fixed` or covariates. The component's density is `1 / V`, with `V`
mclust's `hypvol()`, and it adds two parameters; a supplied `start`
includes `noise_probability`. It is profile `0` in the classification:
`subject_profiles` is `0` there, the assignment table has
`posterior_noise`, and the posterior and profile-probability tables have
profile-0 rows, while the measurement tables describe the Gaussian
profiles only. Inference, the three-step functions,
[`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md),
[`starting_values()`](https://pak.dynasite.org/latents/reference/starting_values.md),
[`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md),
bivariate residuals and the posterior plots refuse a noise fit.

Every refusal above is a classed condition;
[latents-conditions](https://pak.dynasite.org/latents/reference/latents-conditions.md)
lists them.

## See also

[`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
for a model in which every indicator is categorical;
[`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md) for a
single-level model;
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
to compare numbers of profiles and group classes.

## Examples

``` r
# Course sessions nested in students: engagement profiles of sessions, and
# classes of students that differ in how often they are engaged.
activity <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
fit <- multilpa(course_engagement, activity, "student",
                n_profiles = 2, n_group_classes = 2, n_starts = 4,
                seed = 42)
summary(fit)
#> Two-level latent profile analysis: 2 profiles, 2 group classes
#> Individuals: 1422; groups: 106; parameters: 23; converged: TRUE
#> Log likelihood: -8439.816424; AIC: 16925.633
#> BIC (groups): 16986.892; BIC (individuals): 17046.609
#> Best likelihood replicated in 4/4 starts (absolute tolerance 0.00711).
#> 
#> -- profiles --------------------------------------------------------
#>  profile  indicator    mean variance standard_deviation mean_standard_error
#>        1     browse  0.5396   0.5899             0.7680             0.02710
#>        1   lectures  0.4274   0.8447             0.9191             0.03249
#>        1 forum_read  0.6198   0.4815             0.6939             0.02473
#>        1 forum_post  0.5305   0.6888             0.8299             0.02953
#>        1 attendance  0.6490   0.3841             0.6198             0.02226
#>        2     browse -0.7656   0.5285             0.7270             0.03121
#>        2   lectures -0.6065   0.5384             0.7337             0.03093
#>        2 forum_read -0.8792   0.3631             0.6026             0.02622
#>        2 forum_post -0.7527   0.4212             0.6490             0.02768
#>        2 attendance -0.9207   0.3736             0.6113             0.02643
#>  variance_standard_error
#>                  0.02948
#>                  0.04210
#>                  0.02440
#>                  0.03482
#>                  0.01956
#>                  0.03222
#>                  0.03209
#>                  0.02279
#>                  0.02608
#>                  0.02279
#> 
#> -- responses -------------------------------------------------------
#>    (no rows)
#> 
#> -- covariances -----------------------------------------------------
#>  profile  indicator indicator_2 covariance
#>        1     browse      browse     0.5899
#>        1   lectures      browse     0.0000
#>        1 forum_read      browse     0.0000
#>        1 forum_post      browse     0.0000
#>        1 attendance      browse     0.0000
#>        1     browse    lectures     0.0000
#>        1   lectures    lectures     0.8447
#>        1 forum_read    lectures     0.0000
#>        1 forum_post    lectures     0.0000
#>        1 attendance    lectures     0.0000
#>    ... 40 more rows.  get_results(x, what = "covariances")
#> 
#> -- gaussian_diagnostics --------------------------------------------
#>   context   profile gaussian_mass positive_rows weight_effective_n
#>  profiles profile_1         834.1          1422                852
#>  profiles profile_2         587.9          1422                606
#>  min_variance_ratio min_standardized_eigenvalue boundary
#>              0.3930                      0.3930    FALSE
#>              0.3715                      0.3715    FALSE
#> 
#> -- profile_probabilities -------------------------------------------
#>  group_class profile probability group_class_probability
#>            1       1      0.1705                  0.3253
#>            1       2      0.8295                  0.3253
#>            2       1      0.7859                  0.6747
#>            2       2      0.2141                  0.6747
#> 
#> -- counts ----------------------------------------------------------
#>        level class effective_count effective_proportion
#>  individuals     1          834.10               0.5866
#>  individuals     2          587.90               0.4134
#>       groups     1           34.48               0.3253
#>       groups     2           71.52               0.6747
#> 
#> -- posteriors ------------------------------------------------------
#>  row group profile posterior modal
#>    1     1       1 9.996e-01  TRUE
#>    2     1       1 9.851e-01  TRUE
#>    3     1       1 2.421e-06 FALSE
#>    4     1       1 3.863e-06 FALSE
#>    5     1       1 5.749e-06 FALSE
#>    6     1       1 2.278e-05 FALSE
#>    7     1       1 2.301e-05 FALSE
#>    8     1       1 1.107e-05 FALSE
#>    9     1       1 2.173e-06 FALSE
#>   10     1       1 9.510e-06 FALSE
#>    ... 2834 more rows.  get_results(x, what = "posteriors")
#> 
#> -- group_posteriors ------------------------------------------------
#>  group group_size log_likelihood group_class posterior modal
#>      1         13         -61.52           1 1.000e+00  TRUE
#>      2         13         -71.94           1 1.000e+00  TRUE
#>      3         14         -76.69           1 6.005e-10 FALSE
#>      4         13         -67.56           1 1.327e-09 FALSE
#>      5         14         -78.80           1 1.000e+00  TRUE
#>      6         12         -80.90           1 1.984e-02 FALSE
#>      7         15         -89.36           1 1.123e-09 FALSE
#>      8         14         -85.22           1 9.989e-01  TRUE
#>      9         14         -83.80           1 4.309e-10 FALSE
#>     10         15         -82.06           1 1.000e+00  TRUE
#>    ... 202 more rows.  get_results(x, what = "group_posteriors")
#> 
#> -- assignments -----------------------------------------------------
#>  student browse lectures forum_read forum_post attendance profile group_class
#>        1   0.73    -0.08       0.30       0.67       0.72       1           1
#>        1   0.68     0.53      -0.02      -0.55       0.70       1           1
#>        1  -1.51    -0.51      -0.91      -0.80      -0.97       2           1
#>        1   0.64    -1.02      -1.62      -1.20      -1.26       2           1
#>        1  -0.94    -0.37      -0.71      -0.28      -1.52       2           1
#>        1  -0.36    -0.46      -0.77      -0.43      -1.34       2           1
#>        1   0.38    -1.85      -0.84      -1.09      -1.12       2           1
#>        1  -1.44     0.66      -0.67      -1.67      -1.00       2           1
#>        1  -0.76    -0.43      -1.04      -0.69      -1.37       2           1
#>        1   0.26    -0.95      -0.51      -1.89      -1.45       2           1
#>  uncertainty posterior_profile_1 posterior_profile_2
#>    4.436e-04           9.996e-01           0.0004436
#>    1.485e-02           9.851e-01           0.0148535
#>    2.421e-06           2.421e-06           0.9999976
#>    3.863e-06           3.863e-06           0.9999961
#>    5.749e-06           5.749e-06           0.9999943
#>    2.278e-05           2.278e-05           0.9999772
#>    2.301e-05           2.301e-05           0.9999770
#>    1.107e-05           1.107e-05           0.9999889
#>    2.173e-06           2.173e-06           0.9999978
#>    9.510e-06           9.510e-06           0.9999905
#>    ... 1412 more rows.  get_results(x, what = "assignments")
#> 
#> -- classification --------------------------------------------------
#>        level class n_modal proportion_modal estimated_n estimated_proportion
#>  individuals     1     836           0.5879      834.10               0.5866
#>  individuals     2     586           0.4121      587.90               0.4134
#>       groups     1      35           0.3302       34.48               0.3253
#>       groups     2      71           0.6698       71.52               0.6747
#>  average_posterior odds_correct_classification
#>             0.9849                       46.08
#>             0.9818                       76.35
#>             0.9690                       64.80
#>             0.9920                       59.89
#> 
#> -- average_posteriors ----------------------------------------------
#>        level assigned_class class n_assigned average_posterior
#>  individuals              1     1        836          0.984935
#>  individuals              1     2        836          0.015065
#>  individuals              2     1        586          0.018243
#>  individuals              2     2        586          0.981757
#>       groups              1     1         35          0.968985
#>       groups              1     2         35          0.031015
#>       groups              2     1         71          0.007986
#>       groups              2     2         71          0.992014
#> 
#> -- classification_errors -------------------------------------------
#>        level true_class assigned_class probability
#>  individuals          1              1     0.98718
#>  individuals          1              2     0.01282
#>  individuals          2              1     0.02142
#>  individuals          2              2     0.97858
#>       groups          1              1     0.98356
#>       groups          1              2     0.01644
#>       groups          2              1     0.01518
#>       groups          2              2     0.98482
#> 
#> -- bch_weights -----------------------------------------------------
#>        level unit assigned_class class   weight
#>  individuals    1              1     1  1.01327
#>  individuals    1              1     2 -0.01327
#>  individuals    2              1     1  1.01327
#>  individuals    2              1     2 -0.01327
#>  individuals    3              2     1 -0.02218
#>  individuals    3              2     2  1.02218
#>  individuals    4              2     1 -0.02218
#>  individuals    4              2     2  1.02218
#>  individuals    5              2     1 -0.02218
#>  individuals    5              2     2  1.02218
#>    ... 3046 more rows.  get_results(x, what = "bch_weights")
#> 
#> -- entropy ---------------------------------------------------------
#>        level n_classes n_units entropy_sum relative_entropy
#>  individuals         2    1422      60.940           0.9382
#>       groups         2     106       4.435           0.9396
#> 
#> -- replication -----------------------------------------------------
#>  n_starts n_screened n_converged n_failed n_best_replicated share tolerance
#>         4          0           4        0                 4     1   0.00711
#>  best_start log_likelihood
#>           1          -8440
#> 
#> -- residuals -------------------------------------------------------
#>    profile indicator_1 indicator_2     kind observed expected residual
#>  profile_1    lectures  attendance gaussian  0.23393        0  0.23393
#>  profile_2  forum_read  attendance gaussian  0.22869        0  0.22869
#>  profile_1  forum_read  attendance gaussian  0.20880        0  0.20880
#>  profile_1  forum_post  attendance gaussian  0.20314        0  0.20314
#>  profile_2      browse  attendance gaussian  0.19248        0  0.19248
#>  profile_1      browse  attendance gaussian  0.17312        0  0.17312
#>  profile_2    lectures  attendance gaussian  0.12018        0  0.12018
#>  profile_2  forum_post  attendance gaussian  0.11883        0  0.11883
#>  profile_2      browse  forum_read gaussian  0.10925        0  0.10925
#>  profile_2  forum_read  forum_post gaussian  0.06676        0  0.06676
#>  effective_n statistic df   p_value p_adjusted
#>        834.1     6.871 NA 6.374e-12  1.275e-10
#>        587.9     5.630 NA 1.799e-08  8.994e-08
#>        834.1     6.109 NA 1.001e-09  1.001e-08
#>        834.1     5.939 NA 2.867e-09  1.912e-08
#>        587.9     4.714 NA 2.429e-06  8.098e-06
#>        834.1     5.042 NA 4.619e-07  1.847e-06
#>        587.9     2.921 NA 3.493e-03  9.708e-03
#>        587.9     2.888 NA 3.883e-03  9.708e-03
#>        587.9     2.653 NA 7.983e-03  1.774e-02
#>        587.9     1.617 NA 1.059e-01  2.118e-01
#>    ... 10 more rows.  get_results(x, what = "residuals")
#> 
#> -- information_criteria --------------------------------------------
#>  log_likelihood n_parameters   aic   kic bic_groups bic_individual sabic_groups
#>           -8440           23 16926 16952      16987          17047        16914
#>  sabic_individual caic_groups caic_individual awe_groups awe_individual
#>             16974       17010           17070      17292          17412
#>  icl_groups icl_individual clc_groups clc_individual
#>       17116          17176      17009          17009
#> 
#> -- model -----------------------------------------------------------
#>  n_observations n_informative n_groups n_profiles n_group_classes centering
#>            1422          1422      106          2               2      none
#>  scaling covariance_structure n_parameters n_parameters_with_measurement
#>     none                  VVI           23                            23
#>  log_likelihood   aic bic_groups bic_individual converged iterations boundary
#>           -8440 16926      16987          17047      TRUE          9    FALSE
#>  small_classes best_start n_best_replicated weights
#>          FALSE          1                 4    <NA>
#> 
#> -- stages ----------------------------------------------------------
#>  stage group_classes fixed log_likelihood parameters
#>  joint             2  <NA>          -8440         23
#>  parameters_with_measurement converged
#>                           23      TRUE
#> 
#> -- starts ----------------------------------------------------------
#>  start log_likelihood converged iterations error boundary
#>      1          -8440      TRUE          9  <NA>    FALSE
#>      2          -8440      TRUE          9  <NA>    FALSE
#>      3          -8440      TRUE          9  <NA>    FALSE
#>      4          -8440      TRUE          9  <NA>    FALSE
#> 
#> -- data ------------------------------------------------------------
#>  student browse lectures forum_read forum_post attendance
#>        1   0.73    -0.08       0.30       0.67       0.72
#>        1   0.68     0.53      -0.02      -0.55       0.70
#>        1  -1.51    -0.51      -0.91      -0.80      -0.97
#>        1   0.64    -1.02      -1.62      -1.20      -1.26
#>        1  -0.94    -0.37      -0.71      -0.28      -1.52
#>        1  -0.36    -0.46      -0.77      -0.43      -1.34
#>        1   0.38    -1.85      -0.84      -1.09      -1.12
#>        1  -1.44     0.66      -0.67      -1.67      -1.00
#>        1  -0.76    -0.43      -1.04      -0.69      -1.37
#>        1   0.26    -0.95      -0.51      -1.89      -1.45
#>    ... 1412 more rows.  get_results(x, what = "data")
#> 
#> 21 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
as.data.frame(fit)
#>    profile  indicator       mean  variance standard_deviation
#> 1        1     browse  0.5395774 0.5898715          0.7680309
#> 2        1   lectures  0.4274031 0.8446923          0.9190714
#> 3        1 forum_read  0.6198249 0.4814743          0.6938835
#> 4        1 forum_post  0.5304888 0.6887895          0.8299335
#> 5        1 attendance  0.6490498 0.3841263          0.6197793
#> 6        2     browse -0.7655998 0.5284684          0.7269583
#> 7        2   lectures -0.6064852 0.5383733          0.7337392
#> 8        2 forum_read -0.8792137 0.3631130          0.6025886
#> 9        2 forum_post -0.7526883 0.4211624          0.6489703
#> 10       2 attendance -0.9206599 0.3736293          0.6112523
#>    mean_standard_error variance_standard_error
#> 1           0.02710160              0.02948401
#> 2           0.03248649              0.04210054
#> 3           0.02473030              0.02439626
#> 4           0.02953027              0.03481645
#> 5           0.02226482              0.01955858
#> 6           0.03121144              0.03221700
#> 7           0.03093059              0.03209025
#> 8           0.02621510              0.02279295
#> 9           0.02768122              0.02607600
#> 10          0.02642702              0.02278996
get_results(fit, what = "profile_probabilities")
#>   group_class profile probability group_class_probability
#> 1           1       1   0.1704508                0.325297
#> 2           1       2   0.8295492                0.325297
#> 3           2       1   0.7859317                0.674703
#> 4           2       2   0.2140683                0.674703
```
