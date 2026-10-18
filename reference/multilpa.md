# Fit a two-level latent profile model

Fits a latent profile model for observations nested in groups and
returns a fitted `multilpa` object. Each observation belongs to a
Gaussian profile, and each group belongs to a latent group class that
determines the profile proportions. Profile means and variances are
shared across group classes (measurement invariance). By default
indicators are independent within a profile; full residual covariances
and observed-data maximum likelihood for missing indicators are
available.

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
  boundary = c("allow", "exclude"),
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
  cores = 1L
)
```

## Arguments

- data:

  A data frame containing the indicators and the group identifier.

- vars:

  Unique character vector of indicator column names. Columns not named
  in `categorical`, `ordinal` or `count` are continuous indicators.

- id:

  Name of the group identifier column. Character, factor and numeric
  identifiers are supported; missing identifiers are not. `id` has no
  default: omitting it raises `latents_bad_argument`.

  `id = NULL` fits a **single-level** model: each row is its own unit,
  the observations are independent, and `n_group_classes` is one. The
  fit raises a `latents_single_level` message. Every verb of the package
  works on the result. The internal unit column is named `.observation`;
  it is not returned by `get_results(x, "data")`, and a `data` that
  already has a column of that name raises `latents_bad_data`.
  Requesting more than one group class with `id = NULL` raises
  `latents_bad_argument`.

- n_profiles:

  Positive integer number of individual profiles.

- n_group_classes:

  Positive integer number of latent group classes (default 2). With
  `id = NULL` it is one, and any other value is an error.

- profile_covariates, group_covariates:

  Names of numeric columns of `data` that predict individual-profile and
  group-class membership through multinomial logits, with the last class
  as reference. Naming either fits the one-step maximum-likelihood
  covariate model and returns a `multilpa_covariates` object. Profile
  intercepts differ by group class; profile slopes are shared across
  group classes unless `profile_slopes = "group_class"`. A group
  covariate must be constant within each group. Covariates enter in
  their supplied units, without centring or scaling. The covariate model
  refuses `start` and `fixed`; it accepts `missing = "fiml"` for the
  indicators, but the covariates must be complete.
  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
  and [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md)
  are the staged alternatives.

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
  1000). `max_iter = 0` performs no update: with `start`, the model is
  evaluated at exactly those values, for example to score them with
  [`logLik()`](https://rdrr.io/r/stats/logLik.html); without `start`,
  the random initializations are evaluated and the one with the highest
  log likelihood is returned.

- tol:

  Positive relative log-likelihood tolerance (default `1e-8`).
  Convergence requires an absolute change no greater than
  `tol * (1 + abs(previous log likelihood))`.

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

- start:

  Optional list of `means` and `variances` (profiles by indicators),
  `profile_probabilities` (group classes by profiles), and
  `group_probabilities` (vector). Starting probabilities must be
  positive. For full covariance, supply `covariances` (indicators by
  indicators by profiles); `variances` may then be omitted, or must
  match their diagonals. A categorical model also takes
  `response_probabilities`, a list of one profiles-by-categories matrix
  per indicator named in `categorical`. The list is matched by position
  unless it is labelled: when its elements are named after the
  indicators, or its columns after the categories, each block is matched
  to the indicator and category its labels name. A label naming an
  indicator or category the fit does not have raises
  `latents_bad_start`.

- missing:

  `"error"` (the default) rejects missing indicators; `"fiml"` maximizes
  the observed-data likelihood under an ignorable missingness mechanism
  (MAR), integrating missing indicators out rather than filling them in.
  This also holds with membership covariates; the covariates themselves
  must be complete (`latents_bad_data`).

- profile_slopes:

  `"shared"` (the default) gives each profile covariate one slope per
  profile logit, common to every group class. `"group_class"` estimates
  separate slopes in each group class (a cross-level interaction),
  reported with terms such as `z:group_class_1`; it multiplies the
  number of slope parameters by `n_group_classes`, each set estimated
  from the groups in its class. Naming it without `profile_covariates`
  raises `latents_bad_argument`.

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
  as [`sort()`](https://rdrr.io/r/base/sort.html) orders them. Numbers
  stored as text whose alphabetical order is not their numeric order
  (`"2"`, `"10"`) raise `latents_bad_time`; zero-padded numbers and ISO
  dates are accepted.

- fixed:

  Character vector naming measurement blocks to hold at the values in
  `start` instead of estimating them: any of `"means"`, `"variances"`
  and `"response_probabilities"`, or `"measurement"` for every block the
  model has. Under `covariance_model = "full"`, `"variances"` holds the
  residual covariance matrices. A held block stays as supplied in every
  start and is not counted in `n_parameters`, so the result is a
  different model, not a different starting point. `start` must contain
  the named blocks; `starting_values(fit, what = "measurement")`
  produces them, and
  [`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
  performs the two-stage workflow in one call. The mixing parameters
  cannot be held.

- centering:

  How to centre the continuous indicators before fitting. `"none"` (the
  default) fits them as supplied. `"person"` subtracts each group's own
  mean from its rows, so values are deviations from the unit's average,
  the profiles describe within-unit change rather than level, and the
  group classes are types of change pattern rather than types of unit
  (the person-mean-centred within-person design; Quintana, 2021;
  Voelkle, Brose, Schmiedek, & Lindenberger, 2014). `"grand"` subtracts
  one mean per indicator, which moves the origin without changing the
  within-group structure. The offsets are stored on the fit, so
  `get_results(x, "data")` returns the supplied columns and
  row-alignment checks still apply. `"person"` raises `latents_bad_data`
  when it leaves an indicator constant, as happens when a unit has one
  observation of it.

- volume, shape, orientation:

  The covariance structure in its three components. Each profile's
  covariance decomposes as `Sigma_k = lambda_k * D_k * A_k * D_k'`, with
  *volume* `lambda_k = |Sigma_k|^(1/d)`, *orientation* `D_k`
  (eigenvectors) and *shape* `A_k` (diagonal, determinant one).
  Constraining the three across profiles gives the fourteen models
  mclust names with three letters, all of which are available.

  `volume` is `"equal"` or `"varying"`. `shape` is `"equal"`,
  `"varying"` or `"spherical"`; `"spherical"` makes every indicator's
  spread equal within a profile and leaves no orientation to constrain.
  `orientation` is `"axis"` (a diagonal covariance), `"equal"` (one
  orientation shared by every profile) or `"varying"`. Each defaults to
  `NULL`, which follows `variance_model` and `covariance_model`.

  |               |               |               |                    |
  |---------------|---------------|---------------|--------------------|
  | `volume`      | `shape`       | `orientation` | model              |
  | equal         | spherical     | —             | EII                |
  | varying       | spherical     | —             | VII                |
  | equal/varying | equal/varying | axis          | EEI, VEI, EVI, VVI |
  | equal/varying | equal/varying | equal         | EEE, VEE, EVE, VVE |
  | equal/varying | equal/varying | varying       | EEV, VEV, EVV, VVV |

  The estimates are those of Celeux and Govaert (1995); EVE and VVE,
  which combine a shared orientation with a free shape, have no closed
  form and use the minorize-maximize step of Browne and McNicholas
  (2014). Parameter counts match mclust's for all fourteen.

  Structures other than EEI, VVI, EEE and VVV are maximized across all
  profiles jointly and cannot be combined with a held `variances` block.
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  reports Wald standard errors for all fourteen, computed in the
  structure's free coordinates (log volumes, determinant-one log shapes
  and orientations) and transformed to the reported variances and
  covariances by the delta method; `method = "bootstrap"` is also
  available for all fourteen.

- model:

  The covariance structure as an mclust code, one of `"EII"`, `"VII"`,
  `"EEI"`, `"VEI"`, `"EVI"`, `"VVI"`, `"EEE"`, `"VEE"`, `"EVE"`,
  `"VVE"`, `"EEV"`, `"VEV"`, `"EVV"`, `"VVV"`. The letters give the
  volume, shape and orientation, each equal (`E`) or varying (`V`)
  across profiles, with `I` for axis-aligned (no covariances). These are
  the codes
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  takes. Use it instead of `variance_model`/`covariance_model` or
  `volume`/`shape`/`orientation`, not together with them.

- select_start:

  Which start the fit reports. `"likelihood"` (the default) takes the
  highest log likelihood across starts, preferring a converged start
  only among those within a relative `1e-10` of it. `"converged"` takes
  the highest log likelihood among the converged starts when at least
  one converged, and among all starts otherwise. The two differ when an
  unconverged start exceeds a converged one by a negligible amount, as
  happens when a membership logit drifts along a flat ridge of the
  likelihood; only the converged start is a maximum on which
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  can be based. [`summary()`](https://rdrr.io/r/base/summary.html) and
  `get_results(fit, "starts")` show every start.

- boundary:

  Whether a start with a variance at `min_variance` (a covariance
  eigenvalue, under full covariance) may be reported. `"allow"` (the
  default) permits it, as Mplus does, with a `latents_boundary` warning.
  `"exclude"` sets such starts aside and applies `select_start` to the
  rest; for models without covariates, the `boundary` column of the
  starts table identifies them. When every start is at the bound, the
  best is reported with a `latents_all_starts_boundary` warning. A
  Gaussian profile can increase the likelihood without limit by
  shrinking around one observation, so the size of that gain is set by
  `min_variance` rather than by the data, and more starts make such a
  solution more likely to be found. Categorical indicators have no
  variances and are unaffected.

- prior:

  `NULL` (the default) for maximum likelihood, or a
  [`prior_control()`](https://pak.dynasite.org/latents/reference/prior_control.md)
  object for maximum a posteriori estimation of the Gaussian means and
  covariances under the conjugate prior of Fraley and Raftery (2007),
  mclust's `priorControl()`. The hyperparameters default to mclust's
  `defaultPrior()` computed from the indicators, and each M-step is
  mclust's, so a fit from the same start reproduces
  `mclust::Mclust(..., prior = priorControl())`. Available for EII, VII,
  EEI, VEI, EVI, VVI, EEE, EEV, VEV and VVV (mclust has no prior for
  VEE, EVE, VVE and EVV), with complete continuous indicators and
  without `fixed` or membership covariates; otherwise
  `latents_unsupported_prior` is raised. The mixing proportions have no
  prior. As in mclust, the fit's `log_likelihood` is the *unpenalized*
  log likelihood at the posterior mode, and `aic`, `bic` and their
  variants are computed from it with the usual parameter count. EM
  convergence is judged on its relative change; because a MAP step can
  lower that likelihood, the decrease check made under maximum
  likelihood is not applied. `parameter_inference(method = "wald")`
  raises `latents_unsupported_inference` for such a fit, because the
  score is not zero at a posterior mode; `method = "bootstrap"` refits
  every resample with the same prior. The fit stores the request as
  `prior` and the resolved hyperparameters, in the indicators' units
  (centred, under `centering`), as `prior_parameters`.

- noise:

  `FALSE` (the default), or `TRUE` to add a noise component: a mixture
  component with constant density `1 / V` over the hypervolume `V` of
  the data, which absorbs outliers and scatter that no Gaussian profile
  explains (Banfield & Raftery, 1993; mclust's
  `initialization = list(noise = )` and `Vinv`). `V` is mclust's
  `hypvol()`: the smaller of the volumes of the data's axis-aligned and
  principal-component-aligned bounding boxes. The component has no
  measurement parameters; its mixing proportion is estimated with the
  others, and, as in mclust, it adds two parameters (the proportion and
  the hypervolume), so the criteria equal mclust's. A random start puts
  a tenth of the mass on the noise; a `start` supplies
  `noise_probability` alongside `profile_probabilities`, the two summing
  to one. Available for complete continuous indicators with one group
  class (`id = NULL` or `n_group_classes = 1`), under any of the
  fourteen covariance structures, with or without `prior`; otherwise
  `latents_unsupported_noise` is raised. In the fit the noise component
  is profile `0`, as in mclust's classification: `subject_profiles` is
  `0` for an observation assigned to it,
  `get_results(fit, "assignments")` adds a `posterior_noise` column,
  `get_results(fit, "posteriors")` and
  `get_results(fit, "profile_probabilities")` have profile-`0` rows, and
  the classification diagnostics count it as a class. `means`,
  `variances`, `covariances`, `profile_probabilities` and
  `subject_posteriors` describe the Gaussian profiles only, so the last
  two sum to one minus the noise share; the fit adds
  `noise_probability`, `noise_posteriors`, `n_noise` and `hypervolume`.
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md),
  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md),
  [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md),
  [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md),
  [`starting_values()`](https://pak.dynasite.org/latents/reference/starting_values.md),
  bivariate residuals,
  [`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
  and the posterior plots raise `latents_unsupported_noise` for such a
  fit.

- acceleration:

  `"squarem"` (the default) accelerates EM with SQUAREM (Varadhan &
  Roland, 2008): each cycle extrapolates along the last two EM steps and
  ends with an ordinary EM step, and falls back to plain EM when the
  extrapolated point would lower the likelihood, so the iteration stays
  monotone. Convergence is assessed over a whole cycle, so it is never
  looser than plain EM's. `"none"` runs plain EM, which follows the same
  path as mclust's `me()` from the same start. `iterations` and
  `max_iter` count EM steps under either setting. Fits with `prior`, and
  models with membership covariates, always use plain EM. For EVE and
  VVE, whose shared orientation has no closed form, each EM iteration
  takes a fixed number of warm-started orientation steps (a generalized
  EM step) under either setting, which keeps every iteration monotone
  and the EM map a fixed function, as SQUAREM requires.

- family:

  The model family. `"profiles"` (the default) is the model described
  above: individual profiles whose prevalences differ across group
  classes. Three group-class-only families have no individual profiles
  (Houle, Morin & Harvey, 2026): each group has a Gaussian intercept per
  indicator, drawn from its group class's distribution, and its members
  vary around it.

  - `"additive"`: group classes differ in their means (and, optionally,
    in between-group variances); one within-group variance per indicator
    is shared by every class.

  - `"dispersion"`: group classes differ in their within-group
    variances; means and between-group variances are shared.

  - `"additive_dispersion"`: group classes differ in means and in
    within-group variances.

  These families take only `data`, `vars`, `id`, `n_group_classes`,
  `between_variance`, `n_starts`, `max_iter`, `tol`, `min_variance`,
  `seed` and `weights`; any other argument raises
  `latents_bad_argument`. Indicators must be complete and continuous,
  and at least one group must have more than one row
  (`latents_unidentified` otherwise). Between variances estimated at
  zero are reached by a boundary maximization with a Karush-Kuhn-Tucker
  check and flagged (`latents_boundary`); interior fits are finished by
  Newton steps. The result has class `multilpa_additive`; `bic` uses the
  number of groups and `bic_individual` the number of rows. It is read
  with
  [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md);
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html) and
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  work on it,
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  compares families through its `family` argument, and
  [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
  tests one number of group classes against another. These families are
  experimental: `parameter_inference(method = "bootstrap")` raises
  `latents_unsupported_inference`, and missing indicators
  (`latents_bad_data`), `missing`, membership covariates and ordinal or
  count indicators (`latents_bad_argument`) are not available. See the
  [group-class and cross-level models
  article](https://pak.dynasite.org/latents/articles/additive.html).

  Two cross-level families estimate individual profiles and group
  classes from the same indicators, following the manifest-aggregation
  specification of Houle et al. (2026, supplement): each group's means
  of the indicators are between-level indicators of its group class,
  with class-specific means and variances.

  - `"restricted_cross_level"`: profile prevalences are the same in
    every group class, so profiles and group classes are estimated
    separately; the profile composition of each class is reported
    descriptively.

  - `"full_cross_level"`: profile prevalences differ across group
    classes (as in `"profiles"`), and the group means also inform the
    classes.

  These take `data`, `vars`, `id`, `n_profiles`, `n_group_classes`,
  `variance_model` (individual profiles), `between_variance` (group
  classes), `n_starts`, `max_iter`, `tol`, `min_variance`, `seed` and
  `weights`. The group means are computed from the same ratings, so the
  likelihood is the specification's working likelihood: comparable
  between cross-level fits of the same data, not with the other
  families. Standard errors are not available
  (`latents_unsupported_inference`). The result has class
  `multilpa_cross_level`, read with
  [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md).
  Its `get_results(fit, "profiles")` table ends with a `prevalence`
  column: each individual profile's prevalence (estimated and shared for
  the restricted family, the class-weighted marginal for the full
  family) and, on the group-class rows, the class weight.

- between_variance:

  For the additive and additive-dispersion families: `"varying"` (the
  default) estimates between-group variances per group class, `"equal"`
  one set shared by all classes. The dispersion family holds them equal
  and refuses `"varying"`. For the cross-level families, the variances
  of the group means within each group class.

- weights:

  `NULL`, or the name of a numeric column of `data` holding a sampling
  weight for each independent unit: each `id` group of a two-level fit,
  each row of a single-level one. A two-level weight must be constant
  within its group (`latents_bad_weights` otherwise); within-unit
  weights are not supported. The fit maximizes the pseudo log likelihood
  `sum_j w_j log L_j` (Skinner, 1989) with the weights scaled to sum to
  the number of units, as Mplus does, so integer weights give the fit to
  the data with each unit repeated that many times and the criteria stay
  on the sample's scale. Standard errors are sandwich estimates:
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  defaults to `vcov_type = "robust"` and refuses `"observed"` and
  `"opg"` (`latents_unsupported_weights`), and `method = "bootstrap"`
  resamples the units with their weights. Classification tables report
  each unit's own posterior; effective counts and proportions are
  weighted.
  [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md),
  [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md),
  [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md),
  `prior` and `noise` refuse a weighted fit. Weights work with every
  `family`, membership covariates, categorical indicators and
  `missing = "fiml"`.

- ordinal:

  Names of indicators in `vars` to model as ordinal: an (ordered)
  factor, or whole numbers, whose categories are ordered by factor level
  or by value. Each follows an adjacent-category logit with category
  intercepts shared by every profile and one location per profile,
  `log P(k | c) / P(k - 1 | c) = a_k - a_(k-1) + eta_c`, with the last
  profile's location fixed at zero: `(K - 1) + (C - 1)` parameters per
  indicator, against `C (K - 1)` for a `categorical` one. This is Latent
  GOLD's default ordinal model, against which it is checked. The
  category probabilities and locations are returned by
  `get_results(fit, "ordinal")`. A factor level never observed in the
  data raises the `latents_unused_category` warning;
  [`droplevels()`](https://rdrr.io/r/base/droplevels.html) removes it.

- count:

  Names of indicators in `vars` holding non-negative whole numbers,
  modelled as Poisson with one mean per profile (or negative binomial,
  see `count_model`); the means are returned by
  `get_results(fit, "count_means")`. Ordinal and count indicators work
  with `missing = "fiml"`, `weights`, two-level fits, membership
  covariates, Wald and bootstrap inference,
  [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md),
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md). `start`,
  `fixed`, `prior` and `noise` do not accept them
  (`latents_unsupported_indicator`), and the group-class and cross-level
  families refuse them as arguments they do not take
  (`latents_bad_argument`).

- count_model:

  `"poisson"` (the default) or `"negative_binomial"`: NB2, with mean
  `mu` and variance `mu + alpha mu^2`, for counts more variable than a
  Poisson within a profile (Latent GOLD's `poisson overdispersed`). A
  dispersion estimated at zero is the Poisson limit, a boundary fit
  (`latents_boundary` warning; Wald inference is refused there).

- count_dispersion:

  For `count_model = "negative_binomial"`: one dispersion per profile
  (`"varying"`, the default) or one shared by every profile (`"equal"`).

- scaling:

  How to scale the continuous indicators before fitting, at the same
  levels as `centering`. `"none"` (the default) fits them as supplied.
  `"grand"` divides each indicator by its standard deviation over all
  observations; `"person"` divides each group's values by that group's
  own standard deviation. With the matching centring these are z-scores:
  `centering = "grand", scaling = "grand"` standardizes each indicator
  over the sample, and `centering = "person", scaling = "person"`
  standardizes it within each unit. Standard deviations use the observed
  values with denominator `n - 1`. The scales are stored on the fit with
  the centring offsets, so `get_results(x, "data")` returns the supplied
  columns, and [`predict()`](https://rdrr.io/r/stats/predict.html)
  scales new units by their own spread under `"person"` and by the
  fitted spread under `"grand"`. Means, variances and covariances are
  reported in standard-deviation units. The log likelihood and every
  information criterion are those of the scaled values, so only fits
  with the same centring and scaling are comparable. A standard
  deviation that is zero or undefined (under `"person"`, a unit observed
  once on an indicator or always giving the same value) raises
  `latents_bad_data`. Under structures that give each indicator its own
  variance (`EEI`, `VVI`, `EEE`, `VVV`, `EVI`, `VEI`), grand scaling
  changes only the units of the estimates; under the spherical `EII` and
  `VII`, which impose one variance on every indicator, it changes the
  solution.

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
  individual sample size of a `multilpa()` fit, which excludes such rows
  either way.
  [`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)
  and
  [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
  do not compare fits made with different `drop_empty`. The setting is
  stored on the fit and passed on by enumerations and refits; it has no
  effect on data without such rows.

- cores:

  Number of processor cores for the random starts (default 1). The
  starting values are drawn in order in the current process and only the
  EM iterations are distributed, so the result is identical to
  `cores = 1`. Workers are forked processes on macOS and Linux, and a
  socket cluster on Windows and on macOS with Apple's Accelerate BLAS,
  as in
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md);
  a worker that fails as a whole raises `latents_parallel_failure`.
  Models with membership covariates run their starts on one core.

## Value

A `multilpa` object containing `means`, `variances`, optional
`covariances` (indicators by indicators by profiles),
`profile_probabilities`, `group_probabilities`, posterior matrices,
classifications, log likelihood, information criteria, restart
diagnostics and convergence history. Individual rows keep their input
order; groups keep first-occurrence order. `bic` and `bic_groups` use
the observed group count; `bic_individual` uses `n_informative`, the
number of rows with at least one observed indicator. These are
alternative conventions, not interchangeable criteria. `boundary`
identifies variance bounds; `small_classes` flags a profile or group
class whose effective membership (`effective_profile_counts`,
`effective_group_counts`, weighted for a fit with `weights`) is below
one, in which case the fit warns with `latents_small_classes`. The fit
carries no standard errors:
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
computes them, and
[`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
tests nested models. Global optimality is not guaranteed for any
`n_starts`.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the tidy tables, selected by `what`.

## Details

Infinite and constant observed indicators are rejected; each indicator
must have at least two distinct observed values. No rows are silently
dropped. Under FIML, fully missing individuals contribute no direct
measurement likelihood but receive posterior probabilities from their
group's information. Observation counts include these individuals; the
individual-level BIC excludes them through `n_informative`, while the
group-level BIC uses all observed groups. Missing patterns can prevent
parameter identification; no general identification guarantee is made.
Only initialization fills missing values with indicator means; EM uses
conditional Gaussian sufficient statistics and observed marginal
densities. More than one group class requires more than one profile and
at least one group with more than one individual; these checks are
necessary but do not establish identification. By default the highest
finite likelihood across starts is returned, even if that start did not
converge; `select_start = "converged"` prefers a converged start, and
`converged` and `starts` show which was chosen. Profile and group-class
labels are arbitrary.

## References

Vermunt, J. K. (2003). Multilevel latent class models. Sociological
Methodology, 33, 213–239. doi:10.1111/j.0081-1750.2003.t01-1-00131.x.

Banfield, J. D., & Raftery, A. E. (1993). Model-based Gaussian and
non-Gaussian clustering. Biometrics, 49, 803–821. doi:10.2307/2532201.

Houle, S. A., Morin, A. J. S., & Harvey, J.-F. (2026). Multilevel latent
profile analyses: A comprehensive guide. Organizational Research
Methods. doi:10.1177/10944281261469432.

Skinner, C. J. (1989). Domain means, regression and multivariate
analysis. In C. J. Skinner, D. Holt, & T. M. F. Smith (Eds.), Analysis
of complex surveys (pp. 59–87). Wiley.

Asparouhov, T. (2005). Sampling weights in latent variable modeling.
Structural Equation Modeling, 12, 411–434.
doi:10.1207/s15328007sem1203_4.

Fraley, C., & Raftery, A. E. (2007). Bayesian regularization for normal
mixture estimation and model-based clustering. Journal of
Classification, 24, 155–181. doi:10.1007/s00357-007-0004-5.

## See also

[`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
for a model in which every indicator is categorical.

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
#> Best likelihood replicated in 4/4 starts (absolute tolerance 0.00844).
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
#>        834.1     6.871 NA 6.374e-12  6.374e-12
#>        587.9     5.630 NA 1.799e-08  1.799e-08
#>        834.1     6.109 NA 1.001e-09  1.001e-09
#>        834.1     5.939 NA 2.867e-09  2.867e-09
#>        587.9     4.714 NA 2.429e-06  2.429e-06
#>        834.1     5.042 NA 4.619e-07  4.619e-07
#>        587.9     2.921 NA 3.493e-03  3.493e-03
#>        587.9     2.888 NA 3.883e-03  3.883e-03
#>        587.9     2.653 NA 7.983e-03  7.983e-03
#>        587.9     1.617 NA 1.059e-01  1.059e-01
#>    ... 10 more rows.  get_results(x, what = "residuals")
#> 
#> -- information_criteria --------------------------------------------
#>  log_likelihood n_parameters   aic   kic bic_groups bic_individual sabic_groups
#>           -8440           23 16926 16952      16987          17047        16914
#>  sabic_individual caic_groups caic_individual awe_groups awe_individual
#>             16974       17010           17070      17172          17404
#>  icl_groups icl_individual clc_groups clc_individual
#>       16996          17168      16889          17002
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
#>      2          -8440      TRUE         12  <NA>    FALSE
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
#> 19 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
as.data.frame(fit)
#>    profile  indicator       mean  variance standard_deviation
#> 1        1     browse  0.5395773 0.5898715          0.7680309
#> 2        1   lectures  0.4274030 0.8446923          0.9190714
#> 3        1 forum_read  0.6198248 0.4814743          0.6938835
#> 4        1 forum_post  0.5304887 0.6887896          0.8299335
#> 5        1 attendance  0.6490497 0.3841264          0.6197793
#> 6        2     browse -0.7655998 0.5284684          0.7269583
#> 7        2   lectures -0.6064852 0.5383733          0.7337392
#> 8        2 forum_read -0.8792138 0.3631129          0.6025885
#> 9        2 forum_post -0.7526883 0.4211623          0.6489702
#> 10       2 attendance -0.9206600 0.3736293          0.6112522
#>    mean_standard_error variance_standard_error
#> 1           0.02710160              0.02948402
#> 2           0.03248649              0.04210054
#> 3           0.02473030              0.02439627
#> 4           0.02953027              0.03481645
#> 5           0.02226482              0.01955858
#> 6           0.03121144              0.03221700
#> 7           0.03093059              0.03209025
#> 8           0.02621510              0.02279295
#> 9           0.02768122              0.02607600
#> 10          0.02642702              0.02278995
get_results(fit, what = "profile_probabilities")
#>   group_class profile probability group_class_probability
#> 1           1       1   0.1704511               0.3252973
#> 2           1       2   0.8295489               0.3252973
#> 3           2       1   0.7859319               0.6747027
#> 4           2       2   0.2140681               0.6747027
```
