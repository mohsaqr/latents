# Fit a mixture regression model in which the relationship between an outcome and its predictors differs across latent classes

`mixture_regression()` fits a mixture of class-specific regressions, to
ask whether the relationship between an outcome and its predictors
differs across latent classes. Classes can differ in outcome level, in
predictor associations or in trajectory over time. Regression parameters
and class probabilities are estimated jointly by maximum likelihood.

Predictors in `formula` describe variation within classes; covariates in
`membership` describe who belongs to which class. A coefficient may be
shared by every class (`common`). Class membership can belong to each
row or to every row of a unit (`class_level`); with a time predictor and
shared membership the classes are trajectories, and `random` adds
within-class variation around them (a growth mixture model). Outcomes
can be Gaussian, binomial, Poisson, negative binomial or ordinal.

## Usage

``` r
mixture_regression(
  formula,
  data,
  n_classes,
  family = c("gaussian", "binomial", "poisson", "negative_binomial", "ordinal"),
  id = NULL,
  class_level = c("observation", "group"),
  n_group_classes = 1L,
  common = NULL,
  membership = ~1,
  group_membership = ~1,
  variance = c("varying", "equal"),
  n_starts = 10L,
  max_iter = 1000L,
  tol = 1e-08,
  min_variance = 1e-06,
  seed = NULL,
  missing = c("error", "omit"),
  select_start = c("likelihood", "converged"),
  vcov_type = c("observed", "robust", "opg", "none"),
  weights = NULL,
  random = NULL,
  random_covariance = c("varying", "equal", "proportional"),
  random_diagonal = FALSE,
  cluster = NULL,
  residual_correlation = c("none", "ar", "brownian"),
  residual_time = NULL
)
```

## Arguments

- formula:

  A model formula `outcome ~ predictors`. Factors, interactions and
  [`offset()`](https://rdrr.io/r/stats/offset.html) are supported.

- data:

  A data frame.

- n_classes:

  Number of regression classes; a positive integer.

- family:

  The distribution of the outcome within a class:

  - `"gaussian"` (default): a continuous outcome, identity link;

  - `"binomial"`: a binary outcome or counts of successes and failures,
    logit link;

  - `"poisson"`: a count outcome, log link;

  - `"negative_binomial"`: a count outcome more dispersed than a
    Poisson, log link;

  - `"ordinal"`: an ordered categorical outcome, cumulative logit link.

  See *Families* for the accepted outcome codings.

- id:

  `NULL` (default), or the column identifying the units (persons,
  schools) the rows belong to.

- class_level:

  What a class is assigned to:

  - `"observation"` (default): each row has its own class, so two rows
    of the same unit can be in different classes;

  - `"group"`: all rows of a unit share one class, so the classes are
    kinds of units. It needs `id`.

  See *Nesting*.

- n_group_classes:

  Number of second-level group classes (default 1, none). For
  `class_level = "observation"` with `id`, classes of the units: a
  unit's class sets the class probabilities of its rows. For a growth
  mixture with `cluster`, classes of the clusters. More than one needs
  `id`.

- common:

  Predictor terms whose coefficients are shared by every class, as names
  (`"age"`) or a one-sided formula (`~ age`). Each term must appear in
  `formula`; the intercept cannot be shared. `NULL` (default) lets every
  coefficient differ between classes.

- membership:

  Covariates of class membership, a multinomial logit with the first
  class as reference: names or a one-sided formula (`~ 1`, the default,
  for none). Row-level under `"observation"`, constant within units
  under `"group"`.

- group_membership:

  Covariates of the group class, as names or a one-sided formula (`~ 1`,
  the default, for none); constant within the units that carry group
  classes (the `id` units, or the clusters under `cluster`).

- variance:

  For the Gaussian and negative-binomial families (an error when given
  with the others), whether the classes share their residual variation:

  - `"varying"` (default): each class has its own residual standard
    deviation (or negative-binomial dispersion);

  - `"equal"`: one value shared by all classes, so the classes differ
    only in their coefficients.

- n_starts:

  Number of random starts (default 10), besides one start built from the
  residuals of a pooled regression. Each start runs up to 50 EM
  iterations; the better half (at least two) continue to convergence.

- max_iter:

  Maximum EM iterations per start (default 1000).

- tol:

  Convergence tolerance (default `1e-8`), per row: EM stops when a step
  changes the log likelihood by no more than `tol` times the number of
  rows of positive weight.

- min_variance:

  Floor on a Gaussian residual variance, as a fraction of the outcome's
  variance (default `1e-6`). A start in which a class's residual
  variance reaches the floor is excluded from selection. In a growth
  mixture the same floor bounds each random effect's conditional
  variance; a random-effect variance at the floor is kept and warned
  about.

- seed:

  Optional integer seed; when given, the caller's random-number state is
  restored on exit.

- missing:

  What to do with rows that have a missing value in any variable the
  model uses:

  - `"error"` (default): stop, saying how many rows are affected;

  - `"omit"`: drop those rows and fit the complete ones, with a
    `latents_rows_dropped` warning saying how many were dropped.

- select_start:

  Which start is reported:

  - `"likelihood"` (default): the start with the highest log likelihood,
    whether or not it converged;

  - `"converged"`: the start with the highest log likelihood among those
    that converged, or among all starts when none converged.

  Starts in which a class emptied or a residual variance reached its
  floor are set aside under either rule.

- vcov_type:

  The covariance matrix of the estimates stored with the fit, from which
  the standard errors and confidence intervals are taken:

  - `"observed"` (default): the inverse of the observed information;

  - `"robust"`: the sandwich estimator, clustered on the independent
    units;

  - `"opg"`: the outer product of the scores;

  - `"none"`: no covariance matrix is stored, which saves the time of
    computing it during the fit. A table that needs standard errors then
    computes them when it is requested (from the observed information;
    the sandwich for a weighted fit).

  A weighted fit uses `"robust"` by default and accepts only `"robust"`
  or `"none"`.

- weights:

  Name of a numeric column holding a sampling weight per independent
  unit (per row without `id`, per unit with it): pseudo maximum
  likelihood with robust standard errors. A unit of weight zero is
  absent from the fit and every criterion but keeps its row and
  posterior. `NULL` (default) weights every unit equally.

- random:

  `NULL` (default), or the random effects within classes: names such as
  `"time"` (a random intercept and a random slope on time) or
  `"intercept"`, or a one-sided formula such as `~ 1 + time`
  (`~ 0 + time` for a slope alone). Giving `random` fits a growth
  mixture model; it needs the Gaussian family, `id` and
  `class_level = "group"`. See *Growth mixture models*.

- random_covariance:

  With `random`, how the covariance matrix of the random effects differs
  between classes:

  - `"varying"` (default): one matrix per class;

  - `"equal"`: one matrix shared by all classes
    (`lcmm::hlme(nwg = FALSE)`);

  - `"proportional"`: a shared matrix times a class-specific scale, the
    last class's scale being one (`nwg = TRUE`).

- random_diagonal:

  With `random`: `FALSE` (default) estimates the covariances between the
  random effects; `TRUE` sets them to zero, so the random effects are
  uncorrelated (`idiag = TRUE` in lcmm).

- cluster:

  Column grouping the units of `id` into clusters, such as schools, for
  the multilevel growth mixture model (with `random`): each cluster
  belongs to one of `n_group_classes` group classes, which sets the
  trajectory-class probabilities of its persons. The clusters are the
  independent units of the likelihood, the standard errors and the BIC.
  `NULL` (default) has no cluster level. Not combined with `weights`.
  See *Growth mixture models*.

- residual_correlation:

  With `random`, how the residuals of one unit are correlated over time:

  - `"none"` (default): independent residuals;

  - `"ar"`: a continuous-time autoregressive process is added, with
    covariance `sd^2 * exp(-rate * |t - s|)` between times `t` and `s`
    (`lcmm::hlme(cor = AR(time))`);

  - `"brownian"`: Brownian motion is added, with covariance
    `sd^2 * min(t, s)` (`cor = BM(time)`).

  The process is shared by the classes. One the data cannot separate
  from the residual variance or the random effects stays on its bound,
  warns, and has no standard error. These fits refine one class
  placement from the model without the process; `n_best_replicated` is
  `NA`, and `"starts"` reports the correlated fit rather than the
  initialization's replication count.

- residual_time:

  Numeric time column of the residual process; default the first numeric
  variable in `random`, otherwise in the predictors. Brownian motion
  needs non-negative times. Given without `residual_correlation`, it is
  an error.

## Value

A fitted model of class `latents_mixture_regression`, or of class
`latents_growth_mixture` when `random` is given. Printing it shows the
fit, the classes and the coefficients.
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
returns its tables by name, among them:

- `"coefficients"` (also `as.data.frame(fit)`): one row per class and
  term, with the estimate, standard error, confidence limits and
  p-value;

- `"classes"`: the share of each class, the number of rows assigned to
  it, its average posterior probability and its residual standard
  deviation;

- `"membership"`: the coefficients of the membership covariates;

- `"fit"`: the log likelihood, number of parameters, AIC, BIC, SABIC,
  ICL and entropy;

- `"assignments"`: each row's (or unit's) most likely class and
  posterior probabilities;

- `"starts"`: the log likelihood and convergence of every start.

A regression mixture also has `"fitted"`, `"classification"`, `"groups"`
and `"group_classes"`; a growth mixture has `"random"` (the
random-effect variances), `"random_effects"`, `"trajectories"` and
`"individual"`. [`predict()`](https://rdrr.io/r/stats/predict.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) (and
[`plot_trajectories()`](https://pak.dynasite.org/latents/reference/plot_trajectories.md),
[`plot_coefficients()`](https://pak.dynasite.org/latents/reference/plot_coefficients.md)),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html) and
[`nobs()`](https://rdrr.io/r/stats/nobs.html) apply. The within-class
coefficients are conditional associations under the specified model, not
causal effects; check the starts table and the classification
diagnostics before reading substance into the classes.

## Conditions

Errors: `latents_bad_argument` for an argument value or combination that
is not accepted (for example `class_level = "group"` without `id`,
`variance` with a family that has none, or `random` without
`class_level = "group"`); `latents_bad_data` for a variable that `data`
does not have or an outcome the family cannot take;
`latents_missing_data` for missing values under `missing = "error"`;
`latents_not_identified` for the models listed under *Families*;
`latents_no_valid_start` when a class emptied or a variance collapsed in
every start; `latents_unsupported_weights` for a weighted fit asked for
a `vcov_type` other than `"robust"` or `"none"`, or combined with
`cluster`.

Warnings: `latents_unconverged` when the reported start did not
converge; `latents_degenerate_start` when some starts were set aside;
`latents_rows_dropped` for rows left out under `missing = "omit"`;
`latents_random_boundary` for a random-effect variance at its floor.

## Nesting

- `id = NULL`:

  Each row is an independent unit.

- `class_level = "group"`:

  All rows of a unit share a class, the regression mixture for repeated
  measures (`flexmix(y ~ x | id)`); a unit's likelihood is the product
  of its rows' densities.

- `class_level = "observation"` with `id` and `n_group_classes >= 2`:

  Each row has its own class, and a group class of its unit sets the
  class probabilities. With `n_group_classes = 1`, `id` clusters the
  robust standard errors and is the unit of `weights`.

## Growth mixture models

With repeated measurements (`id`, `class_level = "group"`) and time in
`formula`, the classes are trajectories: without random effects a latent
class growth model, with `random` a growth mixture model. For unit `i`
in class `k`,

\$\$y_i = X_i \beta_k + Z_i b_i + e_i, \quad b_i \sim N(0, G_k), \quad
e_i \sim N(0, \sigma_k^2 I),\$\$

with `Z` the design of `random`. Random effects need the Gaussian
family, `id`, `class_level = "group"` and, for a random intercept,
repeated rows in at least some units. With `residual_correlation`, `e_i`
also carries the autoregressive or Brownian process over
`residual_time`.

With `cluster`, persons are nested in clusters and each cluster belongs
to a group class `h` that sets the trajectory-class logits of its
persons (a class-specific intercept with shared covariate slopes):

\$\$L = \prod_j \sum_h \omega_h \prod\_{i \in j} \sum_k \pi\_{k\|h}
f_k(y_i)\$\$

`group_membership` names cluster-level covariates of the group class;
the group classes are in `get_results(fit, "group_classes")` and
`"clusters"`. With `n_group_classes = 1` the fit is comparable by BIC
with fits that add group classes.

## Families

`"gaussian"`: identity link, a residual standard deviation per class or
one shared. `"binomial"`: logit link; the outcome is 0/1, logical, a
two-level factor (second level the success) or
`cbind(successes, failures)`. `"poisson"`: log link,
[`offset()`](https://rdrr.io/r/stats/offset.html) honoured.
`"negative_binomial"`: log link, NB2 (variance `mu + alpha * mu^2`) with
a dispersion per class or one shared; a dispersion estimated at zero is
the Poisson limit, reported with a warning and no standard error.

`"ordinal"`: a cumulative-logit (proportional-odds) regression within
each class, \\P(Y \le c \mid x, k) = \mathrm{plogis}(t\_{kc} -
x'\beta_k)\\, with ordered thresholds per class in place of the
intercept and the sign convention of
[`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html): a positive
slope raises the probability of a higher category. The outcome is an
ordered factor, a factor (levels in order) or whole numbers; the
categories are those observed, at least two. The formula keeps its
intercept, which the thresholds replace. Thresholds appear in the
coefficient table as `"threshold:<lower>|<upper>"`, and `exp_estimate`
of a slope is a cumulative odds ratio. The class mean of an ordinal
outcome (fitted values, trajectories) is the expected category score
with categories scored 1 to C; `predict(type = "probabilities")` gives
the category probabilities.

Identification. Call the thing that carries a class the unit: a row, or
an `id` group under `class_level = "group"`. An ordinal mixture whose
units have one outcome each and no predictors is not identified (a
mixture of categorical distributions is itself one) and is refused. A
binary mixture with one trial per unit, and an ordinal one with one
outcome per unit and two categories, are refused too: whether such a
model is identified depends on the covariate design, which the package
does not assess. Give each unit several outcomes
(`class_level = "group"`) or, for binomial, counts with more than one
trial per row. A single-level ordinal mixture is also refused when it
has more free parameters than the distinct covariate patterns can
identify. Repeated binary or ordinal outcomes without predictors must
also have enough distinct response-count patterns to identify the class
parameters. A growth model is refused when its random-effects and
residual covariance components cannot be separated by the observed
designs. With class-specific residual variances, two common waves do not
identify an unrestricted random intercept and slope covariance. A
diagonal covariance can be identified there, as can proportional class
covariances with a shared residual variance.

Unused factor levels in membership and random-effects formulas are
dropped before fitting. In output tables, an identifier whose name
conflicts with a result column is prefixed with `id_` or `cluster_` (for
example, an `id = "class"` column is returned as `id_class`, beside the
latent `class`). A numeric suffix is added if that name is also
occupied.

Every refusal and warning is a classed condition;
[latents-conditions](https://pak.dynasite.org/latents/reference/latents-conditions.md)
lists them.

## See also

[`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md)
to compare numbers of classes;
[`plot_trajectories()`](https://pak.dynasite.org/latents/reference/plot_trajectories.md),
[`plot_coefficients()`](https://pak.dynasite.org/latents/reference/plot_coefficients.md)
and, for fits that classify persons,
[`plot_classification()`](https://pak.dynasite.org/latents/reference/plot_classification.md).

## Examples

``` r
fit <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2,
                          n_starts = 3, seed = 1)
fit
#> Mixture of gaussian regressions: 2 classes (one class per row)
#> score ~ hours
#> 900 rows | log likelihood -3144.0073 | BIC 6335.63 | entropy 0.443
#> Converged: TRUE | iterations: 53 | best likelihood reached by 2 of 2 completed starts (4 run)
#> 
#> Classes
#> 
#> Class    Share  Expected  Rows assigned  Avg. posterior  Residual SD
#> -------  -----  --------  -------------  --------------  -----------
#> Class 1  0.564    507.71            573           0.789         5.24
#> Class 2  0.436    392.29            327           0.829         6.89
#> 
#> Regression coefficients (95% CI)
#> 
#> Class    Term       Estimate          95% CI      p
#> -------  ---------  --------  --------------  -----
#> Class 1  Intercept     34.36  [32.95, 35.77]  <.001
#> Class 1  hours          4.57  [ 4.36,  4.77]  <.001
#> Class 2  Intercept     54.90  [52.76, 57.04]  <.001
#> Class 2  hours          0.81  [ 0.48,  1.14]  <.001
#> 
#> Every table: get_results(x, what = ), e.g. "classes", "membership", "fit", "assignments".
as.data.frame(fit)
#> Regression coefficients (95% CI)
#> 
#> Class    Term       Estimate          95% CI      p
#> -------  ---------  --------  --------------  -----
#> Class 1  Intercept     34.36  [32.95, 35.77]  <.001
#> Class 1  hours          4.57  [ 4.36,  4.77]  <.001
#> Class 2  Intercept     54.90  [52.76, 57.04]  <.001
#> Class 2  hours          0.81  [ 0.48,  1.14]  <.001
get_results(fit, "classes")
#> Classes
#> 
#> Class    Share  Expected  Rows assigned  Avg. posterior  Residual SD
#> -------  -----  --------  -------------  --------------  -----------
#> Class 1  0.564    507.71            573           0.789         5.24
#> Class 2  0.436    392.29            327           0.829         6.89

# One class per student, several rows per student:
by_student <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2,
                                 id = "student", class_level = "group",
                                 n_starts = 3, seed = 1)
get_results(by_student, "fit")
#> Model fit
#> 
#> Family          gaussian
#> Nesting         group
#> Classes         2
#> Group classes   1
#> Observations    900
#> Groups          150
#> Parameters      7
#> Log likelihood  -3177.76
#> AIC             6369.52
#> BIC             6390.60
#> SABIC           6368.44
#> ICL             6451.31
#> Entropy         0.708
#> Smallest class  46.9%
#> Converged       yes

# \donttest{
# Students nested in schools, each school in one of two group classes:
schools <- mixture_regression(score ~ wave, growth_schools, n_classes = 2,
                              id = "student", class_level = "group",
                              random = "wave", random_covariance = "equal",
                              cluster = "school", n_group_classes = 2,
                              n_starts = 2, seed = 1)
get_results(schools, "group_classes")
#> Group classes: trajectory classes within each
#> 
#> Group class    Class    Probability     SE  Share of clusters  Clusters
#> -------------  -------  -----------  -----  -----------------  --------
#> Group class 1  Class 1        0.842  0.024              0.544        22
#> Group class 1  Class 2        0.158  0.024              0.544        22
#> Group class 2  Class 1        0.233  0.030              0.456        18
#> Group class 2  Class 2        0.767  0.030              0.456        18
# }
```
