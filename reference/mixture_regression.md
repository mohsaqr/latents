# Estimate regression relationships that vary across latent classes

`mixture_regression()` asks whether the relationship between an outcome
and its predictors differs across latent classes. The conditional
outcome distribution is a mixture of class-specific regressions, so
classes can differ in outcome level, in predictor associations or in
trajectory over time. Regression parameters and class probabilities are
estimated jointly by maximum likelihood.

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

  Number of regression classes.

- family:

  `"gaussian"` (default), `"binomial"`, `"poisson"`,
  `"negative_binomial"` or `"ordinal"`. See *Families*.

- id:

  `NULL` (default), or the column identifying the units (persons,
  schools) the rows belong to.

- class_level:

  `"observation"` (default): each row has its own class. `"group"`: all
  rows of a unit share a class (needs `id`). See *Nesting*.

- n_group_classes:

  Number of second-level group classes (default 1, none). For
  `class_level = "observation"` with `id`, classes of the units: a
  unit's class sets the class probabilities of its rows. For a growth
  mixture with `cluster`, classes of the clusters.

- common:

  Predictor terms whose coefficients are shared by every class, as names
  (`"age"`) or a one-sided formula (`~ age`). Each term must appear in
  `formula`; the intercept cannot be shared.

- membership:

  Covariates of class membership, a multinomial logit with the first
  class as reference: names or a one-sided formula (`~ 1`, the default,
  for none). Row-level under `"observation"`, constant within units
  under `"group"`.

- group_membership:

  Covariates of the group class, as names or a one-sided formula;
  constant within the units that carry group classes (the `id` units, or
  the clusters under `cluster`).

- variance:

  `"varying"` (default): a residual standard deviation (or
  negative-binomial dispersion) per class; `"equal"`: one shared.
  Gaussian and negative-binomial families only.

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

  `"error"` (default) rejects rows with a missing value in any variable
  the model uses; `"omit"` drops them with a warning saying how many.

- select_start:

  `"likelihood"` (default): the start with the highest likelihood;
  `"converged"`: the best start that converged.

- vcov_type:

  Covariance of the estimates stored with the fit: `"observed"`
  (default, inverse observed information), `"robust"` (sandwich,
  clustered on the independent units), `"opg"` (outer product of the
  scores) or `"none"`. A weighted fit takes `"robust"` or `"none"` only.

- weights:

  Name of a numeric column holding a sampling weight per independent
  unit (per row without `id`, per unit with it): pseudo maximum
  likelihood with robust standard errors. A unit of weight zero is
  absent from the fit and every criterion but keeps its row and
  posterior.

- random:

  `NULL` (default), or the random effects within classes: names such as
  `"time"` (a random intercept and a random slope on time) or
  `"intercept"`, or a one-sided formula such as `~ 1 + time`
  (`~ 0 + time` for a slope alone). See *Growth mixture models*.

- random_covariance:

  How the random-effect covariance differs across classes: `"varying"`
  (default, one per class), `"equal"` (one shared matrix,
  `lcmm::hlme(nwg = FALSE)`) or `"proportional"` (a shared matrix times
  a class-specific scale, the last class's being one, `nwg = TRUE`).

- random_diagonal:

  `TRUE` for uncorrelated random effects (`idiag = TRUE` in lcmm);
  default `FALSE`.

- cluster:

  Column grouping the units of `id` into clusters, such as schools, for
  the multilevel growth mixture model (with `random`): each cluster
  belongs to one of `n_group_classes` group classes, which sets the
  trajectory-class probabilities of its persons. The clusters are the
  independent units of the likelihood, the standard errors and the BIC.
  See *Growth mixture models*.

- residual_correlation:

  `"none"` (default): independent residuals within a unit. With
  `random`, `"ar"` adds a continuous-time autoregressive process
  (covariance `sd^2 * exp(-rate * |t - s|)`) and `"brownian"` Brownian
  motion (`sd^2 * min(t, s)`) to the residuals, shared by the classes:
  `lcmm::hlme(cor = AR(time))` and `cor = BM(time)`. A process the data
  cannot separate from the residual variance or the random effects stays
  on its bound, warns, and has no standard error. These fits refine one
  class placement from the model without the process;
  `n_best_replicated` is `NA`, and `"starts"` reports the correlated fit
  rather than the initialization's replication count.

- residual_time:

  Numeric time column of the residual process; default the first numeric
  variable in `random`, otherwise in the predictors. Brownian motion
  needs non-negative times.

## Value

A `latents_mixture_regression` object, or `latents_growth_mixture` when
`random` is given.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the coefficient table and
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
every table by name;
[`predict()`](https://rdrr.io/r/stats/predict.html),
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
