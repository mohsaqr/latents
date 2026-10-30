# Estimate regression relationships that vary across latent classes

`mixture_regression()` examines whether the relationship between an
outcome and its predictors differs across latent classes. It represents
the conditional outcome distribution as a mixture of class-specific
regression models, allowing classes to differ in their outcome levels,
predictor associations, or trajectories over time.

Regression parameters and class-membership probabilities are estimated
jointly by maximum likelihood. Coefficients may vary across classes or,
for selected predictor terms, be constrained to a common value.
Predictors in the outcome formula describe variation within classes;
covariates supplied through `membership` describe variation in the
probability of belonging to each class. Posterior membership
probabilities express the uncertainty in assigning observations or
groups to the estimated classes.

The function supports Gaussian, binomial, Poisson, negative-binomial,
and ordinal outcome models, subject to the restrictions described below.
Class membership may be defined for individual observations or shared
across repeated observations from the same unit. With time predictors,
shared membership defines class-specific trajectories; Gaussian random
effects additionally allow individual trajectories to vary within
classes.

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
  `"negative_binomial"` or `"ordinal"` (proportional-odds cumulative
  logit; see *Families*).

- id:

  `NULL` (default), or the name of the column identifying groups
  (persons, schools, ...).

- class_level:

  `"observation"` (default; each row has its own class) or `"group"`
  (all rows of a group share a class; requires `id`).

- n_group_classes:

  Number of second-level group classes (default 1, no latent group
  classes), for `class_level = "observation"` with `id`, or for the
  clusters of a multilevel growth mixture (`cluster`).

- common:

  `NULL`, or the predictor terms whose coefficients are shared by every
  class: variable names such as `"age"`, or a one-sided formula such as
  `~ age`. Each term must appear in `formula`. The intercept cannot be
  shared: `"intercept"` or `~ 1` alone raises `latents_bad_argument`.

- membership:

  Covariates predicting class membership (a multinomial logit, the first
  class as the reference): variable names such as `"age"`, or a
  one-sided formula (`~ 1`, the default, for none). Row-level for
  `"observation"`; constant within groups for `"group"`.

- group_membership:

  Group-level covariates predicting the group class, for the two-level
  model, as names or a one-sided formula. Must be constant within
  groups.

- variance:

  `"varying"` (default; a residual standard deviation, or a
  negative-binomial dispersion, per class) or `"equal"` (one shared).
  Gaussian and negative-binomial families only.

- n_starts:

  Number of random starts (default 10), besides one start built from the
  residuals of a pooled regression. Each start runs for up to 50 EM
  iterations; the better half (at least two) continue to convergence.

- max_iter:

  Maximum EM iterations per start (default 1000).

- tol:

  Convergence tolerance (default `1e-8`): EM has converged when one step
  changes the log likelihood by no more than `tol` times the number of
  positive-weight rows. The tolerance is therefore expressed per
  observation.

- min_variance:

  Floor for a Gaussian residual variance, relative to the outcome's
  variance (default `1e-6`). Starts in which a class reaches the lower
  bound are excluded from selection. In a growth mixture, the same floor
  bounds each random effect's conditional variance, applied in the
  standardized coordinates the model is fitted in (outcome and
  predictors centred and scaled), so it does not depend on the units of
  the time variable.

- seed:

  `NULL` (default) advances the caller's random stream. An integer seed
  restores the caller's `.Random.seed` afterwards. R does not store a
  cached Box-Muller normal draw in `.Random.seed`, so that hidden cached
  draw cannot be preserved.

- missing:

  `"error"` (default) rejects rows with missing values in any variable
  the model uses; `"omit"` drops them with a `latents_rows_dropped`
  warning stating how many.

- select_start:

  `"likelihood"` (default) selects the start with the highest
  likelihood; `"converged"` prefers the best start that converged.

- vcov_type:

  Covariance of the estimates stored with the fit: `"observed"`
  (default; inverse observed information), `"robust"` (sandwich with the
  CR1 factor `G / (G - 1)`, clustered on the top-level unit, or on `id`
  for a single-level fit given `id`) or `"opg"` (outer product of the
  scores). `"none"` skips inference.

- weights:

  `NULL`, or the name of a numeric column of `data` with a sampling
  weight per independent unit: per row without `id`, per `id` group
  otherwise (constant within it). Pseudo maximum likelihood with the
  weights scaled to sum to the units of positive weight; `vcov_type`
  defaults to `"robust"` and rejects `"observed"` and `"opg"`. The
  information criteria use the weighted log likelihood, and the entropy
  term of ICL is weighted by the sampling weights, as the likelihood is.
  A unit of weight zero is absent from the fit: from every start,
  estimate, sample size and criterion; it keeps its row and posterior. A
  row of weight zero whose factor level (or ordinal outcome category) no
  row of positive weight has is left out, with a `latents_rows_dropped`
  warning, since the model has no parameter for that level.

- random:

  `NULL` (default; no random effects), or the random effects within
  classes: variable names such as `"time"` (a random intercept and a
  random slope on `time`) or `"intercept"` (random intercepts only), or
  a one-sided formula such as `~ 1 + time`. A slope without a random
  intercept needs the formula `~ 0 + time`. See *Growth mixture models*.

- random_covariance:

  How the random-effect covariance differs across classes: `"varying"`
  (default; one per class), `"equal"` (one shared matrix, as
  `lcmm::hlme()` with `nwg = FALSE`) or `"proportional"` (a shared
  matrix times a class-specific scale, the last class's being one;
  `nwg = TRUE`).

- random_diagonal:

  `TRUE` for uncorrelated random effects (a diagonal covariance, as
  `idiag = TRUE` in lcmm); default `FALSE`.

- cluster:

  `NULL` (default), or the name of a column grouping the persons (`id`)
  into clusters, such as schools, for the multilevel growth mixture
  model (with `random`): each cluster belongs to one of
  `n_group_classes` group classes, which sets the trajectory-class
  probabilities of its persons. The clusters are the independent units
  of the likelihood, the standard errors and the BIC; with
  `n_group_classes = 1` the model is the single-level growth mixture
  with clusters as those units, comparable by BIC with fits with
  additional group classes. See *Growth mixture models*.

- residual_correlation:

  `"none"` (default) for independent residuals within a person; with
  `random`, `"ar"` adds a continuous-time autoregressive process,
  covariance `sd^2 exp(-rate |t - s|)`, and `"brownian"` Brownian
  motion, `sd^2 min(t, s)`, to the residual variance, shared by the
  classes: `lcmm::hlme(cor = AR(time))` and `cor = BM(time)`. The model
  without the process is fitted first from every start; the process is
  then estimated jointly with all other parameters. A process the data
  cannot separate from the residual variance or the random effects is
  left on its bound and warned about (`latents_random_boundary`); its
  rows then have no standard error. Refused with `latents_bad_argument`
  without `random`.

- residual_time:

  Name of the numeric time column for the residual process. Defaults to
  the first numeric variable in `random`, or otherwise the predictors.
  Brownian motion requires non-negative times (`latents_bad_data`).

## Value

An object of class `latents_mixture_regression`, or
`latents_growth_mixture` when `random` is given.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the coefficient table and
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
every table by name (described there);
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html) and
[`nobs()`](https://rdrr.io/r/stats/nobs.html) also apply.

## Interpretation

The classes describe heterogeneity under the specified regression model.
Their interpretation depends on the outcome distribution, included
predictors, coefficient constraints, and treatment of repeated
observations. Within-class coefficients describe conditional
associations and do not, by themselves, establish causal effects.
Estimation stability and classification uncertainty should be examined
before assigning substantive meaning to the classes.

## Nesting

Arguments `class_level` and `n_group_classes` specify the level of
latent class membership and the likelihood contribution of each
independent unit.

- Single level (`id = NULL`):

  Each row is an independent unit with its own class.

- `class_level = "group"`:

  All rows of a group (`id`) belong to the same class: the regression
  mixture for repeated measures, as `flexmix(y ~ x | id)`. A group's
  likelihood contribution is the product of its rows' densities.

- `class_level = "observation"` with `id` and `n_group_classes >= 2`:

  Each row has its own class, and a second-level group class determines
  the class probabilities within its groups (Vermunt 2003). With
  `n_group_classes = 1`, `id` enters only the `"robust"` standard
  errors, which are then clustered on `id`.

## Growth mixture models

For repeated measurements (`id`, `class_level = "group"`), including
time in `formula` specifies class-specific trajectories. Without random
effects, this is a latent class growth model (Nagin, 2005). Argument
`random` adds within-class variation in individual trajectories,
yielding a growth mixture model (Verbeke and Lesaffre, 1996; Muthen and
Shedden, 1999). For individual i in class k, the model is \$\$y_i = X_i
\beta_k + Z_i b_i + e_i, \quad b_i \sim N(0, G_k), \quad e_i \sim N(0,
\sigma_k^2 I),\$\$ where `Z` is the design matrix specified by `random`.
The Gaussian random effects are integrated analytically. Estimation uses
EM followed by quasi-Newton optimization with analytic scores. The
fitted object has class `latents_growth_mixture`; its tables and plots
are described on
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
and
[plot()](https://pak.dynasite.org/latents/reference/latents-plots.md).
Random effects require the gaussian family, `id` and
`class_level = "group"`.

With `residual_correlation`, the residuals also carry a process shared
by the classes, \\e_i \sim N(0, W_i + \sigma_k^2 I)\\: a continuous-time
autoregressive process (`"ar"`) or Brownian motion (`"brownian"`) over
`residual_time`, corresponding to the correlation structures specified
by `lcmm::hlme(cor = )`.

With `cluster`, persons are nested in clusters (students in schools) and
the multilevel growth mixture model is fitted: each cluster belongs to a
group class h, and the trajectory-class logits of its persons have a
group-class-specific intercept and shared covariate slopes, \$\$L =
\prod_j \sum_h \omega_h \prod\_{i \in j} \sum_k \pi\_{k\|h} f_k(y_i)\$\$
(Asparouhov and Muthen 2008; Vermunt 2003). `group_membership` names
cluster-level covariates of the group class. The group classes are
returned by `get_results(fit, "group_classes")` and `"clusters"`.

## Families

`"gaussian"` (identity link, a residual standard deviation per class, or
one shared under `variance = "equal"`), `"binomial"` (logit link; the
outcome is 0/1, logical, a two-level factor whose second level is the
success, or `cbind(successes, failures)`), `"poisson"` (log link;
[`offset()`](https://rdrr.io/r/stats/offset.html) terms in the formula
are honoured) and `"negative_binomial"` (log link, NB2: variance mu +
alpha mu^2 with a dispersion alpha per class, or one shared under
`variance = "equal"`, for overdispersed counts; a dispersion estimated
at zero is the Poisson limit and raises `latents_boundary`).

`"ordinal"` fits a cumulative-logit (proportional-odds) regression
within each class (McCullagh 1980), \$\$P(Y \le c \mid x, k) =
F(t\_{kc} - x'\beta_k), \quad c = 1, \ldots, C - 1,\$\$ with ordered
thresholds \\t\_{k1} \< \ldots \< t\_{k,C-1}\\ per class in place of the
intercept, and the parametrization of
[`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html): a positive
slope increases the probability of higher categories. The outcome is an
ordered factor, a factor (its levels taken in order) or whole-number
categories; the categories are those observed, at least two. `common`
terms share their slopes across classes, and
[`offset()`](https://rdrr.io/r/stats/offset.html) is honoured; the
formula's intercept is replaced by the thresholds, so it must not be
removed. The thresholds appear in the coefficient table as terms
`"threshold:<lower>|<upper>"`, and `exp_estimate` of a slope is the
cumulative odds ratio of a higher category. The class "mean" of an
ordinal outcome (fitted values,
[`predict()`](https://rdrr.io/r/stats/predict.html), trajectories) is
the expected category score \\\sum_c c P(Y = c)\\, the categories scored
1 to C in order; `predict(type = "probabilities")` gives the category
probabilities. With one outcome per class assignment, an ordinal mixture
with no predictors is not identified. Two-category mixtures follow the
binary support restriction below. Both are rejected with
`latents_not_identified`, unless `class_level = "group"` gives each
class assignment several outcomes. A single-level mixture is also
rejected when it has more free parameters than P (C - 1), the number of
free category probabilities in P distinct covariate patterns (for
example two classes and one binary factor).

Binary mixtures with one trial per class assignment are not supported:
whether such a model is identified depends on the covariate design and
the model constraints, which the package does not assess. A binary
mixture therefore requires `class_level = "group"` (several rows per
class assignment) or binomial counts with more than one trial per row;
other binary specifications are rejected.

## Conditions

`latents_bad_argument` for an inconsistent specification (including a
`common` that names no predictor term); `latents_bad_data` for an
outcome outside the family's support, a covariate that varies within a
group where it may not, or a factor predictor with a single value among
the fitted (positive-weight) rows; `latents_not_identified` for a binary
outcome with one trial per class assignment, or an ordinal mixture with
one outcome per class assignment and two categories, no predictors, or
fewer distinct covariate patterns times C - 1 than free parameters;
`latents_missing_data` for missing values under `missing = "error"`;
`latents_rows_dropped` (warning) under `missing = "omit"`, and for rows
of weight zero carrying a factor level or outcome category that no row
of positive weight has; `latents_no_valid_start` when every start
degenerated; `latents_unconverged` (warning) when the selected start did
not converge; `latents_degenerate_start` (warning) when some starts
degenerated; `latents_separation` (warning) for extreme fitted logits or
log means that may signal separation; `latents_failed_starts` (warning)
when some growth-mixture starts could not be built (they are marked
`"failed"` in the starts table); `latents_boundary` (warning) for a
negative-binomial dispersion at the Poisson limit, or an ordinal
threshold gap at zero (a class that gives no row some middle category;
that gap and the thresholds above it have no standard error);
`latents_random_boundary` (warning) for a degenerate random-effect
covariance; `latents_singular_information` (warning) when some standard
errors cannot be computed.

## References

Asparouhov, T., & Muthen, B. (2008). Multilevel mixture models. In G. R.
Hancock & K. M. Samuelsen (Eds.), *Advances in latent variable mixture
models* (pp. 27–51). Information Age.

McCullagh, P. (1980). Regression models for ordinal data. *Journal of
the Royal Statistical Society B*, 42, 109–142.

Muthen, B., & Shedden, K. (1999). Finite mixture modeling with mixture
outcomes using the EM algorithm. *Biometrics*, 55, 463–469.

Nagin, D. S. (2005). *Group-based modeling of development*. Harvard
University Press.

Proust-Lima, C., Philipps, V., & Liquet, B. (2017). Estimation of
extended mixed models using latent classes and latent processes: the R
package lcmm. *Journal of Statistical Software*, 78(2), 1–56.

Verbeke, G., & Lesaffre, E. (1996). A linear mixed-effects model with
heterogeneity in the random-effects population. *Journal of the American
Statistical Association*, 91, 217–221.

DeSarbo, W. S., & Cron, W. L. (1988). A maximum likelihood methodology
for clusterwise linear regression. *Journal of Classification*, 5,
249–282.

Wedel, M., & DeSarbo, W. S. (1995). A mixture likelihood approach for
generalized linear models. *Journal of Classification*, 12, 21–55.

Follmann, D. A., & Lambert, D. (1991). Identifiability of finite
mixtures of logistic regression models. *Journal of Statistical Planning
and Inference*, 27, 375–381.

Vermunt, J. K. (2003). Multilevel latent class models. *Sociological
Methodology*, 33, 213–239.

Leisch, F. (2004). FlexMix: A general framework for finite mixture
models and latent class regression in R. *Journal of Statistical
Software*, 11(8).

## See also

[`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md)
to compare numbers of classes.

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
