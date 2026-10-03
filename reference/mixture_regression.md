# Fit a finite mixture of regressions

Fits a regression whose coefficients differ across latent classes: the
data are assumed to come from `n_classes` subpopulations, each with its
own regression of the outcome on the predictors, and the model estimates
the regressions, the class sizes and every row's posterior class
membership at once. This is clusterwise, latent class or mixture
regression (DeSarbo and Cron 1988; Wedel and DeSarbo 1995).

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
  cluster = NULL
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

  `"gaussian"`, `"binomial"`, `"poisson"`, `"negative_binomial"` or
  `"ordinal"` (proportional-odds cumulative logit; see *Families*).

- id:

  `NULL`, or the name of the column identifying groups (persons,
  schools, ...).

- class_level:

  `"observation"` (each row has its own class) or `"group"` (all rows of
  a group share a class; needs `id`).

- n_group_classes:

  Number of second-level group classes, for
  `class_level = "observation"` with `id`. `1` fits no second level.

- common:

  `NULL`, or the predictor terms whose coefficients are shared by every
  class: variable names such as `"age"`, or a one-sided formula such as
  `~ age`. Each term must appear in `formula`.

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

  `"varying"` (a residual standard deviation, or a negative-binomial
  dispersion, per class) or `"equal"` (one shared). Gaussian and
  negative-binomial families only.

- n_starts:

  Number of random starts, besides one start built from the residuals of
  a pooled regression. Every start runs 50 EM iterations; the better
  half (at least two) continue to convergence.

- max_iter:

  Maximum EM iterations per start.

- tol:

  Relative convergence tolerance on the log likelihood.

- min_variance:

  Floor for a Gaussian residual variance, relative to the outcome's
  variance. A class that reaches it is fitting too few rows exactly;
  such starts are set aside.

- seed:

  `NULL` or an integer seed. The caller's random state is restored
  afterwards.

- missing:

  `"error"` refuses rows with missing values in any variable the model
  uses; `"omit"` drops them with a `latents_rows_dropped` warning
  stating how many.

- select_start:

  `"likelihood"` keeps the start with the highest likelihood;
  `"converged"` prefers the best start that converged.

- vcov_type:

  Covariance of the estimates stored with the fit: `"observed"` (inverse
  observed information), `"robust"` (sandwich, clustered on the
  top-level unit, or on `id` for a single-level fit given `id`) or
  `"opg"` (outer product of the scores). `"none"` skips inference.

- weights:

  `NULL`, or the name of a numeric column of `data` with a sampling
  weight per independent unit: per row without `id`, per `id` group
  otherwise (constant within it). Pseudo maximum likelihood with the
  weights scaled to sum to the number of units; `vcov_type` defaults to
  `"robust"` and refuses `"observed"` and `"opg"`.

- random:

  `NULL` (no random effects), or the random effects within classes:
  variable names such as `"time"` (a random intercept and a random slope
  on `time`) or `"intercept"` (random intercepts only), or a one-sided
  formula such as `~ 1 + time` (needed for anything names cannot say,
  such as `~ 0 + time`, a slope without a random intercept). See the
  growth mixture section.

- random_covariance:

  How the random-effect covariance differs across classes: `"varying"`
  (one per class), `"equal"` (one shared matrix, as `lcmm::hlme()` with
  `nwg = FALSE`) or `"proportional"` (a shared matrix times a
  class-specific scale, the last class's being one; `nwg = TRUE`).

- random_diagonal:

  `TRUE` for uncorrelated random effects (a diagonal covariance, as
  `idiag = TRUE` in lcmm).

- cluster:

  `NULL`, or the name of a column grouping the persons (`id`) into
  clusters, such as schools, for the multilevel growth mixture model
  (with `random`): each cluster belongs to one of `n_group_classes`
  group classes, which shifts how probable each trajectory class is for
  its persons. The clusters are the independent units of the likelihood,
  the standard errors and the BIC; with `n_group_classes = 1` the model
  is the single-level growth mixture with clusters as those units,
  comparable by BIC with more group classes. See the growth mixture
  section.

## Value

An object of class `latents_mixture_regression`, or
`latents_growth_mixture` when `random` is given. Read it with
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) (the
coefficient table) or
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
(every table, by name), and use
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html) and
[`nobs()`](https://rdrr.io/r/stats/nobs.html) on it. The tables are
described on
[`get_results.latents_mixture_regression()`](https://pak.dynasite.org/latents/reference/get_results.latents_mixture_regression.md).

## Nesting

`class_level` and `n_group_classes` choose one of three likelihoods.

- Single level (`id = NULL`):

  Every row is its own unit and has its own class.

- `class_level = "group"`:

  Every row of a group (`id`) belongs to the same class, so a class is a
  kind of *group*: the regression mixture for repeated measures,
  `flexmix(y ~ x | id)`. Each group's evidence is the product of its
  rows' densities.

- `class_level = "observation"` with `id` and `n_group_classes >= 2`:

  Rows keep their own class, and a second-level group class shifts how
  probable each regression class is within its groups (Vermunt 2003).
  With `n_group_classes = 1` the groups only enter the `"robust"`
  standard errors, which are then clustered on `id`.

## Growth mixture models

With repeated measures of persons (`id`, `class_level = "group"`), a
time predictor in `formula` gives every class its own trajectory: latent
class growth analysis (Nagin 2005). `random` adds random effects
*within* classes, so persons scatter around their class trajectory: the
growth mixture model (Verbeke and Lesaffre 1996; Muthen and Shedden
1999). Within class k, person i's outcomes are \$\$y_i = X_i \beta_k +
Z_i b_i + e_i, \quad b_i \sim N(0, G_k), \quad e_i \sim N(0, \sigma_k^2
I),\$\$ with `Z` built from `random`. The likelihood is exact (Gaussian,
no numerical integration); estimation is EM with the random effects as
missing data, finished by quasi-Newton with analytic scores. Such a fit
has class `latents_growth_mixture`; its tables, plots and inference are
described on
[`get_results.latents_growth_mixture()`](https://pak.dynasite.org/latents/reference/get_results.latents_growth_mixture.md)
and
[`plot.latents_growth_mixture()`](https://pak.dynasite.org/latents/reference/plot.latents_growth_mixture.md).
Random effects need the gaussian family, `id` and
`class_level = "group"`.

With `cluster`, persons are nested in clusters (students in schools) and
the multilevel growth mixture model is fitted: each cluster belongs to a
group class h, and the trajectory-class logits of its persons have a
group-class-specific intercept and shared covariate slopes, \$\$L =
\prod_j \sum_h \omega_h \prod\_{i \in j} \sum_k \pi\_{k\|h} f_k(y_i)\$\$
(Asparouhov and Muthen 2008; Vermunt 2003). `group_membership` names
cluster-level covariates of the group class. The group classes are read
with `get_results(fit, "group_classes")` and `"clusters"`.

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
slope moves the class towards higher categories. The outcome is an
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
probabilities. With one outcome per class assignment the classes are
identified by how continuous predictors shift the category
probabilities: a mixture with no predictors, or with two categories
(then the logistic mixture of `"binomial"`), is refused with
`latents_not_identified` unless `class_level = "group"` gives each class
assignment several outcomes.

A mixture of Bernoulli regressions with one trial per unit is not
identified – any mixture of Bernoullis is again a Bernoulli – so a
binary outcome requires `class_level = "group"` (several rows per class
assignment) or binomial counts with more than one trial per row. The
request is refused otherwise.

## Conditions

`latents_bad_argument` for an inconsistent specification;
`latents_bad_data` for an outcome outside the family's support or a
covariate that varies within a group where it may not;
`latents_not_identified` for a binary outcome with one trial per class
assignment, or an ordinal mixture with one outcome per class assignment
and two categories or no predictors; `latents_missing_data` for missing
values under `missing = "error"`; `latents_rows_dropped` (warning) under
`missing = "omit"`; `latents_no_valid_start` when every start
degenerated; `latents_unconverged` (warning) when the selected start did
not converge; `latents_degenerate_start` (warning) when some starts
degenerated; `latents_separation` (warning) when a coefficient diverged.

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
#> 900 rows | log likelihood -3144.0074 | BIC 6335.63 | entropy 0.443
#> Converged: TRUE | iterations: 42 | best likelihood reached by 3 of 3 completed starts (4 run)
#> 
#> Classes
#> 
#> Class    Share  Expected  Rows assigned  Avg. posterior  Residual SD
#> -------  -----  --------  -------------  --------------  -----------
#> Class 1  0.564    507.55            573           0.788         5.24
#> Class 2  0.436    392.45            327           0.829         6.89
#> 
#> Regression coefficients (95% CI)
#> 
#> Class    Term       Estimate          95% CI      p
#> -------  ---------  --------  --------------  -----
#> Class 1  Intercept     34.35  [32.94, 35.76]  <.001
#> Class 1  hours          4.57  [ 4.36,  4.77]  <.001
#> Class 2  Intercept     54.90  [52.75, 57.04]  <.001
#> Class 2  hours          0.81  [ 0.48,  1.14]  <.001
#> 
#> Every table: get_results(x, what = ), e.g. "classes", "membership", "fit", "assignments".
as.data.frame(fit)
#> Regression coefficients (95% CI)
#> 
#> Class    Term       Estimate          95% CI      p
#> -------  ---------  --------  --------------  -----
#> Class 1  Intercept     34.35  [32.94, 35.76]  <.001
#> Class 1  hours          4.57  [ 4.36,  4.77]  <.001
#> Class 2  Intercept     54.90  [52.75, 57.04]  <.001
#> Class 2  hours          0.81  [ 0.48,  1.14]  <.001
get_results(fit, "classes")
#> Classes
#> 
#> Class    Share  Expected  Rows assigned  Avg. posterior  Residual SD
#> -------  -----  --------  -------------  --------------  -----------
#> Class 1  0.564    507.55            573           0.788         5.24
#> Class 2  0.436    392.45            327           0.829         6.89

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
#> ICL             6451.32
#> Entropy         0.708
#> Smallest class  47.0%
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
