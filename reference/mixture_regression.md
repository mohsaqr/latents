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
  family = c("gaussian", "binomial", "poisson"),
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
  vcov_type = c("observed", "robust", "opg", "none")
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

  `"gaussian"`, `"binomial"` or `"poisson"`.

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

  `NULL`, or a one-sided formula naming predictor terms whose
  coefficients are shared by every class, such as `~ age`. Each term
  must appear in `formula`.

- membership:

  One-sided formula of covariates predicting class membership (a
  multinomial logit, the first class as the reference). Row-level for
  `"observation"`; constant within groups for `"group"`.

- group_membership:

  One-sided formula of group-level covariates predicting the group
  class, for the two-level model. Must be constant within groups.

- variance:

  `"varying"` (a residual standard deviation per class) or `"equal"`.
  Gaussian family only.

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

## Value

An object of class `latents_mixture_regression`. Read it with
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

## Families

`"gaussian"` (identity link, a residual standard deviation per class, or
one shared under `variance = "equal"`), `"binomial"` (logit link; the
outcome is 0/1, logical, a two-level factor whose second level is the
success, or `cbind(successes, failures)`) and `"poisson"` (log link;
[`offset()`](https://rdrr.io/r/stats/offset.html) terms in the formula
are honoured).

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
assignment; `latents_missing_data` for missing values under
`missing = "error"`; `latents_rows_dropped` (warning) under
`missing = "omit"`; `latents_no_valid_start` when every start
degenerated; `latents_unconverged` (warning) when the selected start did
not converge; `latents_degenerate_start` (warning) when some starts
degenerated; `latents_separation` (warning) when a coefficient diverged.

## References

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
#>    class        term estimate std_error   p_value
#>  class_1 (Intercept)  34.3545    0.7193 0.000e+00
#>  class_1       hours   4.5659    0.1047 0.000e+00
#>  class_2 (Intercept)  54.8951    1.0927 0.000e+00
#>  class_2       hours   0.8084    0.1696 1.876e-06
#> 
#> Every table: get_results(x, what = ), e.g. "classes", "membership", "fit", "assignments".
as.data.frame(fit)
#>     class        term   estimate std_error statistic      p_value   conf_low
#> 1 class_1 (Intercept) 34.3545204 0.7192917 47.761595 0.000000e+00 32.9447345
#> 2 class_1       hours  4.5658711 0.1047310 43.596167 0.000000e+00  4.3606021
#> 3 class_2 (Intercept) 54.8951437 1.0927048 50.237855 0.000000e+00 52.7534818
#> 4 class_2       hours  0.8083711 0.1696013  4.766303 1.876367e-06  0.4759587
#>   conf_high   p_adjusted
#> 1 35.764306           NA
#> 2  4.771140 0.000000e+00
#> 3 57.036806           NA
#> 4  1.140784 1.876367e-06
get_results(fit, "classes")
#>     class     share    count n_assigned mean_posterior    sigma sigma_std_error
#> 1 class_1 0.5639412 507.5471        573      0.7884232 5.238049       0.2775751
#> 2 class_2 0.4360588 392.4529        327      0.8294171 6.890695       0.3456934
#>       prior prior_std_error
#> 1 0.5638921      0.03503844
#> 2 0.4361079      0.03503844

# One class per student, several rows per student:
by_student <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2,
                                 id = "student", class_level = "group",
                                 n_starts = 3, seed = 1)
get_results(by_student, "fit")
#>     family nesting n_classes n_group_classes n_observations n_groups
#> 1 gaussian   group         2               1            900      150
#>   log_likelihood n_parameters      aic      bic bic_rows    sabic      icl
#> 1      -3177.761            7 6369.521 6390.596 6403.138 6368.442 6451.323
#>     entropy group_entropy smallest_share converged iterations n_starts
#> 1 0.7079629            NA      0.4696097      TRUE         42        4
#>   n_best_replicated vcov_type
#> 1                 4  observed
```
