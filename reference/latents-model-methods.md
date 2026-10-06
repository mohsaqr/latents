# Coefficients, covariance, intervals and likelihood of latents fits

[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html) and
[`nobs()`](https://rdrr.io/r/stats/nobs.html) for multilevel profile and
class models
([`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)),
covariate models, latent transition models and, for
[`vcov()`](https://rdrr.io/r/stats/vcov.html), pooled imputations.
Mixture regressions, growth mixtures and the additive, cross-level and
general transition families document the same methods with their
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
page.

## Usage

``` r
# S3 method for class 'multilpa_covariates'
vcov(
  object,
  data = NULL,
  step = 1e-04,
  vcov_type = c("observed", "robust", "opg"),
  scale = c("natural", "unconstrained"),
  boundary = c("error", "fix"),
  ...
)

# S3 method for class 'multilpa_covariates'
coef(object, scale = c("natural", "unconstrained"), ...)

# S3 method for class 'multilpa_covariates'
confint(object, parm, level = 0.95, data = NULL, ...)

# S3 method for class 'multilpa_covariates'
logLik(object, ...)

# S3 method for class 'multilpa_covariates'
nobs(object, ...)

# S3 method for class 'multilpa'
coef(object, scale = c("natural", "unconstrained"), ...)

# S3 method for class 'multilpa'
vcov(object, data = NULL, scale = c("natural", "unconstrained"), ...)

# S3 method for class 'multilpa'
confint(object, parm, level = 0.95, data = NULL, ...)

# S3 method for class 'multilpa'
logLik(object, ...)

# S3 method for class 'multilpa'
nobs(object, ...)

# S3 method for class 'latents_pooled'
vcov(object, ...)

# S3 method for class 'multilpa_transitions'
logLik(object, ...)

# S3 method for class 'multilpa_transitions'
nobs(object, ...)

# S3 method for class 'multilpa_transitions'
coef(object, ...)
```

## Arguments

- object:

  A fitted model: from
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  (class `multilpa`), a covariate model (`multilpa_covariates`), a
  latent transition model (`multilpa_transitions`), or pooled
  imputations (`latents_pooled`,
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) only).

- data:

  Optional. The data frame the model was fitted to. When omitted it is
  rebuilt from what the fit stores (indicators, identifiers, occasions
  and designs), which round-trips exactly; supplying it is the stronger
  check that the caller still holds that frame.

- step:

  Finite-difference step for the observed information (covariate
  models).

- vcov_type:

  `"observed"`, `"robust"` or `"opg"`, as for
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  (covariate models).

- scale:

  `"natural"`, the default, or `"unconstrained"`. Natural coefficients
  are in the units [`coef()`](https://rdrr.io/r/stats/coef.html)
  reports: means and variances in the indicators' units, probabilities
  as probabilities, with the covariance carried from the estimation
  scale by the delta method. Unconstrained is the scale the model is
  estimated on: log variances (diagonal), log-Cholesky coordinates (full
  covariance), log volumes, shapes and orientations (a structure that
  constrains them across profiles), and baseline-category logits for
  mixing and response probabilities. In a covariate model, means and
  membership coefficients are the same on both scales.

- boundary:

  `"error"` or `"fix"`, as for
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md):
  `"fix"` holds response probabilities on their bound, whose rows and
  columns are then zero (covariate models).

- ...:

  For [`vcov()`](https://rdrr.io/r/stats/vcov.html) and
  [`confint()`](https://rdrr.io/r/stats/confint.html) of a `multilpa`
  fit, and [`confint()`](https://rdrr.io/r/stats/confint.html) of a
  covariate model, passed to
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  (for example `method = "bootstrap"`, `vcov_type = "robust"` or
  `step`). Ignored otherwise.

- parm:

  Optional coefficient names or indices; defaults to every estimated
  coefficient. Naming a coefficient that `fixed` held raises
  `latents_held_parameter`, because a held value has no interval.

- level:

  Confidence level strictly between zero and one.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html): a named numeric vector.
[`vcov()`](https://rdrr.io/r/stats/vcov.html): a named symmetric matrix
over the same parameters.
[`confint()`](https://rdrr.io/r/stats/confint.html): a matrix with one
row per parameter and the lower and upper bounds as columns.
[`logLik()`](https://rdrr.io/r/stats/logLik.html): a `"logLik"` object
with `df` and `nobs` attributes.
[`nobs()`](https://rdrr.io/r/stats/nobs.html): the number of groups, the
independent units of the likelihood.

## See also

[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
for the tidy table of estimates, standard errors and intervals.

## Examples

``` r
fit <- multilpa(course_engagement, c("browse", "lectures", "forum_read"),
                "student", n_profiles = 2, n_group_classes = 2,
                n_starts = 2, seed = 1)
coef(fit)
#>           measurement.mean.profile_1.browse 
#>                                   0.5388561 
#>         measurement.mean.profile_1.lectures 
#>                                   0.4194852 
#>       measurement.mean.profile_1.forum_read 
#>                                   0.6130248 
#>           measurement.mean.profile_2.browse 
#>                                  -0.7837817 
#>         measurement.mean.profile_2.lectures 
#>                                  -0.6102031 
#>       measurement.mean.profile_2.forum_read 
#>                                  -0.8914110 
#>       measurement.variance.profile_1.browse 
#>                                   0.5896074 
#>     measurement.variance.profile_1.lectures 
#>                                   0.8397745 
#>   measurement.variance.profile_1.forum_read 
#>                                   0.4870705 
#>       measurement.variance.profile_2.browse 
#>                                   0.5052548 
#>     measurement.variance.profile_2.lectures 
#>                                   0.5489694 
#>   measurement.variance.profile_2.forum_read 
#>                                   0.3495553 
#> profile.probability.profile_1.group_class_1 
#>                                   0.7890308 
#> profile.probability.profile_2.group_class_1 
#>                                   0.2109692 
#> profile.probability.profile_1.group_class_2 
#>                                   0.1640865 
#> profile.probability.profile_2.group_class_2 
#>                                   0.8359135 
#>             group.probability.group_class_1 
#>                                   0.6853020 
#>             group.probability.group_class_2 
#>                                   0.3146980 
logLik(fit)
#> 'log Lik.' -5390.953 (df=15)
nobs(fit)
#> [1] 106
```
