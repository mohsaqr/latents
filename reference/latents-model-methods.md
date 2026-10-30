# Extract model coefficients, covariance matrices and likelihood summaries

These methods extract fitted parameters, measures of parameter
uncertainty, and likelihood information for reporting and further
statistical analysis.

[`coef()`](https://rdrr.io/r/stats/coef.html) returns the model's
parameter values. [`vcov()`](https://rdrr.io/r/stats/vcov.html) returns
their estimated sampling variances and covariances, describing the
uncertainty in individual estimates and their joint variation.
[`confint()`](https://rdrr.io/r/stats/confint.html) provides confidence
intervals for selected parameters.

[`logLik()`](https://rdrr.io/r/stats/logLik.html) returns the fitted log
likelihood together with the number of estimated parameters and the
associated sample size. [`nobs()`](https://rdrr.io/r/stats/nobs.html)
reports the number of independent units used in the model's group-based
information criterion: individual observations for a single-level model
and observed groups for a grouped model. The treatment of entirely
missing units follows the fitted model's `drop_empty` setting.

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

# S3 method for class 'multilpa'
coef(object, scale = c("natural", "unconstrained"), ...)

# S3 method for class 'multilpa'
vcov(object, data = NULL, scale = c("natural", "unconstrained"), ...)

# S3 method for class 'multilpa'
confint(object, parm, level = 0.95, data = NULL, ...)
```

## Arguments

- object:

  A fitted model: from
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  (class `multilpa`), a covariate model (`multilpa_covariates`), a
  latent transition model (`multilpa_transitions`), or pooled
  imputations (`latents_pooled`:
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) and
  [`confint()`](https://rdrr.io/r/stats/confint.html)).

- data:

  Optional. The data frame the model was fitted to. When omitted it is
  rebuilt from what the fit stores (indicators, identifiers, occasions
  and designs), to floating-point precision.

- step:

  Finite-difference step for the observed information (covariate
  models).

- vcov_type:

  `"observed"`, `"robust"` or `"opg"`, as for
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  (covariate models).

- scale:

  `"natural"`, the default, or `"unconstrained"`. Natural coefficients
  are expressed on the reported scale: means in indicator units,
  variances in squared indicator units, and probabilities on `[0, 1]`.
  Their covariance is obtained by the delta method. The unconstrained
  scale uses the coordinates in which the model is estimated: log
  variances (diagonal), log-Cholesky coordinates (full covariance), log
  volumes, shapes and orientations (a structure that constrains them
  across profiles), and baseline-category logits for mixing and response
  probabilities. In a covariate model, means and membership coefficients
  are the same on both scales. With `method = "bootstrap"` only the
  natural scale exists, and `scale = "unconstrained"` raises
  `latents_unsupported_inference`.

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
  `step`); an argument
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  does not take raises `latents_bad_argument`.
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) of a covariate model
  warns `latents_unused_argument` for an argument it does not take, and
  refuses `method = "bootstrap"` with `latents_unsupported_inference`;
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) of a transition model
  warns `latents_unused_argument` for any.
  [`confint()`](https://rdrr.io/r/stats/confint.html) of a transition
  model passes them to
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md).
  Ignored otherwise.

- parm:

  Optional coefficient names or indices; defaults to every estimated
  coefficient. A coefficient held by `fixed` has no interval and raises
  `latents_held_parameter`; an unknown name or an out-of-range index
  raises `latents_bad_argument`.

- level:

  Confidence level, strictly between 0 and 1 (default 0.95).

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html): a named numeric vector.
[`vcov()`](https://rdrr.io/r/stats/vcov.html): a named symmetric matrix
over the same parameters.
[`confint()`](https://rdrr.io/r/stats/confint.html): a matrix with one
row per parameter and the lower and upper bounds as columns.
[`logLik()`](https://rdrr.io/r/stats/logLik.html): a `"logLik"` object
with `df` and `nobs` attributes.
[`nobs()`](https://rdrr.io/r/stats/nobs.html): the number of independent
groups counted by the fitted BIC. With `drop_empty = TRUE`, groups
without any observed indicator are excluded; with `drop_empty = FALSE`,
they are counted. Groups with zero sampling weight are excluded under
either setting. For a single-level fit these groups are its individual
observations.

## Details

These methods support profile and class models, models with membership
covariates and homogeneous latent transition models. Pooled
multiple-imputation results support
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html) and
[`confint()`](https://rdrr.io/r/stats/confint.html). Methods for mixture
regressions, growth mixtures, additive models, cross-level models and
general transition models are described in
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md).

Categorical parameter names escape reserved characters within indicator
and category components (`%`, `:`, `,`, `[` and `]`) using percent
codes. This keeps distinct parameters distinct when their original
labels contain delimiters. Ordinary names are unchanged; use
`names(coef(object))` or numeric positions when selecting parameters for
[`confint()`](https://rdrr.io/r/stats/confint.html).

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
#>                                   0.5388558 
#>         measurement.mean.profile_1.lectures 
#>                                   0.4194848 
#>       measurement.mean.profile_1.forum_read 
#>                                   0.6130241 
#>           measurement.mean.profile_2.browse 
#>                                  -0.7837825 
#>         measurement.mean.profile_2.lectures 
#>                                  -0.6102034 
#>       measurement.mean.profile_2.forum_read 
#>                                  -0.8914115 
#>       measurement.variance.profile_1.browse 
#>                                   0.5896076 
#>     measurement.variance.profile_1.lectures 
#>                                   0.8397747 
#>   measurement.variance.profile_1.forum_read 
#>                                   0.4870710 
#>       measurement.variance.profile_2.browse 
#>                                   0.5052541 
#>     measurement.variance.profile_2.lectures 
#>                                   0.5489692 
#>   measurement.variance.profile_2.forum_read 
#>                                   0.3495550 
#> profile.probability.profile_1.group_class_1 
#>                                   0.1640873 
#> profile.probability.profile_2.group_class_1 
#>                                   0.8359127 
#> profile.probability.profile_1.group_class_2 
#>                                   0.7890314 
#> profile.probability.profile_2.group_class_2 
#>                                   0.2109686 
#>             group.probability.group_class_1 
#>                                   0.3146985 
#>             group.probability.group_class_2 
#>                                   0.6853015 
logLik(fit)
#> 'log Lik.' -5390.953 (df=15)
nobs(fit)
#> [1] 106
```
