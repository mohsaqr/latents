# Fit a model to multiply imputed data and pool it with Rubin's rules

Fits
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
to every completed data set, aligns each fit's profile and group-class
labels to the first one, and pools the estimates with Rubin's rules.
This is the route for missing *covariates*, which `missing = "fiml"`
cannot integrate out: the model conditions on covariates, so it has no
distribution for them to integrate over. Impute them first (for example
with the mice package, including the indicators and the group identifier
in the imputation model), then pass the completed data here.

## Usage

``` r
pool_imputations(
  imputed,
  vars,
  id,
  n_profiles,
  ...,
  level = 0.95,
  vcov_type = c("observed", "robust", "opg"),
  adjust = .multilpa_p_adjust_methods
)
```

## Arguments

- imputed:

  The completed data sets: a list of at least two data frames with the
  same columns and rows in the same order, or a `mids` object from the
  mice package (which is then completed here).

- vars, id, n_profiles:

  As in
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).

- ...:

  Further arguments for
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  such as `n_group_classes`, `profile_covariates`, `group_covariates`,
  `categorical`, `n_starts` and `seed`. The same arguments fit every
  imputation.

- level:

  Confidence level of the pooled intervals.

- vcov_type, adjust:

  As in
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md),
  applied to every imputation's fit.

## Value

An object of class `latents_pooled`.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) and
`get_results(x, "estimates")` give one row per parameter, with the
columns of
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
(`level`, `outcome`, `term`, `parameter`, `estimate`, `standard_error`,
`statistic`, `p_value`, `p_adjusted`, `conf_low`, `conf_high`) plus `df`
(Rubin's degrees of freedom), `within` and `between` (the two variance
components), `riv` (the relative increase in variance due to
missingness) and `fmi` (the fraction of missing information).
`statistic` and `p_value` are `NA` where
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
reports no test (variances and probabilities).
`get_results(x, "imputations")` gives each imputation's aligned
estimates, one row per imputation and parameter, and
`get_results(x, "fits")` one row per imputation with its log likelihood,
convergence and the label permutation that aligned it.

## Details

Pooling is on the reported (natural) scale, the scale
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
tests and bounds on: for each parameter the pooled estimate is the mean
over imputations, and its variance is the mean within-imputation
variance plus `(1 + 1/m)` times the between-imputation variance. Degrees
of freedom follow Rubin (1987) with a large complete-data sample, which
is what maximum-likelihood inference assumes.

## Choosing the imputation model

Rubin's rules are only as good as the imputations. The covariate and the
indicators are related through the latent classes, so an imputation
model that is linear in the indicators (mice's defaults) misses that
structure: in `validation/pool-imputations-coverage.R` (a membership
slope with 30% of its covariate missing at random given the indicators,
200 replications, ten imputations) 95% intervals covered 0.925 under
predictive mean matching and 0.935 under Bayesian linear regression,
against 0.995 when the covariate was drawn from its exact conditional
distribution under the generating model. The last shows the pooling
itself is sound; the shortfall is the imputation model's (Meng, 1994).
Give the imputation model what the latent structure implies: every
indicator, and interactions or nonlinear terms among them, or impute
within classes assigned by a covariate-free fit. When missingness
depends only on the indicators, complete-case analysis of the membership
slope is not biased (covered 0.955 in the same study), because selecting
on the outcome side of a logit leaves its slope intact; imputation then
buys precision, not validity.

## Conditions

`latents_bad_data` when `imputed` is not a list of at least two data
frames with identical columns and row counts; `latents_pooling_failed`
when an imputation's fit or inference fails (its message names the
imputation and carries the original reason), since Rubin's rules need
every imputation; `latents_missing_package` for a `mids` object without
mice installed. Warnings from the individual fits are passed on as they
are raised.

## References

Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys*.
Wiley.

Meng, X.-L. (1994). Multiple-imputation inferences with uncongenial
sources of input. *Statistical Science*, 9(4), 538–558.

van Buuren, S. (2018). *Flexible Imputation of Missing Data* (2nd ed.),
section 5.2. Chapman & Hall/CRC.

## Examples

``` r
set.seed(1)
# Two stand-in completed data sets. With real missing covariates, pass the
# `mids` object mice::mice() returns instead.
completed <- lapply(1:2, function(i) {
  within(subset(course_engagement, student <= 40),
         previous_grade <- previous_grade + rnorm(length(previous_grade), sd = 0.1))
})
pooled <- pool_imputations(completed, c("browse", "lectures", "forum_read"),
                           "student", n_profiles = 2, n_group_classes = 1,
                           profile_covariates = "previous_grade",
                           n_starts = 2, seed = 1)
pooled
#> Pooled over 2 imputations (Rubin's rules): multilpa_covariates, 2 profiles, 1 group classes
#> Largest fraction of missing information: 0.002
#>        level   outcome           term   parameter estimate standard_error
#>  measurement profile_1         browse        mean   0.5673        0.04202
#>  measurement profile_1       lectures        mean   0.4578        0.05095
#>  measurement profile_1     forum_read        mean   0.6314        0.04056
#>  measurement profile_2         browse        mean  -0.8397        0.05746
#>  measurement profile_2       lectures        mean  -0.6244        0.05838
#>  measurement profile_2     forum_read        mean  -0.9156        0.04881
#>  measurement profile_1         browse    variance   0.5414        0.04379
#>  measurement profile_1       lectures    variance   0.8138        0.06410
#>  measurement profile_1     forum_read    variance   0.4739        0.03935
#>  measurement profile_2         browse    variance   0.4890        0.05734
#>  measurement profile_2       lectures    variance   0.5656        0.06255
#>  measurement profile_2     forum_read    variance   0.3392        0.04177
#>      profile profile_1    (Intercept) coefficient   0.6142        0.10890
#>      profile profile_1 previous_grade coefficient   0.6516        0.10996
#>  conf_low conf_high       fmi
#>    0.4850    0.6497 8.932e-04
#>    0.3580    0.5577 3.942e-04
#>    0.5519    0.7109 1.837e-03
#>   -0.9523   -0.7271 1.458e-03
#>   -0.7388   -0.5100 7.534e-04
#>   -1.0113   -0.8200 1.307e-03
#>    0.4556    0.6272 5.457e-05
#>    0.6881    0.9394 8.143e-05
#>    0.3968    0.5511 8.404e-04
#>    0.3766    0.6014 1.212e-04
#>    0.4430    0.6882 6.422e-04
#>    0.2573    0.4211 6.997e-05
#>    0.4008    0.8276 3.305e-04
#>    0.4361    0.8671 1.322e-03
get_results(pooled, "fits")
#>   imputation log_likelihood n_parameters converged profile_order group_order
#> 1          1      -2099.880           14      TRUE           1,2           1
#> 2          2      -2099.696           14      TRUE           1,2           1
```
