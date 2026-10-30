# Fit and pool mixture models across multiple imputations

`pool_imputations()` estimates a mixture model from multiply imputed
data and combines the results to account for uncertainty associated with
missing values. It fits the same model specification to each completed
data set, allowing parameter estimates to reflect the alternative
plausible values supplied by the imputation procedure.

Class labels can differ across fitted mixture models even when their
profiles are similar. Before pooling, the function matches individual
profiles and group classes to those in the first imputation's fit. It
then combines the aligned parameter estimates and covariance matrices
using Rubin's rules. The resulting uncertainty incorporates both the
estimated sampling variance within each analysis and variation in
estimates between imputations.

The function returns pooled estimates, standard errors, confidence
intervals, and tests for eligible parameters, together with summaries of
the individual analyses. Completed data sets must be supplied as a list
or a `mids` object from mice; the function does not construct the
imputation model.

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
  mice package, which is completed here.

- vars, id, n_profiles:

  As in
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).

- ...:

  Further arguments for
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  such as `n_group_classes`, `profile_covariates`, `group_covariates`,
  `categorical`, `n_starts` and `seed`, used for every imputation.

- level:

  Confidence level of the pooled intervals (default 0.95).

- vcov_type:

  `"observed"` (default), `"robust"` or `"opg"`, as in
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md),
  applied to every imputation's fit. For fits with sampling weights the
  default is `"robust"`, and the other types are refused with
  `latents_unsupported_weights`; for fits with person centring or
  scaling the default is `"robust"` too, and the other types raise
  `latents_unsupported_inference`.

- adjust:

  Multiplicity correction applied to the pooled p-values to give
  `p_adjusted`, as in
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md):
  `"none"` (default) or a method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). The
  per-imputation tables are computed without correction.

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

Pooling is on the reported (natural) scale: the pooled estimate is the
mean over imputations, and its variance is the mean within-imputation
variance plus `(1 + 1/m)` times the between-imputation variance. Degrees
of freedom follow Rubin (1987) with a large complete-data sample. Pooled
intervals are symmetric t intervals, rather than the transformed Wald
intervals used for bounded parameters by
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md).
A structural constant, such as the probability of one for a single
class, has zero within- and between-imputation variance and zero missing
information. When the within-imputation variance is zero and the
between-imputation variance is positive, `riv` is infinite and `fmi` is
one.

## Interpretation

Pooling requires an imputation model appropriate for the analysis and
missingness assumptions, including the grouping structure where
relevant. The number of classes and other model specifications remain
fixed across imputations. Label matching resolves differences in
numbering but cannot establish that substantially different class
solutions represent the same underlying profiles.

Estimates are pooled on their reported scale. Consequently, symmetric
confidence intervals for probabilities or variances can extend beyond
their permissible ranges.

## Choosing the imputation model

Valid pooling requires an imputation model compatible with the analysis
model and the assumed missingness mechanism. Associations between
covariates and indicators may be nonlinear because they depend on latent
class membership. The imputation model should represent these
associations, account for clustering, include relevant indicators and
predictors of missingness, and propagate uncertainty in its estimated
parameters. Imputing separately within modally assigned classes omits
classification uncertainty. Complete-case analysis can also change the
within-class indicator distributions when missingness depends on the
indicators.

## Conditions

`latents_bad_data` when `imputed` is not a list of at least two data
frames with identical columns and row counts; `latents_pooling_failed`
when any imputation's fit or inference fails (the message names the
imputation and gives the original reason); `latents_missing_package` for
a `mids` object without mice installed. Warnings from the individual
fits are passed on as they are raised.

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
# Illustrative completed data sets for demonstrating the pooling interface.
# For an analysis, supply data generated by an appropriate imputation model.
completed <- lapply(seq_len(2), function(i) {
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
#>  measurement profile_1     forum_read        mean   0.6314        0.04057
#>  measurement profile_2         browse        mean  -0.8397        0.05746
#>  measurement profile_2       lectures        mean  -0.6244        0.05838
#>  measurement profile_2     forum_read        mean  -0.9157        0.04881
#>  measurement profile_1         browse    variance   0.5414        0.04379
#>  measurement profile_1       lectures    variance   0.8138        0.06410
#>  measurement profile_1     forum_read    variance   0.4740        0.03935
#>  measurement profile_2         browse    variance   0.4890        0.05733
#>  measurement profile_2       lectures    variance   0.5656        0.06255
#>  measurement profile_2     forum_read    variance   0.3392        0.04177
#>      profile profile_1    (Intercept) coefficient   0.6143        0.10891
#>      profile profile_1 previous_grade coefficient   0.6516        0.10996
#>  conf_low conf_high       fmi
#>    0.4850    0.6497 9.204e-04
#>    0.3579    0.5577 4.156e-04
#>    0.5519    0.7109 1.910e-03
#>   -0.9523   -0.7271 1.533e-03
#>   -0.7389   -0.5100 7.724e-04
#>   -1.0113   -0.8200 1.362e-03
#>    0.4556    0.6272 5.685e-05
#>    0.6881    0.9394 8.618e-05
#>    0.3968    0.5511 8.793e-04
#>    0.3766    0.6013 1.392e-04
#>    0.4430    0.6882 6.543e-04
#>    0.2573    0.4210 7.721e-05
#>    0.4008    0.8277 3.680e-04
#>    0.4361    0.8672 1.331e-03
get_results(pooled, "fits")
#>   imputation log_likelihood n_parameters converged profile_order group_order
#> 1          1      -2099.880           14      TRUE           1,2           1
#> 2          2      -2099.696           14      TRUE           1,2           1
```
