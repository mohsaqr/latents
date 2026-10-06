# Summarise latents results

Summary methods for fitted models, enumerations, bootstrap tests,
covariate and transition models, mixture regressions and pooled
imputations. A summary collects the model's main tables and diagnostics;
printing it shows each table to a set depth.

## Usage

``` r
# S3 method for class 'multilpa_additive'
summary(object, level = 0.95, vcov_type = c("observed", "robust", "opg"), ...)

# S3 method for class 'multilpa_covariates'
summary(object, ...)

# S3 method for class 'multilpa_enumeration'
summary(object, ...)

# S3 method for class 'multilpa_bootstrap_lrt'
summary(object, ...)

# S3 method for class 'latents_growth_mixture'
summary(object, level = 0.95, vcov_type = NULL, ...)

# S3 method for class 'multilpa'
summary(object, ...)

# S3 method for class 'latents_mixture_regression'
summary(object, level = 0.95, vcov_type = NULL, ...)

# S3 method for class 'latents_pooled'
summary(object, ...)

# S3 method for class 'multilpa_transitions'
summary(object, ...)
```

## Arguments

- object:

  An object returned by a latents function.

- level:

  Confidence level for the reported intervals, strictly between zero and
  one.

- vcov_type:

  `"observed"`, `"robust"` or `"opg"`, as for
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md).
  For mixture regressions and growth mixtures, `NULL` reuses the
  inference stored with the fit, or, when there is none, uses robust
  errors for a weighted fit and observed otherwise.

- ...:

  Ignored; present for compatibility with
  [`summary()`](https://rdrr.io/r/base/summary.html).

## Value

An object of class `summary_<class>` holding the summarised tables.
Print it, or read any table with
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
or [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).

## See also

[latents-print](https://pak.dynasite.org/latents/reference/latents-print.md),
[latents-as-data-frame](https://pak.dynasite.org/latents/reference/latents-as-data-frame.md),
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md).

## Examples

``` r
fit <- multilpa(course_engagement, c("browse", "lectures", "forum_read"),
                "student", n_profiles = 2, n_group_classes = 2,
                n_starts = 2, seed = 1)
print(summary(fit), rows = 3)
#> Multilevel LPA: 2 profiles and 2 group classes
#> Individuals: 1422; groups: 106; parameters: 15; converged: TRUE
#> Log likelihood: -5390.952704; AIC: 10811.905
#> BIC (groups): 10851.857; BIC (individuals): 10890.803
#> Best likelihood replicated in 2/2 starts (absolute tolerance 0.00539).
#> 
#> -- profiles --------------------------------------------------------
#>  profile  indicator   mean variance standard_deviation mean_standard_error
#>        1     browse 0.5389   0.5896             0.7679             0.02838
#>        1   lectures 0.4195   0.8398             0.9164             0.03346
#>        1 forum_read 0.6130   0.4871             0.6979             0.02750
#>  variance_standard_error
#>                  0.03075
#>                  0.04258
#>                  0.02753
#>    ... 3 more rows.  get_results(x, what = "profiles")
#> 
#> -- responses -------------------------------------------------------
#>    (no rows)
#> 
#> -- covariances -----------------------------------------------------
#>  profile  indicator indicator_2 covariance
#>        1     browse      browse     0.5896
#>        1   lectures      browse     0.0000
#>        1 forum_read      browse     0.0000
#>    ... 15 more rows.  get_results(x, what = "covariances")
#> 
#> -- profile_probabilities -------------------------------------------
#>  group_class profile probability group_class_probability
#>            1       1      0.7890                  0.6853
#>            1       2      0.2110                  0.6853
#>            2       1      0.1641                  0.3147
#>    ... 1 more rows.  get_results(x, what = "profile_probabilities")
#> 
#> -- counts ----------------------------------------------------------
#>        level class effective_count effective_proportion
#>  individuals     1          842.63               0.5926
#>  individuals     2          579.37               0.4074
#>       groups     1           72.64               0.6853
#>    ... 1 more rows.  get_results(x, what = "counts")
#> 
#> -- posteriors ------------------------------------------------------
#>  row group profile posterior modal
#>    1     1       1 0.8993756  TRUE
#>    2     1       1 0.8666222  TRUE
#>    3     1       1 0.0003331 FALSE
#>    ... 2841 more rows.  get_results(x, what = "posteriors")
#> 
#> -- group_posteriors ------------------------------------------------
#>  group group_size log_likelihood group_class posterior modal
#>      1         13         -41.08           1 4.167e-05 FALSE
#>      2         13         -44.87           1 2.542e-07 FALSE
#>      3         14         -47.99           1 1.000e+00  TRUE
#>    ... 209 more rows.  get_results(x, what = "group_posteriors")
#> 
#> -- assignments -----------------------------------------------------
#>  student browse lectures forum_read profile group_class uncertainty
#>        1   0.73    -0.08       0.30       1           2   0.1006244
#>        1   0.68     0.53      -0.02       1           2   0.1333778
#>        1  -1.51    -0.51      -0.91       2           2   0.0003331
#>  posterior_profile_1 posterior_profile_2
#>            0.8993756              0.1006
#>            0.8666222              0.1334
#>            0.0003331              0.9997
#>    ... 1419 more rows.  get_results(x, what = "assignments")
#> 
#> -- classification --------------------------------------------------
#>        level class n_modal proportion_modal estimated_n estimated_proportion
#>  individuals     1     837           0.5886      842.63               0.5926
#>  individuals     2     585           0.4114      579.37               0.4074
#>       groups     1      73           0.6887       72.64               0.6853
#>  average_posterior odds_correct_classification
#>             0.9681                       20.86
#>             0.9447                       24.85
#>             0.9853                       30.70
#>    ... 1 more rows.  get_results(x, what = "classification")
#> 
#> -- average_posteriors ----------------------------------------------
#>        level assigned_class class n_assigned average_posterior
#>  individuals              1     1        837           0.96809
#>  individuals              1     2        837           0.03191
#>  individuals              2     1        585           0.05528
#>    ... 5 more rows.  get_results(x, what = "average_posteriors")
#> 
#> -- classification_errors -------------------------------------------
#>        level true_class assigned_class probability
#>  individuals          1              1     0.96162
#>  individuals          1              2     0.03838
#>  individuals          2              1     0.04610
#>    ... 5 more rows.  get_results(x, what = "classification_errors")
#> 
#> -- bch_weights -----------------------------------------------------
#>        level unit assigned_class class   weight
#>  individuals    1              1     1  1.04192
#>  individuals    1              1     2 -0.04192
#>  individuals    2              1     1  1.04192
#>    ... 3053 more rows.  get_results(x, what = "bch_weights")
#> 
#> -- entropy ---------------------------------------------------------
#>        level n_classes n_units entropy_sum relative_entropy
#>  individuals         2    1422     152.457           0.8453
#>       groups         2     106       4.959           0.9325
#> 
#> -- residuals -------------------------------------------------------
#>    profile indicator_1 indicator_2     kind observed expected residual
#>  profile_2      browse  forum_read gaussian  0.09022        0  0.09022
#>  profile_2      browse    lectures gaussian -0.05080        0 -0.05080
#>  profile_1      browse  forum_read gaussian  0.03014        0  0.03014
#>  effective_n statistic df p_value p_adjusted
#>        579.4    2.1718 NA 0.02987    0.02987
#>        579.4   -1.2206 NA 0.22224    0.22224
#>        842.6    0.8737 NA 0.38226    0.38226
#>    ... 3 more rows.  get_results(x, what = "residuals")
#> 
#> -- information_criteria --------------------------------------------
#>  log_likelihood n_parameters   aic   kic bic_groups bic_individual sabic_groups
#>           -5391           15 10812 10830      10852          10891        10804
#>  sabic_individual caic_groups caic_individual awe_groups awe_individual
#>             10843       10867           10906      10977          11350
#>  icl_groups icl_individual clc_groups clc_individual
#>       10862          11196      10792          11087
#> 
#> -- model -----------------------------------------------------------
#>  n_observations n_informative n_groups n_profiles n_group_classes centering
#>            1422          1422      106          2               2      none
#>  covariance_structure n_parameters n_parameters_with_measurement log_likelihood
#>                   VVI           15                            15          -5391
#>    aic bic_groups bic_individual converged iterations boundary small_classes
#>  10812      10852          10891      TRUE         12    FALSE         FALSE
#>  best_start n_best_replicated weights
#>           1                 2    <NA>
#> 
#> -- stages ----------------------------------------------------------
#>  stage group_classes fixed log_likelihood parameters
#>  joint             2  <NA>          -5391         15
#>  parameters_with_measurement converged
#>                           15      TRUE
#> 
#> -- starts ----------------------------------------------------------
#>  start log_likelihood converged iterations error boundary
#>      1          -5391      TRUE         12  <NA>    FALSE
#>      2          -5391      TRUE         15  <NA>    FALSE
#> 
#> -- data ------------------------------------------------------------
#>  student browse lectures forum_read
#>        1   0.73    -0.08       0.30
#>        1   0.68     0.53      -0.02
#>        1  -1.51    -0.51      -0.91
#>    ... 1419 more rows.  get_results(x, what = "data")
#> 
#> 19 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
```
