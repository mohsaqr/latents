# Covariates, outcomes and staged estimation

Variables outside the measurement model can be related to the latent
classes in two ways. The three-step approach fits the profiles first and
then relates them to an external variable, correcting for classification
error, so the profiles stay fixed. The one-step approach adds the
variable to the model, so it can also change the profiles. Staged
estimation applies the same idea of fixed profiles to the group classes.

## Data

`course_engagement` has one row per course enrolment: 1,422 enrolments
from 106 students. `student` identifies the group, and five activity
measures (`browse`, `lectures`, `forum_read`, `forum_post`,
`attendance`) are the indicators. `previous_grade` is the standardized
grade from the course before each enrolment; it is not an indicator, so
it can serve as an external variable.

``` r

library(latents)
```

## The measurement model

To fit the profiles without the external variable, we call
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
with `vars` for the indicators, `id` for the group, `n_profiles`,
`n_group_classes`, `n_starts` and `seed`.

``` r

fit <- multilpa(course_engagement,
                vars = c("browse", "lectures", "forum_read", "forum_post",
                         "attendance"),
                id = "student", n_profiles = 2, n_group_classes = 2,
                n_starts = 10, seed = 1)
fit
#> Two-level latent profile analysis: 2 profiles, 2 group classes
#> 1422 individuals in 106 groups; varying diagonal residual covariance (VVI)
#> Log likelihood: -8439.816434 | AIC: 16925.633 | BIC (groups): 16986.892
#> Converged: TRUE | iterations: 8 | best start: 8/10
#> 
#>  profile browse lectures forum_read forum_post attendance count proportion
#>        1 -0.766   -0.606     -0.879     -0.753     -0.921   588      0.413
#>        2  0.540    0.427      0.620      0.531      0.649   834      0.587
#> 
#> Variances and standard errors: get_results(x, "profiles"). 
#> Every other table: get_results(x, what = ), or get_results(x, "all").
```

## Classification error

To obtain the classification-error matrix, we call
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
with `"classification_errors"` and `level` for the observations or the
groups. Each row gives the probability of each assigned class for cases
in a true class; the three-step methods use it to correct for
misclassification.

``` r

get_results(fit, "classification_errors", level = "individuals")
#>         level true_class assigned_class probability
#> 1 individuals          1              1      0.9786
#> 2 individuals          1              2      0.0214
#> 3 individuals          2              1      0.0128
#> 4 individuals          2              2      0.9872
```

## Class means of an external variable

To estimate the mean of an external variable in each profile, we call
[`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
with `data` for the data frame and `outcome` for the variable. The
default method, BCH, corrects for classification error; standard errors
are clustered on groups.

``` r

three_step(fit, data = course_engagement, outcome = "previous_grade")
#>         level method class estimate standard_error conf_low conf_high effective_n
#> 1 individuals    bch     1   -0.315         0.0359   -0.385    -0.244         564
#> 2 individuals    bch     2    0.222         0.0347    0.154     0.290         810
```

To compare the profiles, we call
[`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
with `contrast = "pairs"`, which returns each difference between profile
means with a Benjamini-Hochberg adjusted p-value.

``` r

three_step(fit, data = course_engagement, outcome = "previous_grade",
           contrast = "pairs")
#>         level method class reference_class estimate standard_error statistic  p_value
#> 1 individuals    bch     2               1    0.536          0.052      10.3 5.72e-25
#>   p_value_adjusted conf_low conf_high
#> 1         5.72e-25    0.435     0.638
```

To see the effect of ignoring classification error, we call
[`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
with `method = "modal"`, which treats each case’s most probable profile
as known.

``` r

three_step(fit, data = course_engagement, outcome = "previous_grade",
           method = "modal")
#>         level method class estimate standard_error conf_low conf_high effective_n
#> 1 individuals  modal     1   -0.305         0.0357   -0.375    -0.235         586
#> 2 individuals  modal     2    0.214         0.0343    0.146     0.281         836
```

## Predicting class membership

To regress profile membership on a covariate, we call
[`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md) with
`data` for the data frame and `covariates` for the predictors. It
returns multinomial logit coefficients, corrected for classification
error. `vcov_type = "robust"` clusters the standard errors on groups.

``` r

r3step(fit, data = course_engagement, covariates = "previous_grade",
       vcov_type = "robust")
#>         level outcome           term estimate standard_error statistic  p_value p_value_adjusted
#> 1 individuals class_1    (Intercept)   -0.376         0.1303     -2.89 3.91e-03               NA
#> 2 individuals class_1 previous_grade   -0.575         0.0593     -9.69 3.35e-22         3.35e-22
#>   conf_low conf_high
#> 1   -0.631    -0.121
#> 2   -0.691    -0.459
```

## One-step estimation

To estimate the covariate together with the profiles, we call
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
with `profile_covariates` for predictors of profile membership.
`group_covariates` does the same for group-class membership and needs a
variable that is constant within groups. Covariate models require
complete data.

``` r

one_step <- multilpa(course_engagement,
                     vars = c("browse", "lectures", "forum_read",
                              "forum_post", "attendance"),
                     id = "student", n_profiles = 2, n_group_classes = 2,
                     profile_covariates = "previous_grade", n_starts = 10,
                     seed = 1)
one_step
#> Multilevel LPA with covariates: 2 profiles, 2 group classes
#> Log likelihood -8422.837634; AIC 16893.675; BIC (groups) 16957.598; converged TRUE
#> 
#>  profile browse lectures forum_read forum_post attendance count proportion
#>        1 -0.766   -0.609     -0.879     -0.753     -0.920   588      0.413
#>        2  0.540    0.429      0.619      0.530      0.648   834      0.587
#> 
#> Variances and standard errors: get_results(x, "profiles"). 
#> Every other table: get_results(x, what = ), or get_results(x, "all").
```

To obtain the membership coefficients, we call
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
with `"coefficients"`, and for their standard errors and intervals,
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md).

``` r

get_results(one_step, "coefficients")
#>     level       outcome           term parameter estimate
#> 1 profile     profile_1  group_class_1     logit    1.538
#> 2 profile     profile_1  group_class_2     logit   -1.268
#> 3 profile     profile_1 previous_grade     logit   -0.422
#> 4   group group_class_1    (Intercept)     logit   -0.762
parameter_inference(one_step)
#>          level       outcome           term   parameter estimate standard_error statistic   p_value
#> 1  measurement     profile_1         browse        mean   -0.766         0.0311    -24.65 3.78e-134
#> 2  measurement     profile_1       lectures        mean   -0.609         0.0309    -19.71  1.67e-86
#> 3  measurement     profile_1     forum_read        mean   -0.879         0.0261    -33.65 3.00e-248
#> 4  measurement     profile_1     forum_post        mean   -0.753         0.0277    -27.14 3.54e-162
#> 5  measurement     profile_1     attendance        mean   -0.920         0.0264    -34.87 2.03e-266
#> 6  measurement     profile_2         browse        mean    0.540         0.0271     19.88  6.40e-88
#> 7  measurement     profile_2       lectures        mean    0.429         0.0324     13.24  4.85e-40
#> 8  measurement     profile_2     forum_read        mean    0.619         0.0248     25.02 3.95e-138
#> 9  measurement     profile_2     forum_post        mean    0.530         0.0294     18.02  1.42e-72
#> 10 measurement     profile_2     attendance        mean    0.648         0.0223     29.12 2.05e-186
#> 11 measurement     profile_1         browse    variance    0.526         0.0321        NA        NA
#> 12 measurement     profile_1       lectures    variance    0.538         0.0320        NA        NA
#> 13 measurement     profile_1     forum_read    variance    0.363         0.0227        NA        NA
#> 14 measurement     profile_1     forum_post    variance    0.422         0.0262        NA        NA
#> 15 measurement     profile_1     attendance    variance    0.375         0.0228        NA        NA
#> 16 measurement     profile_2         browse    variance    0.591         0.0296        NA        NA
#> 17 measurement     profile_2       lectures    variance    0.842         0.0419        NA        NA
#> 18 measurement     profile_2     forum_read    variance    0.483         0.0245        NA        NA
#> 19 measurement     profile_2     forum_post    variance    0.688         0.0347        NA        NA
#> 20 measurement     profile_2     attendance    variance    0.385         0.0196        NA        NA
#> 21     profile     profile_1  group_class_1 coefficient    1.538         0.1676      9.17  4.53e-20
#> 22     profile     profile_1  group_class_2 coefficient   -1.268         0.0943    -13.44  3.44e-41
#> 23     profile     profile_1 previous_grade coefficient   -0.422         0.0736     -5.73  9.79e-09
#> 24       group group_class_1    (Intercept) coefficient   -0.762         0.2212     -3.44  5.73e-04
#>    p_adjusted conf_low conf_high
#> 1   3.78e-134   -0.827    -0.705
#> 2    1.67e-86   -0.669    -0.548
#> 3   3.00e-248   -0.930    -0.828
#> 4   3.54e-162   -0.807    -0.699
#> 5   2.03e-266   -0.972    -0.868
#> 6    6.40e-88    0.486     0.593
#> 7    4.85e-40    0.365     0.492
#> 8   3.95e-138    0.571     0.668
#> 9    1.42e-72    0.473     0.588
#> 10  2.05e-186    0.605     0.692
#> 11         NA    0.463     0.589
#> 12         NA    0.475     0.600
#> 13         NA    0.318     0.407
#> 14         NA    0.371     0.474
#> 15         NA    0.330     0.420
#> 16         NA    0.533     0.649
#> 17         NA    0.760     0.924
#> 18         NA    0.435     0.531
#> 19         NA    0.620     0.756
#> 20         NA    0.347     0.424
#> 21   4.53e-20    1.209     1.866
#> 22   3.44e-41   -1.453    -1.083
#> 23   9.79e-09   -0.567    -0.278
#> 24   5.73e-04   -1.195    -0.328
```

## Staged estimation

To hold the profiles fixed while the group classes are estimated, we
call
[`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
with the same arguments as
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md).
It first fits the profiles with one group class, then estimates the
group classes with the profile means and variances held fixed.

``` r

staged <- fit_staged(course_engagement,
                     vars = c("browse", "lectures", "forum_read",
                              "forum_post", "attendance"),
                     id = "student", n_profiles = 2, n_group_classes = 2,
                     n_starts = 10, seed = 1)
staged
#> Two-level latent profile analysis: 2 profiles, 2 group classes
#> 1422 individuals in 106 groups; varying diagonal residual covariance (VVI)
#> Log likelihood: -8440.081199 | AIC: 16886.162 | BIC (groups): 16894.153
#> Converged: TRUE | iterations: 10 | best start: 3/10
#> Held fixed, not estimated here: means, variances (first stage)
#> Parameters estimated here: 3; with the held measurement: 23
#> 
#>  profile browse lectures forum_read forum_post attendance count proportion
#>        1  0.537    0.430      0.618      0.533      0.653   833      0.586
#>        2 -0.759   -0.608     -0.873     -0.753     -0.923   589      0.414
#> 
#> Variances and standard errors: get_results(x, "profiles"). 
#> Every other table: get_results(x, what = ), or get_results(x, "all").
```

To see what each stage estimated, we call
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
with `"stages"`. `parameters` counts the parameters estimated at a
stage, and `parameters_with_measurement` adds those held fixed.

``` r

get_results(staged, "stages")
#>         stage group_classes            fixed log_likelihood parameters parameters_with_measurement
#> 1 measurement             1             <NA>          -8615         21                          21
#> 2  membership             2 means, variances          -8440          3                          23
#>   converged
#> 1      TRUE
#> 2      TRUE
```

## Reference

The calls below list other functions for external variables and fixed
parameters. Each comment states what the call returns.

``` r

get_results(fit, "classification_errors", level = "groups")  # group-level errors
get_results(fit, "bch_weights", level = "individuals")       # BCH weights per case
three_step(fit, data = course_engagement, outcome = "previous_grade",
           method = "proportional")        # weighting by posterior probabilities
three_step(fit, data = course_engagement, outcome = "previous_grade",
           vcov_type = "independent")      # standard errors without clustering
r3step(fit, data = course_engagement,
       covariates = "previous_grade")      # observed-information standard errors
measurement <- multilpa(course_engagement,
                        vars = c("browse", "lectures", "forum_read",
                                 "forum_post", "attendance"),
                        id = "student", n_profiles = 2,
                        n_group_classes = 1, seed = 1)  # a one-class measurement fit
fit_staged(course_engagement,
           vars = c("browse", "lectures", "forum_read", "forum_post",
                    "attendance"),
           id = "student", n_profiles = 2, n_group_classes = 2,
           measurement = measurement)      # stage two from that fit
multilpa(course_engagement,
         vars = c("browse", "lectures", "forum_read", "forum_post",
                  "attendance"),
         id = "student", n_profiles = 2, n_group_classes = 2,
         start = starting_values(fit, what = "measurement"),
         fixed = "variances")              # hold chosen parameters fixed
parameter_inference(staged)                # standard errors of the estimated stage
```

## Assumptions

- Three-step methods assume classification errors are independent of the
  external variable given the true class, and treat the error matrix as
  known.
- One-step coefficients are conditional on the group class and can
  change the profiles, so they need not equal the three-step
  coefficients.
- Staged estimation does not carry first-stage uncertainty into the
  second stage, and its likelihood is not the joint maximum.
- All associations are observational.
