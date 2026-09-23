# Extract multilevel LPA covariance estimates

Extract multilevel LPA covariance estimates

## Usage

``` r
# S3 method for class 'multilpa'
vcov(object, data = NULL, scale = c("natural", "unconstrained"), ...)
```

## Arguments

- object:

  A fitted `multilpa` model.

- data:

  Optional. The data frame the model was fitted to; when omitted it is
  rebuilt from the indicators, identifiers and occasions the fit stores,
  which round-trip exactly. Supplying it is the stronger check that the
  caller still holds that frame.

- scale:

  Which parameter scale the covariance is on. `"natural"`, the default,
  is the covariance of the estimates
  [`coef()`](https://rdrr.io/r/stats/coef.html) reports – variances in
  their own units and probabilities as probabilities – carried from the
  estimation scale by the delta method with
  `.multilpa_inference_jacobian()`. `"unconstrained"` is the covariance
  on the scale the model is actually estimated on: log variances for a
  diagonal fit, log-Cholesky coordinates for a full-covariance fit, and
  baseline-category logits for the mixing and response probabilities.

- ...:

  Additional arguments passed to
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md),
  including `method = "bootstrap"`, which is how a covariance structure
  the Wald path cannot chart reports one here.

## Value

A square numeric matrix with one row and column per *estimated*
coefficient, named as [`coef()`](https://rdrr.io/r/stats/coef.html)
names them, on the scale `scale` asks for. A fit made with `fixed` held
part of its measurement model at supplied values; those coefficients
were not estimated here, so they carry no row or column. With
`scale = "unconstrained"` the matrix is `n_parameters` square, for any
fit; on the natural scale it is larger by one row and column for each
set of probabilities whose reference category the estimation scale
drops, and singular by construction, because each set of probabilities
sums to one. [`confint()`](https://rdrr.io/r/stats/confint.html) and the
`conf_low`/`conf_high` columns of
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
are built from the natural-scale matrix.

## Examples

``` r
set.seed(3)
example_data <- data.frame(
  school = rep(seq_len(10), each = 8),
  score_a = stats::rnorm(80), score_b = stats::rnorm(80)
)
fit <- multilpa(example_data, c("score_a", "score_b"), "school",
                n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
vcov(fit, data = example_data, scale = "unconstrained")
#>                                            measurement.mean.profile_1.score_a
#> measurement.mean.profile_1.score_a                               0.0153915493
#> measurement.mean.profile_1.score_b                              -0.0010571917
#> measurement.mean.profile_2.score_a                               0.0005337047
#> measurement.mean.profile_2.score_b                              -0.0079173732
#> measurement.log_variance.profile_1.score_a                       0.0083322596
#> measurement.log_variance.profile_1.score_b                      -0.0015439904
#> measurement.log_variance.profile_2.score_a                      -0.0145516133
#> measurement.log_variance.profile_2.score_b                      -0.0048561762
#> profile.logit.profile_1.group_class_1                            0.0240110484
#>                                            measurement.mean.profile_1.score_b
#> measurement.mean.profile_1.score_a                              -0.0010571917
#> measurement.mean.profile_1.score_b                               0.0176709589
#> measurement.mean.profile_2.score_a                              -0.0007020338
#> measurement.mean.profile_2.score_b                              -0.0034551439
#> measurement.log_variance.profile_1.score_a                      -0.0010481218
#> measurement.log_variance.profile_1.score_b                      -0.0005796778
#> measurement.log_variance.profile_2.score_a                       0.0033335425
#> measurement.log_variance.profile_2.score_b                      -0.0030901567
#> profile.logit.profile_1.group_class_1                           -0.0048192257
#>                                            measurement.mean.profile_2.score_a
#> measurement.mean.profile_1.score_a                               5.337047e-04
#> measurement.mean.profile_1.score_b                              -7.020338e-04
#> measurement.mean.profile_2.score_a                               4.714972e-03
#> measurement.mean.profile_2.score_b                               7.550271e-04
#> measurement.log_variance.profile_1.score_a                      -6.091536e-04
#> measurement.log_variance.profile_1.score_b                       6.942619e-05
#> measurement.log_variance.profile_2.score_a                      -5.199690e-03
#> measurement.log_variance.profile_2.score_b                       1.112444e-03
#> profile.logit.profile_1.group_class_1                            4.075934e-03
#>                                            measurement.mean.profile_2.score_b
#> measurement.mean.profile_1.score_a                              -0.0079173732
#> measurement.mean.profile_1.score_b                              -0.0034551439
#> measurement.mean.profile_2.score_a                               0.0007550271
#> measurement.mean.profile_2.score_b                               0.0704434449
#> measurement.log_variance.profile_1.score_a                      -0.0120041870
#> measurement.log_variance.profile_1.score_b                       0.0053541918
#> measurement.log_variance.profile_2.score_a                       0.0180114550
#> measurement.log_variance.profile_2.score_b                       0.0214968486
#> profile.logit.profile_1.group_class_1                           -0.0302364685
#>                                            measurement.log_variance.profile_1.score_a
#> measurement.mean.profile_1.score_a                                       0.0083322596
#> measurement.mean.profile_1.score_b                                      -0.0010481218
#> measurement.mean.profile_2.score_a                                      -0.0006091536
#> measurement.mean.profile_2.score_b                                      -0.0120041870
#> measurement.log_variance.profile_1.score_a                               0.0437274794
#> measurement.log_variance.profile_1.score_b                              -0.0023086443
#> measurement.log_variance.profile_2.score_a                              -0.0190529569
#> measurement.log_variance.profile_2.score_b                              -0.0076756862
#> profile.logit.profile_1.group_class_1                                    0.0320064227
#>                                            measurement.log_variance.profile_1.score_b
#> measurement.mean.profile_1.score_a                                      -1.543990e-03
#> measurement.mean.profile_1.score_b                                      -5.796778e-04
#> measurement.mean.profile_2.score_a                                       6.942619e-05
#> measurement.mean.profile_2.score_b                                       5.354192e-03
#> measurement.log_variance.profile_1.score_a                              -2.308644e-03
#> measurement.log_variance.profile_1.score_b                               3.311631e-02
#> measurement.log_variance.profile_2.score_a                               4.816835e-03
#> measurement.log_variance.profile_2.score_b                              -1.833484e-03
#> profile.logit.profile_1.group_class_1                                   -5.972274e-03
#>                                            measurement.log_variance.profile_2.score_a
#> measurement.mean.profile_1.score_a                                       -0.014551613
#> measurement.mean.profile_1.score_b                                        0.003333543
#> measurement.mean.profile_2.score_a                                       -0.005199690
#> measurement.mean.profile_2.score_b                                        0.018011455
#> measurement.log_variance.profile_1.score_a                               -0.019052957
#> measurement.log_variance.profile_1.score_b                                0.004816835
#> measurement.log_variance.profile_2.score_a                                0.253448574
#> measurement.log_variance.profile_2.score_b                                0.002601722
#> profile.logit.profile_1.group_class_1                                    -0.061986321
#>                                            measurement.log_variance.profile_2.score_b
#> measurement.mean.profile_1.score_a                                       -0.004856176
#> measurement.mean.profile_1.score_b                                       -0.003090157
#> measurement.mean.profile_2.score_a                                        0.001112444
#> measurement.mean.profile_2.score_b                                        0.021496849
#> measurement.log_variance.profile_1.score_a                               -0.007675686
#> measurement.log_variance.profile_1.score_b                               -0.001833484
#> measurement.log_variance.profile_2.score_a                                0.002601722
#> measurement.log_variance.profile_2.score_b                                0.177928294
#> profile.logit.profile_1.group_class_1                                    -0.017913452
#>                                            profile.logit.profile_1.group_class_1
#> measurement.mean.profile_1.score_a                                   0.024011048
#> measurement.mean.profile_1.score_b                                  -0.004819226
#> measurement.mean.profile_2.score_a                                   0.004075934
#> measurement.mean.profile_2.score_b                                  -0.030236469
#> measurement.log_variance.profile_1.score_a                           0.032006423
#> measurement.log_variance.profile_1.score_b                          -0.005972274
#> measurement.log_variance.profile_2.score_a                          -0.061986321
#> measurement.log_variance.profile_2.score_b                          -0.017913452
#> profile.logit.profile_1.group_class_1                                0.176207352
```
