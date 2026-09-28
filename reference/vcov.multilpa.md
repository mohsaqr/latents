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
  diagonal fit, log-Cholesky coordinates for a full-covariance fit, the
  structure's own log volumes, log shapes and orientations for one of
  the ten constrained covariance structures, and baseline-category
  logits for the mixing and response probabilities.

- ...:

  Additional arguments passed to
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md),
  including `method = "bootstrap"`.

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
#> measurement.mean.profile_1.score_a                               0.0153921482
#> measurement.mean.profile_1.score_b                              -0.0010576622
#> measurement.mean.profile_2.score_a                               0.0005338484
#> measurement.mean.profile_2.score_b                              -0.0079230067
#> measurement.log_variance.profile_1.score_a                       0.0083280983
#> measurement.log_variance.profile_1.score_b                      -0.0015442986
#> measurement.log_variance.profile_2.score_a                      -0.0145311011
#> measurement.log_variance.profile_2.score_b                      -0.0048585318
#> profile.logit.profile_1.group_class_1                            0.0240220436
#>                                            measurement.mean.profile_1.score_b
#> measurement.mean.profile_1.score_a                              -0.0010576622
#> measurement.mean.profile_1.score_b                               0.0176688544
#> measurement.mean.profile_2.score_a                              -0.0007025632
#> measurement.mean.profile_2.score_b                              -0.0034556566
#> measurement.log_variance.profile_1.score_a                      -0.0010480052
#> measurement.log_variance.profile_1.score_b                      -0.0005801027
#> measurement.log_variance.profile_2.score_a                       0.0033363737
#> measurement.log_variance.profile_2.score_b                      -0.0030918745
#> profile.logit.profile_1.group_class_1                           -0.0048239988
#>                                            measurement.mean.profile_2.score_a
#> measurement.mean.profile_1.score_a                               0.0005338484
#> measurement.mean.profile_1.score_b                              -0.0007025632
#> measurement.mean.profile_2.score_a                               0.0047158827
#> measurement.mean.profile_2.score_b                               0.0007560018
#> measurement.log_variance.profile_1.score_a                      -0.0006092087
#> measurement.log_variance.profile_1.score_b                       0.0000693490
#> measurement.log_variance.profile_2.score_a                      -0.0052062937
#> measurement.log_variance.profile_2.score_b                       0.0011165477
#> profile.logit.profile_1.group_class_1                            0.0040793551
#>                                            measurement.mean.profile_2.score_b
#> measurement.mean.profile_1.score_a                              -0.0079230067
#> measurement.mean.profile_1.score_b                              -0.0034556566
#> measurement.mean.profile_2.score_a                               0.0007560018
#> measurement.mean.profile_2.score_b                               0.0704854750
#> measurement.log_variance.profile_1.score_a                      -0.0120060866
#> measurement.log_variance.profile_1.score_b                       0.0053598746
#> measurement.log_variance.profile_2.score_a                       0.0179839797
#> measurement.log_variance.profile_2.score_b                       0.0215099790
#> profile.logit.profile_1.group_class_1                           -0.0302743414
#>                                            measurement.log_variance.profile_1.score_a
#> measurement.mean.profile_1.score_a                                       0.0083280983
#> measurement.mean.profile_1.score_b                                      -0.0010480052
#> measurement.mean.profile_2.score_a                                      -0.0006092087
#> measurement.mean.profile_2.score_b                                      -0.0120060866
#> measurement.log_variance.profile_1.score_a                               0.0437111896
#> measurement.log_variance.profile_1.score_b                              -0.0023077078
#> measurement.log_variance.profile_2.score_a                              -0.0190045620
#> measurement.log_variance.profile_2.score_b                              -0.0076773397
#> profile.logit.profile_1.group_class_1                                    0.0320028524
#>                                            measurement.log_variance.profile_1.score_b
#> measurement.mean.profile_1.score_a                                      -0.0015442986
#> measurement.mean.profile_1.score_b                                      -0.0005801027
#> measurement.mean.profile_2.score_a                                       0.0000693490
#> measurement.mean.profile_2.score_b                                       0.0053598746
#> measurement.log_variance.profile_1.score_a                              -0.0023077078
#> measurement.log_variance.profile_1.score_b                               0.0331133909
#> measurement.log_variance.profile_2.score_a                               0.0048108074
#> measurement.log_variance.profile_2.score_b                              -0.0018309343
#> profile.logit.profile_1.group_class_1                                   -0.0059768471
#>                                            measurement.log_variance.profile_2.score_a
#> measurement.mean.profile_1.score_a                                       -0.014531101
#> measurement.mean.profile_1.score_b                                        0.003336374
#> measurement.mean.profile_2.score_a                                       -0.005206294
#> measurement.mean.profile_2.score_b                                        0.017983980
#> measurement.log_variance.profile_1.score_a                               -0.019004562
#> measurement.log_variance.profile_1.score_b                                0.004810807
#> measurement.log_variance.profile_2.score_a                                0.253265379
#> measurement.log_variance.profile_2.score_b                                0.002565869
#> profile.logit.profile_1.group_class_1                                    -0.061945802
#>                                            measurement.log_variance.profile_2.score_b
#> measurement.mean.profile_1.score_a                                       -0.004858532
#> measurement.mean.profile_1.score_b                                       -0.003091875
#> measurement.mean.profile_2.score_a                                        0.001116548
#> measurement.mean.profile_2.score_b                                        0.021509979
#> measurement.log_variance.profile_1.score_a                               -0.007677340
#> measurement.log_variance.profile_1.score_b                               -0.001830934
#> measurement.log_variance.profile_2.score_a                                0.002565869
#> measurement.log_variance.profile_2.score_b                                0.178020487
#> profile.logit.profile_1.group_class_1                                    -0.017928794
#>                                            profile.logit.profile_1.group_class_1
#> measurement.mean.profile_1.score_a                                   0.024022044
#> measurement.mean.profile_1.score_b                                  -0.004823999
#> measurement.mean.profile_2.score_a                                   0.004079355
#> measurement.mean.profile_2.score_b                                  -0.030274341
#> measurement.log_variance.profile_1.score_a                           0.032002852
#> measurement.log_variance.profile_1.score_b                          -0.005976847
#> measurement.log_variance.profile_2.score_a                          -0.061945802
#> measurement.log_variance.profile_2.score_b                          -0.017928794
#> profile.logit.profile_1.group_class_1                                0.176337006
```
