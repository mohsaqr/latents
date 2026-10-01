# Summarize a fitted latent transition model

Summarize a fitted latent transition model

## Usage

``` r
# S3 method for class 'multilpa_transitions'
summary(object, ...)
```

## Arguments

- object:

  A fitted `multilpa_transitions` model.

- ...:

  Reserved for compatibility with
  [`summary()`](https://rdrr.io/r/base/summary.html).

## Value

A `summary_multilpa_transitions` object carrying the measurement
parameters, the initial and transition probabilities and counts, the
effective class counts, the fit statistics and the restart diagnostics.
Its tables are read with
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md);
[`print()`](https://rdrr.io/r/base/print.html) reports the whole model.

## See also

[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
for the transition probabilities as a tidy table.

## Examples

``` r
fit <- lta(subset(course_engagement, student <= 40),
           c("browse", "lectures", "forum_read"), "student",
           n_profiles = 2, time = "sequence", n_starts = 2, seed = 1)
summary(fit)
#> Latent transition model: 2 profiles and 1 group class
#> Observations: 547; groups: 40; up to 15 occasions (unbalanced, observed grid)
#> Parameters: 15; converged: TRUE
#> Log likelihood: -2025.406100; AIC: 4080.812
#> BIC (groups): 4106.145; BIC (individuals): 4145.379
#> Best likelihood replicated in 2/2 starts (absolute tolerance 0.00203).
#> 
#> -- profiles --------------------------------------------------------
#>  profile  indicator    mean variance standard_deviation mean_standard_error
#>        1     browse  0.5552   0.5491             0.7410             0.04015
#>        1   lectures  0.4372   0.8253             0.9085             0.04891
#>        1 forum_read  0.6089   0.4991             0.7065             0.03909
#>        2     browse -0.8706   0.4602             0.6784             0.05194
#>        2   lectures -0.6265   0.5736             0.7574             0.05713
#>        2 forum_read -0.9320   0.3260             0.5709             0.04356
#>  variance_standard_error
#>                  0.04216
#>                  0.06254
#>                  0.03955
#>                  0.04940
#>                  0.06114
#>                  0.03540
#> 
#> -- responses -------------------------------------------------------
#>    (no rows)
#> 
#> -- covariances -----------------------------------------------------
#>  profile  indicator indicator_2 covariance
#>        1     browse      browse     0.5491
#>        1   lectures      browse     0.0000
#>        1 forum_read      browse     0.0000
#>        1     browse    lectures     0.0000
#>        1   lectures    lectures     0.8253
#>        1 forum_read    lectures     0.0000
#>        1     browse  forum_read     0.0000
#>        1   lectures  forum_read     0.0000
#>        1 forum_read  forum_read     0.4991
#>        2     browse      browse     0.4602
#>    ... 8 more rows.  get_results(x, what = "covariances")
#> 
#> -- transitions -----------------------------------------------------
#>  group_class from to probability expected_count stable estimated
#>            1    1  1      0.8906         297.42   TRUE      TRUE
#>            1    1  2      0.1094          36.53  FALSE      TRUE
#>            1    2  1      0.1896          32.81  FALSE      TRUE
#>            1    2  2      0.8104         140.24   TRUE      TRUE
#>  group_class_probability
#>                        1
#>                        1
#>                        1
#>                        1
#> 
#> -- initial ---------------------------------------------------------
#>  group_class profile probability prevalence group_class_probability
#>            1       1      0.7421      0.658                       1
#>            1       2      0.2579      0.342                       1
#> 
#> -- counts ----------------------------------------------------------
#>        level class effective_count effective_proportion
#>  individuals     1           359.9                0.658
#>  individuals     2           187.1                0.342
#>       groups     1            40.0                1.000
#> 
#> -- posteriors ------------------------------------------------------
#>  row group profile posterior modal
#>    1     1       1 0.9989080  TRUE
#>    2     1       1 0.9856001  TRUE
#>    3     1       1 0.0012923 FALSE
#>    4     1       1 0.0010969 FALSE
#>    5     1       1 0.0003721 FALSE
#>    6     1       1 0.0012699 FALSE
#>    7     1       1 0.0021079 FALSE
#>    8     1       1 0.0007051 FALSE
#>    9     1       1 0.0002443 FALSE
#>   10     1       1 0.0099241 FALSE
#>    ... 1084 more rows.  get_results(x, what = "posteriors")
#> 
#> -- group_posteriors ------------------------------------------------
#>  group group_size log_likelihood group_class posterior modal
#>      1         13         -39.79           1         1  TRUE
#>      2         13         -45.17           1         1  TRUE
#>      3         14         -46.96           1         1  TRUE
#>      4         13         -42.33           1         1  TRUE
#>      5         14         -45.73           1         1  TRUE
#>      6         12         -51.34           1         1  TRUE
#>      7         15         -58.06           1         1  TRUE
#>      8         14         -52.27           1         1  TRUE
#>      9         14         -56.47           1         1  TRUE
#>     10         15         -57.90           1         1  TRUE
#>    ... 30 more rows.  get_results(x, what = "group_posteriors")
#> 
#> -- assignments -----------------------------------------------------
#>  student sequence browse lectures forum_read profile group_class uncertainty
#>        1        1   0.73    -0.08       0.30       1           1   0.0010920
#>        1        2   0.68     0.53      -0.02       1           1   0.0143999
#>        1        3  -1.51    -0.51      -0.91       2           1   0.0012923
#>        1        4   0.64    -1.02      -1.62       2           1   0.0010969
#>        1        5  -0.94    -0.37      -0.71       2           1   0.0003721
#>        1        6  -0.36    -0.46      -0.77       2           1   0.0012699
#>        1        7   0.38    -1.85      -0.84       2           1   0.0021079
#>        1        8  -1.44     0.66      -0.67       2           1   0.0007051
#>        1        9  -0.76    -0.43      -1.04       2           1   0.0002443
#>        1       10   0.26    -0.95      -0.51       2           1   0.0099241
#>  posterior_profile_1 posterior_profile_2
#>            0.9989080            0.001092
#>            0.9856001            0.014400
#>            0.0012923            0.998708
#>            0.0010969            0.998903
#>            0.0003721            0.999628
#>            0.0012699            0.998730
#>            0.0021079            0.997892
#>            0.0007051            0.999295
#>            0.0002443            0.999756
#>            0.0099241            0.990076
#>    ... 537 more rows.  get_results(x, what = "assignments")
#> 
#> -- classification --------------------------------------------------
#>        level class n_modal proportion_modal estimated_n estimated_proportion
#>  individuals     1     358           0.6545       359.9                0.658
#>  individuals     2     189           0.3455       187.1                0.342
#>       groups     1      40           1.0000        40.0                1.000
#>  average_posterior odds_correct_classification
#>             0.9884                       44.28
#>             0.9679                       58.01
#>             1.0000                          NA
#> 
#> -- average_posteriors ----------------------------------------------
#>        level assigned_class class n_assigned average_posterior
#>  individuals              1     1        358            0.9884
#>  individuals              1     2        358            0.0116
#>  individuals              2     1        189            0.0321
#>  individuals              2     2        189            0.9679
#>       groups              1     1         40            1.0000
#> 
#> -- classification_errors -------------------------------------------
#>        level true_class assigned_class probability
#>  individuals          1              1     0.98314
#>  individuals          1              2     0.01686
#>  individuals          2              1     0.02220
#>  individuals          2              2     0.97780
#>       groups          1              1     1.00000
#> 
#> -- bch_weights -----------------------------------------------------
#>        level unit assigned_class class   weight
#>  individuals    1              1     1  1.01754
#>  individuals    1              1     2 -0.01754
#>  individuals    2              1     1  1.01754
#>  individuals    2              1     2 -0.01754
#>  individuals    3              2     1 -0.02310
#>  individuals    3              2     2  1.02310
#>  individuals    4              2     1 -0.02310
#>  individuals    4              2     2  1.02310
#>  individuals    5              2     1 -0.02310
#>  individuals    5              2     2  1.02310
#>    ... 1124 more rows.  get_results(x, what = "bch_weights")
#> 
#> -- entropy ---------------------------------------------------------
#>        level n_classes n_units entropy_sum relative_entropy
#>  individuals         2     547       29.14           0.9231
#>       groups         1      40        0.00               NA
#> 
#> -- residuals -------------------------------------------------------
#>    profile indicator_1 indicator_2     kind observed expected residual
#>  profile_2      browse    lectures gaussian -0.15079        0 -0.15079
#>  profile_2      browse  forum_read gaussian  0.12175        0  0.12175
#>  profile_1    lectures  forum_read gaussian  0.07085        0  0.07085
#>  profile_1      browse  forum_read gaussian  0.04309        0  0.04309
#>  profile_2    lectures  forum_read gaussian -0.03844        0 -0.03844
#>  profile_1      browse    lectures gaussian  0.03703        0  0.03703
#>  effective_n statistic df p_value p_adjusted
#>        187.1   -2.0616 NA 0.03925    0.03925
#>        187.1    1.6601 NA 0.09689    0.09689
#>        359.9    1.3407 NA 0.18003    0.18003
#>        359.9    0.8145 NA 0.41537    0.41537
#>        187.1   -0.5218 NA 0.60182    0.60182
#>        359.9    0.6998 NA 0.48402    0.48402
#> 
#> -- information_criteria --------------------------------------------
#>  log_likelihood n_parameters  aic  kic bic_groups bic_individual sabic_groups
#>           -2025           15 4081 4099       4106           4145         4059
#>  sabic_individual caic_groups caic_individual awe_groups awe_individual
#>              4098        4121            4160       4206           4343
#>  icl_groups icl_individual clc_groups clc_individual
#>        4106           4204       4051           4109
#> 
#> -- model -----------------------------------------------------------
#>  n_observations n_informative n_groups n_profiles n_group_classes centering
#>             547           547       40          2               1      none
#>  covariance_structure n_parameters n_parameters_with_measurement log_likelihood
#>                  <NA>           15                            15          -2025
#>   aic bic_groups bic_individual converged iterations boundary small_classes
#>  4081       4106           4145      TRUE          8    FALSE         FALSE
#>  best_start n_best_replicated weights
#>           1                 2    <NA>
#> 
#> -- sequences -------------------------------------------------------
#>  group group_class time profile
#>      1           1    1       1
#>      1           1    2       1
#>      1           1    3       2
#>      1           1    4       2
#>      1           1    5       2
#>      1           1    6       2
#>      1           1    7       2
#>      1           1    8       2
#>      1           1    9       2
#>      1           1   10       2
#>    ... 537 more rows.  get_results(x, what = "sequences")
#> 
#> -- sequence_summary ------------------------------------------------
#>  group_class groups observations mean_length median_length shortest longest
#>            1     40          547       13.68            14       11      15
#>  complete gaps
#>        10    0
#> 
#> -- sequence_lengths ------------------------------------------------
#>  group group_class occasions observations complete
#>      1           1        13           13    FALSE
#>      2           1        13           13    FALSE
#>      3           1        14           14    FALSE
#>      4           1        13           13    FALSE
#>      5           1        14           14    FALSE
#>      6           1        12           12    FALSE
#>      7           1        15           15     TRUE
#>      8           1        14           14    FALSE
#>      9           1        14           14    FALSE
#>     10           1        15           15     TRUE
#>    ... 30 more rows.  get_results(x, what = "sequence_lengths")
#> 
#> -- starts ----------------------------------------------------------
#>  start log_likelihood converged iterations error
#>      1          -2025      TRUE          8  <NA>
#>      2          -2025      TRUE          9  <NA>
#> 
#> -- data ------------------------------------------------------------
#>  student sequence browse lectures forum_read
#>        1        1   0.73    -0.08       0.30
#>        1        2   0.68     0.53      -0.02
#>        1        3  -1.51    -0.51      -0.91
#>        1        4   0.64    -1.02      -1.62
#>        1        5  -0.94    -0.37      -0.71
#>        1        6  -0.36    -0.46      -0.77
#>        1        7   0.38    -1.85      -0.84
#>        1        8  -1.44     0.66      -0.67
#>        1        9  -0.76    -0.43      -1.04
#>        1       10   0.26    -0.95      -0.51
#>    ... 537 more rows.  get_results(x, what = "data")
#> 
#> 22 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
```
