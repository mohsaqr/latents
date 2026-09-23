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
set.seed(7)
example_data <- data.frame(
  person = rep(seq_len(30), each = 5), wave = rep(seq_len(5), times = 30)
)
example_data$score_a <- stats::rnorm(nrow(example_data))
example_data$score_b <- stats::rnorm(nrow(example_data))
fit <- lta(example_data, c("score_a", "score_b"), "person",
                       n_profiles = 2, time = "wave", n_starts = 2, seed = 1)
summary(fit)
#> Latent transition model: 2 profiles and 1 group class
#> Observations: 150; groups: 30; up to 5 occasions (balanced)
#> Parameters: 11; converged: TRUE
#> Log likelihood: -418.211098; AIC: 858.422
#> BIC (groups): 873.835; BIC (individuals): 891.539
#> Best likelihood replicated in 1/2 starts (absolute tolerance 0.000419).
#> 
#> -- profiles --------------------------------------------------------
#>  profile indicator     mean variance standard_deviation
#>        1   score_a  0.09137   0.7804             0.8834
#>        1   score_b -0.09015   1.0865             1.0423
#>        2   score_a  0.76448   1.1599             1.0770
#>        2   score_b  0.89469   0.5083             0.7129
#> 
#> -- responses -------------------------------------------------------
#>    (no rows)
#> 
#> -- covariances -----------------------------------------------------
#>  profile indicator indicator_2 covariance
#>        1   score_a     score_a     0.7804
#>        1   score_b     score_a     0.0000
#>        1   score_a     score_b     0.0000
#>        1   score_b     score_b     1.0865
#>        2   score_a     score_a     1.1599
#>        2   score_b     score_a     0.0000
#>        2   score_a     score_b     0.0000
#>        2   score_b     score_b     0.5083
#> 
#> -- transitions -----------------------------------------------------
#>  group_class from to probability expected_count stable estimated
#>            1    1  1   1.000e+00      1.088e+02   TRUE      TRUE
#>            1    1  2   2.013e-06      1.922e-04  FALSE      TRUE
#>            1    2  1   1.000e-10      7.841e-10  FALSE      TRUE
#>            1    2  2   1.000e+00      1.120e+01   TRUE      TRUE
#>  group_class_probability
#>                        1
#>                        1
#>                        1
#>                        1
#> 
#> -- initial ---------------------------------------------------------
#>  group_class profile probability prevalence group_class_probability
#>            1       1     0.90665    0.90665                       1
#>            1       2     0.09335    0.09335                       1
#> 
#> -- counts ----------------------------------------------------------
#>        level class effective_count effective_proportion
#>  individuals     1             136              0.90665
#>  individuals     2              14              0.09335
#>       groups     1              30              1.00000
#> 
#> -- posteriors ------------------------------------------------------
#>  row group profile posterior modal
#>    1     1       1    0.9979  TRUE
#>    2     1       1    0.9979  TRUE
#>    3     1       1    0.9979  TRUE
#>    4     1       1    0.9979  TRUE
#>    5     1       1    0.9979  TRUE
#>    6     2       1    0.9993  TRUE
#>    7     2       1    0.9993  TRUE
#>    8     2       1    0.9993  TRUE
#>    9     2       1    0.9993  TRUE
#>   10     2       1    0.9993  TRUE
#>    ... 290 more rows.  get_results(x, what = "posteriors")
#> 
#> -- group_posteriors ------------------------------------------------
#>  group group_size log_likelihood group_class posterior modal
#>      1          5         -15.57           1         1  TRUE
#>      2          5         -15.50           1         1  TRUE
#>      3          5         -15.31           1         1  TRUE
#>      4          5         -12.20           1         1  TRUE
#>      5          5         -16.83           1         1  TRUE
#>      6          5         -11.56           1         1  TRUE
#>      7          5         -15.67           1         1  TRUE
#>      8          5         -13.03           1         1  TRUE
#>      9          5         -14.41           1         1  TRUE
#>     10          5         -11.40           1         1  TRUE
#>    ... 20 more rows.  get_results(x, what = "group_posteriors")
#> 
#> -- assignments -----------------------------------------------------
#>  person wave score_a  score_b profile group_class uncertainty
#>       1    1  2.2872 -0.35569       1           1   0.0021375
#>       1    2 -1.1968  1.09730       1           1   0.0021375
#>       1    3 -0.6943 -0.90669       1           1   0.0021376
#>       1    4 -0.4123 -0.20746       1           1   0.0021380
#>       1    5 -0.9707  0.67886       1           1   0.0021397
#>       2    1 -0.9473 -0.79779       1           1   0.0006987
#>       2    2  0.7481 -1.59154       1           1   0.0006990
#>       2    3 -0.1170  1.18035       1           1   0.0007253
#>       2    4  0.1527  1.22257       1           1   0.0007406
#>       2    5  2.1900 -0.01091       1           1   0.0007482
#>  posterior_profile_1 posterior_profile_2
#>               0.9979           0.0021375
#>               0.9979           0.0021375
#>               0.9979           0.0021376
#>               0.9979           0.0021380
#>               0.9979           0.0021397
#>               0.9993           0.0006987
#>               0.9993           0.0006990
#>               0.9993           0.0007253
#>               0.9993           0.0007406
#>               0.9993           0.0007482
#>    ... 140 more rows.  get_results(x, what = "assignments")
#> 
#> -- classification --------------------------------------------------
#>        level class n_modal proportion_modal estimated_n estimated_proportion
#>  individuals     1     135              0.9         136              0.90665
#>  individuals     2      15              0.1          14              0.09335
#>       groups     1      30              1.0          30              1.00000
#>  average_posterior odds_correct_classification
#>             0.9916                       12.18
#>             0.8581                       58.73
#>             1.0000                          NA
#> 
#> -- average_posteriors ----------------------------------------------
#>        level assigned_class class n_assigned average_posterior
#>  individuals              1     1        135          0.991619
#>  individuals              1     2        135          0.008381
#>  individuals              2     1         15          0.141912
#>  individuals              2     2         15          0.858088
#>       groups              1     1         30          1.000000
#> 
#> -- classification_errors -------------------------------------------
#>        level true_class assigned_class probability
#>  individuals          1              1     0.98435
#>  individuals          1              2     0.01565
#>  individuals          2              1     0.08080
#>  individuals          2              2     0.91920
#>       groups          1              1     1.00000
#> 
#> -- bch_weights -----------------------------------------------------
#>        level unit assigned_class class   weight
#>  individuals    1              1     1  1.01732
#>  individuals    1              1     2 -0.01732
#>  individuals    2              1     1  1.01732
#>  individuals    2              1     2 -0.01732
#>  individuals    3              1     1  1.01732
#>  individuals    3              1     2 -0.01732
#>  individuals    4              1     1  1.01732
#>  individuals    4              1     2 -0.01732
#>  individuals    5              1     1  1.01732
#>  individuals    5              1     2 -0.01732
#>    ... 320 more rows.  get_results(x, what = "bch_weights")
#> 
#> -- entropy ---------------------------------------------------------
#>        level n_classes n_units entropy_sum relative_entropy
#>  individuals         2     150       9.271           0.9108
#>       groups         1      30       0.000               NA
#> 
#> -- residuals -------------------------------------------------------
#>    profile indicator_1 indicator_2     kind observed expected residual
#>  profile_2     score_a     score_b gaussian -0.27086        0 -0.27086
#>  profile_1     score_a     score_b gaussian -0.08953        0 -0.08953
#>  effective_n statistic df p_value p_adjusted
#>           14   -0.9214 NA  0.3568     0.3568
#>          136   -1.0353 NA  0.3005     0.3005
#> 
#> -- information_criteria --------------------------------------------
#>  log_likelihood n_parameters   aic   kic bic_groups bic_individual sabic_groups
#>          -418.2           11 858.4 872.4      873.8          891.5        839.6
#>  sabic_individual caic_groups caic_individual awe_groups awe_individual
#>             856.7       884.8           902.5      944.2          998.2
#>  icl_groups icl_individual clc_groups clc_individual
#>       873.8          910.1      836.4            855
#> 
#> -- model -----------------------------------------------------------
#>  n_observations n_informative n_groups n_profiles n_group_classes centering
#>             150           150       30          2               1      none
#>  covariance_structure n_parameters n_parameters_with_measurement log_likelihood
#>                  <NA>           11                            11         -418.2
#>    aic bic_groups bic_individual converged iterations boundary small_classes
#>  858.4      873.8          891.5      TRUE        116    FALSE         FALSE
#>  best_start n_best_replicated
#>           2                 1
#> 
#> -- sequences -------------------------------------------------------
#>  group group_class time profile
#>      1           1    1       1
#>      1           1    2       1
#>      1           1    3       1
#>      1           1    4       1
#>      1           1    5       1
#>      2           1    1       1
#>      2           1    2       1
#>      2           1    3       1
#>      2           1    4       1
#>      2           1    5       1
#>    ... 140 more rows.  get_results(x, what = "sequences")
#> 
#> -- sequence_summary ------------------------------------------------
#>  group_class groups observations mean_length median_length shortest longest
#>            1     30          150           5             5        5       5
#>  complete gaps
#>        30    0
#> 
#> -- sequence_lengths ------------------------------------------------
#>  group group_class occasions observations complete
#>      1           1         5            5     TRUE
#>      2           1         5            5     TRUE
#>      3           1         5            5     TRUE
#>      4           1         5            5     TRUE
#>      5           1         5            5     TRUE
#>      6           1         5            5     TRUE
#>      7           1         5            5     TRUE
#>      8           1         5            5     TRUE
#>      9           1         5            5     TRUE
#>     10           1         5            5     TRUE
#>    ... 20 more rows.  get_results(x, what = "sequence_lengths")
#> 
#> -- starts ----------------------------------------------------------
#>  start log_likelihood converged iterations error
#>      1         -418.3      TRUE        222  <NA>
#>      2         -418.2      TRUE        116  <NA>
#> 
#> -- data ------------------------------------------------------------
#>  person wave score_a  score_b
#>       1    1  2.2872 -0.35569
#>       1    2 -1.1968  1.09730
#>       1    3 -0.6943 -0.90669
#>       1    4 -0.4123 -0.20746
#>       1    5 -0.9707  0.67886
#>       2    1 -0.9473 -0.79779
#>       2    2  0.7481 -1.59154
#>       2    3 -0.1170  1.18035
#>       2    4  0.1527  1.22257
#>       2    5  2.1900 -0.01091
#>    ... 140 more rows.  get_results(x, what = "data")
#> 
#> 22 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
```
