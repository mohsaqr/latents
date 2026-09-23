# Print a latent transition summary

Print a latent transition summary

## Usage

``` r
# S3 method for class 'summary_multilpa_transitions'
print(x, digits = 4L, rows = 10L, ...)
```

## Arguments

- x:

  A `summary_multilpa_transitions` object.

- digits:

  Number of printed significant digits.

- rows:

  How many rows of each table to print. A longer table is shown to that
  depth, with its remaining row count and the
  [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
  call that returns it whole.

- ...:

  Additional arguments passed to matrix printing.

## Value

The summary, invisibly; called for what it prints.

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
print(summary(fit), digits = 3)
#> Latent transition model: 2 profiles and 1 group class
#> Observations: 150; groups: 30; up to 5 occasions (balanced)
#> Parameters: 11; converged: TRUE
#> Log likelihood: -418.211098; AIC: 858.422
#> BIC (groups): 873.835; BIC (individuals): 891.539
#> Best likelihood replicated in 1/2 starts (absolute tolerance 0.000419).
#> 
#> -- profiles --------------------------------------------------------
#>  profile indicator    mean variance standard_deviation
#>        1   score_a  0.0914    0.780              0.883
#>        1   score_b -0.0902    1.086              1.042
#>        2   score_a  0.7645    1.160              1.077
#>        2   score_b  0.8947    0.508              0.713
#> 
#> -- responses -------------------------------------------------------
#>    (no rows)
#> 
#> -- covariances -----------------------------------------------------
#>  profile indicator indicator_2 covariance
#>        1   score_a     score_a      0.780
#>        1   score_b     score_a      0.000
#>        1   score_a     score_b      0.000
#>        1   score_b     score_b      1.086
#>        2   score_a     score_a      1.160
#>        2   score_b     score_a      0.000
#>        2   score_a     score_b      0.000
#>        2   score_b     score_b      0.508
#> 
#> -- transitions -----------------------------------------------------
#>  group_class from to probability expected_count stable estimated
#>            1    1  1    1.00e+00       1.09e+02   TRUE      TRUE
#>            1    1  2    2.01e-06       1.92e-04  FALSE      TRUE
#>            1    2  1    1.00e-10       7.84e-10  FALSE      TRUE
#>            1    2  2    1.00e+00       1.12e+01   TRUE      TRUE
#>  group_class_probability
#>                        1
#>                        1
#>                        1
#>                        1
#> 
#> -- initial ---------------------------------------------------------
#>  group_class profile probability prevalence group_class_probability
#>            1       1      0.9067     0.9066                       1
#>            1       2      0.0933     0.0934                       1
#> 
#> -- counts ----------------------------------------------------------
#>        level class effective_count effective_proportion
#>  individuals     1             136               0.9066
#>  individuals     2              14               0.0934
#>       groups     1              30               1.0000
#> 
#> -- posteriors ------------------------------------------------------
#>  row group profile posterior modal
#>    1     1       1     0.998  TRUE
#>    2     1       1     0.998  TRUE
#>    3     1       1     0.998  TRUE
#>    4     1       1     0.998  TRUE
#>    5     1       1     0.998  TRUE
#>    6     2       1     0.999  TRUE
#>    7     2       1     0.999  TRUE
#>    8     2       1     0.999  TRUE
#>    9     2       1     0.999  TRUE
#>   10     2       1     0.999  TRUE
#>    ... 290 more rows.  get_results(x, what = "posteriors")
#> 
#> -- group_posteriors ------------------------------------------------
#>  group group_size log_likelihood group_class posterior modal
#>      1          5          -15.6           1         1  TRUE
#>      2          5          -15.5           1         1  TRUE
#>      3          5          -15.3           1         1  TRUE
#>      4          5          -12.2           1         1  TRUE
#>      5          5          -16.8           1         1  TRUE
#>      6          5          -11.6           1         1  TRUE
#>      7          5          -15.7           1         1  TRUE
#>      8          5          -13.0           1         1  TRUE
#>      9          5          -14.4           1         1  TRUE
#>     10          5          -11.4           1         1  TRUE
#>    ... 20 more rows.  get_results(x, what = "group_posteriors")
#> 
#> -- assignments -----------------------------------------------------
#>  person wave score_a score_b profile group_class uncertainty
#>       1    1   2.287 -0.3557       1           1    0.002138
#>       1    2  -1.197  1.0973       1           1    0.002138
#>       1    3  -0.694 -0.9067       1           1    0.002138
#>       1    4  -0.412 -0.2075       1           1    0.002138
#>       1    5  -0.971  0.6789       1           1    0.002140
#>       2    1  -0.947 -0.7978       1           1    0.000699
#>       2    2   0.748 -1.5915       1           1    0.000699
#>       2    3  -0.117  1.1803       1           1    0.000725
#>       2    4   0.153  1.2226       1           1    0.000741
#>       2    5   2.190 -0.0109       1           1    0.000748
#>  posterior_profile_1 posterior_profile_2
#>                0.998            0.002138
#>                0.998            0.002138
#>                0.998            0.002138
#>                0.998            0.002138
#>                0.998            0.002140
#>                0.999            0.000699
#>                0.999            0.000699
#>                0.999            0.000725
#>                0.999            0.000741
#>                0.999            0.000748
#>    ... 140 more rows.  get_results(x, what = "assignments")
#> 
#> -- classification --------------------------------------------------
#>        level class n_modal proportion_modal estimated_n estimated_proportion
#>  individuals     1     135              0.9         136               0.9066
#>  individuals     2      15              0.1          14               0.0934
#>       groups     1      30              1.0          30               1.0000
#>  average_posterior odds_correct_classification
#>              0.992                        12.2
#>              0.858                        58.7
#>              1.000                          NA
#> 
#> -- average_posteriors ----------------------------------------------
#>        level assigned_class class n_assigned average_posterior
#>  individuals              1     1        135           0.99162
#>  individuals              1     2        135           0.00838
#>  individuals              2     1         15           0.14191
#>  individuals              2     2         15           0.85809
#>       groups              1     1         30           1.00000
#> 
#> -- classification_errors -------------------------------------------
#>        level true_class assigned_class probability
#>  individuals          1              1      0.9843
#>  individuals          1              2      0.0157
#>  individuals          2              1      0.0808
#>  individuals          2              2      0.9192
#>       groups          1              1      1.0000
#> 
#> -- bch_weights -----------------------------------------------------
#>        level unit assigned_class class  weight
#>  individuals    1              1     1  1.0173
#>  individuals    1              1     2 -0.0173
#>  individuals    2              1     1  1.0173
#>  individuals    2              1     2 -0.0173
#>  individuals    3              1     1  1.0173
#>  individuals    3              1     2 -0.0173
#>  individuals    4              1     1  1.0173
#>  individuals    4              1     2 -0.0173
#>  individuals    5              1     1  1.0173
#>  individuals    5              1     2 -0.0173
#>    ... 320 more rows.  get_results(x, what = "bch_weights")
#> 
#> -- entropy ---------------------------------------------------------
#>        level n_classes n_units entropy_sum relative_entropy
#>  individuals         2     150        9.27            0.911
#>       groups         1      30        0.00               NA
#> 
#> -- residuals -------------------------------------------------------
#>    profile indicator_1 indicator_2     kind observed expected residual
#>  profile_2     score_a     score_b gaussian  -0.2709        0  -0.2709
#>  profile_1     score_a     score_b gaussian  -0.0895        0  -0.0895
#>  effective_n statistic df p_value p_adjusted
#>           14    -0.921 NA   0.357      0.357
#>          136    -1.035 NA   0.301      0.301
#> 
#> -- information_criteria --------------------------------------------
#>  log_likelihood n_parameters aic kic bic_groups bic_individual sabic_groups
#>            -418           11 858 872        874            892          840
#>  sabic_individual caic_groups caic_individual awe_groups awe_individual
#>               857         885             903        944            998
#>  icl_groups icl_individual clc_groups clc_individual
#>         874            910        836            855
#> 
#> -- model -----------------------------------------------------------
#>  n_observations n_informative n_groups n_profiles n_group_classes centering
#>             150           150       30          2               1      none
#>  covariance_structure n_parameters n_parameters_with_measurement log_likelihood
#>                  <NA>           11                            11           -418
#>  aic bic_groups bic_individual converged iterations boundary small_classes
#>  858        874            892      TRUE        116    FALSE         FALSE
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
#>      1           -418      TRUE        222  <NA>
#>      2           -418      TRUE        116  <NA>
#> 
#> -- data ------------------------------------------------------------
#>  person wave score_a score_b
#>       1    1   2.287 -0.3557
#>       1    2  -1.197  1.0973
#>       1    3  -0.694 -0.9067
#>       1    4  -0.412 -0.2075
#>       1    5  -0.971  0.6789
#>       2    1  -0.947 -0.7978
#>       2    2   0.748 -1.5915
#>       2    3  -0.117  1.1803
#>       2    4   0.153  1.2226
#>       2    5   2.190 -0.0109
#>    ... 140 more rows.  get_results(x, what = "data")
#> 
#> 22 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
```
