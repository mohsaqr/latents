# print and summary are stable

    Code
      print(fit)
    Output
      Mixture of gaussian regressions: 2 classes (one class per row)
      score ~ hours
      900 rows | log likelihood -3144.0074 | BIC 6335.63 | entropy 0.443
      Converged: TRUE | iterations: 49 | best likelihood reached by 2 of 2 completed starts (3 run)
      
         class        term estimate std_error   p_value
       class_1 (Intercept)  34.3690    0.7191 0.000e+00
       class_1       hours   4.5635    0.1048 0.000e+00
       class_2 (Intercept)  54.9171    1.0927 0.000e+00
       class_2       hours   0.8046    0.1695 2.067e-06
      
      Every table: get_results(x, what = ), e.g. "classes", "membership", "fit", "assignments".

---

    Code
      print(summary(fit), digits = 3)
    Output
      Model fit
         family     nesting n_classes n_group_classes n_observations n_groups
       gaussian observation         2               1            900       NA
       log_likelihood n_parameters  aic  bic bic_rows sabic  icl entropy
                -3144            7 6302 6336     6336  6313 7030   0.443
       group_entropy smallest_share converged iterations n_starts n_best_replicated
                  NA          0.435      TRUE         49        3                 2
       vcov_type
        observed
      
      Regression coefficients
         class        term estimate std_error statistic  p_value conf_low conf_high
       class_1 (Intercept)   34.369     0.719     47.79 0.00e+00   32.960     35.78
       class_1       hours    4.564     0.105     43.56 0.00e+00    4.358      4.77
       class_2 (Intercept)   54.917     1.093     50.26 0.00e+00   52.776     57.06
       class_2       hours    0.805     0.169      4.75 2.07e-06    0.472      1.14
       p_adjusted
               NA
         0.00e+00
               NA
         2.07e-06
      
      Classes
         class share count n_assigned mean_posterior sigma sigma_std_error prior
       class_1 0.565   508        573          0.789  5.24           0.278 0.565
       class_2 0.435   392        327          0.828  6.89           0.346 0.435
       prior_std_error
                 0.035
                 0.035
      
      Class membership (multinomial logit)
       model   class        term reference estimate std_error statistic p_value
       class class_2 (Intercept)   class_1   -0.261     0.142     -1.83  0.0674
       conf_low conf_high odds_ratio
          -0.54    0.0186      0.771
      

