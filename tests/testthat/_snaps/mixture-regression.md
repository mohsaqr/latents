# print and summary are stable

    Code
      print(fit)
    Output
      Mixture of gaussian regressions: 2 classes (one class per row)
      score ~ hours
      900 rows | log likelihood -3144.0074 | BIC 6335.63 | entropy 0.443
      Converged: TRUE | iterations: 49 | best likelihood reached by 2 of 2 completed starts (3 run)
      
      Classes
      
      Class    Share  Expected  Rows assigned  Avg. posterior  Residual SD
      -------  -----  --------  -------------  --------------  -----------
      Class 1  0.565    508.26            573           0.789         5.24
      Class 2  0.435    391.74            327           0.828         6.89
      
      Regression coefficients (95% CI)
      
      Class    Term       Estimate          95% CI      p
      -------  ---------  --------  --------------  -----
      Class 1  Intercept     34.37  [32.96, 35.78]  <.001
      Class 1  hours          4.56  [ 4.36,  4.77]  <.001
      Class 2  Intercept     54.92  [52.78, 57.06]  <.001
      Class 2  hours          0.80  [ 0.47,  1.14]  <.001
      
      Every table: get_results(x, what = ), e.g. "classes", "membership", "fit", "assignments".

---

    Code
      print(summary(fit), digits = 3)
    Output
      Model fit
      
      Family          gaussian
      Nesting         observation
      Classes         2
      Group classes   1
      Observations    900
      Parameters      7
      Log likelihood  -3144.01
      AIC             6302.01
      BIC             6335.63
      SABIC           6313.40
      ICL             7030.37
      Entropy         0.443
      Smallest class  43.5%
      Converged       yes
      
      Regression coefficients (95% CI)
      
      Class    Term       Estimate          95% CI      p
      -------  ---------  --------  --------------  -----
      Class 1  Intercept     34.37  [32.96, 35.78]  <.001
      Class 1  hours          4.56  [ 4.36,  4.77]  <.001
      Class 2  Intercept     54.92  [52.78, 57.06]  <.001
      Class 2  hours          0.80  [ 0.47,  1.14]  <.001
      
      Classes
      
      Class    Share  Expected  Rows assigned  Avg. posterior  Residual SD
      -------  -----  --------  -------------  --------------  -----------
      Class 1  0.565    508.26            573           0.789         5.24
      Class 2  0.435    391.74            327           0.828         6.89
      
      Class membership (log odds, 95% CI)
      
      Class    Term       Estimate         95% CI     p  Odds ratio
      -------  ---------  --------  -------------  ----  ----------
      Class 2  Intercept     -0.26  [-0.54, 0.02]  .067        0.77
      

