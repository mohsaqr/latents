# print and summary are stable

    Code
      print(fit)
    Output
      Mixture of gaussian regressions: 2 classes (one class per row)
      score ~ hours
      900 rows | log likelihood -3144.0073 | BIC 6335.63 | entropy 0.443
      Converged: TRUE | iterations: 53 | best likelihood reached by 2 of 2 completed starts (3 run)
      
      Classes
      
      Class    Share  Expected  Rows assigned  Avg. posterior  Residual SD
      -------  -----  --------  -------------  --------------  -----------
      Class 1  0.564    507.71            573           0.789         5.24
      Class 2  0.436    392.29            327           0.829         6.89
      
      Regression coefficients (95% CI)
      
      Class    Term       Estimate          95% CI      p
      -------  ---------  --------  --------------  -----
      Class 1  Intercept     34.36  [32.95, 35.77]  <.001
      Class 1  hours          4.57  [ 4.36,  4.77]  <.001
      Class 2  Intercept     54.90  [52.76, 57.04]  <.001
      Class 2  hours          0.81  [ 0.48,  1.14]  <.001
      
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
      ICL             7030.54
      Entropy         0.443
      Smallest class  43.6%
      Converged       yes
      
      Regression coefficients (95% CI)
      
      Class    Term       Estimate          95% CI      p
      -------  ---------  --------  --------------  -----
      Class 1  Intercept     34.36  [32.95, 35.77]  <.001
      Class 1  hours          4.57  [ 4.36,  4.77]  <.001
      Class 2  Intercept     54.90  [52.76, 57.04]  <.001
      Class 2  hours          0.81  [ 0.48,  1.14]  <.001
      
      Classes
      
      Class    Share  Expected  Rows assigned  Avg. posterior  Residual SD
      -------  -----  --------  -------------  --------------  -----------
      Class 1  0.564    507.71            573           0.789         5.24
      Class 2  0.436    392.29            327           0.829         6.89
      
      Class membership (log odds, 95% CI)
      
      Class    Term       Estimate         95% CI     p  Odds ratio
      -------  ---------  --------  -------------  ----  ----------
      Class 2  Intercept     -0.26  [-0.54, 0.02]  .070        0.77
      

