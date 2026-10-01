# Latent transition extensions: external comparison

Run 2026-10-01; R version 4.5.2 (2025-10-31). **30 of 30 quantities agree.**

| Case | Quantity | latents | Reference | Difference | Tolerance | Agree |
|---|---|---|---|---|---|---|
| depmix_covariate | log likelihood | -2152.25 | -2152.25 | 1.12e-07 | 1e-04 | yes |
| depmix_covariate | parameter count | 13 | 13 | 0.00e+00 | 0e+00 | yes |
| depmix_covariate | max |conditional transition diff| | 7.07646e-06 | 0 | 7.08e-06 | 1e-03 | yes |
| depmix_covariate | max |abs slope diff| | 4.01227e-05 | 0 | 4.01e-05 | 1e-03 | yes |
| depmix_covariate | max relative slope SE diff | 0.00160459 | 0 | 1.60e-03 | 5e-02 | yes |
| ex8.13 | log likelihood (+ known-class term) | -14596.2 | -14596.2 | 3.57e-04 | 5e-03 | yes |
| ex8.13 | parameter count (+ known-class proportion) | 32 | 32 | 0.00e+00 | 0e+00 | yes |
| ex8.13 | max |threshold diff| | 0.000546869 | 0 | 5.47e-04 | 2e-03 | yes |
| ex8.13 | max |C1 ON CG#1 diff| | 0.000664825 | 0 | 6.65e-04 | 2e-03 | yes |
| ex8.14 | log likelihood | -12902.1 | -12902.1 | 2.70e-04 | 5e-03 | yes |
| ex8.14 | parameter count | 31 | 31 | 0.00e+00 | 0e+00 | yes |
| ex8.14 | max |threshold diff| | 0.000421762 | 0 | 4.22e-04 | 2e-03 | yes |
| ex8.14 | max |initial logit on x diff| | 0.00039306 | 0 | 3.93e-04 | 2e-03 | yes |
| ex8.14 | max |log odds ratio diff| | 0.000275804 | 0 | 2.76e-04 | 5e-03 | yes |
| lmest_A | log likelihood | -3282.79 | -3282.79 | 1.23e-11 | 1e-04 | yes |
| lmest_A | parameter count | 17 | 17 | 0.00e+00 | 0e+00 | yes |
| lmest_A | max |transition diff| per occasion | 1.706e-08 | 0 | 1.71e-08 | 1e-04 | yes |
| lmest_A | max |initial diff| | 4.75197e-09 | 0 | 4.75e-09 | 1e-04 | yes |
| lmest_B | log likelihood | -2958.85 | -2958.85 | 4.06e-08 | 1e-04 | yes |
| lmest_B | parameter count | 26 | 26 | 0.00e+00 | 0e+00 | yes |
| lmest_B | max |transition diff| per occasion | 8.87936e-07 | 0 | 8.88e-07 | 1e-04 | yes |
| lmest_B | max |initial diff| | 2.53119e-07 | 0 | 2.53e-07 | 1e-04 | yes |
| lmest_C | log likelihood | -3229.12 | -3229.12 | 5.14e-11 | 1e-04 | yes |
| lmest_C | parameter count | 19 | 19 | 0.00e+00 | 0e+00 | yes |
| lmest_C | max |transition diff| per occasion | 6.04552e-08 | 0 | 6.05e-08 | 1e-04 | yes |
| lmest_C | max |initial diff| | 1.42769e-08 | 0 | 1.43e-08 | 1e-04 | yes |
| mplus_measurement | log likelihood | -4015.05 | -4015.05 | 3.33e-05 | 5e-03 | yes |
| mplus_measurement | parameter count | 15 | 15 | 0.00e+00 | 0e+00 | yes |
| mplus_measurement | max |threshold diff| | 0.000452575 | 0 | 4.53e-04 | 2e-03 | yes |
| mplus_measurement | max |transition diff| | 0.000359838 | 0 | 3.60e-04 | 2e-03 | yes |
