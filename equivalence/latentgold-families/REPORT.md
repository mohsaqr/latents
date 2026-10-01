# Latent GOLD comparison: group-class families — results

Run 2026-10-01; R version 4.5.2 (2025-10-31); latents 0.8.8. Tolerances declared in README.md before running.

**16 of 24 quantities agree.**

| Case | Family | Between | Quantity | latents | Latent GOLD | Difference | Tolerance | Agree |
|---|---|---|---|---|---|---|---|---|
| g01 | additive | varying | parameter count |     6.000000 |     6.0000 | 0.00e+00 | 0e+00 | yes |
| g01 | additive | varying | log likelihood | -1524.521874 | -1524.5242 | 2.33e-03 | 1e-04 | **no** |
| g01 | additive | varying | group posteriors |  |  | 1.85e-02 | 1e-03 | **no** |
| g02 | additive | equal | parameter count |     5.000000 |     5.0000 | 0.00e+00 | 0e+00 | yes |
| g02 | additive | equal | log likelihood | -1510.249882 | -1510.2499 | 1.80e-05 | 1e-04 | yes |
| g02 | additive | equal | group posteriors |  |  | 2.79e-06 | 1e-03 | yes |
| g03 | dispersion | equal | parameter count |     5.000000 |     5.0000 | 0.00e+00 | 0e+00 | yes |
| g03 | dispersion | equal | log likelihood | -1618.323096 | -1618.3232 | 1.04e-04 | 1e-04 | **no** |
| g03 | dispersion | equal | group posteriors |  |  | 1.17e-05 | 1e-03 | yes |
| g04 | additive_dispersion | varying | parameter count |     7.000000 |     7.0000 | 0.00e+00 | 0e+00 | yes |
| g04 | additive_dispersion | varying | log likelihood | -1587.897773 | -1587.8978 | 2.67e-05 | 1e-04 | yes |
| g04 | additive_dispersion | varying | group posteriors |  |  | 4.97e-07 | 1e-03 | yes |
| g05 | additive_dispersion | equal | parameter count |     6.000000 |     6.0000 | 0.00e+00 | 0e+00 | yes |
| g05 | additive_dispersion | equal | log likelihood | -1580.951227 | -1580.9512 | 2.70e-05 | 1e-04 | yes |
| g05 | additive_dispersion | equal | group posteriors |  |  | 3.78e-06 | 1e-03 | yes |
| g06 | additive | varying | parameter count |     9.000000 |     9.0000 | 0.00e+00 | 0e+00 | yes |
| g06 | additive | varying | log likelihood | -1673.755468 | -1673.7556 | 1.32e-04 | 1e-04 | **no** |
| g06 | additive | varying | group posteriors |  |  | 2.24e-04 | 1e-03 | yes |
| g07 | additive | varying | parameter count |    11.000000 |    11.0000 | 0.00e+00 | 0e+00 | yes |
| g07 | additive | varying | log likelihood | -3118.092970 | -3118.1614 | 6.84e-02 | 1e-04 | **no** |
| g07 | additive | varying | group posteriors |  |  | 1.91e-02 | 1e-03 | **no** |
| g08 | dispersion | equal | parameter count |     9.000000 |     9.0000 | 0.00e+00 | 0e+00 | yes |
| g08 | dispersion | equal | log likelihood | -3080.964240 | -3080.7368 | 2.27e-01 | 1e-04 | **no** |
| g08 | dispersion | equal | group posteriors |  |  | 2.84e-02 | 1e-03 | **no** |
