# Growth mixture models: simulation of estimates, standard errors and R3STEP

Predeclared study of `mixture_regression(random = ~ 1 + time)` (phases G2 and
G4 of `todo.MD`). Script [growth-simulation.R](growth-simulation.R); registry
[growth-simulation-registry.csv](growth-simulation-registry.csv) written before
the first replicate (the script refuses to run if it changes). 1000 datasets
per cell, 5 starts each, 2 classes, a random intercept and slope with
class-specific covariances and residual variances. Class membership depends
on a person covariate (logit 0.2 + 0.8 x); the growth mixture is fitted
without it, so its class shares estimate E[P(class | x)], and
`r3step(fit, data, "x")` estimates the slope. Fitted classes are named by the
true class whose trajectory they match. Run 2026-10-17 on latents 0.9.37
(1.2 hours, 8 cores). Outputs: `growth-simulation.txt`,
`growth-simulation-summary.csv`, `growth-simulation-replicates.csv`,
`growth-simulation-r3step.csv`.

Gates as in ADDITIVE_SIMULATION.md: availability >= 98%, |bias| <= 0.1
empirical SD, observed SE ratio (mean SE / empirical SD) 0.9-1.1, observed
95% coverage 0.925-0.975. 18 parameters per cell: per class two coefficients,
the residual SD, two random-effect variances and SDs, their correlation, and
the class share.

| Cell | Design | Pass | max abs(bias)/SD | SE ratio (observed) | Coverage observed | Coverage robust |
|---|---|---|---|---|---|---|
| balanced | 300 persons, 6 occasions, separated trajectories | 17/18 | 0.101 | 0.97-1.03 | 0.940-0.963 | 0.939-0.965 |
| unbalanced | 300 persons, times jittered by +-0.3, 25% of occasions missing | 15/18 | 0.145 | 0.96-1.05 | 0.939-0.965 | 0.939-0.965 |
| small_variance | as balanced, random-slope SD 0.1 | 16/18 | 0.129 | 0.92-1.02 | 0.936-0.964 | 0.936-0.962 |
| weak | 500 persons, intercepts 2 vs 3.2, slopes 0.6 vs 0.1 | 3/18 | 0.139 | 0.70-0.94 | 0.885-0.943 | 0.903-0.946 |

Availability was at least 99.5% in every cell; the 7 unavailable fits (of
4000) had a degenerate random-effect covariance, which the fit reports and
inference refuses. Gates were not changed after the results.

With separated classes, 48 of 54 parameters meet every gate. The six misses
fail on bias alone (0.10-0.15 empirical SD) and all on random-effect SDs and
correlations: the small-sample bias of maximum-likelihood variance
components. Their coverage is 0.949-0.965.

Under weak separation the standard errors are too small: the SE ratio is
0.70-0.94 and coverage 0.885-0.943, worst for the trajectory coefficients and
class shares. The sandwich (robust) errors are only slightly better
(coverage 0.903-0.946). Weakly separated classes trade persons from one
replicate to the next, a source of variation that the information at one
maximum does not see; the person bootstrap (`parameter_inference(method =
"bootstrap")`) is the better choice there. On a correctly specified,
separated case, bootstrap and Wald errors agreed (ratios 0.85-1.29 over 60
replicates; `tests/testthat/test-growth-parameter-inference.R` checks the
machinery).

## R3STEP

`r3step()` slope on x (true 0.8 in absolute value), sandwich errors over persons:

| Cell | Available | bias/SD | SE ratio | Coverage |
|---|---|---|---|---|
| balanced | 1.000 | 0.088 | 0.963 | 0.944 |
| unbalanced | 0.998 | 0.104 | 0.947 | 0.946 |
| small_variance | 1.000 | 0.093 | 0.980 | 0.954 |
| weak | 0.984 | 0.014 | 0.206 | 0.920 |

The weak cell's SE ratio comes from 2 of 984 replicates with
near-separated estimates (largest 29.6). Without them the ratio is 0.87
(empirical SD 0.211, mean SE 0.183) and coverage 0.920: the three-step
correction treats the classification errors as known, and with weakly
separated classes they are not.

Not covered: correlated residuals (`residual_correlation`), multilevel
growth mixtures, membership covariates in the growth model itself, and the
bootstrap intervals' coverage.
