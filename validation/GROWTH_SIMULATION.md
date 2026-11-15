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
| small_variance | as balanced (see the correction below) | 16/18 | 0.129 | 0.92-1.02 | 0.936-0.964 | 0.936-0.962 |
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

**Correction (2026-11-15).** The `small_variance` cell was declared as the
balanced design with a random-slope SD of 0.1, but the script built it with
`utils::modifyList()`, which leaves an unnamed list (the per-class
variances) unchanged. The cell as declared in the registry and run drew the
base variances (slope variances 0.09 and 0.05, as its summary's truth
column shows): it is a second balanced cell with its own seeds, and a
random-slope variance near zero has not been studied. The script now builds
the cell as it ran, with a comment; the registry row is unchanged.

Not covered: a random-slope variance near zero (above), multilevel growth
mixtures, and membership covariates in the growth model itself.

## Correlated residuals (added 2026-11-15)

Two cells, `ar` and `brownian`, declared in the registry (commit "Growth
simulation: declare correlated-residual and bootstrap coverage cells")
before their first replicate: 300 persons, 8 occasions (0-7), the trajectory
coefficients and class shares of the balanced cell, a random intercept
(variances 1 and 0.6), a class-specific white noise (SD 0.5 and 0.6), and a
process shared by the classes, fitted with the matching
`residual_correlation`:

- `ar`: `sd^2 exp(-rate |t - s|)` with rate 0.7 and sd 0.6;
- `brownian`: `sd^2 min(t, s)` with sd 0.3.

The generator draws the process from the covariance the package fits (and
`lcmm::hlme(cor = AR(time))`/`BM(time)` define);
`tests/testthat/test-growth-simulation-coverage.R` checks the generator's
empirical covariance against `.growth_process_covariance()`. The likelihood
itself is checked against `hlme` in
`equivalence/test-growth-correlated-lcmm.R`.

Criteria (the verification article's Layer 4, fixed before the first
replicate): availability >= 98%, |mean estimate - truth| within 3.5 Monte
Carlo SE, observed-information 95% coverage inside the central 99.9%
binomial range for the replicates available. The gates of the four cells
above are also reported (`pass_gates`); they are not the criterion here,
because the reduced run below cannot resolve 0.925-0.975.

**Reduced run.** 150 datasets per cell of the 500 declared (45.8 minutes on
7 socket workers of a shared Windows machine, 2026-11-15, latents 1.3.4).
Full run: `GROWTH_SIM_STUDY=correlated Rscript validation/growth-simulation.R`.
Outputs: `growth-simulation-correlated-{replicates,summary}.csv`,
`growth-simulation-correlated.txt`.

| Cell | Parameters | Available | max abs(bias)/MCSE | SE ratio (observed) | Coverage observed | Coverage robust | Pass |
|---|---|---|---|---|---|---|---|
| ar | 14 | 142/150 (0.947) | 2.57 | 0.94-1.19 | 0.908-1.000 | 0.908-0.979 | 0/14 (availability) |
| brownian | 13 | 150/150 | 1.68 | 0.90-1.09 | 0.920-0.987 | 0.927-0.973 | 13/13 |

Binomial range for 95% coverage: 0.880-1.000 at 142 replicates, 0.880-0.993
at 150.

Brownian motion meets every criterion on every parameter, its process SD
included (bias -0.002, SE ratio 0.95, coverage 0.947).

In the autoregressive cell every parameter is unbiased within 3.5 Monte
Carlo SE and covered within its binomial range, the rate and the process
SD included (rate: bias 0.005, coverage 142/142, SE ratio 1.09: its Wald
interval, formed on the log scale, is conservative; process SD: coverage
0.965). The cell misses the availability criterion: 8 of 150 fits (5.3%)
put one class's white-noise SD on its floor while the process took a
higher rate (refits of replicates 3, 10 and 32: rate 1.1-1.6 against 0.7,
residual SD of class 2 at 0.0017), and `parameter_inference()` refuses a
fit with a coordinate on its bound. That is the documented behaviour (an
autoregressive process with a high rate is nearly white noise, so the data
separate it from the residual variance only weakly), and the replicates
record the reason. A pilot of 6 replicates with a random intercept and slope
over 6 occasions, before the design was declared, left the process on its
bound in 1 of 6 with its rate barely identified (standard errors as large as
the estimate), which is why the declared cells use a random intercept over
8 occasions.


## Bootstrap interval coverage (declared 2026-11-15, not yet run)

The registry declares a cell for the coverage of the person-bootstrap
intervals of `parameter_inference(method = "bootstrap")`. Its first run was
interrupted on 2026-11-15 before any replicate was retained, so no result is
reported and no claim rests on it. Run it with
`GROWTH_SIM_STUDY=bootstrap Rscript validation/growth-simulation.R`.
