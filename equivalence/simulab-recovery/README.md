# Recovery of known parameters from simulab data

Every public function of latents is fitted to data that simulab drew from
known parameters, many times over, and its estimates are compared with the
parameters that generated the data. This is the fourth reference every
function needs, beside an R package, Mplus and Latent GOLD: those check that
latents computes what other programs compute; recovery checks that what it
computes is the truth.

| Case | Function | Design (simulab) | Replicates | Compared |
|---|---|---|---|---|
| `lca_k3` | `lca()` | `simulate_lca()`: 3 classes (.5, .3, .2), 6 binary items, n = 1000 | 100 | class shares and P(1) of every item and class |
| `lca_r3step` | `r3step()` | the same, with a covariate normal within each class (means 0, 1, -1; variance 1) | 100 | the multinomial-logit intercepts and slopes, known in closed form |
| `lca_three_step_bch`, `_ml` | `three_step()` | the same, with a distal outcome normal within each class (means 0, 1.5, 3) | 100 | each class's mean, BCH and ML |
| `lca_enumeration` | `enumerate_lca()` | the same | 30 | BIC selects 3 classes |
| `lpa_eei_k3` | `lpa()` | `simulate_lpa()`: 3 profiles (.4, .35, .25), 4 indicators, unit variances, n = 800 | 100 | shares, means, variances |
| `lpa_enumeration` | `enumerate_lpa()` | the same | 30 | BIC selects 3 profiles |
| `lta_hmm_k2` | `lta()` | `simulate_hmm()`: 2 states, 3-category observation, 5 occasions, 1000 chains | 100 | initial, transition and response probabilities |
| `mixture_regression_k2` | `mixture_regression()` | `simulate_lpa()` draws class (.6, .4) and covariate; y = class line + N(0, 0.5²) | 100 | intercepts, slopes, residual sd, membership logit |
| `growth_k2` | `mixture_regression(random =)` | `simulate_growth()` per class (binomial sizes, .6): random intercept and slope sd .6 and .2, residual sd .5 | 50 | fixed effects, residual and random-effect sds, membership logit |
| `multilpa_k3_g2` | `multilpa()` | `simulate_ml_lpa()`: 60 clusters of 20, 3 profiles, 2 cluster classes | 40 | cluster-class shares, profile shares within each, means, variances |
| `multilca_k3_g2` | `multilca()` | `simulate_lca()` per cluster, its shares those of the cluster's class: 60 clusters of 25 | 40 | cluster-class shares, class shares within each, P(1) of every item |
| `bootstrap_lrt_size`, `_power` | `bootstrap_lrt()` | `simulate_lca()` with 2 and with 3 classes, n = 500; 2 against 3 classes, 19 draws | 20, 10 | rejection rate at 0.05 |
| `pool_imputations_lpa_k2` | `pool_imputations()` | `simulate_lpa()` (2 profiles), `inject_missingness()` 15% MCAR, mice `norm` with m = 5 | 30 | pooled shares and means |

```
Rscript equivalence/simulab-recovery/compare.R   # comparison.csv, RUN.txt; about an hour
```

`equivalence/test-simulab-recovery.R` runs the same comparison under testthat;
ci runs it as its own job, `mini-ci equivalence-r`.

## Criteria, fixed before the first run

- **Bias.** The mean of the estimates (or of estimate minus truth) is within
  3.5 Monte Carlo standard errors (3.5 sd / sqrt(R)) of the true value: a
  two-sided test at about 0.0005 per parameter.
- **Coverage.** The share of the 95% intervals that contain the true value
  lies in the central 99.9% binomial range for R draws at 0.95. Intervals are
  latents' own (`parameter_inference()`, `three_step()`, `r3step()`,
  `pool_imputations()`); for `mixture_regression()` they are Wald intervals
  from `vcov()` by the delta method on the reported scale.
- **Selection and tests.** BIC picks the generating number of classes in at
  least 90% of replicates (well separated designs); the bootstrap LRT rejects
  in the central 99.9% binomial range at 0.05 under two classes and in at
  least 90% under three.
- **Every replicate fits.** A replicate whose fit or inference fails is
  counted, and its design's "replicates fitted" row then disagrees.

Fitted classes are matched to the generating ones by the order that
minimises the squared distance of their class-specific parameters.

## Where simulab draws only part of the data

simulab draws the latent classes and their indicators; the kit adds from the
replicate's seed, in base R, what simulab has no generator for: the
covariate and distal outcome of the r3step() and three_step() designs, the
outcome of the regression mixture, the cluster classes of the two-level LCA,
and the binomial class sizes of the growth mixture. The pool_imputations()
design imputes with the generating profile as a predictor, so the
imputations are proper and the check isolates Rubin's pooling.

## Known limitation found

`parameter_inference(boundary = "fix")` is documented for categorical
response probabilities but transition fits (`lta()`) do not apply it: a fit
with a response probability on its bound is refused
(`R/transition-inference.R`). The lta design keeps every probability at 0.15
or more so that such a maximum is rare; one that occurs is counted as a
replicate not fitted.
