# Latent transition extensions: simulation of standard errors

Predeclared study of the general transition model (`lta()` with covariate,
occasion-varying or second-order transitions, categorical measurement).
Script [lta-simulation.R](lta-simulation.R); registry
[lta-simulation-registry.csv](lta-simulation-registry.csv) written before the
first replicate (the script refuses to run if it changes). 1000 datasets per
cell, 5 starts each, 2 profiles. Output: `lta-simulation.txt`,
`lta-simulation-summary.csv`, `lta-simulation-replicates.csv`; the final code
reproduces stored replicates to 4.7e-15. Gates as in ADDITIVE_SIMULATION.md:
availability >= 98%, |bias| <= 0.1 empirical SD, observed SE ratio 0.9-1.1,
observed 95% coverage 0.925-0.975.

| Cell | Design | Pass | max abs(bias)/SD | SE ratio (observed) | Coverage observed | Coverage robust |
|---|---|---|---|---|---|---|
| covariate | 300 groups, 5 occasions, initial logit on w, transitions on time-varying z | 6/6 | 0.088 | 0.93-1.03 | 0.935-0.959 | 0.936-0.963 |
| occasion | 300 groups, 4 occasions, a transition matrix per move | 6/7 | 0.123 | 0.95-1.04 | 0.952-0.958 | 0.953-0.962 |
| second_order | 300 groups, 5 occasions, second-order logits | 5/7 | 0.118 | 0.94-1.01 | 0.942-0.967 | 0.947-0.968 |
| categorical | 400 groups, 3 occasions, 5 binary items, transitions on z | 6/6 | 0.073 | 0.97-1.02 | 0.933-0.958 | 0.934-0.956 |

No fit failed, every fit converged and none was on a boundary. 23 of 26
coefficients meet every gate. The three misses fail on bias alone:
occasion-2 move 2->1 (-0.123 SD), the first-order move 1->2 under second
order (-0.118 SD) and the second-order logit of origin pair (2,2) (-0.113
SD), all estimated slightly too negative, the small-sample bias of a
maximum-likelihood logit with few informative moves; their coverage is
0.948-0.956. Gates were not changed after the results.

Not covered: the mover-stayer class (validated against the exact path sum
only), missing data and covariance structures in the extended model.

## Constrained covariance structures in the homogeneous `lta()` (2026-11-15)

Predeclared study of the Wald intervals that the homogeneous `lta()`
(`multilpa_transitions`) reports for the constrained covariance structures.
Script [lta-structure-coverage.R](lta-structure-coverage.R); registry
[lta-structure-coverage-registry.csv](lta-structure-coverage-registry.csv)
written before the first replicate (the script refuses to run if it changes).
Four cells, one per chart form: EII (spherical), VEI (diagonal), EVE
(rotation), VEE (determinant-one Cholesky); 250 sequences of five occasions,
two profiles, two Gaussian indicators, initial (0.6, 0.4), transitions
(0.85, 0.15 / 0.25, 0.75); three starts per fit. Gates fixed with the
registry: availability >= 98%, |bias| <= 0.1 SD, SE ratio 0.9-1.1, observed
95% coverage inside the central 99.9% binomial range (0.916-0.978 at 500).

**Reduced run: 500 datasets per cell** (2000 fits, about 16 minutes on 8
workers), not the registry's default 1000. Each replicate draws from its own
seed, so the full run reproduces these 500 and adds 500 more:
`LTA_STRUCT_N=1000 LTA_STRUCT_CORES=8 Rscript validation/lta-structure-coverage.R`.
Outputs: `lta-structure-coverage-{replicates,summary}.csv`,
`lta-structure-coverage.txt`.

| Cell | Parameters | Pass | max abs(bias)/SD | SE ratio | Coverage |
|---|---|---|---|---|---|
| EII | initial, two stays, two means, shared variance | 5/6 | 0.162 | 0.957-1.021 | 0.932-0.970 |
| VEI | initial, two stays, two means, two variances | 7/7 | 0.099 | 0.945-1.037 | 0.930-0.964 |
| EVE | initial, two stays, two means, three covariances | 8/8 | 0.052 | 0.918-1.076 | 0.930-0.966 |
| VEE | initial, two stays, two means, three covariances | 7/8 | 0.120 | 0.949-1.090 | 0.926-0.972 |

Every fit converged, none was on a boundary, no warning was raised. All 29
coverages lie in the binomial range and every SE ratio in 0.9-1.1. 27 of 29
parameters meet every gate; the two misses fail on bias alone and are
variances estimated slightly too small, the downward small-sample bias of a
maximum-likelihood variance: EII's shared variance (-0.162 SD) and VEE's
profile-2 variance of `b` (-0.120 SD), with coverage 0.932 and 0.940. Gates
were not changed after the results.
