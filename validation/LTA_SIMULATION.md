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
