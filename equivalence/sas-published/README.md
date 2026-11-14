# lca(), lpa(), enumerate_lpa() and mixture_regression() against SAS output

SAS output printed in tutorials and in SAS's own documentation, with the data
those pages use (`data/`, `DATA.md5`, `PROVENANCE.txt`), is reproduced:

| Case | Source | Compared |
|---|---|---|
| `proclca_drinking_k3` | PROC LCA, UCLA OARC data analysis example "Latent Class Analysis" (`lca1.dat`: 1000 persons, nine yes/no items, three classes) | log likelihood, G2, AIC, BIC, CAIC, adjusted BIC, entropy, degrees of freedom; every gamma and rho with its error |
| `mccutcheon_p16_k2`, `_p33_k1`, `_p33_k2`, `_p47_k1` | PROC LCA, UCLA OARC textbook examples of McCutcheon (1987), *Latent Class Analysis*, pages 16, 33 (Table 3.2) and 47 (Table 4.2) | log likelihood, G2, AIC, BIC, degrees of freedom; every gamma and rho |
| `proclca_help_k3` | PROC LCA, Kleinman and Horton, *SAS and R*, Example 8.21 (HELP study: homeless, CESD > 20, substance abuse treatment, linkage; 431 persons, three classes) | log likelihood, G2, AIC, BIC, CAIC, adjusted BIC, entropy, degrees of freedom; every gamma and rho |
| `fmm_galaxies_varying_k3`, `fmm_galaxies_equal_k4` | PROC FMM, SAS/STAT 14.2 User's Guide ch. 40, Getting Started: galaxy velocities (Roeder 1990), Figures 40.14 and 40.17 | -2 log likelihood, AIC, AICC, BIC, parameters; every mean and variance with its error; mixing probabilities, GLogits and their errors |
| `fmm_galaxies_equal_k3` ... `_k7` | the same, the equal-variance grid of Figure 40.17 (`enumerate_lpa()`) | -2 log likelihood and parameters of each model |
| `fmm_yeast_k2` | PROC FMM, the same chapter: "Student"'s yeast cell counts, a two-component binomial(5) mixture, Figures 40.2 and 40.3 | intercepts with errors, success probabilities, mixing probability, logit with error, the posteriors of every cell |

```
Rscript equivalence/sas-published/compare.R     # comparison.csv, RUN.txt
Rscript equivalence/sas-published/mad-stop.R    # mad-stop.csv (about 6 minutes)
```

`equivalence/test-sas-published.R` runs the comparison under testthat.
`mad-stop.csv` is the macOS run of `mad-stop.R`; on ci (Linux) the script
gives the same figures to two significant digits.

Tolerances, fixed before the comparison: half a unit of the last printed
digit plus 1e-6 of the value; parameter counts and degrees of freedom exactly.

## Conventions

- **PROC LCA's criteria** are on the deviance scale: G2 = 2 (L_saturated - L),
  AIC = G2 + 2p, BIC = G2 + p log n, CAIC = G2 + p (log n + 1), adjusted BIC
  = G2 + p log((n + 2) / 24), degrees of freedom = cells - p - 1. latents has
  no whole-table G2, so the kit computes it from latents' log likelihood and
  the observed response patterns; the degrees of freedom check latents'
  parameter count exactly.
- **PROC LCA's errors** invert the Hessian of the log likelihood (users'
  guide 1.3.2, section 4.3): latents' observed-information errors, unscaled.
- **PROC LCA's codes.** Items are coded 1, 2, ... as PROC LCA requires; the
  kit uses the pages' own recodes (UCLA: 2 - item; SAS and R: + 1, CESD > 20).
- **PROC FMM's mixing logits** are GLogit = log(p_k / p_K), the last printed
  component the reference; errors by the delta method on the observed
  covariance. AICC = -2L + 2pn / (n - p - 1).
- **Frequency tables** (McCutcheon's, the yeast cells) are expanded to one
  row per person, the same likelihood as SAS's FREQ statement.
- **Galaxies.** The manual's DATA step, v = velocity / 1000. Its 78th value is
  26960; `MASS::galaxies` has 26690, a typo its help page notes.

## Known differences (`known-differences.csv`)

- **PROC LCA's stopping rule.** PROC LCA stops EM when no parameter moved more
  than 1e-6 in an iteration. On a flat likelihood (drinking) or one whose
  maximum lies on a bound (McCutcheon page 16, HELP) that is short of the
  maximum by up to 1.3e-3. `mad-stop.R` runs a plain EM with that rule from
  50 random starts: in every case its nearest stopped point lies within 5e-5
  of every printed estimate, and for the drinking example the inverted
  Hessian there gives every printed error within 6.4e-5 (latents' errors at
  the maximum: within 2.4e-4). Every log likelihood, G2, criterion, entropy
  and degree of freedom agrees. As a check of the script, its Hessian errors
  at latents' maximum equal latents' own observed-information errors to 1e-7.
- **PROC FMM's equal-variance grid.** At three, five, six and seven
  components PROC FMM stopped at lower maxima (478.74; its four-component
  416.49 for five to seven), which latents' random starts also reach; latents'
  best are 425.36, 410.68, 394.58 and 388.86. Equal variances keep the
  likelihood bounded. The selected four-component model agrees in full.
- **Yeast, third cell's posterior.** 1.4e-5 off at the fifth decimal: PROC
  FMM's value needs a mixing probability 6e-6 from the maximum (both print as
  0.8799).

## Not compared

PROC LCA's raw entropy; McCutcheon's models with restrictions (pages 33 to
70: fixed or equated probabilities, which latents does not fit); the PROC LCA
users' guide examples (both printed outputs use the stabilizing prior
`RHO PRIOR=1`, a posterior mode); PROC FMM's Pearson statistic and maximum
gradients; the unequal-variance grid above three components (that
likelihood is unbounded: latents puts a component on one or two velocities,
PROC FMM printed its three-component value); and the chapter's examples with
distributions latents does not fit (binomial cluster and beta-binomial,
Weibull, multinomial cluster, a zero-inflated Poisson built from a constant
component).
