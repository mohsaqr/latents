# Latent transition extensions: external references

`lta()` with occasion-varying or covariate-dependent transitions, a covariate
initial distribution, occasion-specific measurement or second-order
transitions (the general engine, class `multilpa_lta`) compared with external
programs. The references and data come from the JStats (Carm) repository,
which must sit beside this one (`../JStats`). Run
`Rscript equivalence/lta-extensions/run-all.R`; it writes `comparison.csv` and
`REPORT.md` (30 of 30 quantities agree, 1 October 2026).

| Script | Reference | Extension | Compared |
|---|---|---|---|
| compare-ex814.R | Mplus 8.8, User's Guide ex. 8.14 (2000 persons, 2 occasions, 5 binary items, 3 classes) | covariate on initial classes and origin-specific covariate effects on transitions | log-likelihood, parameter count, thresholds, initial logits on x, transition odds ratios on x (against staying) |
| compare-ex813.R | Mplus 8.8, User's Guide ex. 8.13 | known group on initial classes and group-specific transitions (as a covariate) | log-likelihood (plus the known-class term), parameter count, thresholds, `C1 ON CG#1` |
| compare-mplus-measurement.R | Mplus 9 demo, occasion-specific thresholds (1200 persons) | measurement not invariant across occasions | log-likelihood, parameter count, thresholds at both occasions, transition matrix |
| compare-depmix-covariate.R | depmixS4 1.5-1 `transition = ~ x` (cross-checked with LMest in JStats) | covariate-dependent transitions | log-likelihood, parameter count, conditional transitions at x = -1.5, 0, 1.5, slopes, slope SEs |
| compare-lmest-occasion.R | LMest 3.2.8 `modBasic = 0`, three scenarios | occasion-varying transitions | log-likelihood, parameter count, transition matrix per occasion, initial probabilities |

Tolerances: Mplus prints three decimals, so likelihoods are compared within
5e-3 and estimates within 2e-3; depmixS4 and LMest are compared within 1e-4
(likelihood) to 1e-3 (estimates), and depmixS4's numerical-Hessian standard
errors within 5% relative.

## Not covered externally

- **Second-order transitions:** no reference program or fixture fits them.
  They are checked against an exact sum over every profile path
  (`tests/testthat/test-lta-general.R`) and by nesting the first-order model.
- **Covariance structures in `lta()`** (`model =`): no external LTA reference;
  checked by nesting and by `model = "VVI"` reproducing the default fit.
- **Occasion-specific measurement with continuous indicators:** the external
  reference has binary items; the continuous case is checked against the exact
  path sum.

The bundled test fixtures (`tests/testthat/fixtures/lta-references.rds`) carry
the depmixS4, LMest (scenario A) and Mplus-measurement data and reference
values, so the package tests do not need JStats.
