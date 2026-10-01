# Latent GOLD comparison: group-class families

External reference for `multilpa(family = "additive" | "dispersion" |
"additive_dispersion")`. Latent GOLD 6.1 (Syntax module, free academic
licence) fits each family as a two-level model with a group-level nominal
latent `GClass` and one group-level continuous factor per indicator, loading
fixed at 1 (the random intercept). Written 1 October 2026, before any case
below was run; the tolerances are fixed here, not after seeing results.

## Latent GOLD specification per family

| Family | Group factor variance | Indicator mean | Within variance |
|---|---|---|---|
| additive, varying | `F1 \| GClass;` | `y1 <- 1 + GClass + (1) F1;` | `y1;` |
| additive, equal | `F1;` | `y1 <- 1 + GClass + (1) F1;` | `y1;` |
| dispersion | `F1;` | `y1 <- 1 + (1) F1;` | `y1 \| GClass;` |
| additive-dispersion, varying | `F1 \| GClass;` | `y1 <- 1 + GClass + (1) F1;` | `y1 \| GClass;` |
| additive-dispersion, equal | `F1;` | `y1 <- 1 + GClass + (1) F1;` | `y1 \| GClass;` |

Every model: pure ML (`bayes ... = 0`), `tolerance=1e-10`, 20 random start
sets, Gauss-Hermite quadrature (50 nodes with one indicator, 15 per factor with
two). latents: `tol = 1e-14`, 20 starts.

## Cases

Simulated by `make-kit.R` (seed per case) from each family's own model, 150
groups of sizes 4, 6, 8, 10 cycled. `g07`/`g08` have two indicators.

| Case | Family | Classes | Indicators |
|---|---|---|---|
| g01 | additive, varying | 2 | 1 |
| g02 | additive, equal | 2 | 1 |
| g03 | dispersion | 2 | 1 |
| g04 | additive-dispersion, varying | 2 | 1 |
| g05 | additive-dispersion, equal | 2 | 1 |
| g06 | additive, varying | 3 | 1 |
| g07 | additive, varying | 2 | 2 |
| g08 | dispersion | 2 | 2 |

## Agreement, declared before running

- **Parameter count:** exact.
- **Log likelihood:** latents' maximum within 1e-4 of Latent GOLD's printed
  value (printed rounding 5e-5 plus quadrature error).
- **Group posteriors:** after aligning Latent GOLD's class labels to latents'
  by the permutation that brings them closest, the largest absolute
  difference over groups and classes at most 1e-3 (6-decimal output plus
  quadrature error; the pilot on g01-like data gave 6.2e-5).

## Running

1. `Rscript equivalence/latentgold-families/make-kit.R` writes `kit/` and
   copies it to `~/.wine-latentgold/drive_c/lgfamilies/`.
2. `sh equivalence/latentgold-families/run-kit.sh` runs Latent GOLD under Wine
   on each case and copies listings and posterior files to `returned/`.
3. `Rscript equivalence/latentgold-families/compare.R` writes
   `comparison.csv` and `REPORT.md` and exits 1 if any quantity disagrees.

## Results (1 October 2026)

**Against the declared criteria, base run:** 17 of 24 quantities agree
(`comparison.csv`, `REPORT.md`). Parameter counts agree in all 8 cases.
Cases g02, g04 and g05 agree on every quantity. The log-likelihood criterion
fails in g01, g03, g06, g07 and g08, and the posterior criterion in g01, g07
and g08. The 1e-4 likelihood tolerance did not allow for Latent GOLD's own
quadrature error, which the README underestimated.

**Diagnosis** (`diagnose.R`, `diagnose.csv`): the exact latents likelihood
evaluated at Latent GOLD's printed estimates. For g02-g06 it equals latents'
maximum to 1.4e-6 or better, so both programs found the same optimum and the
reported-likelihood gaps (1.7e-5 to 1.3e-4) are Latent GOLD's quadrature
error. In g01, g07 and g08, Latent GOLD maximized its quadrature
approximation and stopped below the exact maximum (by 5.7e-4, 6.6e-3 and
1.1e-2).

**Node refinement** (`refine.R`, `refine.csv`, `returned/refine*`): rerunning
Latent GOLD with finer quadrature moves it to latents in every case.

| Case | Nodes | Log-lik gap | Posterior gap |
|---|---|---|---|
| g01 | 50 -> 100 | 2.3e-3 -> 2.6e-5 | 1.9e-2 -> 1.4e-6 |
| g03 | 50 -> 100 | 1.0e-4 -> 4.4e-6 | 1.2e-5 -> 4.9e-7 |
| g06 | 50 -> 100 | 1.3e-4 -> 3.2e-5 | 2.2e-4 -> 6.7e-6 |
| g07 | 15 -> 30 per factor | 6.8e-2 -> 1.3e-4 | 1.9e-2 -> 2.0e-4 |
| g08 | 15 -> 30 -> 60 per factor | 2.3e-1 -> 5.4e-2 -> 1.6e-4 | 2.8e-2 -> 2.6e-3 -> 8.4e-6 |

Conclusion: every disagreement traces to Latent GOLD's numerical integration
and vanishes as its nodes are refined; latents' exact likelihood, which needs
no integration, is the more precise of the two. Parameter counts and the
family parameterizations agree exactly. After refinement the remaining
log-likelihood gaps (at most 1.6e-4) still exceed the 1e-4 declared at the
outset for g07 and g08; posteriors agree within 2.0e-4 everywhere.
