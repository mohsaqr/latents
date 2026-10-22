# Direct effects: equivalence and property tests

`lpa(..., direct_effects = "y3 ~ z")`: a covariate moves an indicator's mean by
the same amount in every profile, `y = mean_k + b z + e` (Latent GOLD
`y3 <- 1 | Cluster + z`; Mplus `y3 ON z;` in `%OVERALL%`). Tested against
genuine **Latent GOLD 6.1** (Wine) and the **Mplus 9 Demo (Mac)**, run on
2026-10-20, and against an independent implementation. Findings:
[FINDINGS.md](FINDINGS.md).

## Results (2026-10-20)

- `Rscript equivalence/run.R features-direct-effects`: 6,125 expectations
  pass, 2 fail (both heavy_02, finding 1).
- `tests/testthat/test-features-direct-effects.R` (NOT_CRAN=true): 79 pass,
  3 fail (findings 2 and 3).

## Files

| File | What it is |
|---|---|
| `make-kit.R` | Builds `kit/` deterministically: 67 cases from simulab, each stored (data, spec, truth, seed) in `kit/cases.rds`, plus `<case>.lgs`, `<case>-lg.dat`, `<case>.inp`, `<case>.dat` |
| `returned/` | Retained outputs: LG `.lst`, `_write.txt` (full-precision parameters and SEs), `_posteriors.txt`; Mplus `.out`, `-results.dat`, `-posteriors.dat` |
| `helpers.R` | Readers (LG write file, Mplus TECH1-mapped results file), latents fits (3 seeds x 20 starts, tol 1e-10), alignment, `de_compare()` with every tolerance and its justification |
| `reference.R` | Independent implementation (no latents code): `dnorm` mixture likelihood, EM with a joint WLS M-step for intercepts and effects, numDeriv Hessian SEs, and the 84-dataset battery |
| `compare.R` | Writes `comparison.csv` (5,641 rows), `maxima.csv`, `reference.csv` |
| `sabotage.R` | In-memory defects against the unit tests; `sabotage.csv` |
| `repro/` | Reproductions of the findings |
| `../../test-features-direct-effects.R` | The equivalence assertions |
| `../../../tests/testthat/test-features-direct-effects.R` | Unit and property tests |

## Design (external kit, 67 cases)

Four indicators; z ~ N(0.5, 1.5^2) unless stated; b = 0.8 on y3 (two
covariates: 0.8 and -0.5 on y3; two indicators: 0.8 on y3 and -0.6 on y1).

| Configuration | Model | n | Effects | Cases |
|---|---|---|---|---|
| vvi2, eei2 | VVI / EEI, 2 profiles | 500 | y3 ~ z | 5 + 5 |
| vvi3, eei3 | VVI / EEI, 3 profiles | 500 | y3 ~ z | 5 + 5 |
| vvi2cov2 | VVI, 2 | 500 | y3 ~ z1 + z2 | 5 |
| vvi2ind2 | VVI, 2 | 500 | y3 ~ z, y1 ~ z | 5 |
| vvi2small, eei2small | VVI / EEI, 2 | 40, 50, 60, 70, 80 | y3 ~ z | 5 + 5 |
| vvi2large | VVI, 2 | 5000 | y3 ~ z | 5 |
| eei3large | EEI, 3 | 5000 | y3 ~ z1 + z2 | 5 |
| vvi2huge | VVI, 2 | 20000 | y3 ~ z | 1 |
| farmean | VVI, 2 | 500 | z mean 1000 | 2 |
| unitsbig / unitssmall | VVI, 2 | 500 | z in units 1e4 / 1e-4 | 1 + 1 |
| zero | VVI, 2 | 500 | true b = 0 | 2 |
| collinear | VVI, 2 | 500 | z = 2 [class 2] + N(0, 0.3^2) | 2 |
| binary | VVI, 2 | 500 | z ~ Bernoulli(0.4) | 2 |
| heavy | VVI, 2 | 500 | t(3) noise | 2 |
| k1, k1cov2 | 1 profile (= linear regression) | 300 | two indicators / two covariates | 1 + 1 |
| fewrows | EEI, 2 | 400 | a 3% profile (about 12 rows) | 2 |

LG: 100 sets x 250 iterations, NR tolerance 1e-12. Mplus: STARTS 400 80,
observed information. Mplus needs `y1 WITH y3@0` (it frees that residual
covariance by default) and `DEFINE` rescaling for the unit cases (it refuses
variances above 1e6 or near zero); the reader converts back.

## Comparisons and largest differences (from `comparison.csv`)

Quantities are compared where a program reached latents' maximum: 66 cases
against both programs except farmean_01 (Latent GOLD stopped at a degenerate
solution; Mplus only) and heavy_02 (latents lower; finding 1). Classes are
aligned once from the posteriors (`montecarlo_relabel`). d is the per-case
largest posterior difference.

| Quantity | Tolerance | Latent GOLD max diff | Mplus max diff |
|---|---|---|---|
| log likelihood | 1e-6 relative | 8.1e-8 | 4.3e-3 (3.5e-8 relative, n = 20000) |
| parameters | exact | 0 | 0 |
| BIC | 2e-6 relative (+5e-5 printed) | 5.0e-5 | 8.6e-3 |
| posteriors, predict() (row max) | 5e-4 | 1.04e-4 | 1.04e-4 |
| means | (1e-4 + 10 d) x max(1, \|value\|) | 5.8e-5 | 6.0e-5 |
| variances | same | 5.0e-5 | 5.2e-5 |
| effects | same | 2.9e-4 (1.5e-6 relative, unitssmall) | 7.5e-5 |
| class proportions | same | 1.7e-5 | 1.7e-5 |
| SE of means | (1e-3 + 10 d) relative | 5.7e-6 (5.4e-5 rel.) | 5.6e-6 |
| SE of variances | same | 9.7e-6 (5.8e-5 rel.) | 9.2e-6 |
| SE of effects | same | 9.3e-5 (1.7e-6 rel.) | 4.7e-5 (1.4e-6 rel.) |
| SE of proportions | same | 2.6e-6 | not printed by Mplus |

Also asserted: latents reaches the best maximum of either program (fails on
heavy_02 only), every seed reaches the same maximum, k = 1 equals `lm()`
(coefficients and ML variance to 1e-9, log likelihood to 1e-12), and the
estimated effects lie within 4.5 SE of the generating value (collinear
excluded: the model is misspecified there by construction, and all three
programs agree on b = 1.36 and 1.48 against a truth of 0.8, the bias of
ignoring that the covariate predicts class membership).

## Independent implementation (84 datasets, `reference.csv`)

n in {60, 200, 600, 2000}, k in {1, 2, 3}, VVI/EEI, the three effect sets,
normal and binary covariates.

| Check | Tolerance | Worst |
|---|---|---|
| (a) independent likelihood at latents' estimates | 1e-12 relative | 3.6e-15 relative |
| (b) latents vs independent maximum | latents not below by 1e-9 relative | latents below by 5.4e-8 (2.7e-11 relative); latents higher in 5 datasets where the independent EM's 12 starts missed the basin |
| (c) SEs vs numDeriv Hessian | 1e-4 relative | 2.3e-5 |
| (d) score per SE at latents' estimate | 1e-3 | 7.2e-5 |

Two n = 60, k = 3 VVI fits end at an active variance bound; (c) and (d) are
skipped there (inference is not defined at the bound).

## Unit and property tests

`tests/testthat/test-features-direct-effects.R` (18 s with NOT_CRAN=true):
row permutation, column/term order and unused columns, group relabelling,
covariate shift (means move by -b c; b, variances, likelihood, SEs fixed),
covariate scaling (b by 1/s, SE by 1/|s|), indicator scaling (mean, effect, SD
by a; likelihood by -n log a), integer weights = repeated rows, predict() on
the training data, k = 1 = `lm()`, bootstrap and sensitivity refits carry the
effects (traced), every documented error class. Data are drawn in base R:
simulab is not a declared dependency of the package's tests.

Sabotage (`sabotage.R`, `sabotage.csv`): the control changes nothing; a
half-size shift in the likelihood, an M-step without the class variances,
means left at the covariate centre, and a score without its covariate are
each detected. The last was first missed (NA standard errors compared equal
to NA); the shift test now requires finite positive SEs.

## Rerunning

From the project root:

```
Rscript equivalence/features/direct-effects/make-kit.R
JOBS=2 sh equivalence/features/run.sh direct-effects all     # about 30 min
Rscript equivalence/features/direct-effects/compare.R         # CSVs, about 8 min
Rscript equivalence/run.R features-direct-effects             # assertions, about 8 min
NOT_CRAN=true Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-features-direct-effects.R")'
NOT_CRAN=true Rscript equivalence/features/direct-effects/sabotage.R
```
