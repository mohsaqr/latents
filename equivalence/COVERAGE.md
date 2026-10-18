# External-equivalence coverage

Audit of 2026-10-17, package version 0.95.0; rows updated on 2026-10-18
(0.95.1) for the comparisons Equivalence 2.0 added. Each row names an exported
function or a user-facing option and the external reference that checks it.
The rows were filled in by reading the test files and comparison scripts.
File names, report names or READMEs were not taken as evidence of what is
compared.

## How to read it

- **Reference**: the program or source the package is compared with: Mplus,
  Latent GOLD (LG), mclust, lcmm, flexmix, depmixS4, LMest, tidySEM,
  tidyLPA, poLCA, glca, multilevLCA, sandwich, MASS/stats (one-class limits),
  published tables, JStats-pinned fixtures, or the in-house
  exhaustive-enumeration likelihood (`tests/equivalence/helper-independent-likelihood.R`).
- **Where**: `E/` = `equivalence/`, `T/` = `tests/equivalence/`. Test names
  are quoted. "script" means a comparison script that writes a report
  (`comparison.csv`/`REPORT.md`) rather than a testthat assertion.
  "harness" means `E/harness/suites/*.R`.
- **Verdict**:
  - **covered**: an external reference compares this quantity in at least one
    configuration.
  - **partial**: an external check exists but misses a named part, or the only
    independent check is an in-house formula (dense or path-sum likelihood,
    `numDeriv`) in `T/`.
  - **none**: only package-internal tests in `tests/testthat/` (internal
    only), or no test at all.
  - **n/a**: an optimizer control, plotting or reporting function with no
    estimand to compare. Not counted.
- Tolerances are the ones in the code. Absolute unless marked "rel".
- Status of the evidence: the 2026-10-17 audit read the tests without running
  them. The 2026-10-18 update ran every file it cites as new evidence
  (`E/test-montecarlo*.R`, `E/test-enumeration-mplus.R`, `E/test-mplus-*.R`,
  `E/test-latentgold*.R`, `E/test-external-scripts.R`, `E/test-claims.R`).
  The harness was rerun on 0.95.1: 1,018 of 1,018 quantities agree.
  The documentation's claims of agreement are listed in `E/CLAIMS.csv`, each
  with its test; `E/test-claims.R` fails when a claim or its test changes.

## 1. Fitting functions and their options

### Exported functions

| Function | Reference | Where | Compared | Tol. | Verdict |
|---|---|---|---|---|---|
| `multilpa()` | Mplus, LG, mclust, glca, multilevLCA, poLCA, in-house enumeration | see rows below | LL, parameters, posteriors, SEs, criteria | see below | covered |
| `lpa()` | Mplus 9 Demo, LG 6.1, tidyLPA formulas | `E/test-missing-patterns.R`, `E/test-tidylpa-classification.R` | LL, N, BIC, posteriors, entropy, AvePP | 5e-4 (printed) | covered |
| `lca()` | none calls the wrapper; engine covered by poLCA (harness `polca`) through `multilpa(categorical=)` | harness | item probabilities, LL, class shares | 1e-4 | partial (wrapper never called externally) |
| `multilca()` | as `lca()`; engine covered by glca, multilevLCA, Mplus ex10.7 | harness | as above | 1e-4 | partial (wrapper never called externally) |
| `lta()` | depmixS4, Mplus LTA, LMest, LG, Muthén & Asparouhov (2022) | `E/test-jstats-lta.R`, `E/lta-extensions/`, `T/test-lta-general.R`, LG c10/c11 | see `lta()` options | see below | covered |
| `mixture_regression()` | flexmix, lcmm, MASS, stats, sandwich | `E/test-mixture-regression-flexmix.R`, `E/test-growth-*-lcmm.R`, `T/` | see options | see below | covered |
| `enumerate_classes()`, `enumerate_lpa()`, `enumerate_lca()`, `enumerate_regressions()` | Mplus 9 Demo, one run per k | `E/test-enumeration-mplus.R` (`lpa_vvi_01` k = 1–5, `lca_binary_01` k = 1–4) | LL, n_par, BIC per candidate; BIC minimum | printed precision; only where Mplus replicated its best | covered (`enumerate_regressions()` internal only) |
| `three_step()` | tidySEM, LG `step3` | `E/test-tidysem-bch.R`, `E/three-step/compare-tidysem.R`, `E/latentgold/` c09 | BCH/modal/proportional estimates | exact / 5e-5 | partial (SEs not comparable) |
| `r3step()` | LG `step3 modal ml` | `E/latentgold/` c09, `E/test-latentgold-step3-markov.R` | LL, logits, SEs | 5e-5 / 5.1e-5 | covered |
| `bootstrap_lrt()` | Mplus 9 Demo TECH14 | `E/test-mplus-tech14.R` | LR statistic, parameter difference, null LL; p-value where the null is true | printed precision; Monte Carlo error | covered |
| `compare_models()` | Mplus BIC per k | `E/test-enumeration-mplus.R` | BIC and BIC weights | printed precision | covered |
| `parameter_inference()` | Mplus ML/MLR, LG, depmixS4, in-house structure reference | `E/test-mplus-inference.R`, `-robust.R`, `-covariates.R`, `-categorical.R`, `E/test-structure-inference.R`, `E/test-depmixs4-standard-errors.R` | SEs, estimates | 1e-6 to 5% | covered (Wald); bootstrap none |
| `sensitivity()` | none (harness `stability` refits the package itself across seeds) | internal only: `test-sensitivity.R` | — | — | none |
| `pool_imputations()` | none | internal only: `test-pool-imputations.R`; `T/test-audit-services-fixes.R` "Rubin pooling handles both zero-variance limits" | — | — | none |
| `fit_staged()` | none (Houle et al. supplement stages were not reproduced; no data) | internal only: `test-staged-identity.R` | — | — | none |
| `prior_control()` | mclust `priorControl()` | `E/test-mclust-prior.R` | MAP M-step, LL, BIC | 1e-8 / 1e-10 | covered |
| `starting_values()` | used only to start at Mplus estimates | `E/test-mplus-robust.R` | (input, not compared) | — | partial (only as a plumbing input) |
| `get_results()` | see section 3, one row per table | — | — | — | n/a |
| `candidate_fit()`, `get_candidate()` | accessors of enumeration grids | internal only | — | — | none |
| `descriptives()`, `diagnostics()`, `report()` | none | internal only (`test-diagnostics.R`, `test-descriptives-safety.R`) | — | — | none |
| `get_tna()`, `get_group_tna()` | none | internal only | — | — | none |
| `plot_*()` (9 functions) | — | — | — | — | n/a |

### `multilpa()` options

| Option | Reference | Where (test name) | Compared | Tol. | Verdict |
|---|---|---|---|---|---|
| `n_profiles`, `n_group_classes` (two-level core) | Mplus 9 Demo; LG 6.1; in-house enumeration | `E/test-mplus-twolevel-synthetic.R` "native multilevel LPA reproduces Mplus9…"; `E/test-mplus-twolevel-public.R` "matched multilevel fits…"; LG c01, c02, c12; `E/test-numerical-equivalence.R` | LL, means, variances, profile/group probabilities, both posteriors, AIC, BIC | 1e-4 LL, 1e-6/1e-4 params; LG 6e-5 LL, 1e-4 posteriors; enumeration 1e-10 | covered |
| `variance_model` varying / equal | Mplus, LG, mclust (VVI/EEI) | as above; `E/test-numerical-equivalence.R` "single group class reproduces mclust diagonal…" | LL, means, variances, posteriors | 1e-8 LL, 1e-6 params | covered |
| `covariance_model = "full"` | mclust VVV/EEE; Mplus FIML full; LG c03 | `E/test-mclust.R` "full covariance mixture limits…"; `E/test-mplus-missing-full.R`; LG c03 | LL, means, covariances, posteriors, n_par | 1e-7 LL, 3e-5 cov; Mplus 1e-5 | covered |
| `volume`, `shape`, `orientation`, `model` (14 structures) | mclust | `E/test-mclust.R` "every M-step reproduces mclust's own…", "the parameter counts…", "the likelihood is the mixture likelihood…" | M-step covariances, parameter counts, LL at mclust estimate (EVI) | 1e-8 (EVE/VVE: objective ≤ +1e-8) ; counts exact; LL 1e-10 | covered (single level). Two-level: SEs only, via in-house reference (`E/test-structure-inference.R`) |
| `profile_covariates`, `group_covariates` | Mplus (two-level covariates); LG c07 | `E/test-mplus-covariates.R` "one-step covariate model agrees…", "standard errors reproduce…", "the group-class contrast matches Mplus through vcov()" | LL, means, variances, both coefficient sets, posteriors, AIC, BIC, SEs | 5e-5 LL, 1e-5 params, 5e-3 SE; LG 5e-5 | covered |
| `profile_slopes = "group_class"` | none | internal only: `test-covariate-slopes.R` | — | — | none |
| `missing = "fiml"` | Mplus, LG, glca | `E/test-missing-patterns.R` (5 patterns); `E/test-mplus-missing-full.R`; LG c06; harness `glca` | LL, N, BIC, posteriors row by row, entropy, AvePP, counts | 5e-4 (printed); Mplus full 5e-5 LL, 1e-5 params | covered |
| `categorical` (binary/nominal) | Mplus categorical two-level; poLCA; glca; multilevLCA; Mplus ex10.7; LG c04, c05 | `E/test-mplus-categorical.R` "thresholds and their standard errors…", "categorical fits reproduce genuine Mplus two-level…"; harness | thresholds, SEs, probabilities, LL | 1e-4 thresholds, 2e-3 SE, 1e-3 LL; harness 1e-4 | covered |
| `ordinal` | LG 6.1 adjacent-category | `T/test-extra-indicators.R` "fits agree with Latent GOLD 6.1"; `E/latentgold-ordinal/compare.R` (script) | LL, n_par; category probabilities (script only) | 1e-4 LL; 5e-5 params (script) | covered |
| `count`, `count_model` poisson / negative_binomial, `count_dispersion` | LG `poisson` / `poisson overdispersed` | same as `ordinal` (o01–o07) | LL, n_par; count means (script) | 1e-4 | covered |
| `family = "additive"`, `"dispersion"`, `"additive_dispersion"`, `between_variance` | LG 6.1 (two-level + random intercept) | `E/test-latentgold-families.R` (finest quadrature: refine60, else refine, else base) | n_par, LL, group posteriors | 2e-4 LL, 1e-3 posteriors | covered (no parameter or SE comparison) |
| `family = "restricted_cross_level"`, `"full_cross_level"` | in-house published-formula likelihood | `T/test-cross-level.R` "full cross-level likelihood equals the published formula" | LL | 1e-8 | partial (in-house formula only) |
| `prior` | mclust | `E/test-mclust-prior.R` "each MAP M-step…", "the whole MAP EM reproduces mclust::me…" | covariances, means, LL, BIC, posteriors | 1e-8 to 1e-10 | covered |
| `noise` | mclust | `E/test-mclust-noise.R` (14 structures, with prior) | hypervolume, LL, means, covariances, noise prob., posteriors, n_par, BIC | 1e-8 LL, 1e-6 params, 1e-7 BIC | covered |
| `acceleration` | mclust (`"none"` in noise test), all others default squarem | `E/test-mclust-noise.R` | via LL | — | covered (both values used) |
| `start` | mclust | `E/test-mclust.R` (`max_iter = 0` at mclust's estimate) | LL at given values | 1e-10 | covered |
| `fixed` | in-house enumeration | `T/test-inference-fixed-measurement.R` "held-measurement standard errors match an independent numerical information" | LL, SEs | 1e-10 / 1e-6 | covered |
| `weights` | Mplus 9 `WEIGHT` (MLR) | `E/test-mplus-mixture-weights.R` (three-profile mixture); `E/weights/` via `E/test-external-scripts.R` | pseudo LL, posteriors, estimates, robust SEs (CR1 factor documented) | 1e-4 | covered (no two-level weights) |
| `centering` | none | internal only: `tests/testthat/test-centering.R` | — | — | none |
| `scaling` | Mplus DEFINE: STANDARDIZE; Mplus on person-standardized data | `E/test-mplus-standardize.R`, `E/test-mplus-person-scaling.R` | LL, estimates (up to the n vs n − 1 constant), posteriors | printed precision | covered |
| `drop_empty` | Mplus (TRUE); LG `missing includeall` (FALSE) | `E/test-missing-patterns.R`; `E/test-montecarlo.R` (`lpa_missing`, 5% empty rows) | N, counts, group N, entropy/AvePP; FALSE counts against LG | 0.05 counts; 5e-4 | covered |
| `min_variance`, `min_probability` | none | internal only | — | — | none |
| `select_start` | none | internal only | — | — | none |
| `time` (stored only, for sequence tables) | — | — | — | — | n/a |
| `n_starts`, `max_iter`, `tol`, `seed` | — | — | — | — | n/a |

### `lta()` options

| Option | Reference | Where (test name) | Compared | Tol. | Verdict |
|---|---|---|---|---|---|
| homogeneous first order (`time`, `n_profiles`) continuous | depmixS4 (JStats); LG c10 | `E/test-jstats-lta.R` "Gaussian LTA reproduces depmixS4…"; harness `depmixs4`; LG c10 | LL, n_par, AIC, BIC, means, SDs, initial, transitions; LG LL, posteriors, means, variances | 1e-4 LL, 1e-3 params; harness 1e-5 | covered |
| `categorical` | depmixS4; Mplus LTA (2 runs); Muthén & Asparouhov (2022) Table 5 | `E/test-jstats-lta.R` "binary LTA reproduces depmixS4…", "…the two Mplus LTA runs…", "plain LTA reproduces Table 5…" | LL, n_par, AIC, BIC (groups), SABIC, initial, transitions, item probabilities | 5e-4 LL, 1.5e-3 criteria, 6e-4 transitions | covered |
| `n_group_classes > 1` (mixture LTA) | LG c11 | `E/latentgold/compare.R` | LL, both posteriors, means, variances, group proportions | 6e-5 LL, 1e-4 posteriors | covered |
| `variance_model = "equal"` | LG 6.1 Markov | `E/test-montecarlo-extended.R` (`lta_equal`, 10 datasets) | n_par, LL, posteriors of every row | 2.5e-4 LL, 1e-3 | covered |
| `covariance_model = "full"`, `model` | none | `E/lta-extensions/README.md`: "no external LTA reference" | — | — | none |
| `missing = "fiml"` | LG 6.1 Markov, `missing includeall`; in-house path sum | `E/test-montecarlo-extended.R` (`lta_fiml`, 10 datasets, occasions entirely unobserved); `T/test-lta-general.R` | n_par, LL, posteriors | 2.5e-4 LL, 1e-3; 1e-9 | covered |
| `occasions` | depmixS4 (`"observed"`, unbalanced case) | harness `depmixs4` gaussian-unbalanced | LL, means, transitions | 1e-5 | partial (`"grid"` none) |
| `transitions = "occasion"` | LMest 3.2.8 | `T/test-lta-general.R` "occasion-varying transitions agree with LMest"; `E/lta-extensions/compare-lmest-occasion.R` | LL, n_par; transitions per occasion, initial (script) | 1e-8; 1e-4 | covered |
| `transition_covariates` | depmixS4; Mplus ex8.13, ex8.14 | `T/test-lta-general.R` "covariate transitions agree with depmixS4"; `E/lta-extensions/compare-depmix-covariate.R`, `compare-ex81[34].R` | LL, n_par, slopes, slope SEs, odds ratios | 1e-7 LL, 1e-3 slopes, 5% rel SE; Mplus 5e-3 LL, 2e-3/5e-3 | covered |
| `initial_covariates` | Mplus ex8.13, ex8.14 | `E/lta-extensions/compare-ex81[34].R` (script) | initial logits on x, `C1 ON CG#1` | 2e-3 | covered (script only) |
| `measurement = "occasion"` | Mplus 9 Demo (binary) | `T/test-lta-general.R` "occasion-specific measurement agrees with Mplus"; `E/lta-extensions/compare-mplus-measurement.R` | LL, n_par; thresholds, transitions (script) | 5e-3; 2e-3 | partial (continuous indicators: in-house path sum only) |
| `order = 2` | in-house path sum | `T/test-lta-general.R` "second-order transitions match the exact path sum…" | LL | 1e-9 | partial (in-house only) |
| `mover_stayer` | in-house path sum | `T/test-lta-general.R` "mover-stayer…" | LL | 1e-9 | partial (in-house only) |
| `ordinal`, `count`, `count_model`, `count_dispersion` | none | — | — | — | none |
| `weights` | none | internal only: `test-weights.R` (duplication) | — | — | none |
| `drop_empty`, `min_probability`, `select_start` | none | internal only | — | — | none |

### `mixture_regression()` options

| Option | Reference | Where (test name) | Compared | Tol. | Verdict |
|---|---|---|---|---|---|
| `family = "gaussian"` | flexmix; `lm()` | `E/test-mixture-regression-flexmix.R` "Gaussian rows"; `T/test-mixture-regression.R` "one class reduces to lm() and glm()" | LL at flexmix's estimates; optimum ≥ flexmix | 1e-9; −1e-7 | covered (estimates not compared: flexmix's sigma is not ML) |
| `family = "binomial"` | flexmix; `glm()` | "binomial, one class per group" (± concomitant); `T/test-post-submission-review.R` | LL, coefficients | 1e-9, 1e-5 | covered |
| `family = "poisson"` | flexmix; `glm()` | "Poisson rows and Poisson per group" | LL, coefficients; one-class SEs | 1e-9, 1e-5 | covered |
| `family = "negative_binomial"` | `MASS::glm.nb` (one class) | `T/test-negative-binomial-regression.R` "one class is MASS::glm.nb" | LL, coefficients, dispersion, SEs | 1e-8, 1e-6, 1e-5, 5% | partial (mixture: in-house dense likelihood only) |
| `family = "ordinal"` | `MASS::polr` (one class) | `T/test-ordinal-regression.R` "one class is MASS::polr…" | LL, thresholds, coefficients, SEs | 1e-8, 1e-5, 1e-3 | partial (mixture: in-house dense likelihood only) |
| `class_level = "group"`, `id` | flexmix (`y ~ x \| id`); lcmm | flexmix tests "…one class per group" | LL, coefficients | 1e-9, 1e-5 | covered |
| `n_group_classes > 1` (two-level regression) | in-house direct likelihood | `T/test-mixture-regression.R` | LL | 1e-10 | partial (in-house only) |
| `common` | flexmix `FLXMRglmfix`; lcmm | "Poisson with a coefficient shared across classes"; lcmm "common covariate and membership covariate" | LL, shared coefficient | 1e-7, 1e-5 | covered |
| `membership` | flexmix `FLXPmultinom`; lcmm `classmb` | "Gaussian rows with a concomitant model", "binomial… with a concomitant model"; lcmm | LL, coefficients | 1e-9 | covered |
| `group_membership` | none | — | — | — | none |
| `variance` varying / equal | flexmix (`varFix`); lcmm (equal) | "Gaussian with equal variances"; lcmm configurations | optimum ≥ flexmix; LL at hlme's estimates | −1e-7; 1e-8 | partial (flexmix: inequality only; equal checked fully via lcmm) |
| `random`, `random_covariance` equal / proportional, `random_diagonal` | lcmm `hlme` (`nwg`, `idiag`) | `E/test-growth-mixture-lcmm.R` (6 configurations) | LL at hlme's estimates; optimum ≥ hlme | 1e-8; −1e-6 | covered |
| `random_covariance = "varying"` | in-house dense likelihood (no hlme counterpart) | `T/test-growth-mixture.R` | LL | — | partial (in-house only) |
| `residual_correlation` ar / brownian, `residual_time` | lcmm `cor = AR()/BM()` | `E/test-growth-correlated-lcmm.R` | LL at hlme's estimates, optimum, n_par | 1e-8, 1e-4 | covered (`residual_time` only at its default) |
| `cluster` (multilevel growth) | in-house dense likelihood | `T/test-multilevel-growth.R` "the likelihood equals the dense multilevel likelihood" | LL | 1e-10 | partial (in-house only) |
| `vcov_type = "observed"` | `glm`, `MASS` one-class SEs | `T/test-mixture-regression.R`, `T/test-post-submission-review.R` | SEs | 1e-5 to 5% | partial (one-class limit only; mixture SEs: numDeriv in `T/test-growth-mixture.R`) |
| `vcov_type = "robust"` | `sandwich::vcovCL` | `T/test-cluster-robust-cr1.R` "a one-class robust regression reproduces sandwich::vcovCL()" | SEs | 1e-4 | partial (one class only) |
| `vcov_type = "opg"` | none | — | — | — | none |
| `weights` | none | internal only: `test-weights.R` | — | — | none |
| `missing = "omit"`, `select_start` | none | internal only | — | — | none |

## 2. Model families

| Family | Reference | Where | Compared | Tol. | Verdict |
|---|---|---|---|---|---|
| Single-level LPA, diagonal (EII, VII, EEI, VEI, EVI, VVI) | mclust; Mplus ex7.9/7.10; LG c09 | `E/test-mclust.R`, `E/test-mplus-public.R`, harness `mclust` | M-step, n_par, LL, means, variances, posteriors | 1e-8 | covered |
| Single-level LPA, full (EEE, VEE, EVE, VVE, EEV, VEV, EVV, VVV) | mclust | `E/test-mclust.R`, `E/test-mclust-noise.R`, `T/test-predict.R` | M-step, LL, posteriors | 1e-8 (EVE/VVE objective) | covered |
| Two-level LPA, diagonal | Mplus, LG | section 1 | LL, params, posteriors, SEs | 1e-4 to 1e-6 | covered |
| Two-level LPA, full covariance | Mplus FIML; LG c03 | `E/test-mplus-missing-full.R` | LL, covariances, posteriors | 1e-5 | covered |
| Two-level LPA, the other 12 mclust structures | Mplus 9 Demo and LG 6.1 (EEE, VVV); in-house structure reference (the rest) | `E/test-montecarlo-extended.R` (`ml_lpa_eee`, `ml_lpa_vvv`, 10 each); `E/test-structure-inference.R` | LL, n_par, profile and group posteriors; SEs | 5e-4, 1e-3; 1e-5 rel | partial (EEE and VVV covered; 10 structures in-house only) |
| LCA binary | poLCA, Mplus categorical, depmixS4/Mplus LTA, published tables (Agresti & Lang, Dayton, Linzer & Lewis, McCutcheon, Zhou & Lange) | harness `polca`, `published`; `E/test-mplus-categorical.R` | LL, item probabilities, shares, G², χ² | 1e-4; printed digits | covered |
| LCA polytomous / nominal | poLCA (values, election), glca, LG c04 | harness; LG | LL, item probabilities, posteriors | 1e-4 | covered |
| LCA ordinal | LG adjacent-category | `T/test-extra-indicators.R` | LL, n_par (+ category probabilities in script) | 1e-4 | covered |
| LCA count (Poisson, NB) | LG | same | LL, n_par | 1e-4 | covered |
| Mixed continuous + categorical | LG c05 | `E/latentgold/compare.R` | LL, posteriors, means, variances | 6e-5, 1e-4 | covered |
| Two-level LCA | Mplus ex10.7, glca, multilevLCA, Mplus categorical, LG c04 | harness `mplus-twolevel-lca`, `glca`, `multilevlca` | LL, item/prevalence probabilities, AIC, BIC both N | 1e-6 rel; 1e-4 | covered |
| Group-class additive / dispersion / additive-dispersion | LG 6.1 | `E/test-latentgold-families.R` | LL, n_par, group posteriors | 2e-4, 1e-3 | covered (no parameters or SEs) |
| Restricted / full cross-level | in-house formula | `T/test-cross-level.R`, `T/test-review-families.R` | LL | 1e-8 | partial (in-house only) |
| Transition model, homogeneous | depmixS4, Mplus, LG | section 1 `lta()` | — | — | covered |
| Transition model, general (occasion transitions, covariates, occasion measurement) | LMest, depmixS4, Mplus ex8.13/8.14, Mplus 9 | section 1 `lta()` | — | — | covered (second order, mover-stayer: partial) |
| Growth mixture (random effects, correlated residuals) | lcmm | `E/test-growth-*-lcmm.R` | LL, optimum, predictions | 1e-8 | covered |
| Multilevel growth mixture (`cluster`) | in-house only | `T/test-multilevel-growth.R` | LL | 1e-10 | partial |
| Mixture regression (row / group level) | flexmix | `E/test-mixture-regression-flexmix.R` | LL, coefficients, n_par | 1e-9, 1e-5, exact | covered |
| Two-level mixture regression | in-house only | `T/test-mixture-regression.R` | LL | 1e-10 | partial |
| One-profile random intercept (deferred, not exported) | Mplus | `E/test-mplus-random-intercept.R` | params, LL, AIC, BIC | 1e-5 | n/a (not exported) |

## 3. Reported quantities

| Quantity | Reference | Where | Tol. | Verdict |
|---|---|---|---|---|
| Log likelihood | every reference above | — | 1e-10 (mclust) to 5e-4 (Mplus printed) | covered |
| Parameter estimates (`profiles`, `responses`, `covariances`, `profile_probabilities`, `transitions`, `initial`, `coefficients`) | mclust, Mplus, LG, poLCA, glca, depmixS4, flexmix | sections 1–2 | 1e-8 to 1e-3 | covered |
| SEs, observed (Wald) | Mplus ML (diag and full FIML two-level), Mplus covariates/categorical, LG c01/c02/c07, depmixS4, LG step-3 | `E/test-mplus-inference.R` "ordinary ML Hessian…", "full-covariance FIML Hessian…" | 1e-6 (Mplus), 5e-5 (LG), 5% rel (depmixS4) | covered |
| SEs, robust (MLR sandwich) and scaling correction | Mplus MLR (CR0, rescaled by √(G/(G−1))); sandwich | `E/test-mplus-robust.R` "robust standard errors reproduce genuine Mplus MLR results" | 1e-6 | covered |
| SEs, `opg` | none | — | — | none |
| Individual posteriors (`posteriors`) | mclust, Mplus, LG, missing patterns | many | 1e-10 to 5e-4 | covered |
| Group posteriors (`group_posteriors`) | Mplus, LG, missing pattern m5 | `E/test-mplus-twolevel-*.R`, LG | 1e-6 to 5e-4 | covered |
| Modal assignments (`assignments`) | Mplus public refit | `E/test-mplus-twolevel-public.R` | identical | covered |
| `information_criteria`: AIC | Mplus, LG, depmixS4, glca, poLCA, multilevLCA | many | 5e-4 / 1.7e-4 / 1e-4 | covered |
| BIC, individuals convention | Mplus, LG, depmixS4, multilevLCA | many | 5e-4 / 1.7e-4 | covered |
| BIC, groups convention | LG, Mplus LTA (subjects), Table 5, multilevLCA | `E/test-jstats-lta.R`; LG | 1.5e-3 / 1.7e-4 / 1e-4 | covered |
| SABIC, both conventions | Mplus MLR (individuals), Mplus LTA (groups), LG (both) | `E/test-mplus-robust.R`, `E/test-jstats-lta.R`, LG | 1e-3 / 1.7e-4 | covered |
| CAIC, both conventions | LG | `E/latentgold/compare.R` | 1.7e-4 | covered |
| BIC on mclust's sign (`fit$bic`) | mclust | `E/test-mclust-prior.R`, `E/test-mclust-noise.R` | 1e-10 / 1e-7 | covered |
| Criteria arithmetic vs. published table | Houle et al. (2026) Table S1 | `E/houle-2026/check-fit-indices.R` (script) | 1e-9 | covered (AIC, CAIC, BIC, SABIC) |
| AWE, ICL, CLC, KIC, deviance | Mplus quantities (LL, n_par, N, entropy) | `E/test-montecarlo.R` (every criterion recomputed from Mplus's own quantities) | printed precision | covered |
| Entropy, individuals | Mplus (ex7.9/7.10, two-level FIML, 5 missing patterns), tidyLPA formula | `E/test-mplus-classification.R`, `E/test-missing-patterns.R`, `E/test-tidylpa-classification.R` | 5e-4 (printed); 1e-10 | covered |
| Entropy, groups | Mplus `CB` | `E/test-mplus-classification.R` "two-level FIML entropy and AvePP…"; m5 through Mplus's row averaging | 5e-4 | covered |
| AvePP (diagonal) | Mplus, tidyLPA | as entropy | 5e-4; 1e-12 | covered |
| `average_posteriors` off-diagonal cells | Mplus (full matrix) | `E/test-montecarlo.R` (every cell) | printed precision + modal slack | covered |
| OCC (`odds_correct_classification`) | Mplus AvePP and class proportions | `E/test-montecarlo.R` | 5% rel | covered |
| Class counts, estimated (`counts`, `estimated_n`) | Mplus posterior counts | `E/test-mplus-public.R`; `E/test-missing-patterns.R` "class counts reproduce Mplus's counts…" | 1e-3; 0.05 | covered |
| Class counts, modal (`n_modal`) | Mplus ex7.9/7.10 | `E/test-mplus-public.R` (`compare_public_mplus`) | exact | covered |
| `classification` table | Mplus (AvePP, counts, OCC) | `E/test-montecarlo.R`, `E/test-missing-patterns.R` | — | covered |
| `classification_errors` | tidySEM | `E/test-tidysem-bch.R` "the correction agrees with an independent implementation" | exact (default `expect_equal`) | covered |
| `bch_weights` | tidySEM | `E/three-step/compare-tidysem.R` (script) | 1.11e-15 recorded | covered (script only) |
| `residuals` (bivariate residuals) | LG c08 | `E/latentgold/compare.R` | ranking of pairs only | partial (values not compared) |
| Sample sizes (`n_informative`, `n_groups`, `n_informative_groups`) | Mplus, LG | `E/test-missing-patterns.R` | exact | covered |
| `sequences`, `sequence_summary`, `sequence_lengths`, `stages`, `starts`, `model`, `data` | — | — | — | n/a |

## 4. Inference and tests

| Item | Reference | Where | Compared | Tol. | Verdict |
|---|---|---|---|---|---|
| `bootstrap_lrt()` | Mplus 9 Demo TECH14 | `E/test-mplus-tech14.R` | LR statistic, parameter difference, null LL, p-value where the null is true | printed precision | covered |
| Lo-Mendell-Rubin | Mplus TECH11 | `E/test-mplus-lmr.R` | the adjustment formula applied to the fixture | 5e-3 | n/a (the package has no LMR function; the test does not call package code) |
| `parameter_inference(method = "wald")` | see SE rows in section 3 | — | — | — | covered |
| `parameter_inference(method = "bootstrap")` | none | internal only: `test-bootstrap-inference.R` | — | — | none |
| `parameter_inference(boundary = "fix")` | none | internal only | — | — | none |
| `parameter_inference` `level` (CI), `adjust` | none | — | — | — | none |
| `vcov()` contrasts | Mplus | `E/test-mplus-covariates.R` "the group-class contrast matches Mplus through vcov()" | contrast estimate, SE | 1e-3, 5e-3 | covered |
| `three_step(method = "bch")` | tidySEM; LG `step3 bch` | `E/test-tidysem-bch.R`; LG c09 | estimates | exact; 5e-5 | partial (SEs: LG gives model-based, the package cluster-robust; marked "not comparable") |
| `three_step(method = "modal")`, `"proportional"` | LG `step3` | LG c09 (script) | estimates | 5e-5 | partial (SEs not comparable) |
| `three_step(level = "groups")` | none | — | — | — | none |
| `three_step(contrast = "pairs")`, `adjust`, `vcov_type` | none | internal only: `test-three-step.R`, `test-three-step-clusters.R` | — | — | none |
| `r3step()` (individual level, `vcov_type = "observed"`) | LG `step3 modal ml` | LG c09 `covariate_ml` | logits, SEs | 5e-5 | covered |
| `r3step(initial_covariates, transition_covariates)` | LG 6.1 dynamic step 3 | `E/test-latentgold-step3-markov.R` "three-step latent Markov model matches Latent GOLD: …" | LL, every logit, every SE | 5.1e-5 | covered |
| `r3step(vcov_type = "robust")` (default), `level = "groups"`, `by_group_class` | none | internal only: `test-growth-r3step.R` | — | — | none |
| `sensitivity()` | none | internal only | — | — | none |

## 5. Data handling

| Item | Reference | Where | Tol. | Verdict |
|---|---|---|---|---|
| `missing = "fiml"`, single and two level | Mplus, LG, glca | `E/test-missing-patterns.R` (MCAR, empty rows, MAR, 70% one indicator, two-level with empty groups), `E/test-mplus-missing-full.R`, LG c06 | 5e-4 / 1e-5 | covered |
| `missing = "fiml"` in `lta()` | LG 6.1 Markov; in-house path sum | `E/test-montecarlo-extended.R`; `T/test-lta-general.R` | 2.5e-4; 1e-9 | covered |
| `drop_empty = TRUE` | Mplus | `E/test-missing-patterns.R` | exact N; 0.05 counts | covered |
| `drop_empty = FALSE` | LG `missing includeall` | `E/test-missing-patterns.R`; `E/test-montecarlo.R` (counts) | 5e-4 | covered |
| `weights` | Mplus `WEIGHT` (mixture, MLR) | `E/test-mplus-mixture-weights.R`; `E/weights/` (via `E/test-external-scripts.R`) | 1e-4 | partial (no two-level, no `lta()`/regression) |
| `centering` | none | internal only: `test-centering.R` | — | none |
| `scaling` | Mplus STANDARDIZE; Mplus on person-standardized data | `E/test-mplus-standardize.R`, `E/test-mplus-person-scaling.R` | printed precision | covered |
| Categorical binary / nominal | Mplus, poLCA, glca, LG | sections 1–2 | 1e-4 | covered |
| Ordinal (adjacent-category) | LG | `T/test-extra-indicators.R` | 1e-4 | covered |
| Count Poisson / NB | LG | same | 1e-4 | covered |
| `time` / occasions (`lta()`) | depmixS4 (unbalanced, `"observed"`) | harness `depmixs4` | 1e-5 | partial (`occasions = "grid"` none) |

## 6. Infrastructure

| Item | Reference | Where | Compared | Tol. | Verdict |
|---|---|---|---|---|---|
| Enumeration grids (`enumerate_*`, `candidate_fit`, `get_candidate`) | Mplus per k | `E/test-enumeration-mplus.R` | — | — | covered (`candidate_fit`, `get_candidate` through the grid) |
| `cores` (parallel) | equality on an externally checked grid | `E/test-enumeration-mplus.R` (`cores = 2` identical to `cores = 1`) | every candidate | identical | covered |
| `predict.multilpa` `type = "posterior"`, `"density"` | mclust `estep` | `T/test-predict.R` "predictions match mclust's E-step at the same parameters" | posteriors, total log density | 1e-10 | covered |
| `predict.multilpa` `type = "class"` | Mplus posteriors | `E/test-montecarlo.R` | posteriors and modal classes | 1e-3 | covered |
| `predict.latents_growth_mixture` (all six types) | lcmm `predictY`, `predictClass`, `predictRE` | `E/test-growth-predict-lcmm.R` | fitted values, posteriors, random effects | 1e-8 / 1e-6 | covered |
| `predict.latents_mixture_regression` (`response`, `class_response`, `posterior`, `probabilities`) | none | internal only | — | — | none |
| `predict.multilpa_transitions`, `predict.multilpa_covariates` | none | internal only | — | — | none |
| `compare_models()` | Mplus BIC per k | `E/test-enumeration-mplus.R` | — | — | covered |
| Start replication (`n_best_replicated`) | Mplus ("best LL replicated" printed) | `E/test-mplus-public.R` | count ≥ 2 | — | n/a |

## Verdict counts

Counted from the table rows above, n/a excluded: **covered 104, partial 26,
none 29** (8 rows n/a); at the 2026-10-17 audit, covered 82, partial 37,
none 40. A row whose verdict cell names two parts, such as
`parameter_inference()` ("covered (Wald); bootstrap none"), is counted by its
first word; the uncounted part has its own row in section 4.

## Gaps

Every row marked partial or none, each with a reference that could close it.
The available tools are Mplus 9 Demo (at most 6 dependent, 2 independent and
2 between variables), Latent GOLD 6.1 under Wine, and simulab for data with
known truth.

### Fitting functions and options
- `lca()` / `multilca()` wrappers (partial): run one poLCA and one glca harness case through the wrappers rather than `multilpa()`.
- `parameter_inference(method = "bootstrap")` (none): LG `bootstrap` SEs or Mplus `BOOTSTRAP =` on the c01 case; compare within Monte Carlo error.
- `parameter_inference(boundary = "fix")`, `level`, `adjust` (none): Mplus CIs (`CINTERVAL`) for `level`; `adjust` is `p.adjust` (assert against it directly).
- `sensitivity()` (none): simulab truth plus multiple-seed LG runs (`seed=` sets); no program reports this directly.
- `pool_imputations()` (none): Mplus `TYPE = IMPUTATION` on simulab data with 5 imputed sets (Rubin pooling of estimates and SEs).
- `fit_staged()` (none): Mplus with measurement fixed (`@` values) in a second run, the Houle et al. two-stage set-up.
- `starting_values()` (partial): assert that its output reproduces an Mplus `SVALUES` block.
- `candidate_fit()`/`get_candidate()`, `descriptives()`/`diagnostics()`/`report()`, `get_tna()`/`get_group_tna()` (none): accessors and summaries; check `descriptives()` against Mplus `SAMPSTAT` and leave the rest internal.
- `profile_slopes = "group_class"` (none): Mplus 9 Demo two-level with `%CB#1% CW ON x` (class-specific slopes), in the covariates kit.
- `family = "*_cross_level"` (partial): LG 6.1 with group-level indicators and a group-class `GClass` (cross-level), or the Mplus Demo restricted cross-level syntax from Houle et al. on simulab data.
- `centering` (none): LG/Mplus on the pre-centred or pre-scaled data (the transform done outside), compared with the package's in-call option.
- `min_variance`, `min_probability`, `select_start` (none): optimizer safeguards; an LG/Mplus case on a boundary solution (Heywood) would show whether the floors change the reported optimum.
- `lta()` `covariance_model = "full"`, `model` (none): LG dynamic model with `y1 <-> y2 | State`.
- `lta()` `occasions = "grid"` (partial): LG dynamic model with explicit missing occasions (rows for skipped waves, indicators missing).
- `lta()` `measurement = "occasion"`, continuous (partial): LG dynamic model with `y1 <- 1 + State + State*time` measurement per occasion, or Mplus Demo with occasion-specific means.
- `lta()` `order = 2` (partial): LG `State[-2]` in the transition equation.
- `lta()` `mover_stayer` (partial): LG mover-stayer (a `Class` with a fixed-identity transition matrix) or Mplus `ex8.15`-type syntax.
- `lta()` `ordinal`, `count`, `count_model`, `count_dispersion` (none): LG dynamic model with `ordinal` / `poisson` dependents.
- `lta()` `weights`, `drop_empty`, `min_probability`, `select_start` (none): LG `caseweight` for weights; the rest as above.
- `mixture_regression(family = "negative_binomial")` mixture (partial): flexmix has no NB driver; LG regression model with `poisson overdispersed` dependent, 2 classes.
- `mixture_regression(family = "ordinal")` mixture (partial): LG regression with an `ordinal` dependent, which needs a cumulative logit (LG offers `ordinal` as adjacent-category, so check the link first) or Mplus Demo `CATEGORICAL` regression mixture (cumulative logit).
- `mixture_regression(n_group_classes > 1)` (partial): LG two-level regression mixture (`GClass` plus `Class <- 1 + GClass`).
- `group_membership` (none): the same LG model with `GClass <- 1 + w`.
- `variance = "varying"/"equal"` against flexmix (partial): add an estimate comparison against LG regression (ML variances), where flexmix's (n − rank) divisor does not apply.
- `random_covariance = "varying"` (partial): Mplus Demo GMM with class-specific `i s` variances.
- `cluster` (multilevel growth) (partial): Mplus Demo `TYPE = TWOLEVEL MIXTURE` GMM with a between-level class.
- `vcov_type = "observed"` in mixtures (partial): lcmm `hlme` and flexmix `refit()` SEs at the shared estimates.
- `vcov_type = "robust"` in mixtures (partial): Mplus MLR on a two-class regression mixture.
- `vcov_type = "opg"` (none): no program reports OPG directly; numDeriv-based outer product of the scores.
- `mixture_regression` `weights`, `missing = "omit"`, `select_start` (none): flexmix `weights=` for case weights; the rest needs no external reference.

### Model families
- Two-level LPA, 12 non-diagonal mclust structures (partial): LG 6.1 two-level with the matching restrictions on `<->` (e.g. EEE = `y1 <-> y2;` without `| Cluster`).
- Cross-level families (partial): as above.
- Multilevel growth mixture (partial): as above.
- Two-level mixture regression (partial): as above.

### Reported quantities
- SEs `opg` (none): as above.
- `residuals` bivariate residual values (partial): compare LG's printed BVR values, not only their ranking (scale differences documented first).

### Inference and tests
- `parameter_inference(method = "bootstrap")` (none), `boundary`, `level`, `adjust` (none): as above.
- `three_step(method = "bch"/"modal"/"proportional")` SEs (partial): rerun LG c09 with `standarderrors=robust` to match the package's cluster-robust errors.
- `three_step(level = "groups")` (none): LG two-level step 3 with a group-level distal (`step3` on `GClass`).
- `three_step(contrast = "pairs")`, `adjust`, `vcov_type` (none): LG `step3` prints paired Wald tests; compare those.
- `r3step(vcov_type = "robust")`, `level = "groups"`, `by_group_class` (none): LG step 3 with robust SEs; LG two-level step 3 with `GClass` covariates.
- `sensitivity()` (none): as above.

### Data handling
- `weights` beyond a single-level mixture (partial), `centering` (none), `occasions = "grid"` (partial): as above.

### Infrastructure
- `predict.latents_mixture_regression` (none): flexmix `predict()` / `posterior()` at the shared estimates.
- `predict.multilpa_transitions`, `predict.multilpa_covariates` (none): depmixS4 `posterior()` (forward-backward) at the shared estimates; Mplus `SAVEDATA` posteriors for covariate models.

### Process gaps
- The harness suites (`E/harness/`) write a report rather than run as
  testthat assertions. Rerun on 0.95.1 (2026-10-18, source unchanged during
  the run): 1,018 of 1,018 quantities agree across 17 sources.
