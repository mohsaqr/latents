# latents roadmap

Rewritten 2026-09-28. **Models come first**: new model families and extensions
of the existing ones, ahead of inference refinements, workflow reach and
infrastructure. Every item keeps the package rule that whatever ships is
defensible: an independent likelihood check, parameter recovery by simulation,
interval coverage where standard errors are reported, and classed refusals for
what is not supported.

Size is a rough guide: **S** under half a day, **M** one to two days,
**L** several days or dependent on outside material.

Evidence for existing claims: `validation/MATH_AUDIT.md`,
`equivalence/houle-2026/REPORT.md`, `validation/*.R` and `equivalence/`.

---

## Priority 1 — Models

### 1.1 Multilevel model families (Houle et al.)

The core model has the dispersion-heterogeneity structure. The other five
families need explicit identification constraints, a simulation recovery study
each, and ideally Mplus runs as the oracle; the published supplement is
ambiguous on syntax and parameter counts, so it cannot be the oracle alone.

- [x] **Additive**, **dispersion**, **additive-dispersion** — shipped in 0.9.0
  as `multilpa(family = )` (experimental): exact likelihood, tables,
  inference, weak-class diagnostic, enumeration, bootstrap LRT; simulations
  in `validation/ADDITIVE_SIMULATION.md`, Latent GOLD comparison in
  `equivalence/latentgold-families/`.
- [x] **Restricted / full cross-level** — shipped in 0.9.0 (manifest
  group-mean specification); no standard errors.
- [ ] Group-class follow-ups: Mplus references; REML-type correction of the
  between-variance ML bias; a J = 100 cell for the weak-class threshold. **S–M**
- [ ] Cross-level follow-ups: enumeration over profiles x classes, recovery
  simulation, latent-intercept (doubly latent) variant. **M–L**
- [ ] Extend fixed-measurement (`fit_staged()`, `fixed =`) and the
  predictor/outcome workflows to the new families; keep their likelihoods
  distinct from BCH and three-step. **M**
- [ ] Latent aggregation and multilevel factor-score workflows, with the
  measurement-error assumptions stated. **M**

### 1.2 Latent transition models (`lta()`)

Standard errors are done (2026-09-27). The model itself is still first-order,
homogeneous and measurement-invariant.

- [x] **Covariate-dependent transitions** (fixed and time-varying covariates)
  and covariate initial distribution — done 2026-10-01; Mplus ex8.13/8.14
  and depmixS4 agree.
- [x] **Occasion-varying transitions** — done; LMest modBasic = 0 agrees.
- [x] **Measurement invariance across occasions** (`measurement = "occasion"`)
  — done; Mplus occasion-specific thresholds agree.
- [x] **Second-order transitions** (`order = 2`) — done; exact path sum.
- [x] The fourteen covariance structures in `lta()` (`model =`) — done for the
  homogeneous model; SEs only for EEI/VVI/EEE/VVV.
- [x] Mover-stayer class; FIML and covariance structures in the general engine (2026-10-01).
- [ ] Gaps vs JStats: random-intercept LTA (needs the exact Mplus RI-LTA
  specification), absorbing profiles
  (seqHMM-style). Sampling weights done (0.9.3). **M-L**
- [ ] Follow-ups: covariance structures and FIML in the general engine;
  structure-aware standard errors in lta(); likelihood-ratio and bootstrap
  tests between nested transition models; a simulation study of the new
  standard errors. **M**

### 1.3 Membership-covariate models

- [ ] The fourteen covariance structures in covariate models (today only the
  four `variance_model` × `covariance_model` combinations). **S–M**
- [ ] `start` and `fixed` for covariate models (both refused today). **M**
- [ ] Missing covariates beyond multiple imputation: integrate out a normal
  covariate by quadrature, for one or two continuous covariates only. Decide
  first whether the assumption is worth offering. **M**

### 1.4 Random effects

`fit_random_intercept()` is parked in `future/`; restoring it is a
precondition for everything here.

- [ ] Standard errors and intervals for random-intercept fits, handling the
  variance boundary and integration accuracy. **M**
- [ ] Validate random-intercept fits beyond the externally checked one-profile
  limit, including multiple-profile parameter recovery. **M**
- [ ] Random slopes and more general random-effects structures beyond a single
  Gaussian group intercept with unit loadings. **L**

### 1.5 Mixture regression (`mixture_regression()`)

- [ ] Multiple imputation for missing predictors: extend `pool_imputations()`
  (it pools `multilpa()` fits only). **S**
- [ ] Non-canonical links (probit, cloglog, identity-Poisson). **M**
- [ ] Random effects within classes (a random intercept per group inside each
  regression class). **L**

### 1.6 Across all models

- [x] **Sampling weights** (2026-10-01, 0.9.3): unit-level pseudo ML in every
  engine, sandwich SEs, integer weights = duplication, Mplus check.
- [ ] Weights follow-ups: within-unit (level-1) weights with WTSCALE-style
  scaling; stratification/PSU design-based variance; weighted `three_step()`
  and `r3step()`; weighted descriptive tables (`descriptives()`, `get_tna()`,
  classification counts are per-unit and unweighted); a weighted-mixture
  external reference (Mplus `TYPE = MIXTURE COMPLEX` with `WEIGHT`). **M**
- [ ] `fit_staged()` with a constrained covariance structure (refused today,
  because a structure maximized across profiles cannot be held one profile at
  a time). **M**

---

## Priority 2 — Choosing among models

- [x] `enumerate_classes(time = )` and `bootstrap_lrt()` for transition models (done 2026-10-01).
- [x] (done 0.9.2) `enumerate_classes()` and `bootstrap_lrt()` for transition models, with
  profile count and group-class count as a joint grid. **M**
- [ ] `bootstrap_lrt()` for random-intercept models, and nested comparisons
  beyond the one-class increment (for example diagonal against full
  covariance at equal class counts). **M**
- [ ] VLMR/LMR test with an independently validated reference distribution
  (`lmr_lrt()` is parked in `future/`; a chi-square approximation is not
  acceptable). **L**
- [ ] Calibrate the local-dependence diagnostics (bivariate residuals) for
  fitted measurement, classification uncertainty and clustering, for example
  by a parametric bootstrap; current p-values are approximate. **M**
- [ ] `compare_models(a, b, labels = )`, or `information_criteria()` over
  several fits, returning one labelled table. **S**
- [ ] Document the supported comparison families and failed-replicate
  behaviour; claim Mplus TECH14 equivalence only if separately established. **S**

---

## Priority 3 — Workflow reach

- [ ] **Public label alignment** between two fits of any family. Internal
  versions exist (`.multilpa_align_labels()` for `multilpa`,
  `.multilpa_pool_align()`/`.multilpa_cov_permute()` for covariate fits); an
  LTA version must reorder `transition_probabilities`, `transition_counts`,
  `initial_probabilities` and `empty_transition_rows` on both margins, and a
  two-profile test cannot catch a margin error because `c(2, 1)` is its own
  inverse. Prerequisite for the next item. **M**
- [ ] `sensitivity()` for `lta()` and covariate fits (refused today with
  `latents_unsupported_sensitivity`; needs refit paths and the alignment
  above). **M**
- [ ] Transition-aware `three_step()`, BCH, `r3step()` and classification
  errors: they read an LTA fit's posteriors but treat occasions as
  independent; the covariate model does not accept an LTA fit. **M**
- [ ] Error bars in `plot(what = "bars")` for transition fits (standard errors
  now exist; the plot does not draw them). **S**
- [x] `plot(what = "bars")` on an unconverged or bound-active fit now draws
  without whiskers and says why (2026-09-28).
- [ ] Scalar accessors so prose can quote one number with one verb:
  `entropy(fit, level = )`, `worst_residual(fit)`. **S**

---

## Naming consistency

- [ ] `three_step()` reports a `class` column where every other verb says
  `profile`/`group_class`; `mixture_regression()` uses `class_k` throughout
  its own family (consistent within it). Decide one vocabulary.

## Priority 4 — Decisions for the owner

- [ ] Cluster-robust variance is CR0. Apply the `G/(G-1)` correction? Moves
  every cluster-robust standard error by about 1% and the calibrated bounds in
  `tests/testthat/test-tidy-three-step.R`.
- [ ] Should a separated covariate fit be refused `parameter_inference()` /
  `vcov()`? It is today, because it is reported unconverged.
- [ ] `group_probabilities` on a transition fit is one M-step stale, so
  "class-weighted prevalence = overall profile share" holds only to EM
  tolerance (4.7e-6 at `max_iter = 60`).
- [ ] Missing covariates: keep refusing under FIML and point to
  `pool_imputations()` (current), or add the quadrature option in 1.3.

---

## Priority 5 — Release and infrastructure

- [x] Bump `DESCRIPTION` (0.8.4) before any push. **S**
- [x] (2026-10-01, 0.9.4: CI repaired) CI has never run: GitHub Actions on release, devel and oldrel across
  Linux, macOS and Windows; only macOS with R 4.5.2 so far. **S**
- [ ] Re-audit the README feature matrix. Its residual-covariance row still
  says four of fourteen structures; all fourteen are fitted with Wald standard
  errors, and transition models, pooled imputations and the extended bootstrap
  are new since the last audit. **S**
- [ ] Vignette coverage: exercised 14 of 28 exports when audited on
  2026-09-20; add the verbs added since (`mixture_regression()`,
  `pool_imputations()`, transition inference, `predict()`). Re-count against
  `NAMESPACE` before claiming coverage. **M**
- [ ] Reproduce the Houle et al. fitted analyses once participant data, fitted
  input/output and preprocessing are available. **L**
- [ ] The six rendered plots in `Tutorial_docs/two-level-lpa-engagement.Rmd`
  have not been looked at by a person.

---

## Done

| Date | Item | Evidence |
|---|---|---|
| 2026-09-28 | `bootstrap_lrt()` for FIML and membership-covariate fits | size 0.040 / 0.053 at 5% (`validation/bootstrap-lrt-size.R`) |
| 2026-09-27 | Standard errors for `lta()` fits | Richardson Hessian 1e-6, depmixS4; coverage 0.947–0.960 |
| 2026-09-27 | `pool_imputations()` for missing covariates | `mice::pool.scalar` to 1e-12; coverage by imputation model |
| 2026-09-27 | Two-stage bootstrap for `fit_staged()` | coverage 0.97 against 0.86–0.89 conditional |
| 2026-09-27 | Profile slopes by group class (`profile_slopes`) | coverage 0.933–0.970 |
| 2026-09-27 | FIML in membership-covariate models | likelihood 1e-12; coverage 0.942–0.954 |
| 2026-09-27 | Standard errors for categorical covariate fits (found already done) | score 4.4e-08 |
| 2026-09-27 | Faster EVE/VVE, Ward first start, `predict()`, single-level enumeration | 13/14 structures faster than mclust |
| 2026-09-27 | SQUAREM EM acceleration, vectorised E-step | E-step 14× faster, identical output |
| 2026-09-27 | `mixture_regression()`; mclust prior, noise, SEs for all 14 structures | flexmix 22 checks; coverage 0.915–0.965 |
| 2026-09-22 | `sensitivity()` for `multilpa()`; transition-model `plot()` | fixed a permutation defect |
| 2026-09-21 | Centering (`person`, `grand`); `lmr_lrt()` and random intercept moved to `future/` | |
| 2026-09-20 | External review remediated; tutorial rebuilt; feature matrix audited | `NEWS.md` 0.11.0 |
| 2026-09-19 | `fixed =`, `fit_staged()`, covariate models with categorical, mixed and full covariance, `lta()`, CLC and KIC, estimates with SEs from one verb | `validation/MATH_AUDIT.md` |

- [x] Ordinal (adjacent-category) and Poisson count indicators in the core
  engine (2026-10-01, 0.9.5; Latent GOLD-checked).
- [ ] Ordinal/count follow-ups: covariate models, `lta()`, `start`/`fixed`,
  negative binomial and zero-inflated counts, a cumulative-logit ordinal
  option, ordinal/count in the bivariate-residual table and plots. **M**
