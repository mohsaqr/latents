# latents roadmap

Rewritten 2026-09-28. **Models come first**: new model families and
extensions of the existing ones, ahead of inference refinements,
workflow reach and infrastructure. Every item keeps the package rule
that whatever ships is defensible: an independent likelihood check,
parameter recovery by simulation, interval coverage where standard
errors are reported, and classed refusals for what is not supported.

Size is a rough guide: **S** under half a day, **M** one to two days,
**L** several days or dependent on outside material.

Evidence for existing claims: `validation/MATH_AUDIT.md`,
`equivalence/houle-2026/REPORT.md`, `validation/*.R` and `equivalence/`.

------------------------------------------------------------------------

## Priority 0 — One engine (plan)

Consolidate the eight estimation engines into measurement blocks ×
latent structures on one kernel with shared services. Plan, phases,
acceptance and risks:
[validation/ENGINE_DESIGN.md](https://pak.dynasite.org/latents/validation/ENGINE_DESIGN.md).
Public API unchanged; every phase is a release that reproduces golden
fixtures.

**E0. Golden baseline** of every engine and option. **M** (251 cases,
bit-identical reruns)

**E1. Shared inference service.** **M**

**E2. Shared alignment and simulation services.** **M**

**E3. Kernel; port mixture regression and growth mixture.** **L**

**E4. Port profile models.** **L**

**E5. Port transition models.** **L**

**E6. Port Houle families.** **M-L** (closed-form additive likelihood
kept)

\[~\] **E7. New combinations**: multilevel growth and negative-binomial
regression mixtures done; growth with categorical/ordinal items,
transitions of regression/growth blocks, ordinal regression mixtures
open. **M each**

------------------------------------------------------------------------

## Priority 1 — Models

### 1.1 Multilevel model families (Houle et al.)

The core model has the dispersion-heterogeneity structure. The other
five families need explicit identification constraints, a simulation
recovery study each, and ideally Mplus runs as the oracle; the published
supplement is ambiguous on syntax and parameter counts, so it cannot be
the oracle alone.

**Additive**, **dispersion**, **additive-dispersion** — shipped in 0.9.0
as `multilpa(family = )` (experimental): exact likelihood, tables,
inference, weak-class diagnostic, enumeration, bootstrap LRT;
simulations in `validation/ADDITIVE_SIMULATION.md`, Latent GOLD
comparison in `equivalence/latentgold-families/`.

**Restricted / full cross-level** — shipped in 0.9.0 (manifest
group-mean specification); no standard errors.

Group-class follow-ups: Mplus references; REML-type correction of the
between-variance ML bias; a J = 100 cell for the weak-class threshold.
**S–M**

Cross-level follow-ups: enumeration over profiles x classes, recovery
simulation, latent-intercept (doubly latent) variant. **M–L**

Extend fixed-measurement
([`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md),
`fixed =`) and the predictor/outcome workflows to the new families; keep
their likelihoods distinct from BCH and three-step. **M**

Latent aggregation and multilevel factor-score workflows, with the
measurement-error assumptions stated. **M**

### 1.2 Latent transition models (`lta()`)

Standard errors are done (2026-09-27). The model itself is still
first-order, homogeneous and measurement-invariant.

**Covariate-dependent transitions** (fixed and time-varying covariates)
and covariate initial distribution — done 2026-10-01; Mplus ex8.13/8.14
and depmixS4 agree.

**Occasion-varying transitions** — done; LMest modBasic = 0 agrees.

**Measurement invariance across occasions** (`measurement = "occasion"`)
— done; Mplus occasion-specific thresholds agree.

**Second-order transitions** (`order = 2`) — done; exact path sum.

The fourteen covariance structures in
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) (`model =`)
— done for the homogeneous model; SEs only for EEI/VVI/EEE/VVV.

Mover-stayer class; FIML and covariance structures in the general engine
(2026-10-01).

Gaps vs JStats: random-intercept LTA (needs the exact Mplus RI-LTA
specification), absorbing profiles (seqHMM-style). Sampling weights done
(0.9.3). **M-L**

Covariance structures and FIML in the general engine; bootstrap LRT
between nested transition models; simulation study of the new standard
errors (all 0.9.2; audited 2026-10-01).

Structure-aware standard errors in the extended
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) beyond the
diagonal EEI and VVI (refused today). **M**

### 1.3 Membership-covariate models

The fourteen covariance structures in covariate models (today only EEI,
VVI, EEE, VVV). Until 0.9.6 any other `model =` was silently fitted as
one of those four; it is now refused (`latents_unsupported_structure`).
**S–M**

`start` and `fixed` for covariate models (both refused today). **M**

~~Missing covariates beyond multiple imputation by quadrature~~ —
decided against (2026-10-15);
[`pool_imputations()`](https://pak.dynasite.org/latents/reference/pool_imputations.md)
is the route.

### 1.4 Random effects

`fit_random_intercept()` is parked in `future/`; restoring it is a
precondition for everything here.

- \[~\] Parked: `fit_random_intercept()` was removed in 0.8.0 and lives
  in `future/`. Bringing it back needs standard errors (variance
  boundary, integration accuracy), multiple-profile validation, and then
  random slopes. Only if the continuous group-effect family is wanted.
  **M–L**

### 1.5 Mixture regression (`mixture_regression()`)

Multiple imputation for missing predictors: extend
[`pool_imputations()`](https://pak.dynasite.org/latents/reference/pool_imputations.md)
to
[`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
(it calls
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
only; confirmed 2026-10-01). **S**

Non-canonical links (probit, cloglog, identity-Poisson). **M**

Random effects within classes: planned as growth mixture models, see
1.6.

### 1.6 Growth models (latent class growth and growth mixture models)

Status (2026-10-02). Latent class growth analysis (group-based
trajectories, Nagin) already works through
`mixture_regression(y ~ time, data, n_classes, id = "id", class_level = "group")`
for Gaussian, binomial and Poisson outcomes. It matches `lcmm::hlme`
without random effects: log likelihood -2857.281151 vs -2857.281176,
same estimates up to class order (300 persons x 6 occasions, 3 classes).
Growth mixture models (random intercepts/slopes *within* classes,
`hlme(random = ~ time)`) are not supported.

**Interface (decided before coding).** One verb, named arguments: extend
[`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
rather than add a function. `random = ~ 1 + time` (within-class random
effects; default `NULL` keeps today’s model),
`random_covariance = c("varying", "equal", "proportional")`
(class-specific, common, or common times a class scale, the `lcmm`
`nwg`), `random_diagonal = FALSE` (`lcmm` `idiag`), and later
`residual_correlation = c("none", "ar1", "brownian")`. A growth mixture
model needs `id` and `class_level = "group"`; anything else is refused
by class.

**References.** `lcmm::hlme` 2.2.1 (installed; exact for Gaussian, no
integration) and Latent GOLD 6.1 Syntax (installed under Wine; latent
growth and mixture growth models with CFactors, Syntax Guide sections
5.3.4-5.3.5 and 5.4.3; quadrature, so evaluate our exact likelihood at
its estimates before calling a gap a disagreement).

**Phases** (each ships only with its evidence; sizes as at the top):

**G0. Document what exists.** (2026-10-02: vignette `growth-mixture`,
dataset `growth_scores`, trajectory views for models without random
effects.) Example and vignette section “trajectory classes” with
[`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md);
a tidy `trajectories` table (predicted mean per class and time) from
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md);
add the `lcmm::hlme` match as an `equivalence/` fixture. **S** (half a
day)

**G1. Gaussian growth mixture model, core.** (2026-10-02;
`equivalence/test-growth-mixture-lcmm.R`, six configurations.) Within
class k, y_i ~ N(X_i b_k, Z_i G_k Z_i’ + s_k^2 I), unbalanced times and
missing occasions allowed. Closed-form marginal likelihood (Cholesky per
person); EM over classes with a per-class weighted linear-mixed-model
M-step (inner EM on the random effects, or Newton on the log-Cholesky of
G_k); starts from the latent class growth fit plus random starts;
SQUAREM. `random_covariance` (`varying`/`equal`/`proportional`) and
`random_diagonal`. Classed refusals for boundary variances. Evidence:
independent dense-matrix likelihood and analytic score checks; `hlme`
matches (log likelihood within 1e-6, estimates within 1e-4) for at least
6 configurations (1-4 classes; random intercept / intercept+slope; each
covariance option; unbalanced data). **M** (about 1.5 days)

\[~\] **G2. Inference.** (Observed and sandwich errors done and checked
against a numerical Hessian; bootstrap and the coverage study remain.)
Observed-information and sandwich (weighted) standard errors on the
estimation scale, delta method to variances/correlations; boundary
handling for random-effect variances; person bootstrap via the existing
machinery. Evidence: numerical-Hessian agreement 1e-6; coverage
simulation (\>= 200 datasets, gates as in
`validation/LTA_SIMULATION.md`). **M** (about 1 day)

\[~\] **G3. Integration.** (Tables, print, summary, plots, weights,
[`simulate()`](https://rdrr.io/r/stats/simulate.html), enumeration with
the bootstrap LRT and distal outcomes in
[`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
done; [`predict()`](https://rdrr.io/r/stats/predict.html) for new
persons and
[`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md)
remain, both refused with classed conditions.) `print`/`summary`/`plot`
(class trajectories with individual-curve spread),
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
[`predict()`](https://rdrr.io/r/stats/predict.html) (marginal and
conditional on the random effects, i.e. empirical Bayes), posterior
class table, `enumerate` over classes,
[`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
(needs simulation from the fitted model), label alignment, sampling
weights,
[`sensitivity()`](https://pak.dynasite.org/latents/reference/sensitivity.md).
Evidence: tests per verb; bootstrap LRT size check. **M** (about 1 day)

**G4. Validation report.** Recovery and coverage study (balanced and
unbalanced designs, weak and strong separation, random-effect variance
near zero); Latent GOLD cross-check on two Gaussian configurations;
write-up in `validation/`. **S-M** (half a day, mostly compute)

**G5. Correlated residuals** (`residual_correlation = "ar1"`,
`"brownian"`). Changes only the residual covariance in G1;
`hlme(cor = ...)` as reference. **M** (about 1 day)

**G6. Binary and count outcomes with random effects.** Adaptive
Gauss-Hermite quadrature (Laplace as a fast option) inside the class EM;
references Latent GOLD (binary, Poisson growth mixtures) and `lme4` /
`GLMMadaptive` for one class. Accuracy check against increasing nodes.
**L** (1.5-2 weeks)

**G7. Out of scope unless asked:** latent-process link functions
(`lcmm::lcmm`), multivariate outcomes (`multlcmm`), joint latent class
and survival models (`Jointlcmm`). Each is a separate model family
(weeks).

Order: G0 -\> G1 -\> G2 -\> G3 -\> G4 gives a releasable Gaussian growth
mixture model in about 4 working days; G5 next; G6 only with a concrete
need. Risk: convergence with random-effect variances on the boundary or
weakly separated classes may add a day to G1-G2.

### 1.7 Across all models

**Sampling weights** (2026-10-01, 0.9.3): unit-level pseudo ML in every
engine, sandwich SEs, integer weights = duplication, Mplus check.

Weights follow-ups: within-unit (level-1) weights with WTSCALE-style
scaling; stratification/PSU design-based variance; weighted
[`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
and [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md);
weighted descriptive tables
([`descriptives()`](https://pak.dynasite.org/latents/reference/descriptives.md),
[`get_tna()`](https://pak.dynasite.org/latents/reference/get_tna.md),
classification counts are per-unit and unweighted); a weighted-mixture
external reference (Mplus `TYPE = MIXTURE COMPLEX` with `WEIGHT`). **M**

[`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
with a constrained covariance structure (refused today, because a
structure maximized across profiles cannot be held one profile at a
time). **M**

------------------------------------------------------------------------

## Priority 2 — Choosing among models

`enumerate_classes(time = )` and
[`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
for transition models (done 2026-10-01).

(done 0.9.2)
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
and
[`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
for transition models, with profile count and group-class count as a
joint grid. **M**

[`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
for nested comparisons beyond the one-class increment in
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
(diagonal against full covariance at equal class counts is refused today
with `latents_incomparable_models`; transition and family fits already
have a nesting rule). **M**

\[~\] Parked: VLMR/LMR test (`lmr_lrt()` in `future/`) until its
reference distribution is independently validated; a chi-square
approximation is not acceptable. **L**

Calibrate the local-dependence diagnostics (bivariate residuals) for
fitted measurement, classification uncertainty and clustering, for
example by a parametric bootstrap; current p-values are approximate.
**M**

`compare_models(a, b, labels = )`, or `information_criteria()` over
several fits, returning one labelled table. **S**

Document the supported comparison families and failed-replicate
behaviour; claim Mplus TECH14 equivalence only if separately
established. **S**

------------------------------------------------------------------------

## Priority 3 — Workflow reach

**Public label alignment** between two fits of any family. Internal
versions exist (`.multilpa_align_labels()` for `multilpa`,
`.multilpa_pool_align()`/`.multilpa_cov_permute()` for covariate fits);
an LTA version must reorder `transition_probabilities`,
`transition_counts`, `initial_probabilities` and `empty_transition_rows`
on both margins, and a two-profile test cannot catch a margin error
because `c(2, 1)` is its own inverse. Prerequisite for the next item.
**M**

[`sensitivity()`](https://pak.dynasite.org/latents/reference/sensitivity.md)
for [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) and
covariate fits (refused today with `latents_unsupported_sensitivity`;
needs refit paths and the alignment above). **M**

Transition-aware
[`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md),
BCH, [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md)
and classification errors: they read an LTA fit’s posteriors but treat
occasions as independent; the covariate model does not accept an LTA
fit. **M**

Error bars in `plot(what = "bars")` for transition fits (standard errors
now exist; the plot does not draw them). **S**

`plot(what = "bars")` on an unconverged or bound-active fit now draws
without whiskers and says why (2026-09-28).

Scalar accessors so prose can quote one number with one verb:
`entropy(fit, level = )`, `worst_residual(fit)`. **S**

Named plot functions for fitted models (2026-10-07):
[`plot_profiles()`](https://pak.dynasite.org/latents/reference/latents-plots.md):
the `plot(fit, what = )` sweep hides each view’s arguments from editor
completion.
[`plot_enumeration()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
(2026-10-07) is the pattern: the generic delegates to the named
function, which has no `...`. Owner’s request. **M**

------------------------------------------------------------------------

## Naming consistency

[`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
reports a `class` column where every other verb says
`profile`/`group_class`;
[`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
uses `class_k` throughout its own family (consistent within it). Decide
one vocabulary.

## Priority 4 — Decisions for the owner

All four decided 2026-10-15 (0.9.31).

Cluster-robust variance: **CR1** (`G/(G-1)`) everywhere, matching
`sandwich::vcovCL()`; one helper, `.inference_cluster_factor()`, and a
required `n_units` on `.inference_sandwich()` so no robust path can skip
it.

Separated covariate fit: **refuse** Wald inference
(`latents_boundary_fit`) on a unit-free check (coefficient per RMS unit
of its covariate \> 20), shared with the `latents_extreme_coefficients`
warning. The old refusal relied on a singular information and failed on
some seeds.

Transition `group_probabilities` vs posteriors: **tighter tolerance**,
`lta(tol = 1e-10)`; gap 2e-4 -\> 2e-5 for ~1.5x iterations.

Missing covariates: **keep refusing** under FIML (`latents_bad_data`),
pointing to
[`pool_imputations()`](https://pak.dynasite.org/latents/reference/pool_imputations.md);
quadrature dropped from 1.3.

------------------------------------------------------------------------

## Priority 5 — Release and infrastructure

Bump `DESCRIPTION` (0.8.4) before any push. **S**

(2026-10-01, 0.9.4: CI repaired) CI has never run: GitHub Actions on
release, devel and oldrel across Linux, macOS and Windows; only macOS
with R 4.5.2 so far. **S**

README update: it is now a tutorial (the feature matrix with the
four-of-fourteen claim is gone) and does not mention sampling weights,
ordinal/count indicators, the Houle families, the
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) extensions
or
[`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md).
**S**

Vignette coverage: exercised 14 of 28 exports when audited on
2026-09-20; add the verbs added since
([`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md),
[`pool_imputations()`](https://pak.dynasite.org/latents/reference/pool_imputations.md),
transition inference,
[`predict()`](https://rdrr.io/r/stats/predict.html)). Re-count against
`NAMESPACE` before claiming coverage. **M**

Reproduce the Houle et al. fitted analyses once participant data, fitted
input/output and preprocessing are available. **L**

The six rendered plots in `Tutorial_docs/two-level-lpa-engagement.Rmd`
have not been looked at by a person.

------------------------------------------------------------------------

## Done

| Date | Item | Evidence |
|----|----|----|
| 2026-10-02 | 0.9.11 (pending push): false non-convergence of transition fits (L-BFGS-B code 52) | 5/120 test datasets flagged at the maximum, 0 after; restart gain \< 1e-12 |
| 2026-10-02 | 0.9.9: full-covariance transition simulation; transition bootstrap inference | coverage 0.924-0.928, bootstrap LRT size 0.065 (200 datasets) |
| 2026-10-02 | 0.9.8: review of 0.9.0-0.9.7 plus six follow-up fixes | 5,332 assertions; `validation/review-2026-10-0*` |
| 2026-09-28 | [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md) for FIML and membership-covariate fits | size 0.040 / 0.053 at 5% (`validation/bootstrap-lrt-size.R`) |
| 2026-09-27 | Standard errors for [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fits | Richardson Hessian 1e-6, depmixS4; coverage 0.947–0.960 |
| 2026-09-27 | [`pool_imputations()`](https://pak.dynasite.org/latents/reference/pool_imputations.md) for missing covariates | [`mice::pool.scalar`](https://amices.org/mice/reference/pool.scalar.html) to 1e-12; coverage by imputation model |
| 2026-09-27 | Two-stage bootstrap for [`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md) | coverage 0.97 against 0.86–0.89 conditional |
| 2026-09-27 | Profile slopes by group class (`profile_slopes`) | coverage 0.933–0.970 |
| 2026-09-27 | FIML in membership-covariate models | likelihood 1e-12; coverage 0.942–0.954 |
| 2026-09-27 | Standard errors for categorical covariate fits (found already done) | score 4.4e-08 |
| 2026-09-27 | Faster EVE/VVE, Ward first start, [`predict()`](https://rdrr.io/r/stats/predict.html), single-level enumeration | 13/14 structures faster than mclust |
| 2026-09-27 | SQUAREM EM acceleration, vectorised E-step | E-step 14× faster, identical output |
| 2026-09-27 | [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md); mclust prior, noise, SEs for all 14 structures | flexmix 22 checks; coverage 0.915–0.965 |
| 2026-09-22 | [`sensitivity()`](https://pak.dynasite.org/latents/reference/sensitivity.md) for [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md); transition-model [`plot()`](https://rdrr.io/r/graphics/plot.default.html) | fixed a permutation defect |
| 2026-09-21 | Centering (`person`, `grand`); `lmr_lrt()` and random intercept moved to `future/` |  |
| 2026-09-20 | External review remediated; tutorial rebuilt; feature matrix audited | `NEWS.md` 0.11.0 |
| 2026-09-19 | `fixed =`, [`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md), covariate models with categorical, mixed and full covariance, [`lta()`](https://pak.dynasite.org/latents/reference/lta.md), CLC and KIC, estimates with SEs from one verb | `validation/MATH_AUDIT.md` |

Ordinal (adjacent-category) and Poisson count indicators in the core
engine (2026-10-01, 0.9.5; Latent GOLD-checked).

Ordinal/count in covariate models and
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md);
negative-binomial counts (2026-10-02, 0.9.7; Latent GOLD-checked).

Ordinal/count follow-ups: `start`/`fixed`, zero-inflated counts, a
cumulative-logit ordinal option, ordinal/count in the bivariate-residual
table and plots. **M**
