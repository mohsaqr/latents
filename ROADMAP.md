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

## Next (set 2026-10-17, after 0.9.37)

Owner priority remains growth, transitions and
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md). Each item
ships with its evidence (an external reference, a simulation, or both).

| \# | Item | Reference / evidence | Size | Target |
|----|----|----|----|----|
| 1 | Coverage cell for correlated residuals (`residual_correlation`), and coverage of the growth bootstrap intervals | add cells to `validation/growth-simulation.R` (new registry rows) | S | 2026-10-19 |
| 2 | Transition-aware [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md) (distal outcomes of an [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) profile sequence) and BCH | Latent GOLD `step3` on a dynamic latent variable | M | 2026-10-21 |
| 3 | Group classes in the three-step latent Markov chain ([`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md) on `lta(n_group_classes > 1)`) | Latent GOLD step3 with a group latent variable | M | 2026-10-23 |
| 4 | Standard errors for the ten constrained structures in the homogeneous [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) (delta method from the general engine’s logits) | numerical Hessian; agreement with the extended engine | S-M | 2026-10-24 |
| 5 | Multilevel growth bootstrap (match group classes too) | simulation | M | 2026-10-27 |
| 6 | Vignette sections: growth [`predict()`](https://rdrr.io/r/stats/predict.html)/[`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md)/correlated residuals; three-step latent Markov | rendered articles, claims tested in tests/equivalence | M | 2026-10-29 |
| 7 | G6: binary and count growth mixtures (adaptive quadrature) | Latent GOLD, GLMMadaptive | L | November |

## Priority E — Equivalence 2.0 (plan, set 2026-10-17)

Owner priority: equivalence is to be very strong. The suite has breadth
(Mplus, Latent GOLD, mclust, lcmm, flexmix, depmixS4, tidySEM, tidyLPA;
700 expectations at 0.95.0) but the work of 2026-10-17 exposed three
weaknesses:

- **Tests that could not fail.** Every group-level statistic replaced by
  a perfect 1 still passed 26 of 26 expectations, because the two-level
  references are so well separated that everything prints as 1.000.
- **Tests that passed for the wrong reason.** AvePP matched by sorting
  apart from the class alignment; reference ids unchecked; `run.sh`
  reusing stale outputs; a two-level count parser reading the
  joint-pattern table.
- **Gaps no test covered.** Rows with nothing observed were counted in
  entropy for several releases: no reference dataset had such a row, and
  the Mplus helper computed entropy itself instead of reading the
  package’s value.

Resources on this machine: genuine **Mplus 9 Demo**
(`/Applications/MplusDemo/mpdemo`) and **Latent GOLD 6.1** under Wine
(`~/.wine-latentgold`), both run in batch; **simulab** 0.4.10
(`../simulab`), whose generators store the true parameters with every
dataset (`simulate_lpa()`, `simulate_ml_lpa()`, `simulate_lca()`,
`simulate_hmm()`, `simulate_markov()`, `simulate_growth()`,
`inject_missingness()` with MCAR/MAR/MNAR, `apply_batch()`); and
**JStats** (`../JStats/validation/r-reference/`), 95 genuine Mplus
mixture outputs, 62 with their data file, of which 4 transition fixtures
are used so far.

Constraint to verify first: the Mplus **Demo** limits model size (number
of dependent variables and covariates); configurations beyond it need a
full licence or Latent GOLD alone.

| \# | Component | Deliverable | Acceptance | Size |
|----|----|----|----|----|
| E2.1 | **Coverage audit** | `equivalence/COVERAGE.md`: every exported verb and every option of the fitting functions against the external reference that checks it (program, file, tolerance) or “none” | Every row filled; every “none” either gets an item below or a stated reason | S |
| E2.2 | **Moderately separated references** | Two-level LPA and LCA cases with entropy 0.5–0.8 at both levels, through Mplus and Latent GOLD | A deliberately perfect group-level statistic fails the suite (the regression the 2026-10-17 review found) | S |
| E2.3 | **Monte Carlo equivalence** | `equivalence/montecarlo/`: simulab generates N datasets per configuration (single and two level; LPA over the structures Mplus parameterizes; LCA binary, polytomous, nominal; missing values MCAR/MAR; covariates); batch Mplus and Latent GOLD runs; a test compares log likelihood, parameters, posteriors (row by row, one relabelling), N, criteria, entropy, AvePP, class counts and standard errors | Per configuration, agreement within printed precision on every dataset; results and the generating seeds retained so the test runs offline | M–L |
| E2.4 | **Global-maximum audit** | Across many seeds per dataset, the package’s log likelihood is never below the reference’s beyond tolerance | No dataset where a reference program finds a higher maximum; any exception investigated, not tolerated | S |
| E2.5 | **Recovery against the truth** | On the simulab datasets, estimates and classifications against the generating parameters (bias, RMSE, coverage where standard errors are reported) | Catches errors shared by every program (parameterization, conventions); reported per configuration in `validation/` | M |
| E2.6 | **JStats Mplus runs** | Refit the in-scope runs of the 62 (LPA EEI/VVI, LCA binary, polytomous, nominal, count, covariates, R3STEP/DU3STEP/DE3STEP, TECH14 statistics) with model sections parsed, not names | Each in-scope run compared, or its exclusion stated (factor and censored mixtures are out of scope) | M |
| E2.7 | **Sabotage battery** | `equivalence/sabotage/`: a catalogue of deliberate defects (wrong denominator, swapped labels, wrong rows or groups counted, dropped weights, shifted entropy, wrong parameter count), each injected in turn; the suite must fail on every one | 100% of the catalogue detected; a defect that survives is a missing test | M |
| E2.8 | **Claims registry** | Every statement of agreement in README, NEWS, vignettes and equivalence READMEs mapped to the test that establishes it; a check fails on a claim without one | No unbacked claim (the 2026-10-17 README overclaim is the case it prevents) | S–M |
| E2.9 | **Known gaps with genuine runs** | Houle families against Mplus; sampling weights against `TYPE = MIXTURE COMPLEX` with `WEIGHT`; the bootstrap LRT statistic against TECH14; ordinal and count items against Latent GOLD; `scaling` against Mplus on standardized data | Each gap closed or its limit documented | M each |

Order: E2.1 (where to aim), then E2.2 and E2.7 (tests that can fail),
then E2.3 with E2.4 and E2.5 (the main body), then E2.6, E2.8 and E2.9.

### Status (2026-10-18)

**E2.1** `equivalence/COVERAGE.md`: 82 rows covered, 37 partial, 40 none
at the audit; its gap list drives E2.9.

**E2.2** Two-level LPA and LCA at entropy 0.55 to 0.87 (main kit) and
two-level EEE/VVV (extended kit), group level included.

**E2.3** `equivalence/montecarlo/`: 160 simulated datasets through Mplus
and Latent GOLD; 3,461 + 403 + 49 expectations.

**E2.4** 1,200 of 1,200 fits reach the reference maximum at the true
number of classes (retained `audit-maxima.csv`).

**E2.5** `validation/MONTECARLO_RECOVERY.md`: mean bias at most 0.015,
item probability bias 0.003, logits within one Monte Carlo SE.

**E2.6** `test-jstats-mplus.R`: 42 genuine Mplus runs, 450 expectations;
19 runs out of scope (factor, censored, growth with occasion variances,
constraints, KNOWNCLASS, TRAINING, direct effects).

**E2.7** `equivalence/sabotage/run.R`: 11 injected defects, all
detected, with an identity control that must stay undetected.

**E2.8** `equivalence/CLAIMS.csv` with a test that every claim names a
passing test or is marked unbacked.

\[~\] **E2.9** Closed: `scaling` (grand and person) against Mplus,
weights in a mixture against Mplus WEIGHT with robust SEs,
[`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
against TECH14, families at Latent GOLD’s finest quadrature, ordinal and
count against Latent GOLD (now a test), `drop_empty = FALSE` against
Latent GOLD, criteria/AvePP
matrix/OCC/[`predict()`](https://rdrr.io/r/stats/predict.html),
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) with
missing values and with shared variances. Open (each M): Houle families
against Mplus;
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) order 2,
mover-stayer, `occasions = "grid"`, full covariance, ordinal and count;
[`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
negative binomial, ordinal and two-level beyond one class;
[`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
and [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md)
robust SEs and group level; `parameter_inference(method = "bootstrap")`;
[`pool_imputations()`](https://pak.dynasite.org/latents/reference/pool_imputations.md);
[`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md);
[`sensitivity()`](https://pak.dynasite.org/latents/reference/sensitivity.md).

### Found and fixed by Equivalence 2.0

- [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) general
  engine: the quasi-Newton finish stopped early near a boundary (factr
  tied to `tol` since 0.9.13); LMest occasion model off by 1e-4 in
  transitions. Fixed (factr = 10, as in 0.9.1).
- [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md):
  BFGS from zero settled on a separation ridge (Mplus ex7.3);
  `latents_no_converge`. Fixed (second start from the proportional
  estimate).
- Comparisons kept as scripts could not fail the suite; the Latent GOLD
  kit compared against saved targets. Fixed (tests, live targets).
- `enumerate_classes(cores = )` on macOS with Accelerate BLAS: forked
  workers segfaulted once the session had used BLAS, and a worker killed
  by a signal came back as an empty candidate. Fixed (socket cluster
  where forking is unsafe; lost workers raise
  `latents_parallel_failure`).
- The documentation claimed plain EM follows Mplus’s path from the same
  start; untested, so narrowed to mclust’s, which is tested.

### New item from Equivalence 2.0

**Screened random starts.** On over-extracted models (4 profiles fitted
to 3-profile data) the random starts, converged k-means partitions,
explore the likelihood less widely than Mplus’s screened random starts:
Mplus reaches -3106.613 (twice) and a degenerate -3105.390, this package
-3108.709 from 30 to 300 starts and -3106.051 from 1000. A two-stage
search (many short starts, the best run to convergence, as Mplus’s
STARTS = n m) would close it. Changes results for existing seeds, so it
needs an opt-in or a deliberate release. **M**

Rules for every item: the reference output is genuine and retained with
its inputs and version; the comparison reads the package’s reported
values, never a formula written in the test; class alignment uses one
relabelling chosen from the posteriors and reused for every class-level
quantity; ids are checked unique and complete; and each new test is
shown to fail on a deliberate defect before it is accepted.

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
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md)
(2026-10-17, 0.9.37): every structure and `covariance_model = "full"`,
by chart scores (analytic Sigma gradients through the coordinates’
Jacobian); scores equal the likelihood’s numerical gradient to ~2e-7.

The same for the homogeneous
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md)
(multilpa_transitions), which still refuses the ten constrained
structures: its table reports natural probabilities, so it needs the
delta method from the general engine’s logits. **S-M**

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

**G2. Inference.** (2026-10-17, 0.9.37:
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
method, Wald and person bootstrap with replicates matched by posterior
membership; coverage in `validation/GROWTH_SIMULATION.md`.)
Observed-information and sandwich (weighted) standard errors on the
estimation scale, delta method to variances/correlations; boundary
handling for random-effect variances; person bootstrap via the existing
machinery. Evidence: numerical-Hessian agreement 1e-6; coverage
simulation (\>= 200 datasets, gates as in
`validation/LTA_SIMULATION.md`). **M** (about 1 day)

**G3. Integration.** (2026-10-17, 0.9.37:
[`predict()`](https://rdrr.io/r/stats/predict.html) for new persons,
matching lcmm’s predictY/predictClass/predictRE at hlme’s estimates;
[`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md) for
trajectory models;
[`sensitivity()`](https://pak.dynasite.org/latents/reference/sensitivity.md).)
`print`/`summary`/`plot` (class trajectories with individual-curve
spread), [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
[`predict()`](https://rdrr.io/r/stats/predict.html) (marginal and
conditional on the random effects, i.e. empirical Bayes), posterior
class table, `enumerate` over classes,
[`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
(needs simulation from the fitted model), label alignment, sampling
weights,
[`sensitivity()`](https://pak.dynasite.org/latents/reference/sensitivity.md).
Evidence: tests per verb; bootstrap LRT size check. **M** (about 1 day)

**G4. Validation report.** (2026-10-17:
`validation/GROWTH_SIMULATION.md`, 4 cells x 1000; separated designs
48/54 pass, misses bias-only on variance components; weak separation
under-covers (0.885-0.943); R3STEP coverage 0.944-0.954. Latent GOLD
cross-check not run: lcmm::hlme is exact and already matched.) Recovery
and coverage study (balanced and unbalanced designs, weak and strong
separation, random-effect variance near zero); Latent GOLD cross-check
on two Gaussian configurations; write-up in `validation/`. **S-M** (half
a day, mostly compute)

**G5. Correlated residuals** (2026-10-17, 0.9.37:
`residual_correlation = "ar"`/`"brownian"`, `residual_time`; dense
per-pattern engine, L-BFGS-B within variance floors, Newton finish;
hlme(cor = AR/BM) likelihood equal at its estimates to 1e-8.) Coverage
study of the process parameters not done.

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
growth mixtures (2026-10-17, 0.9.37; profiles matched on the
measurement, classes by posterior membership). Covariate fits still
refused. **M**

Transition-aware
[`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md)
(2026-10-17, 0.9.37):
`r3step(initial_covariates =, transition_covariates =)`, the three-step
latent Markov model; Latent GOLD `step3 modal ml` agrees to its printed
precision on three designs. Still open: a transition-aware
[`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
(distal outcomes) and BCH, and group classes in the step-three chain.
**M**

Error bars in `plot(what = "bars")`, `"profiles"` and `"responses"` for
transition fits (2026-10-17, 0.9.37; `intervals =`).

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
[`sandwich::vcovCL()`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html);
one helper, `.inference_cluster_factor()`, and a required `n_units` on
`.inference_sandwich()` so no robust path can skip it.

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
