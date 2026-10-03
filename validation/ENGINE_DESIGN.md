# One engine: measurement blocks × latent structures

Plan, 3 October 2026, against latents 0.9.12; implemented the same day (see
**Status** at the end). The public API (every exported verb, argument and
result table) does not change; this is an internal consolidation with new
model combinations as its payoff.

## Why

latents has grown eight estimation engines, each complete in itself:

| Engine | Lines (R/) | What it fits |
|---|---|---|
| `multilpa` (+ structures, extra indicators) | ~2,000 | profiles, LCA, mixed indicators, group classes |
| covariates | ~2,100 | membership covariates |
| `lta` | ~2,900 | general latent transitions |
| transitions | ~2,000 | homogeneous latent transitions |
| additive / cross-level | ~2,900 | Houle families |
| mixture regression | ~3,100 | regression mixtures, trajectory classes |
| growth mixture | ~2,000 | random effects within classes |

Every engine re-implements the same machinery: E-step and M-step (6-8
copies), EM loop (7), parameter packing (4), analytic scores (10), label
alignment and permutation (4-5), simulation (6), plus its own starts,
convergence rules, weight handling, inference and tables.

The duplication is where defects came from. In the 0.9.x reviews alone:
sampling weights were handled differently per engine (classification tables
and cross-level composition ignored them); ordinal relabelling was fixed in
label matching but needed separate rebasing in each place it is used; the
L-BFGS-B line-search stop (code 52) was judged correctly in the growth
engine's finish but falsely reported non-convergence in the transition
engine's; `bootstrap_lrt()` swallowed a refusal in one engine's replicate
loop. Each fix had to be found and made per engine.

The package already proved the alternative once: ordinal, Poisson and
negative-binomial indicators were added as one *extra block* (E-step term,
M-step, scores, coordinates, labels, draws) and then plugged into the
profile, covariate and transition engines. This plan generalizes that.

## Target architecture

A model is **measurement blocks** combined with one **latent structure**,
estimated by one **kernel**, and served by shared **services**.

```text
           measurement blocks                 latent structure
   gaussian profiles (14 structures)     independent classes
   categorical                            membership logit (covariates)
   ordinal / Poisson / negative binomial  two-level group classes
   GLM regression (gaussian/binomial/     Markov chain over occasions
     Poisson)                              (first/second order, mover-stayer)
   gaussian random effects (growth,
     additive group intercepts)
                   \                        /
                    kernel: EM + SQUAREM + quasi-Newton finish,
                    starts, convergence, boundaries, weights, missing data
                                   |
     services: inference (observed/robust/OPG, delta, Wald), bootstrap,
     BLRT, label alignment, simulation, compare_models, three-step/BCH,
     report tables, plots
```

### Likelihood contract

For independent units u (persons or groups) with latent class sequence c,

```text
log L = sum_u w_u log sum_c P_structure(c | u; gamma) prod_b f_b(y_ub | c; theta_b)
```

`P_structure` is a class prior (independent, membership logit, two-level) or
a Markov path probability; the sum over paths is done by forward-backward.
Each block's density `f_b` is evaluated at the level its latent variable
lives: a row for a profile indicator, a person for a growth trajectory
(random effects integrated out exactly), a group for an additive group
intercept. A block declares that level; the structure declares which
latent variable governs each level.

### Block interface

A block is a list of functions over its own parameter list `theta`:

| Function | Contract |
|---|---|
| `level()` | `"row"`, `"person"` or `"group"`. |
| `initialize(data, posterior)` | Starting `theta` from soft or hard class weights. |
| `log_density(theta, data)` | Units x classes matrix of log f_b. Must be finite; `-Inf` only for impossible data. |
| `maximize(theta, data, weights)` | Weighted M-step: the `theta` that maximizes `sum weights * log f_b`, or an ascent step (ECM) with a non-decreasing guarantee. |
| `scores(theta, data, posterior)` | Units x coordinates matrix of d log L_u / d coordinate, posterior-weighted. |
| `pack(theta)` / `unpack(vector, theta)` | Estimation coordinates with stable names; round trip to 1e-12. |
| `natural(theta)` | Reported parameters and their names (for tables and delta-method errors). |
| `permute(theta, order)` | Relabel classes, preserving the likelihood exactly (reference rebasing included). |
| `signature(theta, reference)` | Classes x features matrix for label matching, on a scale-free footing. |
| `simulate(theta, data, classes)` | Draw the block's data given classes, on the fitted design. |
| `boundary(theta)` | Which coordinates sit on a bound, for refusals and warnings. |
| `describe(theta, inference)` | Tidy tables with display plans (`latents_table`). |
| `count()` | Number of free parameters. |

### Structure interface

| Function | Contract |
|---|---|
| `prior(gamma, design)` | Log class priors per unit, or (Markov) initial and transition log probabilities. |
| `posterior(log_prior, log_density)` | Unit posteriors and the log likelihood; forward-backward for Markov structures, including pair posteriors. |
| `maximize(gamma, posterior, design, weights)` | Weighted M-step of the structure (multinomial logit by Newton, Markov counts). |
| `scores`, `pack`/`unpack`, `permute`, `count`, `describe` | As for blocks. |

### Kernel and services

- **Kernel:** EM driver over blocks and structure, SQUAREM acceleration,
  the quasi-Newton finish with the code-52 restart rule, screened multistart,
  one convergence rule, one boundary and degeneracy report, sampling weights
  (multiply posteriors once), missing data (blocks integrate their own
  missing entries).
- **Services:** observed, robust and OPG covariance from packed scores; the
  delta method and Wald tables; nonparametric and parametric bootstrap;
  BLRT; label alignment from block signatures; `simulate()`;
  `compare_models()`; three-step and BCH from posteriors; report-style
  tables; the plot views.

## Phases

Each phase is its own release. Public results must not move: every ported
model reproduces its golden fixture (below) to 1e-10 in the log likelihood
and 1e-8 in estimates and standard errors, the full suite and the external
equivalence harness pass unchanged, `R CMD check` is clean, and no fit
becomes more than 20% slower on the performance fixtures.

| Phase | Work | Acceptance | Size |
|---|---|---|---|
| **0. Golden baseline** | Record, for every engine and its documented options (structures, families, weights, missing data, covariates, orders), the fitted log likelihood, estimates, standard errors (observed and robust), class tables and timing, as `equivalence/engine-golden/`. A comparison script reruns them. | Fixtures reproduce on the current code to machine precision; timing baseline recorded. | M |
| **1. Inference service** | One implementation of covariance from packed scores (observed by Jacobian of analytic scores, robust sandwich with weights, OPG), delta method and Wald tables; engines call it. | Golden standard errors to 1e-8 for every engine; the ten score implementations stay, behind one interface. | M |
| **2. Alignment and simulation services** | Label matching from signatures with exact permutation (ordinal rebasing, Markov rebasing, proportional covariance rescaling) and one simulation entry point. Bootstrap, BLRT, imputation pooling and sensitivity use them. | Permutation invariance tests per block; golden bootstrap results with fixed seeds. | M |
| **3. Kernel, first ports** | Block and structure interfaces; the EM kernel. Port mixture regression (GLM block; independent, membership and two-level structures) and growth mixture (random-effects block). | Golden fixtures; the `lcmm` and `flexmix` harnesses unchanged. | L |
| **4. Profile models** | Gaussian profile block (14 structures), categorical and extra blocks; membership-covariate and two-level structures. `lpa()`, `lca()`, `multilpa()`, `multilca()` and covariate models run on the kernel. | Golden fixtures; Mplus, mclust and Latent GOLD harnesses. | L |
| **5. Transitions** | Markov structure (first and second order, occasion-varying, covariates, mover-stayer); `lta()` and homogeneous transitions on the kernel. | Golden fixtures; depmixS4, Mplus and JStats transition fixtures. | L |
| **6. Houle families** | Additive group intercepts as a group-level random-effects block (the growth random intercept at another level); dispersion and cross-level families. | Golden fixtures; the Houle equivalence report. | M-L |
| **7. New combinations** | Each is a pairing, documented and validated as a model: multilevel growth (random-effects block × two-level structure), growth with categorical or ordinal items, transitions of regression or growth blocks, NB and ordinal outcomes in regression mixtures. | Per model: an independent likelihood, parameter recovery, and an external reference where one exists (Latent GOLD, Mplus). | per model, M |

Sizes as in `todo.MD`: S under half a day, M one to two days, L several days.
Phases 0-2 deliver most of the defect reduction on their own; 3-6 complete
the consolidation; 7 is the payoff and can proceed model by model.

## Order and checkpoints

1. Phase 0 before anything else: without golden fixtures a consolidation
   cannot show that results did not move.
2. Phases 1 and 2 touch only services, so they can ship while every engine
   still estimates as it does today.
3. Ports go newest first (growth, regression), where the code is smallest and
   the external references are exact; the profile engine, the most tuned and
   most used, goes after the kernel has proved itself.
4. After each phase: a review pass on the ported engine (silent failures,
   refusals, weights, labels), then the release.

## Risks

| Risk | Mitigation |
|---|---|
| Behaviour drift (a result moves by more than rounding) | Golden fixtures at 1e-10/1e-8; any intended change is a separate, announced fix with its own test. |
| Speed loss from generic code | Blocks keep their vectorized internals (pattern sharing, closed forms); only orchestration is shared; a 20% budget per fixture. |
| A port stalls half done | Phases ship independently; an unported engine keeps working unchanged. |
| Interface churn | The interfaces are fixed after phase 3; later phases add blocks and structures, not new functions in the interface. |
| Review load | One engine per release; each release's diff is a port plus its fixtures. |

## Decisions for the owner

1. Keep the public API exactly as it is through phase 6 (recommended), with
   new combinations exposed through existing verbs and arguments.
2. Accept golden-fixture tolerance of 1e-10 (log likelihood) and 1e-8
   (estimates, errors), with any larger change treated as a fix.
3. Which new combination (phase 7) comes first; multilevel growth is the
   recommendation.

## Status (3 October 2026)

All seven phases are implemented in the working tree (uncommitted). The
acceptance bar held throughout: after every step the golden baseline
reproduced **bit for bit** (`identical()` fingerprints), not merely within
tolerance, except for the intended changes listed under E0 and E7.

| Phase | Where | What was done |
|---|---|---|
| E0 | `equivalence/engine-golden/` | 251 cases (every engine and documented option, 38 recorded refusals, seeded services) with `record.R` / `compare.R`; full record ~3 min. Defects found while recording were fixed first and their cases re-recorded: noise fits (`summary()`, BCH and error tables), `sensitivity()` on noise fits, the LTA table ignoring `step`, `get_tna()` on general LTA fits, an additive negative-definite information, unclassed refusals. |
| E1 | `R/inference-service.R` | One implementation of score-difference information (central / five-point; plus1 / max1 / absolute steps), positive-definiteness rules, Cholesky inversion, rank refusal, sandwich, delta covariance and standard errors, numeric Jacobian, simplex Jacobian, critical value and Wald columns; the profile family's optimHess, OPG, cross-product, cluster refusal, Wald bounds and p adjustment moved there. Every engine calls it, each with the variant that reproduces its arithmetic. |
| E2 | `R/alignment-service.R`, `R/simulation-service.R` | One matcher, one measurement signature (with a `mixing` column per family), shared group signature, logit rebasing, one align driver; one seed scope (`.latents_local_seed()`), one inverse-CDF draw, shared cluster-bootstrap and BLRT replicate loops (each engine's policy kept as flags). |
| E3 | `R/kernel-em.R`, `R/kernel-structures.R` | One EM driver for all eight engines (decrease guard, condition, non-finite stop, SQUAREM, inner-solve stall), one start runner, one quasi-Newton finish. The regression and growth models are blocks (GLM densities; person densities with random effects integrated out) under shared structures (observation, group, two-level E-step, weighting, mixing M-step, mixing scores, packing, class order). |
| E4 | `R/kernel-blocks.R` | One measurement block for profiles, covariates, cross-level and both transition engines; the profile E-steps read as block + named structure (`.multilpa_nested_structure()`, `.multilpa_logit_structure()`). |
| E5 | `R/kernel-structures.R` | One forward-backward pass for both transition engines (the homogeneous chain is the general one with one matrix repeated). |
| E6 | — | The Houle families run on the shared EM driver, start runner, inference service and (cross-level) measurement block. |
| E7 | growth and regression files | Multilevel growth mixture (growth block x two-level structure, `cluster =`), negative-binomial regression mixtures (`family = "negative_binomial"`) and ordinal regression mixtures (`family = "ordinal"`, `R/mixture-regression-ordinal.R`: cumulative-logit density, a joint Newton M-step over thresholds and slopes, thresholds packed as first threshold + log gaps as in `MASS::polr()`). |

### Decisions taken during the work

- **Bit identity over uniformity.** Where engines differed only in floating
  point (critical value `qnorm((1 + level) / 2)` vs `qnorm(1 - (1 - level) /
  2)`, `J V J'` vs `rowSums((J V) * J)`, `solve` vs `chol2inv`, the
  covariate engine's `log(2 pi v)` vs `log(2 pi) + log(v)`), the shared
  function takes the variant as an argument. Merging the covariate
  normalizer moved its robust errors by up to 1.6e-8 relative, beyond the
  bar, so it stays a named variant.
- **Structures kept distinct where their parametrizations are.** The
  profile models' fixed-probability nesting and the covariate models' logit
  nesting are separate structure functions: one could be written as the
  other, but not with the same arithmetic.
- **The additive families keep their closed-form likelihood** on group
  sufficient statistics rather than becoming the growth random-intercept
  block at the group level: it is exact, faster, and a port would not be bit
  identical.
- **Growth starts.** A start from each person's own least-squares trajectory
  (Ward-grouped, deterministic) was added because the level-only starts of
  the model without random effects trapped growth fits 9-11 log-likelihood
  units below the maximum. It is appended after the existing starts, so
  their random streams are unchanged; two golden cases improved and every
  growth case's starts table gained a row.
- **Multilevel growth with one group class** is allowed: it is the
  single-level model with clusters as the independent units, so the number
  of group classes can be enumerated from one.
- **New models refuse what they do not yet support** with classed
  conditions: weights and three-step with `cluster`.

### Not done

- The other E7 pairings: growth with categorical or ordinal items,
  transitions of regression or growth blocks. Each is now a block x
  structure pairing on the shared kernel.
- Multilevel latent class growth analysis (`cluster` without `random`) and
  sampling weights for multilevel growth.
- External references for the new models beyond `MASS::glm.nb()` and
  `MASS::polr()` (one class): a Latent GOLD or Mplus fit of a multilevel
  growth mixture or of an ordinal regression mixture.
