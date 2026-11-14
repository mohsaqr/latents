# Argument sweep

The comparison kits beside this folder check each function at a few settings
against one reference each. The sweep covers the other axis: every argument of
the fitting and post-fit functions, at each of its options, on several
datasets, each time against a reference that does not use the package's
likelihood code.

```
mini-ci sweep                                        # every module, on the ci VM
Rscript equivalence/argument-sweep/run.R             # the same, from the project root
Rscript equivalence/argument-sweep/run.R lta growth  # modules whose name matches
```

`run.R` writes `tmp/argument-sweep/results.csv` (one row per check: module,
case, quantity, difference, tolerance, status, the conditions the case
signalled, seconds, and the function, reference and kind of reference the
check is evidence for) with its provenance sidecar, and `SUMMARY.md`, and
exits non-zero when a check disagrees or a case stops with an error. The
function and reference come from the verification matrix's
`sweep-functions.csv` and `sweep-references.csv` (a regex on
"<case>: <quantity>" within the module; the first match wins): a check of
another package's EM or search is that package, the likelihood or posteriors
rebuilt from the result tables and the forward/backward recursions are
independent code, a one-class reduction (lm, glm.nb, lme4, closed forms at
one class) is a component check, and invariances, nesting, refusals and plots
are internal. A check no rule names stops the run
(`latents_sweep_unattributed`), so new evidence cannot drop out unattributed.
The CVS report carries the same reference on each check, and the
verification matrix counts the sweep as evidence. Reference packages are installed by
`Rscript equivalence/install-reference-packages.R` (`mini-ci deps`).

## References

| Kind | Where it is used |
|---|---|
| The mixture likelihood rebuilt from the fit's public result tables (`helpers.R`): Gaussian, categorical, ordinal, Poisson and negative-binomial densities over observed entries, mixed within a profile, profiles mixed within group classes (Vermunt 2003), sampling weights | every `multilpa()`-family fit |
| A local-maximum check on that likelihood: no small move of a mean, the common covariance scale, a share, a response probability, a count mean or dispersion, an ordinal location or intercept may raise it | every interior `multilpa()`-family fit |
| The hidden Markov likelihood by the forward recursion, and smoothed posteriors by the backward recursion, from the transition tables: first and second order, several chain classes, movers and stayers, occasion-specific transitions and measurement | `lta()` |
| The mixture-of-regressions likelihood, and the marginal Gaussian likelihood of each person's trajectory with random effects and a residual process, from the regression tables | `mixture_regression()` |
| mclust: its EM continued from the package's posteriors, its one-component fits, parameter counts and own search | the fourteen covariance structures |
| poLCA: its EM continued from the package's probabilities, its search, its latent class regression | `lca()`, class covariates |
| flexmix: its EM continued from the package's posteriors | Poisson indicators; Gaussian, binomial and Poisson regressions, with a concomitant model |
| `lm()`, `glm()`, `MASS::glm.nb()`, `MASS::polr()`, `lme4::lmer()`, `stats::cov.wt()` | one-class limits |
| `stats::optim()` on the step-three likelihoods; `numDeriv` Hessians and Jacobians of the rebuilt likelihood | `three_step(method = "ml")`, `r3step()`, `parameter_inference()` |
| Published formulas recomputed from the posteriors | entropy, AvePP, classification errors, BCH weights, criteria, BIC weights |
| Invariances | row and column order, identifier and time types, seeds, `cores`, `acceleration`, `search`, weights of one, integer weights against repeated units, rescaled indicators and covariates, FIML on complete data, argument spellings, wrappers against `multilpa()` |
| Classed refusals | each argument given an invalid value must raise a `latents_` condition |

`known-differences.csv` registers a disagreement that has an explanation
(module, case, quantity, reason); such a row is reported as `known
difference` and does not fail the run, and an entry that no longer explains
anything fails it. It holds one row: ten k-means starts miss the maximum
poLCA reaches on five binary items with missing values.

A reference that cannot judge a case is recorded with status `no reference`
and its reason, never dropped: a fit on a boundary is a constrained maximum,
mclust's VVE M-step loses likelihood (`diagnostics/vve.R`), and flexmix's
Gaussian M-step is not the maximum likelihood one.

## Modules

| Module | Functions and arguments | Datasets |
|---|---|---|
| `gaussian` | `lpa()`: `model`, `volume`/`shape`/`orientation`, `variance_model`, `covariance_model`, `missing`, `weights`, `centering`, `scaling`, `search`, `n_screen`, `select_start`, `boundary`, `acceleration`, `cores`, `seed`, `start`, `drop_empty`, `prior`, `noise` | `faithful`, `iris`, `srl`, simulated blobs |
| `categorical` | `lca()`: `n_classes`, `missing`, `min_probability`, `local_dependence`, `weights`, input types | poLCA `carcinoma` and `values`, simulated nominal items, `student_esm` |
| `discrete` | `count`, `count_model`, `count_dispersion`, `ordinal`, mixed indicators with FIML and weights | simulated counts and ordered items |
| `twolevel` | `multilpa()`, `multilca()`: `n_group_classes`, both variance and covariance models, `missing`, group `weights`, person `centering`/`scaling`, `time`, identifier types | simulated two-level sets, `course_engagement` |
| `lta` | `lta()`: `transitions`, `order`, `mover_stayer`, `n_group_classes`, `measurement`, `occasions`, `variance_model`, `covariance_model`, `missing`, `weights`, mixed indicators | four simulated panels (balanced, three-state, ragged, binary, mixed), `study_hours` |
| `regression` | `mixture_regression()`: five families, `class_level`, `common`, `membership`, `variance`, `vcov_type`, `missing`, `weights`, predictor types; `enumerate_regressions()` | `study_hours`, flexmix `NPreg`, simulated nested set |
| `growth` | `mixture_regression()` with `random`, `random_covariance`, `random_diagonal`, `residual_correlation`, `residual_time`, `membership` | `growth_scores`, `study_hours` |
| `covariates` | `profile_covariates`, `group_covariates`, `profile_slopes`, `direct_effects`, `start` with `fixed`, `fit_staged()` | simulated two-level set, simulated nominal items |
| `families` | `multilpa(family = )`: `additive`, `dispersion`, `additive_dispersion`, `restricted_cross_level`, `full_cross_level`; `n_group_classes`, `between_variance`, `variance_model`, `weights`, `search`, identifier types | simulated group-intercept set, simulated two-level set, `course_engagement` |
| `pooling` | `pool_imputations()`: `imputed` (list or mice `mids`), number of imputations, `vcov_type`, `level`, `adjust`, single-level and two-level models | simulated blobs, simulated two-level set |
| `plots` | `plot()` with every `what` of every method, each `plot_*()` verb on a fit of every family, an enumeration grid and a pooled result; the arguments of `plot_profiles()`, `plot_diagnostics()`, `plot_responses()`, `plot_enumeration()`, `plot_selection()` | simulated sets, `study_hours`, `growth_scores` |
| `postfit` | classification tables, `three_step()` (methods, contrasts, adjustments, levels), `r3step()`, `parameter_inference()` (three covariance types), `compare_models()`, `enumerate_lpa()`, `bootstrap_lrt()`, `sensitivity()` | `srl`, `faithful`, simulated blobs |

Further references of these three modules: each family's likelihood written
from its definition (a group's ratings exchangeable Gaussian in the group-class
families; Houle, Morin and Harvey's manifest-aggregation likelihood in the
cross-level ones) and evaluated at the result tables, `lme4::lmer()` at one
group class, Rubin's rules written out and `mitools::MIcombine()`. A picture
has no reference package: a view must draw to a device or refuse by class,
and the numbers behind two views are compared with the fit's tables.

Not swept here: the numerical values of `prior` and `noise`
(`test-mclust-prior.R`, `test-mclust-noise.R`); the sweep checks only their
relations to the maximum likelihood fit.

`diagnostics/` keeps the scripts that settled a flagged row: which side holds
the higher likelihood under VVE, how the categorical searches compare, and the
ordinal regression that reached `eigen()` with non-finite information.
