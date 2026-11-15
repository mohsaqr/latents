# Standard errors for the constrained covariance structures in the homogeneous `lta()`

The homogeneous `lta()` (class `multilpa_transitions`) reports natural
parameters: means, variances or covariances, and initial and transition
probabilities. For the ten covariance structures that tie volume, shape or
orientation across profiles (EII, VII, EVI, VEI, EEV, VEV, EVE, VVE, VEE,
EVV) its observed information is taken in the structure's chart (log volumes,
log shapes, Cayley orientations or a determinant-one log-Cholesky factor,
`R/structure-inference.R`) with initial and transition logits, and carried to
the natural table by the delta method. Until 1.3.4 these fits were refused
(`latents_unsupported_inference`); EEI, VVI, EEE and VVV had standard errors
already.

## Evidence

| Oracle | What it is | Cases | Criterion (fixed before running) | Where |
|---|---|---|---|---|
| Independent likelihood, numerical Hessian | Own forward recursion and Gaussian density, Helmert/Givens/matrix-log charts, first-category logits (`tests/equivalence/helper-lta-structure-reference.R`); central second differences, numDeriv Jacobian | all ten structures (one group class); VEI, EVE, VEE with two group classes; EVI, VVE, EVV with 10% missing values and ragged sequences (FIML); VEV, EVE, VEE with three indicators | SE within 1e-5 relative; log likelihood within 1e-8; the reference gradient at latents' estimate below 1e-3 | `tests/equivalence/test-lta-structure-inference.R` |
| Extended engine | `multilpa_lta` (`R/lta-engine.R`, `R/lta-inference.R`): its own likelihood and analytic chart scores; its covariance of its logits and chart coordinates is carried to the same natural parameters by the delta method (numDeriv Jacobian) | all ten structures, plus one case each of two group classes, FIML and three indicators | SE within 1e-4 relative; same log likelihood (1e-9) and parameter count | `tests/equivalence/test-lta-structure-inference.R` |
| depmixS4 1.5.4 | An independent HMM implementation; EII and VII are its Gaussian responses with linear equality constraints on the standard deviations. Its likelihood at latents' estimates, differentiated by numDeriv (Richardson) | EII, VII | likelihood within 1e-8; SE within 1e-4 relative (its `standardError()` within 5%) | `equivalence/test-lta-structure-depmixs4.R` |
| Simulation coverage | Predeclared registry; EII, VEI, EVE, VEE; 250 sequences of 5 occasions | see `validation/LTA_SIMULATION.md` | availability >= 98%, abs(bias) <= 0.1 SD, SE ratio 0.9-1.1, coverage within the central 99.9% binomial range | `validation/lta-structure-coverage.R` |

The retained fixture `tests/equivalence/fixtures/lta-structure-references.rds`
holds the four datasets, latents' estimates (the point at which the
reference is evaluated), and the reference log likelihoods, gradients,
natural estimates and standard errors, with the R, numDeriv and latents
versions. It is built by `make-references.R` (run from the project root).

depmixS4's own `standardError()` is compared as well, within the 5% that
`equivalence/test-depmixs4-standard-errors.R` allows it (it agrees within
0.1% here). It is not the primary reference because it differentiates by
simple finite differences: on an exploratory dataset (not retained) it was
8.7% high on one EII mean, while numDeriv's Hessian of depmixS4's own
likelihood reproduced latents to 1e-7.

## Pending: commercial programs

Neither Mplus nor Latent GOLD is installed on the machine that prepared this
evidence, so **no commercial reference exists yet** for these standard
errors, and none is claimed.

* **Mplus.** `mplus/export-mplus.R` writes the first three waves of the
  `one_class` dataset (`one_class_w3.dat`, six dependent variables, inside the
  Mplus 9 demo's limit) and inputs for EII, VII (equality labels) and EVI,
  VEI (`MODEL CONSTRAINT`) with homogeneous transitions. They have never been
  run through Mplus and may need syntax repairs. To close the gap: run
  them with the Mplus 9 demo, keep the `.out` files here, fit `lta()` to the
  same three waves with the same structure, and compare the log likelihood,
  parameter count, means, variances and their standard errors (Mplus prints
  three decimals; criterion half a unit of the last printed digit) and the
  transition logits' standard errors (Mplus parameterises the move as
  `[c2#1]` and `c2#1 ON c1#1`; latents' table is natural, so its unconstrained
  covariance, or the extended engine's logits, must be mapped to Mplus's
  logits first). The full-covariance structures need `MODEL CONSTRAINT`s on
  eigen-decompositions that Mplus cannot express directly; they are left to
  Latent GOLD or to the references above.
* **Latent GOLD.** Its Markov models restrict error variances and covariances
  by class-independence and equality, which gives EEI, VVI, EEE and VVV and,
  through equality restrictions across indicators, possibly EII and VII; it
  has no volume-shape-orientation restrictions. No input is prepared.
