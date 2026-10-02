# Transition/LTA review — 2026-10-02

Scope: `R/transitions.R`, `R/transition-inference.R`, `R/lta-engine.R`,
`R/lta-inference.R`, `R/lta-methods.R`, `R/lta-compare.R` and their tests.
No git mutations were performed. Root agent owns the consolidated handoff.

## Confirmed defects fixed

- General LTA coefficient packing previously used diagonal coordinates for
  every covariance structure. It now reuses the core constrained covariance
  charts and log-Cholesky encoders/decoders. `coef()` contains exactly the free
  parameter count for all 14 structures; decode reconstructs all fitted means,
  variances/covariances, and the likelihood. Both invariant Gaussian and
  occasion-specific categorical/ordinal/negative-binomial measurement blocks
  are tested across all 14 structures. EEI/VVI names and order are preserved.

- Second-order occasion-specific transitions previously counted and optimized
  first-order intercepts for occasions 3 onward and second-order intercepts for
  occasion 2. These coordinates never enter the likelihood. First-order and
  second-order free coordinates now use only their applicable occasions in
  names, pack/unpack, scores and parameter counts. Storage retains the original
  rectangular arrays; omitted coordinates decode to zero. Designs' column
  labels survive initialization, EM, and stored-parameter reconstruction.
  With two profiles, four occasions and a single continuous indicator the
  resulting count is 15, rather than 23. A real fitted model now gives finite
  Wald errors and identical coefficient/covariance row names.
- Second-order probability tables now show exact identity transitions for
  stayers. First-order transition tables and boundary checks read only the
  move into occasion 2 when fitting second-order transitions. Stayer coefficient
  labels agree with the `stayers` class label.
- Zero-Gaussian LTA's profiles table now returns a correctly typed empty frame,
  rather than erroring on a scalar occasion column. Public all-table, coefficient,
  refit and simulation surfaces are checked for pure ordinal, count, categorical,
  and mixed non-Gaussian measurement, both invariant and by occasion.
- Homogeneous transition inference checks original category labels and actual
  time values in addition to indicator codes and group ordering.
- `parameter_inference.multilpa_lta()` now includes second-order coefficients
  and their previous-state labels. It validates supplied fitting data, confidence
  level, and differentiation step rather than ignoring or propagating invalid
  controls.
- The quasi-Newton finish enforces existing count-mean and negative-binomial
  dispersion floors. LTA Wald inference refuses count means at zero and
  dispersions at the Poisson limit.
- The general engine now rejects second-order models with fewer than three
  occasions and impossible profile/class counts, including the added stayer.
- Transition nesting now checks categorical indicator roles, observed/grid
  occasion semantics, missing-data handling, covariance-structure restrictions,
  and ordinal levels. Equal negative-binomial dispersion nests inside varying
  dispersion. Exact `[["count"]]` access prevents R partial matching of a missing
  count field to `count_dispersion`.
- Transition BLRT validates both fits' original indicators, groups, row order,
  actual time values, complete sequence layouts, categorical levels, ordinal
  levels/codes, count data, and initial/transition designs. Previously it refitted
  the null at zero iterations and checked only row count. General fits now retain
  their original `time_values` and `categorical_levels` for these checks.
- Transition BLRT validates replicate/start/iteration/tolerance controls before
  dispatching; positive-integer checks avoid overflow-producing coercion.
- Covariates beginning `occasion_` are initialized as covariates, rather than
  receiving a transition intercept; homogeneous second-order models keep such
  covariates free.

## Independent evidence and coverage

`tests/testthat/test-review-lta.R` adds an exact enumeration of every profile
path. It uses direct `dnorm`, `dpois`, and `dnbinom`, independent softmax and
ordinal probabilities, and explicit first/second-order moves and stayer paths.
It does not call package density, softmax, forward-backward, or score helpers
for the reference. The fitted coordinate decoder is used to expose the point
being evaluated; a separate pack/unpack round-trip check verifies that map.

The reference evaluates nonstationary finite coordinates so an optimizer cannot
hide missing score blocks. Central differences compare **each group's entire
score vector**, including weighted scores, with the independently enumerated
likelihood at tolerance `2e-7`. The likelihood comparison uses `1e-11`.

| Feature | Review evidence |
| --- | --- |
| Three profiles and multiple group classes | Independent exact path sums and per-group derivatives |
| Invariant and occasion measurement | Both represented in mixed-indicator exact checks |
| Homogeneous and occasion transitions | Exact likelihood, gradients, free-coordinate counts |
| Initial and time-varying covariates | Nonzero designs in exact checks |
| Second order | Exact paths; free-coordinate round trip; fitted finite Wald inference |
| Mover/stayer | Exact class/path sums and identity probability table |
| FIML | Missing continuous, ordinal and count cells in exact checks |
| Unequal spans/grid gaps | Missing inside-span occasions and shorter trailing sequences in exact checks |
| Categorical and ordinal | Direct probability calculations, labels/types in data/simulation checks |
| Poisson/NB, equal/varying dispersion | Exact path likelihood/scores and actual fitted NB stationary checks |
| Sampling weights | Each group's finite-difference score scaled by normalized sampling weight |
| All 14 covariance structures | Direct observed Gaussian marginal calculation using determinant and solve, including missing indicators, tolerance `1e-10`; complete coefficient count/covariance/basepoint likelihood reconstruction in invariant Gaussian and occasion mixed categorical/ordinal/NB measurement |
| Simulation/refit/comparison | Existing BLRT and refit tests; new missingness/type/count preservation; both-fit data validation; strengthened nesting |
| Homogeneous numerical/inference path | Existing transition-numerics and transition-inference suites |

A separate synthetic negative-binomial fit reached log likelihood
`-1725.828`, dispersions approximately `0.53016/0.47267`, and maximum absolute
analytic score `9.96e-7`; the likelihood history never decreased. The test repeats
this optimization with equal and varying dispersion.

## Validation commands and results

From the repository root:

```r
pkgload::load_all(".", quiet = TRUE)
Sys.setenv(NOT_CRAN = "true")
testthat::test_file("tests/testthat/test-review-lta.R")
testthat::test_file("tests/testthat/test-lta-general.R")
testthat::test_file("tests/testthat/test-lta-compare.R")
testthat::test_file("tests/testthat/test-transitions-numerics.R")
testthat::test_file("tests/testthat/test-transition-inference.R")
```

The existing suites passed 58, 31, 47 and 25 assertions respectively, with no
failures/warnings/skips. `NOT_CRAN=true` enabled bundled external cases:
`depmixS4` covariate transitions, `LMest` occasion transitions, and Mplus
occasion measurement. These validate stored external results, not newly
rerun third-party estimators.

The focused review passed 225 assertions with zero failures/warnings/skips.
A final rerun after hardening generated-intercept counts and adding two more
regression assertions is recorded below. Consolidated full-package checks are
the root agent's responsibility.

`lintr::lint()` reported no lints for all five modified LTA/inference source files and
`test-review-lta.R`. Expected warnings in deliberately unfitted/boundary fixtures
are filtered by class through `quietly()`; unexpected warnings remain visible.

## Limits and next-agent checks

- General LTA Wald inference remains explicitly supported only for diagonal
  EEI/VVI measurement; other structures are refused with
  `latents_unsupported_inference`. Full-covariance simulation is also an explicit
  existing refusal. The review validates their estimation densities, not
  unimplemented inference/simulation.
- Full package tests, R CMD check, and regenerated documentation are the root
  agent's release validation responsibility.
- Preserve design dimnames if changing the rectangular transition storage:
  second-order occasion coefficient packing relies on those names.
- Freshly fitted general objects retain time/category metadata. Saved historical
  objects without those fields cannot prove exact fitting-data identity and
  should be refitted before BLRT/data-dependent inference.
- Covariance orientations use the core local chart anchored at the fitted
  covariance eigenvectors. Zero rotation coordinates represent that anchor;
  fitted coefficient reconstruction is verified. Re-encoding an arbitrarily
  rotated decoded covariance chooses a new local anchor, so arbitrary-point
  coordinate-vector equality is not a global invariant. The original fit must
  accompany its chart coefficients when decoding.

## Final focused result

The final review source was loaded fresh and `test-review-lta.R` passed **325
assertions**, zero failures, warnings, errors or skips. This includes all eight
non-Gaussian public-surface combinations and the final covariance/intercept
metadata and homogeneous-data validation changes. Structured test results were
saved to `/tmp/lta-review-final-results.rds` and console output to
`/tmp/lta-review-complete-final.log`; these temporary files are local validation
artifacts. Root's consolidated validation directory contains the durable release
checks.
