# Findings: ML three-step (`three_step(method = "ml")`)

## 1. Step three returned a local maximum from its one start (fixed; one case left)

- Tests: `equivalence/test-features-ml-three-step.R`, "the independent
  reference agrees on the sweep of 120 datasets" and "on the sweep's
  multimodal cases latents reaches the best known maximum"; unit test "the
  step-three model reaches the best of many random starts".
- Repro: `Rscript equivalence/features/ml-three-step/repro-local-maximum.R`
- EM ran from Latent GOLD's proportional start only, and Latent GOLD's
  default stops at the same point. With skewed or heavy-tailed outcomes, or
  four classes, a higher regular maximum existed in 12 of 239 sweep fits,
  0.7 to 203 log-likelihood units higher (sweep_105: -6616.44 against
  -6589.19, class means off by up to 3 SE). Latent GOLD's 100 random start
  sets reach the higher value. On 200 synthetic heavy-tailed seeds (3
  classes) the single start missed by more than 0.5 in 75.
- Fix (2026-10-21, `.multilpa_step_ml_fit()` in `R/three-step-ml.R`): EM
  runs from a fixed start set and keeps the highest non-floor run; a later
  start replaces the first only when it is higher by more than 1e-10
  (relative), so a likelihood with one maximum keeps Latent GOLD's
  solution exactly. The set: the proportional start (first), the modal-class
  means, each class in turn at the mean of the outcome's top or bottom 2%,
  the class means at the outcome's quantiles in every order (cyclic orders
  above five classes) with unit and proportional variances, and 20 random
  quantile / standard-deviation starts drawn from a fixed seed in a local
  stream (the caller's seed neither read nor moved; unit test). A start that
  empties a class or gives a non-finite likelihood is dropped.
- After the fix: 1 of 239 sweep fits misses. sweep_108 (equal variances, 4
  classes, n = 1000) stops 0.037 below the reference; the outlier class
  (effective n about 5) has means 1287 against 1382 (1.7 SE). EM started at
  the reference's estimates stays there, so it is a regular maximum, but EM
  from the reference's own 25 starts and from 100 further random starts
  never reaches it (only the reference's BFGS does). Both failing assertions
  are kept.
- Cost: the 240 sweep fits take 247 s instead of 158 s; the slowest single
  fit 10.3 s instead of 8.4 s.

## 2. A class carried by one unit got a standard error of 0 (fixed)

- Tests: unit tests "a class carried by one unit gets no standard error,
  with a warning" and "a class carried by two units keeps its standard
  error"; the sweep comparison checks sweep_041.
- Repro: `Rscript equivalence/features/ml-three-step/repro-zero-se.R`
- sweep_041 (four classes on 50 units, one unit alone in class 3) reported
  SE 0.000 and a zero-width interval under `variance = "equal"`; the
  model-based SE is 66. The unit's score for its class mean is zero at the
  estimate, so the sandwich is zero.
- Fix: a class whose Kish effective number of units, (sum r)^2 / sum r^2,
  is below 1.5 gets NA standard errors, intervals and pair tests, with a
  `latents_small_classes` warning. A first rule (effective n < 2) also caught
  sweep_095's two-unit class (effective n 1.997), whose sandwich SE (17.3)
  is close to its model SE (17.1); the Kish count keeps it.

## 3. A start that emptied a class stopped the fit (introduced by fix 1, fixed)

- Test: unit test "a start that empties a class is dropped, not carried
  into the fit" (sweep_041's data, rounded).
- With varying variances, the single unit's class has its proportional
  spread at the variance floor; a tail start at that spread underflowed all
  of the class's responsibilities, and EM stopped on "missing value where
  TRUE/FALSE needed" instead of the classed `latents_singular_information`
  refusal. Each start now stops and is dropped on an empty class or a
  non-finite likelihood; if every start fails, the refusal is raised. With
  both guards removed in memory, the unit test fails.
