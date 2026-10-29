# Growth, latent class, profile, transition and enumeration fits against Latent GOLD 6.1

Before this kit no Latent GOLD run checked growth mixtures, unweighted
single-level `lca()` and `lpa()`, plain or covariate-driven `lta()`, or class
enumeration. Each case is fitted by both programs. Latent GOLD 6.1 (free
academic licence) ran in batch under Wine on this machine on 2026-10-28
(`~/.wine-latentgold`, `lg61.exe <case>.lgs /b /o <case>.lst`). Its listings
and posterior files are kept in `returned/`, so the comparison runs offline:

```
Rscript equivalence/latentgold-coverage/make-kit.R    # data and syntax (needs no Latent GOLD)
sh equivalence/latentgold-coverage/run-kit.sh         # Latent GOLD under Wine
Rscript equivalence/latentgold-coverage/compare.R     # comparison.csv, search.csv, RUN.txt
```

`equivalence/test-latentgold-coverage.R` runs the same comparison under
testthat. `returned/INPUTS.md5` and `OUTPUTS.md5` fingerprint the data, syntax
and listings, and the test refuses a kit that has changed.

**Result (2026-10-29): 440 quantities compared, 434 agree, 6 known differences
in Latent GOLD's standard errors (below), 0 unexplained.** Every estimate and
standard error Latent GOLD prints, each case's log likelihood and parameter
count, and every row's posterior probabilities.

| case | function | model |
|---|---|---|
| g01 | growth | latent class growth, 2 classes, 400 persons x 5 waves, class-specific residual variance (long format) |
| g02 | growth | growth mixture, random intercept and slope with class-specific covariance (`random = ~time`, `random_covariance = "varying"`), written exactly as a factor model |
| g03 | growth | growth mixture, random intercept with one variance for both classes (`random = ~1`, `random_covariance = "equal"`) |
| c01 | lca | 3 classes, 6 binary items, n 800 |
| c02 | lca | 2 classes, 4 three-category nominal items, n 600 |
| c03 | lca | 2 classes, 5 binary items, membership on a covariate (`profile_covariates`) |
| p01, p02, p03 | lpa | VVI 3 profiles; EEI 3 profiles; VVV 2 profiles with covariances |
| t01 | lta | 2 states, 2 Gaussian indicators, 400 persons x 4 occasions |
| t02 | lta | 2 states, 3 binary indicators |
| t03 | lta | 2 states, Gaussian indicators, transitions on a person covariate (`transition_covariates`) |
| r06 | mixture_regression | negative-binomial (NB2) mixture regression, 2 classes, dispersion by class (Latent GOLD `poisson overdispersed`), n 900 |
| r07 | mixture_regression | cumulative-logit ordinal mixture regression, 4 categories, 2 classes with opposite slopes (Latent GOLD `cumlogit`), n 1500 |
| e1k1-e1k4 | enumeration | `lca()` at 1 to 4 classes (c01's data): likelihood, parameters, BIC |
| e2k1-e2k4 | enumeration | `lpa()` VVI at 1 to 4 profiles (p01's data) |

## How the two are put on one scale

`compare-functions.R` reads every parameter table of a listing (one to five
term columns; the width of each row is its tab count less the header's, plus
one), computes Latent GOLD's reference-coded quantities from `coef()` and
`vcov()` of the latents fit with delta-method errors, and matches the
arbitrary class labels by the ordering that brings the estimates closest.
Category logits are against the first category, class logits against Latent
GOLD's first class, transition terms against the first state given the
previous one. Growth: the class means of the intercept and slope factors are
the latents coefficients, their (co)variances the products of the random-effect
Cholesky factors, each wave's residual variance `exp(2 log_sigma)`.

- **Exact growth mixtures.** A Gaussian growth mixture in long format with
  continuous latent variables makes Latent GOLD integrate numerically. The same
  model written as a factor model on the five waves (loadings fixed at 1 and
  0..4, intercepts absent, one residual variance per class shared by the waves
  through a common label) is computed exactly, so g02 and g03 are in that form.
- **Transition terms.** In Latent GOLD `State <- 1 + z | State[-1]` conditions
  only `z` on the previous state. Each term needs its own condition,
  `State <- 1 | State[-1] + z | State[-1]`, for the model latents fits.
- **Enumeration** is compared at Latent GOLD's solution: latents is started
  there (`start =`) and must reach the same stationary point. Beyond the true
  number of classes the likelihood has several maxima, and for VVI unbounded
  spikes at the variance floor, so "the" maximum is not a well-defined target.
  latents' own 50-start search is recorded beside it in `search.csv`.
- **Regressions.** Latent GOLD's `poisson overdispersed` reports the NB2
  dispersion d (variance mu + d mu^2), latents' exp(log_dispersion). Its
  `cumlogit` intercepts model "category c or higher", so each is minus latents'
  threshold below c (the first threshold plus the exponentiated gaps); the
  slopes are the same.
- **Convergence.** Both sides run to 1e-12; the covariate latent class model to
  1e-15, where at 1e-12 one item logit stopped 1e-5 short of the maximum; the
  ordinal regression to 1e-16, where at 1e-14 one threshold stopped 9e-6 short
  (Latent GOLD at a tighter tolerance prints the same values).

## Known differences (`known-differences.csv`)

- **g03 standard errors.** In the factor form a fixed slope is a factor whose
  variance is fixed at zero. For that model Latent GOLD's errors for the
  random-intercept variance and the second class's residual variance differ
  from the observed information by up to 0.8%. latents agrees with an
  independent numDeriv Hessian of the likelihood written from its definition to
  1e-9 (`check-g03-errors.R`), and with Latent GOLD's own long-format fit of the
  same model (residual error 0.0284; random-intercept variance error
  2 x 0.6612 x 0.0270 = 0.0357 from its Cholesky factor). The estimates,
  likelihood and posteriors agree.

## What the comparison found in latents

- **Slow ordinal convergence.** EM moves slowly along the ordinal thresholds of
  a mixture regression: `tol = 1e-14` left a threshold 9e-6 from the maximum
  (the log likelihood 5e-13 lower). The default `tol` (1e-8) leaves such
  estimates further away. Recorded, not changed.
- **A degenerate first design.** With similar slopes in the two classes the
  first r07 data had no interior maximum: latents stopped at a zero threshold
  gap (flagged), Latent GOLD at a class with thresholds near +-1000. The
  replacement design separates the classes by their slopes.

- **Search at an overfit VVI model.** At four profiles on three-profile data
  (e2k4) latents' default 50-start search reached a likelihood of -2489.105
  against Latent GOLD's -2486.604; started from Latent GOLD's solution it
  reaches that solution exactly. With 200 starts latents finds -2481.739, a
  variance-floor spike that Latent GOLD did not report. Recorded, not a
  defect: the start search is frozen, and the likelihood here has no
  well-defined global maximum.
