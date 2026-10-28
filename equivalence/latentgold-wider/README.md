# Mixture regression and sampling weights against Latent GOLD 6.1

Before this kit, no Latent GOLD run checked `mixture_regression()`, or sampling
weights in any family. Each case below is fitted by both programs: Latent GOLD 6.1
(free academic licence) run in batch under Wine on this machine
(`~/.wine-latentgold`, `lg61.exe <case>.lgs /b /o <case>.lst`) on 2026-10-28.
Its listings are kept in `returned/`, so the comparison runs offline
(`equivalence/test-latentgold-wider.R`, or `Rscript
equivalence/latentgold-wider/compare.R` from the project root).
`returned/INPUTS.md5` fingerprints the data and syntax Latent GOLD was given,
and the test refuses a kit that has changed since.

**Result: 172 quantities compared, 172 agree.** These are every estimate and
standard error Latent GOLD prints, plus each case's log likelihood and
parameter count.

| case | model | standard errors | LL (Latent GOLD) | LL (latents) | Npar | estimates: worst gap | s.e.: worst gap |
|---|---|---|---|---|---|---|---|
| r01 | Gaussian mixture regression, 2 classes, `y ~ x`, `membership = ~z`, class-specific variances, n 800 | model-based (`observed`) | -1107.9381 | -1107.938116 | 8 | 3.8e-05 | 4.9e-05 |
| r02 | r01's data and model | robust (sandwich) | -1107.9381 | -1107.938116 | 8 | 3.8e-05 | 4.9e-05 |
| r03 | Poisson mixture regression, as r01 | model-based | -1326.5220 | -1326.522034 | 6 | 2.7e-05 | 4.8e-05 |
| r04 | logistic mixture regression, 5 binary outcomes per person (`id`, `class_level = "group"`), 300 persons | model-based | -882.7874 | -882.787417 | 6 | 4.9e-05 | 5.0e-05 |
| r05 | Gaussian mixture regression with sampling weights | robust | -1156.9806 | -1156.980591 | 8 | 2.9e-05 | 4.6e-05 |
| w01 | `lpa()`, 3 profiles, 3 indicators, class-specific variances, sampling weights, n 900 | robust | -3765.4152 | -3765.415200 | 20 | 5.0e-05 | 4.8e-05 |
| w02 | `lta()`, 2 states, 4 occasions, 2 indicators, a sampling weight per person, 400 persons | robust | -4321.0221 | -4321.022105 | 11 | 4.9e-05 | 3.5e-05 |
| w03 | two-level `multilpa()`, 2 profiles x 2 group classes, a sampling weight per group, 120 groups | robust | -2813.8222 | -2813.822161 | 11 | 4.5e-05 | 3.7e-05 |

Every gap is within the rounding of Latent GOLD's four printed decimals.

## How the two are put on one scale

Latent GOLD prints reference-coded parameters (`parameters=first`, class 1
the reference): class logits, an intercept and contrasts for each indicator
mean, transition logits against staying, and variances. `compare-functions.R`
computes the same quantities from `coef()` and `vcov()` of the latents fit,
with delta-method standard errors, and matches Latent GOLD's arbitrary class
labels by the ordering that brings the estimates closest. Mixture-regression
variances come from `log_sigma` (variance = exp(2 log_sigma)).

- **Weights.** Latent GOLD's `samplingweight w rescale` scales the weights to
  the number of cases (groups, in w03), as latents does, so the pseudo log
  likelihoods agree. With a sampling weight, Latent GOLD reports model-based
  errors unless asked for `standarderrors=robust`. latents refuses
  model-based errors for a weighted fit, so the weighted cases are compared
  with robust errors.
- **Logistic mixtures.** latents refuses a mixture of logistic regressions
  with one binary outcome per unit, which the covariate design alone barely
  identifies. Latent GOLD fits it. r04 therefore uses the design latents
  supports: several outcomes per person, with the class at the person level
  (Latent GOLD's `caseid`).
- **Convergence.** Mixture regressions are compared at `tol = 1e-15`. At
  `1e-12`, r03's membership intercept, which has a standard error of 0.69,
  stopped 1.8e-4 short of the maximum, where the log likelihood is 4e-8 lower.

## What the comparison found in latents

- **`coef()` of a general `lta()` fit** returned its means centred on each
  indicator's grand mean, the engine's internal coordinates. The profile
  tables were right. Fixed on 2026-10-28: `coef()` now reports the data's
  units. The affected fits are the ones the general engine runs, among them
  every `lta()` fit with sampling weights, transition or initial covariates,
  occasion-varying transitions, or second order.
- **`mixture_regression()`'s stored variances** kept the names of their
  pre-relabelling positions (`fit$params$sigma2`); no reported number used
  the names. Fixed on the same day.

## Files

- `make-kit.R`: simulates the eight datasets and writes `kit/<case>.dat`,
  `kit/<case>.lgs`, and `kit/cases.rds`.
- `run-kit.sh`: runs Latent GOLD on every case under Wine and copies the
  listings to `returned/`.
- `compare-functions.R`, `compare.R`: the comparison; `comparison.csv` holds
  one row per compared quantity.
