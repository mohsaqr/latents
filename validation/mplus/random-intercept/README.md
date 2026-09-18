# Scalar Gaussian group random intercept

`fit_random_intercept()` fits continuous Gaussian group heterogeneity by
maximum likelihood. One scalar group intercept adds with loading one to every
indicator. Residual variances may vary across profiles or be shared; profile
probabilities are constant across groups. The class is intentionally separate
from `multilpa`: inference methods for a discrete group mixture do not apply.

Mixtures use fixed Gaussian quadrature and compare the fitted-parameter
likelihood with a higher-order rule. A failed check requires more nodes and a
new fit. One-profile models have an exact Gaussian marginal likelihood and
posterior intercept distribution. This is not general random-effects support;
there are no random slopes, missing indicators, covariates, or SEs in this API.

## Genuine Mplus comparison

`univariate.inp/.out/-results.dat` are a successful Mplus 9 Demo ML run on the
first indicator of `one-profile.dat` (300 people in 60 groups). This is a
one-profile, one-indicator random-intercept model. At the independent native R
optimum, largest parameter difference is 3.55e-8 and likelihood difference is
4.91e-6, within the eight-significant-digit Mplus saved-results precision.
`comparison.csv`, `comparison-log.txt`, and the offline test fixture preserve
the comparison. `compare.R` recreates them from genuine Mplus saved results.

The deterministic simulated dataset was generated in R with seed 983:
60 group intercepts `rnorm(60, sd=.65)`, each repeated five times;
`y1 = 1.2 + intercept + rnorm(300, sd=.5)` and
`y2 = -.6 + intercept + rnorm(300, sd=.7)`.

Reproduction, from the package root:

```r
source("validation/mplus/random-intercept/compare.R")
```

To rerun Mplus, use `/Applications/MplusDemo/mpdemo univariate.inp univariate.out`
with the shell working directory set to this directory. The executable is not
needed for package tests.

## Additional checks and unsuccessful external specifications

- Package tests verify dense multivariate Gaussian likelihoods, exact Gaussian
  posterior moments, the zero-random-variance ordinary-mixture limit, independent
  adaptive integration for a genuine two-profile mixture, and synthetic recovery.
- `one-profile.inp` and `one-profile-mixture.inp` fit the two-indicator rank-one
  between covariance in Mplus. Despite specifying zero between residuals, their
  Mplus likelihood differs from the exact rank-one model by 0.00532. Adding
  1e-4 to each indicator-specific between residual variance reproduces Mplus's
  likelihood to 4.1e-6. This is an empirical numerical diagnosis, consistent with
  the minimum-variance setting printed by Mplus, not a verified general claim
  about Mplus internals. Raw successful outputs are retained, but these models
  are not counted as exact parity validations.
- `mixture.inp/.out` retains an unsuccessful attempt to express two within-level
  profiles plus a continuous group intercept. Mplus rejected the class-specific
  within means/variances in that specification. It provides no usable results
  and is not an oracle. `mixture.dat` is a separate 150-person synthetic dataset.
  The mixture branch is validated by independent numerical integration and
  synthetic recovery, not by a matching Mplus mixture run.
