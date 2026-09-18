# Native R method extensions — 2026-09-17, extended 2026-09-18

Version 0.2.0 implements the six capability areas identified after the original
Mplus comparison. Version 0.3.0 adds robust sandwich inference, the full set of
information criteria, classification diagnostics, tidy accessors and the
Lo-Mendell-Rubin adjusted statistic. Runtime dependencies remain base R and stats. These APIs
define specific models; they are not a general reproduction of Mplus.

| Method | Implementation | Independent validation |
|---|---|---|
| Missing indicators | Gaussian observed-data ML, EM conditional moments; diagonal/full, equal/varying residual covariance | Analytic observed-count Gaussian fit, direct optimization, genuine Mplus with missing indicators |
| Residual covariance | Shared/profile-specific full matrices; constrained eigenvalue floor | mclust EEE/VVV fixed points, conditional-normal moments, genuine Mplus |
| SEs and CIs | Observed likelihood Hessian from analytic scores; log-variance/log-Cholesky/logit coordinates; delta-method natural scale | Analytic Gaussian SEs, numerical derivatives, genuine Mplus ML SEs |
| Membership covariates | One-step nested EM with weighted multinomial-logit M-steps at both levels | No-covariate equivalence, binomial regression limit, genuine Mplus |
| Continuous group intercept | Scalar normal intercept, unit loading on all indicators; analytic K=1, quadrature for mixtures | Dense Gaussian algebra, zero-variance limit, independent adaptive integration, Mplus K=1 |
| Class enumeration | Grid fits with diagnostics, entropy and both BIC conventions; optional separate parametric bootstrap LRT | Simulation/refit tests, nesting/data checks and failed-replicate handling |
| Robust sandwich SEs | Per-group analytic scores, cross-product matrix, `A^-1 B A^-1`, MLR scaling correction | Score-decomposition invariant across all seven model families, genuine Mplus `ESTIMATOR = MLR` |
| Information criteria | AIC, BIC, SABIC, CAIC, AWE, ICL under both sample-size conventions | Hand-computed formulas, genuine Mplus AIC/BIC/SABIC |
| Classification diagnostics | Modal and estimated sizes, average posteriors, odds of correct classification, entropy | Internal consistency invariants, degenerate and single-class limits |
| LMR adjusted statistic | `2(L1 - L0)` and the `1 + 1/(df log n)` adjustment | Two genuine Mplus `TECH11` runs, end to end from native fits |
| Categorical indicators | Unrestricted profile-specific response probabilities; binary, ordinal and unordered; mixed with Gaussian blocks by conditional independence | Genuine Mplus `CATEGORICAL =` two-level run, `poLCA` single-level limit, closed-form M-step and bounded-multinomial invariants |

## Added genuine Mplus runs

All successful references are actual Mplus VERSION 9 DEMO (Mac) outputs,
not simulated reference values. Raw specifications, data, output, saved results
and reproduction scripts are retained in the linked directories.

| Matching specification | Largest parameter difference | Largest posterior difference | Absolute likelihood difference |
|---|---:|---:|---:|
| Full varying covariance + missing indicators | 7.11e-8 | 1.60e-7 | 2.58e-5 |
| Full shared covariance + missing indicators | 3.61e-8 | 5.85e-8 | 2.16e-6 |
| Individual and group membership covariates | 3.43e-6 | 5.90e-6 | 3.63e-5 |
| One-profile, one-indicator continuous group intercept | 3.55e-8 | Not compared | 4.91e-6 |

Version 0.3.0 adds three further genuine Mplus comparisons.

| Matching specification | Quantity compared | Largest difference |
|---|---|---:|
| `ESTIMATOR = MLR`, two discrete levels, varying variances | robust standard errors | 1.82e-7 |
| `ESTIMATOR = MLR`, two discrete levels, varying variances | MLR scaling correction factor | 7.03e-7 |
| `ESTIMATOR = MLR`, two discrete levels, varying variances | AIC, BIC, sample-size adjusted BIC | 4.00e-5 |
| `TECH11` on published Example 7.9, 1 versus 2 classes | LMR adjusted statistic | 5.53e-4 |
| `TECH11` on published Example 7.9, 2 versus 3 classes | LMR adjusted statistic | 4.57e-4 |
| `CATEGORICAL =`, two-level, five binary indicators | thresholds | 5.75e-6 |
| `CATEGORICAL =`, two-level, five binary indicators | profile and group probabilities | 1.06e-6 |
| `CATEGORICAL =`, two-level, five binary indicators | log likelihood | 2.75e-6 |

The MLR artifacts are in [robust/](robust/), the TECH11 artifacts in
[lmr/](lmr/), and the categorical artifacts in [categorical/](categorical/);
each has its own `compare.R`. The categorical comparison needs an explicit
`MODEL cw:` block: without it Mplus frees thresholds across both latent class
variables and estimates 23 parameters rather than 13, which is a different
model. Its README explains the check. TECH11 is available only under MLR,
so the robust estimator had to be implemented before the LMR reference could be
generated at all.

The original six comparisons remain in [COMPARISON.md](COMPARISON.md).
The additional specifications are documented in
[extensions-missing-full/README.md](extensions-missing-full/README.md),
[covariates/README.md](covariates/README.md), and
[random-intercept/README.md](random-intercept/README.md).
They add four successful external model comparisons. Mplus saved likelihoods
retain eight significant digits; test tolerances account for this precision.

Observed-information SE comparisons use the genuine Mplus ML standard errors
for the original diagonal models and both new full-covariance/FIML models.
The full/FIML maximum SE difference is 1.18e-7 (varying) and 8.18e-8 (shared);
diagonal SE differences are below 1e-6. Tests explicitly transform the R logits
to Mplus's group intercept/profile intercept/group-class effect parameterization.
ML standard errors are compared against Mplus ML, and robust standard errors
against Mplus MLR; the two are never compared across estimators.

The MLR scaling correction factor is computed as `tr(A^-1 B) / q`. That
definition was confirmed empirically against the value Mplus reports, rather
than assumed; several variants appear in the literature.

## Scope limits

- Missing indicators use an ignorable missingness likelihood. This does not
  estimate a missingness process or establish that a dataset is missing at random.
  Fully unobserved individuals are retained for posterior prediction; they do
  not directly contribute indicator information. The Mplus tests exclude such
  rows to avoid differing row-exclusion policies.
- Full covariance and FIML can be combined with SEs in the discrete-group API.
  Covariate and continuous-intercept APIs currently require complete indicators
  and diagonal residuals; neither currently provides SEs/CIs.
- Membership predictors must be numeric (categorical predictors can be explicitly
  dummy-coded). Group predictors must be constant within group. Profile slopes
  are shared across group classes. There are no direct covariate effects on
  indicator means in this API.
- The continuous intercept does not model random slopes, separate correlated
  intercepts per indicator, or simultaneous discrete group types. Its mixture
  likelihood is checked by independent adaptive integration. A matching Mplus
  two-profile specification was not obtained, and it is not claimed externally
  validated. The unsuccessful specification and a two-indicator numerical
  discrepancy are preserved and explained in its README.
- The native parametric bootstrap is not Mplus TECH14. It preserves observed
  group sizes and simulates complete data under a smaller discrete model.
  Invalid or nonconverged replicates make the p-value NA rather than being
  discarded. Mixture boundaries and local optima remain substantive concerns;
  inspect all start/boundary/replication diagnostics.
- Robust inference covers the discrete-group API only, and treats groups as the
  independent units. It relaxes the within-group distributional assumption, not
  the assumption that groups are independent. Covariate and continuous-intercept
  APIs still provide no standard errors of any kind.
- The Vuong-Lo-Mendell-Rubin reference distribution is **not** implemented. Mplus
  reports a mean, a standard deviation and a p-value beside the LMR statistic;
  none of those is reproduced, and no substitute p-value is returned, because a
  chi-square tail probability is not a valid approximation to it. On the second
  retained TECH11 run the chi-square tail gives 0.153 where Mplus reports 0.582,
  which is why the naive value is withheld rather than shipped with a caveat.
  Use `bootstrap_lrt_ml_lpa()` for a calibrated class-count p-value.
- Categorical indicators are supported for estimation, classification, tidy
  output and plotting, but **not** for standard errors or the parametric
  bootstrap. Both refuse a categorical fit with a classed condition rather than
  returning a number from the Gaussian score functions. Supplying `start` is
  likewise not yet supported for them.
- Count, censored and nominal-with-covariate measurement models are still
  absent; a nominal indicator can be fitted as an unordered categorical one,
  which is the same unrestricted likelihood, but no nominal-specific
  parameterization is reported.
- No unrestricted cross-level measurement heterogeneity, survey weights,
  longitudinal transitions, three-step and BCH classification methods, or
  general Mplus syntax interface.

The official [Mplus mixture examples](https://www.statmodel.com/HTML_UG/chapter10V8.htm)
describe several richer model families. The
[Mplus OUTPUT documentation](https://www.statmodel.com/download/usersguide/Chapter18.pdf)
also limits TECH14 availability for models with multiple categorical latent
variables; the native bootstrap here is deliberately identified separately.

## Reproduction

From the project root, run each directory's `compare.R` to compare retained
Mplus results against fresh R fits. Tests under `tests/testthat/test-mplus-*.R`
use compact offline fixtures, requiring no Mplus installation or network.
`test-inference.R` includes the SE comparisons. `test-random-intercept.R`
includes the continuous-intercept comparison and independent integration tests.

```sh
Rscript -e 'pkgload::load_all("."); testthat::test_dir("tests/testthat")'
R CMD build .
R CMD check --no-manual mllpa_0.3.0.tar.gz
```
