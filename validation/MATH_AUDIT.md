# Mathematical and numerical audit — 2026-09-19

The review covered every file in `R/`, the public result interfaces and
documentation, existing test coverage, and the retained external validation
framework. The original test suite passed; additional independent checks found
defects in cases that suite did not exercise.

## Corrections

| Area | Defect and correction | Independent evidence |
|---|---|---|
| Nested EM | Very large negative Gaussian log densities erased mixing priors and categorical evidence, yielding posterior rows whose sums exceeded one. Remove common density and group-evidence offsets before adding priors, then normalize shifted exponentials. | Identical measurement profiles must retain their prior class probabilities, including observations of magnitude `1e8`; exhaustive summation of all latent assignments for mixed full-covariance, incomplete data. |
| Covariate EM | Converting a tiny prior probability to a log discarded finite log probabilities after underflow. Work directly with log-softmax probabilities and center evidence before mixing. | A component with log prior `-1000` must still dominate a sufficiently unlikely alternative. |
| Measurement inference | Full covariance helpers counted categorical indicators as Gaussian dimensions, producing incorrect coefficients and failed information calculations. Zero-dimensional Gaussian blocks now remain empty. | Closed-form normal and binomial standard errors, plus a full covariance matrix calculated from independent cluster influence contributions. |
| Covariate inference | Empty membership blocks failed, changes of measurement or predictor units could produce false singularity, and altered data/nonregular fits were accepted. Preserve empty shapes, standardize curvature calculations, and validate fitting data and regularity. | Normal-model standard errors, cluster sandwich identities, and invariance under changes of units. |
| Inference boundaries | Constrained categorical response estimates require a boundary check; natural variance and probability parameters should not receive regular Wald tests against zero. | Bound-active categorical fits are rejected; interior estimates retain independently verified standard errors. |
| Information curvature | An entirely negative Hessian spectrum passed a minimum/maximum eigenvalue-ratio check, allowing invalid zero standard errors. Require a strictly positive smallest eigenvalue. | A concave quadratic must be rejected as negative observed information. |
| Residual diagnostics | Pooled Gaussian residuals included association explained by profile means; categorical profile diagnostics mixed profiles again. Missing pairs also had incorrect effective sample sizes. | Constructed independent within-profile distributions, exact categorical tables, and pairwise complete-case calculations. |
| Information criteria | The criteria table and enumeration used total rows while the fitted BIC used rows with observed indicators. | One-normal-population likelihood and BIC calculated directly, compared across result interfaces. |
| Random intercepts | A single indicator with multiple profiles produced transposed initialization arrays; generic diagnostics incorrectly required discrete group classes. | Univariate mixture fitting, independent integration, and diagnostics that retain ordinary criteria but leave undefined group classification criteria unavailable. |
| Model comparisons | Bootstrap comparisons lost categorical probability bounds and admitted incomplete or relabeled categorical data; LMR comparison checks did not establish data identity. | Regression checks for fitting constraints, changed categories, missing responses, and changed data/group layouts. |
| Results and plots | One-indicator warm starts, unused factor levels, distinct numeric group IDs with identical printed labels, and degenerate plot dimensions could fail or report incorrect output. | Analytic mixture likelihood, exact sequence counts and graphics-call coordinates. |

The new regression tests are in `tests/testthat/test-math-*-audit.R`. Plot
examples are rendered separately in `tmp/math-plot-audit.html`.

## Validation

- Full package suite: **1,381 assertions passed**, zero failures, warnings,
  skips or errors. The five new audit files contribute **174 assertions**.
- An additional independent review reproduced the updated full-covariance mixed
  FIML likelihood within `1.78e-15` and posteriors within `1.11e-16`.
- `R CMD build` and `R CMD check --no-manual --no-build-vignettes` completed
  with **Status: OK**. Installation, examples, documentation and packaged tests
  passed. Repository-index access emitted network messages; installed
  dependencies were checked successfully.
- All **974 external comparisons agreed within their specified tolerances**:
  glca, mclust, Mplus two-level categorical output, multilevLCA, poLCA and
  published reference results. Comparisons were rerun from the existing harness
  and written to `tmp/math-audit-equivalence.csv` and
  `tmp/math-audit-equivalence.log`, preserving the retained report. Printed
  reference precision and machine-computed differences remain distinct metrics.

## Interpretation limits

These checks establish correctness for the tested models and numerical cases;
they do not prove global maximum-likelihood optimization for every dataset.
Mixtures can have local maxima, unidentified parameters and active constraints.
Inspect convergence, replicated starts and boundary diagnostics.

The bivariate diagnostic p-values remain approximate and unadjusted. The LMR
function still withholds a p-value because its mixture reference distribution
is not implemented. The bootstrap is the package's own parametric procedure;
it is not a claim of Mplus TECH14 equivalence. The [Mplus output
documentation](https://www.statmodel.com/HTML_UG/chapter18V8.htm) also distinguishes
TECH14's simulated likelihood-ratio distribution from ordinary chi-square
testing and specifies narrower supported model families.

Random-intercept mixtures still require a quadrature-sensitivity check. External
Mplus evidence for that family covers the one-profile limit; no new Mplus or
Latent GOLD executable runs were performed in this audit.

## Comparison with Houle et al. (2026)

The subsequent [article and supplement comparison](../equivalence/houle-2026/REPORT.md)
identifies the core package model with the article's dispersion-heterogeneity
family, distinguishes joint estimation from fixed-measurement stages, and
checks all 24 published fit-table rows plus a synthetic exhaustive likelihood.
It does not claim reproduction of the published fitted models: the supplied
materials lack participant data and contain illustrative syntax inconsistencies.
