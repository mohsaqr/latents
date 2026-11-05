# Independent review of inference repairs

Reviewed the current 1.1.17 working source and `../inference/REPORT.md`, covering `R/three-step.R`, `R/lta-three-step.R`, `R/inference.R`, `R/covariate-inference.R`, `R/bootstrap-inference.R`, `R/inference-service.R` and the categorical contrast helper.

No high-impact defect was found in this independent pass. No production source or existing test was modified. The checks here supplement the inference agent's tests and root's installed package suite; they do not repeat either suite.

## Source findings

- The estimated-step-one refusal checks staged fits and any held measurement block before constructing an incomplete covariance. The fixed conditional path remains available. Modal output on staged fits is conservatively refused under `step_one = "estimated"` too, as its new test explicitly specifies.
- Staged bootstrap provenance compares the multiset of complete indicator rows, retaining duplicates, missing-value positions and category labels. Reordering factor declarations or rows does not alone reject equivalent observations. Supplied measurement with weights or transformations is already excluded by `fit_staged()`, so ignoring those attributes in this provenance predicate does not admit those incompatible estimators.
- OPG covariance is retained while transform sensitivity uses observed information in both profile and covariate paths. The weighted population mean/SD nuisance equations use the same row weights and denominator as the fitting transform.
- The R3step affine mapping follows `X = Q R`, with natural coefficients obtained by triangular inversion. The group-class transformation applies the intercept shift separately within each group class and non-reference logit. The boundary check uses the same coordinates. The transition guard checks both initial and outgoing blocks before Wald covariance construction.
- Binary categorical multiplicity adjustment counts one complementary pair per class pair; larger category tables retain the explicitly documented family across category differences. Invalid tests remain unavailable. The changed profile/covariate/three-step statistic paths guard zero standard errors.

## Fresh independent checks

`inference-review.R` uses one core and synthetic data (seed 251105); output is retained in `inference-review.log`.

1. Three correlated covariates with independent affine shifts as large as `1e8`, plus rescaling and a sign reversal: the orthonormal Gram errors are `8.04e-16` and `4.02e-9`. Fitted probabilities at corresponding coefficient vectors differ by at most `9.51e-10`. Both ordinary and transformed interior checks accept the regular coefficients.
2. Three group classes and three profile classes, with the same affine covariate transformation: the independently transformed natural/search coordinates preserve the group-model objective within `2.03e-9`. This checks the group-class intercept mapping beyond the new single-covariate fixture.
3. Weighted FIML grand centering and population scaling on 270 rows in 45 groups, including zero-weight groups and different missingness patterns in two indicators: all four summed nuisance equations are below `1.45e-13` in absolute value. Replaying the stored grand transform reproduces the fitted indicator matrix within `4.45e-16`.

## Narrow scope caution

The shared `.inference_wald()` helper at `R/inference-service.R:357` still computes `estimate / std_error` directly. A direct call with estimate 1 and SE 0 returns `Inf` and p 0; regression and growth methods call that helper. This was not changed by the repairs, and no live regression/growth counterexample with zero estimated SE was established here. The inference repair report's zero-variance claim should therefore remain scoped to the corrected profile/covariate/three-step paths rather than asserting a package-wide guarantee. The historical NEWS wording already scopes its analogous claim to `three_step()` and `r3step()`.

No external fitting, new coverage study or full package check was run in this slice.
