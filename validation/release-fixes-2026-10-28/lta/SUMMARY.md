# REL-2 — general LTA inference in different measurement units

## Changed

- `R/lta-inference.R`: general LTA now differentiates displacements in standardized parameter coordinates. Inverse column score norms condition the information, robust score cross-product, and OPG inversion. Returned covariance, information, and scores retain the original public/internal coordinate conventions.
- A separate derivative-width cap prevents huge inverse score norms in weakly identified directions from moving the numerical derivative outside local curvature. `.lta_inference_scales()` uses measurement SDs for means and Cholesky entries, existing scales for constrained covariance charts, and unit scale for logits/log parameters. The cap therefore remains meaningful when outcomes change units.
- Existing identified-coordinate masks, inactive covariance rows, boundary checks, public arguments, and covariance names are preserved.
- `tests/testthat/test-release-audit-lta.R`: exact-fit unit transformations test all covariance entries and public transition/initial SEs for observed, robust, and OPG covariance. Transformations include 1e-4, 1e-5, 1e6, and two indicators scaled differently by 1e-5 and 1e6. A separate public weighted LTA refit compares original and tiny-unit data, adjusted likelihood, full covariance, and reported SEs.

## Verification

- Original code: dedicated initial regression test failed with 5 failures, including the false singular-information error (`before-tests.log`). Audit synthetic fitting/inference evidence retained in `before-probe.log`.
- Final dedicated tests: **41 passed**, no failures/warnings/skips. Runtime **4.904 seconds** including fitting (`weak-curvature-probe.log`). No CRAN skip added by this agent; parent decides the check-time policy.
- Relevant six-file run before the derivative-width-cap refinement: **588 expectations passed**, no failures (`relevant-tests.log`). Includes general LTA score checks, empty-origin inference, empty trailing-occasion masks, covariate unit/origin checks, and bootstrap-related tests. Parent's final full suite validates the frozen cap refinement across all these paths again.
- Weak-curvature numerical probe: a Bernoulli initial/transition logit at 16 has inverse score norm 8,886,112. Uncapped standardized steps give relative curvature error 5,832.334; capped steps recover analytic information within **3.82e-12 relative**. This directly verifies the need for different conditioning and derivative scales.
- Lint: no lints in either changed production file or new test file (`lint.log`).

No shared inference-service file, release metadata, handoff, or learning log was edited. No commit or version bump was performed. Other defects remain owned by the parent and other repair agents.

## Independent covariance-chart followup

After the source freeze, `covariance-charts.R` fitted actual general-LTA EEE, VVV, VEE, and EVE models. Their fitted objects were exactly rescaled, including covariance matrices, by 1e-5 and 1e6. All 24 observed/robust/OPG comparisons passed. Maximum covariance-entry error after mapping back and dividing by marginal SE products was 2.65e-7 (VEE robust); unrestricted/shared full covariance EEE/VVV had maximum 5.86e-9. Every `.lta_inference_scales()` vector had exactly `coef()` names and length, all finite positive entries, and no added coordinates. This exercises full Cholesky, constrained shape-Cholesky, and constrained orientation chart branches. Logs and full numerical results are `covariance-charts.log` and `covariance-charts.csv`. No source or test file changed in this followup.

## Step-forwarding test review

Independently reviewed the parent's edit to `test-engine-fixes.R`: recording the exact `.lta_inference()` call sequence tests the caller-step and table-reuse contract without requiring numerical error to be detectably large. An isolated copy of that test passed 9/9 expectations. In-memory mutation checks confirmed that (1) a registered S3 method dropping the caller's step fails the sequence assertion and (2) the coefficient-table helper discarding supplied inference and recomputing defaults fails both the extra-call assertion and covariance-to-table agreement. Mutants existed only in separate R processes. Scripts and logs: `step-contract-mutations.R`, `step-baseline.log`, `step-ignored.log`, `step-recomputed.log`. No production/test source changed during this review.
