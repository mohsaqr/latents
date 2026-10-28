# Release audit repairs — 2026-10-28

All six correctness defects from the [release audit](../release-audit-2026-10-28/REPORT.md) are repaired in the working tree. Package version remains 0.98.0; no commit, push, publication or CRAN submission was made. The existing uncommitted release-engine test correction was preserved.

## Repairs

| Finding | Change | Regression evidence |
| --- | --- | --- |
| REL-1: missing shared regression dispersions | Extract unnamed scalar variances when packing shared Gaussian/NB parameters, preserving their intended public names. | One/two-class Gaussian and negative-binomial fits, including shared slopes, retain finite dispersion summaries and aligned coefficient/covariance names. |
| REL-2: unit-dependent general LTA covariance | Differentiate and invert in score-scaled coordinates, use the same coordinates for observed/robust/OPG covariance, and map back. Bound derivative widths separately using intrinsic parameter scales so weakly identified coordinates do not cause oversized steps. | Exact fitted-object transformations at small/large/mixed units and independent weighted public refits; 41 assertions. Additional EEE/VVV/VEE/EVE probes pass 24 covariance comparisons. |
| REL-3: zero-weight rows distort regression inference | Run inference on the same positive-weight specification used for fitting, retaining excluded rows in the fitted object's data and posteriors. | Row-level, clustered, group-class and two-level fits retain the covariance from explicit row/group deletion, including membership designs and a zero-weight group that changes group indices. |
| REL-4: additive fit aborts on excluded sentinel groups | Normalize row-max-shifted exponentials directly, independently of the log-likelihood's normalizing calculation. | Adding zero-weight groups at 1e6 preserves likelihood, means, variance components, mixing weights and contributing posteriors; all posterior rows normalize. |
| REL-5: wrappers silently drop comparison models | Expand forwarded dots with a dots-only call signature before reconstructing model order. | Generic/nested wrappers, wrappers adding fits, all three reserved labels, and multiple reserved labels retain the complete comparison and BIC weights. |
| REL-6: noise classification described as trivial | Count actual posterior categories, including noise, when deciding whether a level has anything to classify. | One Gaussian profile plus noise receives entropy, class-size and average-posterior checks; an ordinary one-profile fit remains appropriately trivial. |

Source files: `R/mixture-regression-inference.R`, `R/lta-inference.R`, `R/additive.R`, `R/compare-models.R`, `R/diagnostics-checks.R`. Four new `tests/testthat/test-release-audit-*.R` files add **132 assertions**. Their before-fix runs reproduced the old failures; targeted repaired runs pass. `NEWS.md` describes the corrected behavior.

The broader external tests exposed a pre-existing stale assertion in `tests/equivalence/test-negative-binomial-regression.R`: the package's class-named dispersion was compared to an unnamed `glm.nb` scalar. Its numeric values already agreed. The test now compares numeric values explicitly with the **same 1e-5 tolerance**; the new package tests separately enforce the parameter-name contract. The original failure and the six-assertion passing rerun are retained in `release/nb-reference-*.log`.

The workflow documentation check also exposed an older false monotonicity claim. Original and repaired packages produce **identical** enumeration tables (`release/workflow-grid-before.csv`, `workflow-grid-after.csv`): EEI BIC decreases through six profiles, but VVI reaches its minimum at a bound-active five-profile fit and rises at six. `vignettes/articles/workflow-lpa.Rmd` and its test now describe and verify that actual sequence. No fitting settings, numerical tolerances or plotting code were changed for this correction.

## Numerical review

- The ordinary, unpatched flexmix harness now passes **433 quantities / 29,933 individual values**, with no uncovered output. The audit's workaround is no longer needed.
- Independent inspection of the LTA mapping found no coordinate inconsistency. Reconstructing covariance from returned raw information and scores agrees with the public covariance within 1.3e-15 for observed, sandwich and OPG inference.
- Weak-logit probe: allowing an uncapped score-derived step produced relative curvature error 5,832; the bounded derivative step agrees with analytic curvature within 3.82e-12. Conditioning and derivative widths therefore deliberately use separate scales.
- General-LTA covariance-structure probes at outcome multipliers 1e-5 and 1e6 pass for EEE, VVV, VEE and EVE under all three covariance estimators. Worst error, standardized by marginal SEs, is 2.65e-7.
- Independent mutation checks confirm the corrected step-routing test fails both when the caller's step is ignored and when coefficient tables recompute inference with the default step; the unchanged implementation passes. Evidence is retained in `lta/step-*.log`.
- New regression file runtimes are approximately 1.5 seconds (mixture regression), 1.4 seconds (core), and 4.9 seconds (LTA); these checks remain active on CRAN.

Details: [regression repairs](regression/REPORT.md), [LTA repairs](lta/SUMMARY.md), [core repairs](core/README.md), and `surface/regression-before.log` / `surface/regression-after.log`.

## Final validation

- Full package tests: **13,472 passed, 0 failures/errors/warnings**, across 224 files and 1,961 blocks; 16 fork tests skip because Apple Accelerate cannot safely run in forked workers. The original full run had one stale assertion requiring a coarse derivative to be less accurate by at least 1e-8. The repaired derivatives differed by only 2e-12. `test-engine-fixes.R` now records the actual inference calls to verify the caller's step and table reuse directly, while retaining real covariance calculations and table-to-covariance equality. All 26 assertions in that file pass. `release/unit-summary.log` and `unit-consolidated.csv` replace that file's original result with its rerun without double-counting; `unit.log` preserves the original failure.
- Strict equivalence harness: **24/24 suites ran; 12,468/12,468 quantities agree, covering 3,331,999 values across 37 sources**. No unexplained disagreement, uncovered output or stale known-open entry. Four diagnostic warnings arise in the deliberately difficult `hard-mclust` cases; they are recorded, not suppressed. Worst fraction of the allowed tolerance is 0.994 for computed references and approximately 1 for published rounded references. See [generated equivalence report](equivalence/REPORT.md).

- External numerical-reference tests: **412 passed, 0 failures/errors/skips**, after replacing the two corrected test files with their reruns; six existing example warnings remain (variance boundary, thin profiles and LCA score precision). The original test failures, corrections, reruns and non-double-counted consolidation are retained in `release/reference-summary.log` and `external-tests-consolidated.csv`. No warnings were suppressed to obtain these results.
- Retained Mplus/Latent GOLD/missing-pattern/documentation reference tests: **441 passed, 0 failures/errors/skips/warnings** (`release/references.log`). This includes the four Latent GOLD test files and all 172 quantities in the wider kit.
- Package-wide lint: **0 lints** (`release/lint.log`).
- Source build with vignettes: **passed** (`release/build.log`).
- Built-tarball `R CMD check --as-cran --no-tests --no-manual`: **0 errors, 0 warnings, 1 NOTE** for network-dependent incoming/URL checks unavailable in this environment (`release/check.log`). Examples, donttest examples, documentation and vignette rebuilding pass. Unit tests run separately; PDF-manual compilation is not checked.
- Full installed-package tests run with socket access against a fresh installation at `/private/tmp/latents-release-fixes-20261028/lib`; the runner records its package location and session information.
- `release/source-manifest.json` fingerprints 337 source/test/metadata files at the start of combined validation; `source-manifest-final.json` verifies that production source and the pre-existing user edit remain unchanged. Only the subsequently corrected and rerun `test-engine-fixes.R` differs. The checked tarball predates that test-only correction; its production source is identical to the final working tree. External-test/article corrections have their own manifest and are excluded from the package artifact.

## Remaining release work

The six audited defects are closed with the validation above. The separate lower-priority weighted diagnostic wording and already-known sensitivity best-flag issue were outside the six requested repairs and remain unchanged. The broader external-reference coverage gaps recorded in the audit also remain.

The previously inspected remote CI runs test the old commit, not these working-tree changes. Cross-platform CI needs to run on the repaired candidate before submission. No inference of Windows/Linux success is made from local macOS checks. Versioning, committing/pushing, and submission remain separate actions.
