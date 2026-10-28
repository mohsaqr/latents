# Astra core corrections — 2026-10-27

Implemented A02, A03, A06, A07 and A08, plus the confirmed extreme-growth starting-value limitation. Existing uncommitted work was retained; no commit was created. Production sources and the regression file are ready for the parent's clean install and full validation.

## Corrections

- **A02 — LTA covariate charts:** compute the chart from the same active initial/move masks used by the rank checks. Fitting and inference both exclude zero-weight persons and transitions beyond the last observed response. Coefficients are still mapped back to the original input units. Original and padded fixtures now have identical likelihood and coefficients, and matching covariance. A trailing empty occasion with an extreme covariate also leaves the fit unchanged.
- **A03 — constant raw responses:** the Gaussian guard now accounts for the offset and intercept. A constant raw response with a varying offset, or a no-intercept design, remains estimable. A constant adjusted response that the intercept can fit exactly is refused. The variance reference stays positive and relative to the appropriate outcome scale; weighted fits recompute it on carrying rows. Both reported valid examples match independent `lm()` coefficients, ML residual variances and likelihoods, weighted and unweighted.
- **A06 — LTA boundary masks:** boundary checks consider only initial/move probabilities contributing to the observed likelihood. Excluded extreme covariate values no longer disable inference; active zero probabilities still trigger the guard.
- **A07 — categorical support:** ordinary categories are discovered on positive-weight support, so excluded-only values cannot enlarge the parameter count or alter BIC. Excluded rows receive NA measurement codes for these unsupported categories. Their original typed values are saved sparsely outside the parameter dictionary and restored by data accessors, preserving the original-data contract. Training prediction tolerates unsupported categories only on excluded training rows; explicit new-data prediction still refuses an unseen category. Inference, covariate and transition consumers use the same encoding mask. Declared joined-item grids remain full Cartesian grids, as required for previously unobserved joint combinations.
- **A08 — bootstrap covariance:** coefficients marked `estimated=FALSE` now have NA covariance rows and columns as well as NA table standard errors and intervals. Identified diagonal entries equal the reported bootstrap standard errors squared.
- **Extreme growth starts:** all uncorrelated growth starts and optimization use carrying persons only; the returned object then restores expectations for all original persons. This removes absent outcomes from the pooled regression, trajectory clustering and random-covariance initialization together. The reported ±1e100 padding now gives identical likelihood, parameters and start table to deleting those persons, while retaining all 34 persons' posteriors.

## Modified implementation files

`R/lta-engine.R`, `R/lta-inference.R`, `R/lta-bootstrap.R`, `R/mixture-regression.R`, `R/categorical.R`, `R/growth-mixture.R`, `R/multilpa.R`, `R/inference.R`, `R/covariate-inference.R`, `R/covariates.R`, `R/predict.R`, `R/accessors.R`, `R/lta-compare.R`, `R/transition-inference.R`.

The edits in accessors/covariates restore sparse excluded originals; the edits in inference/transition consumers carry the new support mask. No report/get-data code was edited. No exported signature changed and these helper documentation changes are `@noRd`.

## Regression proof

`tests/testthat/test-astra-core.R` contains seven deterministic test blocks:

1. LTA excluded covariates and trailing empty occasions.
2. LTA boundary/inference exclusion, with an active-boundary control.
3. Valid constant-response regression models, both weighted/unweighted, and genuine-degeneracy controls.
4. Ordinary character/factor/numeric category support, covariance, supplied-data checking, raw reconstruction, training/new-data prediction, reconstructed refitting, and a full joined-grid control.
5. Actual person-bootstrap refits with unidentified trailing-occasion coefficients.
6. Extreme zero-weight growth-start contamination.
7. The category-support convention in membership-covariate and general transition models, including inference and data reproduction.

Run against the pre-fix installed library `/private/tmp/astra-review/lib`: **every one of the seven blocks fails or errors**; 11 failed assertions, six errors, five passing guards, zero warnings. The summary reporter truncates its displayed failures at ten, but `before-results.rds` and the printed per-test table contain all seven results. An explicit assertion verified `all(failed > 0 | error)` across the seven rows.

Run against the final fixed tree: **seven blocks, 78 passing assertions, zero failures, errors or warnings**. An explicit assertion verified all seven results pass.

Commands:

```sh
OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NOT_CRAN=true Rscript validation/astra-fixes-2026-10-27/core/run-tests.R /private/tmp/astra-review/lib
OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NOT_CRAN=true Rscript validation/astra-fixes-2026-10-27/core/run-tests.R
```

Logs: `before.log`, `after.log`; serialized results: `before-results.rds`, `after-results.rds`.

## Focused compatibility checks

`run-focused.R` exercised 17 existing/new test files with one BLAS thread: **180 test blocks, 1,149 passing assertions, zero failures, errors or warnings**. This run included the initial six-block/66-assertion version of the new regression file. The final updated regression file was subsequently rerun with all 78 assertions passing; the full focused run was not needlessly repeated. The sparse-original restoration was introduced while that process was running, so the parent's fresh full-suite snapshot is the final whole-tree verification. The final new regression explicitly verifies that restoration in ordinary and covariate models.

Coverage: categorical, missing covariance, local dependence, covariate inference/missingness, prediction, get-data, mixture regression, growth mixture, correlated growth, growth prediction, prior Codex interaction fixes, round-two growth, round-two LTA, open LTA regression fixes and LTA bootstrap. `focused-summary.csv` gives exact per-file results; `focused.log` and per-file RDS files retain the evidence. A requested filename `test-lta-inference.R` does not exist and was not counted; general LTA inference is covered by the round-two/open files and the new tests.

Additional exploratory public-family checks in `category-families.log` reproduced exact deletion-reference LL and covariance for the membership-covariate fit and covariance within 2.8e-17 for general LTA. An exploratory `vcov(lta_fit, data=...)` emitted the documented ignored-argument warning; the permanent regression uses `.lta_check_data()` to verify explicit-data encoding instead.

No broad full suite, network reference campaign or new Monte Carlo calibration was run by this agent. The parent owns clean package checks and final integration validation.
