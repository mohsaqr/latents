# Extra indicator and covariate review handoff

Review completed 2026-10-02 in the existing worktree. Coverage: `R/extra-indicators.R`, `R/covariates.R`, and `R/covariate-inference.R`; root integration changes are reviewed separately.

## Fixed findings

1. **Invalid extra indicator column shapes and ordinal infinity.** Numeric `Inf`/`-Inf` ordinal values were admitted as categories (`Inf == round(Inf)`), and matrix-valued ordinal/count columns could be flattened into observations. Fitting now requires vector columns and finite observed ordinal numbers. Prediction refuses matrix-valued columns for both types. Missing values remain accepted under FIML.
2. **Poisson density cancellation.** The hand-expanded log density `y log(mu) - mu - lgamma(y + 1)` lost substantial precision for large counts. At `y = mu = 1e15`, it returned `-20` instead of R's stable `dpois()` result of approximately `-18.1883`. The E-step now uses `stats::dpois(..., log = TRUE)`.
3. **Negative-binomial score/curvature cancellation at the Poisson limit.** Digamma/trigamma subtraction at size `1e8` produced a maximum score error `4.236887e-7` in a `y = 0:15, mu = 3, alpha = 1e-8` probe, despite the true scores being at most `6.449999e-7`. Stable fourth-order series in alpha now replaces both dispersion derivatives when `alpha * max(y, mu) < 0.001`. At the same probe, the error dropped to `3.660529e-16`. Outside that region, `log1p()` and an algebraically stable reciprocal difference also avoid unnecessary cancellation. Both the analytic score and the NB M-step Hessian consume the same stable derivatives.
4. **Ordinal inference accepted different categories.** Covariate inference compared integer category codes without comparing the categories those codes represented. Multiplying all numeric category values by ten was accepted because ranks remained unchanged. Inference now compares stored category labels too and raises `latents_bad_inference_data` on changed identities.
5. **Single-level covariate inference rejected original data.** `vcov(lpa(..., profile_covariates = ...), original_data)` demanded the internally synthesized `.observation` column. The covariate inference path now rebuilds that identifier for flagged single-level fits. The root agent also changed the covariate dispatch in `R/multilpa.R` to preserve the `single_level` flag; both changes are required.
6. **One-profile NB dispersion standard error lost through inconsistent labels.** Root reproduced `get_results(fit, "count_means")$dispersion_standard_error` as `NA` for a one-profile NB fit with `count_dispersion = "varying"`. `.latents_extra_coordinates()` named any one-row dispersion `shared`, whereas `.latents_extra_labels()` and the count table used `profile_1` for varying dispersion. Root-owned repair: determine the `shared` outcome from the explicit `count_dispersion = "equal"` setting, preserving the profile outcome for varying dispersion even with one profile. Root owns this source correction and the regression in `tests/testthat/test-review-integration.R`; the final combined source suite passed this regression for both dispersion settings (see `validation/review-2026-10-01/source-verification.log`).

The covariate `confint()` roxygen and Rd return descriptions were also corrected: positive parameters use log-scale intervals and response probabilities use logit-scale intervals, matching the existing implementation.

No final change to `R/covariates.R` was needed. A suspected empty-slope dimnames issue with `profile_slopes = "group_class"` and no profile covariates is unreachable through the public API, which correctly refuses that combination; the provisional defensive edit was reverted.

## Numerical review and validation

- `tests/testthat/test-review-extra.R`: **29 assertions passed**, no failures, warnings, or skips. Includes input-shape/infinity regressions; independent finite-sum NB score and curvature references for counts 0..30, means 0.1/3/20, and alpha 1e-8/1e-6/1e-5/1e-4/0.01/0.2; unchanged-category and singleton inference checks; and stable large-count Poisson density.
- New full-model derivative checks include weighted two-level covariate fits with Gaussian full covariance, missing Gaussian/ordinal/count indicators, group-specific profile slopes, group covariates, and both equal/varying NB dispersion. Every analytic likelihood score is checked against central finite differences away from the maximum; natural-coordinate Jacobians are checked against an independent numerical map. Both variants passed.
- `NOT_CRAN=true` existing `test-extra-indicators.R`: **69 assertions passed**, no failures, warnings, or skips. This includes all seven retained Latent GOLD 6.1 fixtures, their likelihood and parameter-count checks, ordinal M-step BFGS comparison, count recovery, missingness, weights/row duplication, simulation, prediction, inference, bootstrap, and LTA wiring. This run preceded the final isolated Poisson-density replacement; root full-suite validation should use the final worktree.
- `NOT_CRAN=true` all seven `covariate*` test files: passed without failures, warnings, or skips, including existing Mplus likelihood/score/inference fixtures and predictor-scale tests.
- `lintr` reported no lints in the two changed R files; `git diff --check` passed. Root should repeat the repository-wide gates after integrating all agents' changes.

Reproduce:

```sh
NOT_CRAN=true Rscript -e 'pkgload::load_all(".", quiet=TRUE); testthat::test_dir("tests/testthat", filter="review-extra|extra-indicators|covariate", reporter="summary", stop_on_failure=TRUE)'
```

Transient probe and test logs are under `tmp/review-extra/` (ignored). Durable regression tests and this record are the handoff artifacts.

## Next-agent double check

- Verify both the `single_level` dispatch flag and normalization, including original explicit singleton-ID validation and changed indicator rows.
- Check the root-owned one-profile NB labeling regression for both dispersion settings: coordinate names, tidy inference outcomes, and the count table's dispersion standard errors must agree.
- Independently differentiate the NB-to-Poisson log-density expansion used in the new stable helper; confirm score and second derivative at the switching threshold remain continuous within numerical precision. The finite-sum regression intentionally shares no digamma/trigamma formula with the implementation.
- Re-run final repository tests/check/lint, including Latent GOLD fixtures with `NOT_CRAN=true`, against the combined changes rather than an earlier worktree snapshot.
- The Poisson-limit dispersion boundary still refuses Wald inference; this review repaired numerical derivatives and did not change that inferential contract. Existing weak-separation/near-boundary information checks remain necessary.
