# Post-submission code correctness review

Baseline: `34f2ae6` (version 0.8.4), the first baseline option retained after the user answered “1”. Reviewed through `1f11861` (0.9.14), plus the local 0.9.15 corrections. Nothing was submitted, pushed or committed.

The audit read the added or modified production code in all 72 changed R files, the related shared helpers, and the four new dataset-generation scripts. It also checked the executable documentation, ran the full non-CRAN suite, compared saved reference cases, and used independent numerical and API probes. The pre-fix suite passed all 5,824 expectations; the defects below demonstrate why those tests alone were insufficient.

## Findings corrected locally

| # | Defect and consequence | Correction and evidence |
| --- | --- | --- |
| 1 | A one-class Bernoulli regression was rejected as an unidentified mixture. | Apply the identification guard only to multiple classes; coefficients, likelihood and SEs match ordinary `glm()`. |
| 2 | Expanded Poisson/NB log densities lost precision. At a Poisson count of 1e15 the log density was wrong by 1.81. | Use the stable R density functions; check large counts and the supported NB dispersion floor. |
| 3 | Robust regression/growth covariance accepted too few or rank-deficient independent score units. One cluster could produce SEs around 1e-16. | Require sufficient independent units and full rank after boundary coordinates are removed. |
| 4 | Multilevel growth `nobs()` and `logLik()` counted persons instead of clusters, so generic BIC disagreed with the fit table. | Report clusters consistently; direct BIC/table equality check. |
| 5 | Extracting enumerated growth models used regression inference; weighted models received observed rather than robust inference. | Resolve inference by model family and weighting policy; verify extracted coefficients and covariance. |
| 6 | A regression class with no modal assignments broke the classification table. | Preserve every class, with NA averages for an empty assigned class. |
| 7 | Trajectory bands omitted formula offsets. | Include the evaluated offset in continuous, count, binary and ordinal trajectories. |
| 8 | Non-Gaussian trajectory estimates were on the response scale but their SEs were on the link scale. | Apply the inverse-link derivative; compare logistic trajectory SEs with `glm()`. |
| 9 | Growth residual cross-products cancelled after large response translations: residual sums 0.07 and 0.09 became zero. | Keep raw design/outcome references and recompute residuals where subtraction is unstable; verify against a dense Gaussian likelihood. A fitted model shifted by 1e8 differs in log likelihood by only 1.39e-7. |
| 10 | An NA covariance for a boundary parameter erased SEs of unrelated transformed parameters. | Use only nonzero Jacobian coordinates for each transformed output. |
| 11 | Explicit alternate time columns broke growth table printing and both trajectory plot types. Random-effect spread could also vary an unrelated predictor. | Pass the chosen column through individual tables and plotting; hold other random-effect predictors at typical values. |
| 12 | Regression BLRT included nonconverged/warned fits or silently dropped failed replicates, sometimes reporting p=1 after every refit failed. | Reject invalid/nonfinite/reversed replicates and withhold the p-value unless all requested replicates validate. Test failed, nonconverged and fully valid refits. Correct the old real-data test that classified failed refits as successes. |
| 13 | Regression model comparisons accepted different trials, weights, cluster layout and retained rows. | Include those inputs in the likelihood/data signature; refuse incompatible fits. |
| 14 | Inf outcomes or nonfinite evaluated predictor/offset/membership matrices could reach the engines. | Reject them explicitly in fitting and prediction, including ordinal outcomes. |
| 15 | Matrix-valued continuous prediction columns could be flattened/recycled. | Require numeric vectors and reject matrix-valued indicators. |
| 16 | Weighted LTA group-class counts ignored sampling weights. | Sum weighted group posteriors; independently check known counts. |
| 17 | Noise profile percentages were renormalized over Gaussian profiles, overstating their shares. | Include noise mass in counts/percentages, including summaries; retain noise as class zero in the counts table. |
| 18 | Noise assignments disappeared from wide sequence tables and crashed sequence plots. | Preserve factor level zero, label noise cells explicitly and verify that every observation is plotted. |
| 19 | Observed growth plot means ignored person sampling weights. | Include the person weights in the posterior-weighted observed means. |
| 20 | NB log-dispersion derivatives lost precision just above the Poisson-series switch. | Use a stable large-size digamma expansion, combining cancelling likelihood terms first. Across 80 independent 256-bit finite-sum comparisons (including means of 1e20), maximum relative error decreases from 2.69e-7/1.38e-7 to 1.06e-13/5.30e-13 for score/curvature. |

Targeted regression tests: `tests/testthat/test-post-submission-review.R`. Probe programs/logs are in ignored `tmp/review-2026-10-05/`. The stable digamma expansion follows [NIST DLMF §5.11](https://dlmf.nist.gov/5.11); the high-precision reference uses finite harmonic sums instead of the implementation expansion.

## Review coverage

“Read” means all new code and the relevant context were inspected, rather than merely searching for selected patterns. New files were read in full; existing files were reviewed through their complete code diffs and dependent helpers. Tests and reference checks cover interactions between families.

| Production file | Status |
| --- | --- |
| `R/acceleration.R` | Read |
| `R/accessors.R` | Read |
| `R/additive-compare.R` | Read |
| `R/additive-inference.R` | Read |
| `R/additive-methods.R` | Read |
| `R/additive.R` | Read |
| `R/alignment-service.R` | Read |
| `R/bivariate-residuals.R` | Read |
| `R/bootstrap-inference.R` | Read |
| `R/categorical.R` | Read |
| `R/centering.R` | Read |
| `R/compare-models.R` | Read |
| `R/conditions.R` | Read |
| `R/covariance-structure.R` | Read |
| `R/covariate-inference.R` | Read |
| `R/covariates.R` | Read |
| `R/cross-level.R` | Read |
| `R/data.R` | Read |
| `R/descriptives.R` | Read |
| `R/diagnostics.R` | Read |
| `R/enumeration.R` | Read |
| `R/extra-indicators.R` | Read |
| `R/fixed-measurement.R` | Read |
| `R/gaussian-moments.R` | Read |
| `R/get-data.R` | Read |
| `R/growth-mixture-engine.R` | Read |
| `R/growth-mixture-inference.R` | Read |
| `R/growth-mixture-methods.R` | Read |
| `R/growth-mixture-plot.R` | Read |
| `R/growth-mixture-simulate.R` | Read |
| `R/growth-mixture.R` | Read |
| `R/inference-service.R` | Read |
| `R/inference.R` | Read |
| `R/kernel-blocks.R` | Read |
| `R/kernel-em.R` | Read |
| `R/kernel-structures.R` | Read |
| `R/latents-table.R` | Read |
| `R/lpa.R` | Read |
| `R/lta-bootstrap.R` | Read |
| `R/lta-compare.R` | Read |
| `R/lta-engine.R` | Read |
| `R/lta-inference.R` | Read |
| `R/lta-methods.R` | Read |
| `R/methods.R` | Read |
| `R/mixture-regression-engine.R` | Read |
| `R/mixture-regression-enumerate.R` | Read |
| `R/mixture-regression-inference.R` | Read |
| `R/mixture-regression-methods.R` | Read |
| `R/mixture-regression-ordinal.R` | Read |
| `R/mixture-regression-predict.R` | Read |
| `R/mixture-regression-starts.R` | Read |
| `R/mixture-regression.R` | Read |
| `R/multilca.R` | Read |
| `R/multilpa.R` | Read |
| `R/noise.R` | Read |
| `R/plot-gg.R` | Read |
| `R/plot-style.R` | Read |
| `R/plot.R` | Read |
| `R/pool-imputations.R` | Read |
| `R/predict.R` | Read |
| `R/prior.R` | Read |
| `R/report.R` | Read |
| `R/robust.R` | Read |
| `R/simulation-service.R` | Read |
| `R/structure-inference.R` | Read |
| `R/three-step.R` | Read |
| `R/tna.R` | Read |
| `R/trajectory-views.R` | Read |
| `R/transition-inference.R` | Read |
| `R/transitions.R` | Read |
| `R/utils.R` | Read |
| `R/weights.R` | Read |

`R/sequences.R` was also reviewed and corrected as a dependent helper. Dataset scripts: `data-raw/growth-schools.R`, `growth-scores.R`, `study-hours.R`, and `srl.R`.

## Verification

- Pre-fix baseline: 908 cases, 5,824 passing expectations, zero failures/errors/warnings/skips.
- Final frozen-source suite: **927 cases, 5,891 passing expectations, zero failures/errors/warnings/skips**. Includes 19 targeted defect cases with 65 expectations. Source hashes match after the run (`SOURCE.csv`); results are in `tests-final.csv` and `tests-final.log`.
- Final `R CMD check --as-cran --no-tests --no-manual`: **zero errors, zero warnings, one NOTE** (maintainer/new submission/archive size), including ordinary/extended examples and all 11 vignette rebuilds. Tests were checked separately in the full non-CRAN run. Log: `R-CMD-check.log`. The command set `_R_CHECK_FORCE_SUGGESTS_=true` and `_R_CHECK_SYSTEM_CLOCK_=false`: the external time service was unavailable, while local UTC date matched 2026-10-05. CRAN incoming network checks otherwise ran with network access.
- Lint: zero lints under the repository configuration.
- Pre-fix golden comparison: 260/260 identical.
- Post-fix golden comparison: 248/260 identical, three within tolerance, nine reviewed differences and zero harness errors. No golden files were regenerated.

### Saved-reference differences

| Cases | Explanation |
| --- | --- |
| `profiles_noise_single`, `profiles_noise_VVV`, `profiles_prior_noise_EII` | Corrected percentages and the added noise row in counts; fitted likelihoods/parameters are unchanged. |
| `growth_multilevel_one_class` | Independent-unit count changes from 600 persons to 40 schools. |
| `reg_enumerate_bootstrap` | Failed refits now emit a validation warning and withhold the p-value. |
| `profiles_count_nb_varying`, `lta_general_count_nb_varying` | Stable NB derivatives change optimizer/inference rounding. Likelihood differences are below 9e-9; the largest table estimate/SE difference is 0.00037. |
| `reg_negative_binomial` | Stable density/derivatives change near-Poisson estimates by at most 4.51e-6; the expanded old likelihood differs by 1.62e-5. The new density is checked directly against R and derivatives against high-precision recurrence sums. |
| `growth_multilevel` | Stable residual calculations change coefficient rounding by at most 2.96e-6; likelihood difference is 8.10e-11. |

Detailed reference differences: `golden-differences.csv`. The comparison intentionally retains the original saved references; a difference is recorded even when it represents a verified correctness fix.

## Documentation and release boundary

All 18 documents that changed display options now restore them. Eleven vignettes are rebuilt into `inst/doc`; package plotting functions restore `par()` on exit and no shipped working-directory changes were found. This follows the [CRAN Cookbook restoration guidance](https://contributor.r-project.org/cran-cookbook/code_issues.html#change-of-options-graphical-parameters-and-working-directory).

Review evidence is local. The review does not establish correctness for every possible dataset or R platform; randomized fits, boundary identification and numerical conditioning remain governed by the package’s existing checks. No live external statistical-software service was contacted, and existing reference packages/fixtures were exercised through the tests.

## Artifact identity

The final archive matches 334 source/data files (see `ARCHIVE-SOURCE.csv`), with build-added DESCRIPTION fields and whitespace normalized. Its 11 installed R scripts restore options; no generated `Rplots.pdf` is included. All 18 document restoration pairs were independently checked using digits=11 and width=141.

Archive: `tmp/review-2026-10-05/latents_0.9.15.tar.gz`; size 5,484,693 bytes; SHA256 `bb2d2b90ae08c3862690eaeb239abdc6d28fc7f4206cc73c4a17bebea1a6809b`. The earlier archive in `tmp/cran-resubmit-0.9.15/` predates the model corrections.
