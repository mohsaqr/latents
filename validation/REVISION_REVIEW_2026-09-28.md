# Revision review — 2026-09-28

Scope: the ten commits after `34f2ae6`, through `bffa930` (latents 0.8.5),
plus the working-tree corrections documented below. The pre-existing untracked
`todo.MD` was left unchanged. This is a review of implementation, mathematical
identities, public results, documentation and reproducibility of the evidence;
it is not a proof of correctness for every data set.

## Findings and corrections

| Priority | Finding | Correction and evidence |
|---|---|---|
| High | Regression predictions rebuilt membership factors from new data. Reversing factor levels reversed fitted responses; a single character-category row failed. Transformed membership terms were also recomputed on the prediction sample. | Retain training terms, factor levels and contrasts for both membership levels. Tests cover reversed levels, one-row predictions, changed global contrast options, unknown/missing categories and a polynomial membership term. |
| High | A one-column group-membership design was assumed to be an intercept, including `~0 + w`. Its M-step replaced the slope with a class-share logit. | Recognize an intercept by its column name as well as its width; support an empty design without Newton iteration. The slope update is checked against an independent binomial GLM. |
| High | `bootstrap_lrt()` accepted posterior-mode fits made with `prior`, then refitted bootstrap samples without that prior. The observed and bootstrap statistics used different estimators. | Refuse either prior-fitted model with `latents_unsupported_prior`, before simulation. Regression tests cover a prior on the null, alternative and both. |
| Medium | `enumerate_lpa()` always forwarded `model`, turning an omitted default into an explicit specification. Documented legacy covariance switches were rejected, as were all-categorical inputs through this wrapper. | Preserve omission when forwarding. Tests check the resulting model and likelihood against `lpa()`, categorical inputs, and continued rejection of explicitly conflicting specifications. |
| Medium | Regression posterior prediction retained training group IDs for an observation-level fit with an `id`, even on a different subset. | Rebuild the group map from prediction rows. Check a reordered two-row subset; reject missing prediction group IDs. |
| Medium | Regression plots lacked membership/grouping columns and could not reconstruct transformed predictors such as `log(x)`. | Build plotting data from retained raw columns, identify numeric predictors from the formula, and retain the other required variables. A grouped model with a transformed predictor and membership factor is exercised. |
| Medium | `r3step(by_group_class = TRUE)` failed when there was only one group class: an empty group-logit block received a nonempty label. | Preserve the empty block. Compare the one-group-class limit with the ordinary pooled regression for both observed and robust SEs. |
| High, evidence | The report harness checked for the obsolete `multilpa` namespace, so it stopped before checking any comparison. Its published-mode experiment repeatedly used the new deterministic first start. | Check the `latents` namespace. Use one deterministic and one random start under each seed, retaining both likelihoods and using plain EM to investigate local basins. Existing comparison thresholds are unchanged. |
| Medium, evidence | Four harness parameter comparisons missed their thresholds despite matching likelihoods. | Tighten fitting tolerances in the mclust and latents fits and refine poLCA's winning solution with poLCA itself. In carcinoma, refinement changes a probability by 0.000108 while gaining only 1.63e-8 in log likelihood. Agreement thresholds are unchanged; before/after tables are retained. |
| Medium, evidence | Regenerating Latent GOLD targets failed when it tried to muffle the single-level message with a warning restart. | Use the message restart. Regenerate targets in a scratch copy; preserve the retained original kit and external output. |
| Medium, evidence | Latent GOLD's retained Step-3 listing uses observed-information SEs, but regenerated targets inherited the package's newer robust default and disagreed. | Request `vcov_type = "observed"` explicitly in that comparison, preserving the external reference, estimator definition and tolerance. |
| High, claims | The imputation documentation treated known-parameter oracle coverage as proof of valid pooling, called a finite importance-resampling approximation exact, and generalized a logistic case-control result to latent outcomes with estimated measurement. | Remove those conclusions. Explain parameter uncertainty, classification uncertainty, the limits of a single simulation design, and natural-scale pooled intervals. Preserve the historical CSV and explain its legacy `exact` label. |
| Medium, claims | FIML bootstrap documentation did not clearly distinguish MAR estimation from the stronger assumptions of a fixed-mask bootstrap. README still said complete indicators were required. | State that a fixed mask does not reproduce general indicator-dependent MAR or MNAR. The retained size experiment covers MCAR. Update source documentation, help and README together. |
| Low | Internal roxygen links targeted deliberately undocumented functions; a formula beginning with `r ` was parsed as inline R. A test expected only the first of several legitimate convergence warnings. | Render internal names as code, fix the formula markup, and capture/assert the convergence-warning class across report panels. Documentation generation and linting are clean. |

## Mathematical checks and commit coverage

| Commit | Main change | Evidence exercised |
|---|---|---|
| `b876c25` | Mixture regression, prior/noise, all-structure inference | Direct mixture likelihoods, analytic scores, flexmix comparisons, mclust prior/noise comparisons, independent covariance charts and numerical Hessians for all 14 structures. Additional prediction/design and group-slope regressions above. |
| `483c8c0` | Vectorized E-step and SQUAREM | Exhaustive likelihood identities, posterior normalization, monotonicity and accelerated/plain-EM comparisons. Acceleration can change attraction basins; identical starting seeds need not visit identical local modes. |
| `6639bd0` | Covariance optimization, initialization, prediction | mclust structures, numerical inference, prediction identities, seed stability. Update the external mode experiment for the deterministic first start. |
| `7c7f41b` | Covariate FIML | Observed-data likelihood and score checks for missing indicators, including full covariance and categorical blocks. Missing covariates remain outside FIML. |
| `d33b572` | Slopes by group class | Numerical score checks, nesting, relabeling and covariate inference. |
| `74af7cc` | Staged bootstrap | Both stages refitted per resample, constraints preserved, label alignment and reproducibility. |
| `8372e3c` | Multiple imputation | Rubin formula agreement with `mice::pool.scalar`, covariance/table consistency and label reparameterization. Formula agreement does not validate an imputation model. |
| `bead246` | Transition standard errors | Analytic score versus numerical differentiation, Hessian-based SEs, depmixS4 and retained transition references. Added a three-profile/two-group-class score check to exercise distinct transition margins. |
| `5841a40` | Extended bootstrap LRT | Fixed missingness mask, conditional covariate simulation, nesting/refit checks and failed-replicate behavior; added prior-fit refusals. |
| `bffa930` | Single-level verbs, grids, result intervals, vignettes | Wrapper equivalence, grid construction, tidy result classes/columns, interval/accessor consistency and executable assertions for vignette prose. |

## Validation results

| Check | Result |
|---|---|
| Full suite, `NOT_CRAN=true` | **3,919 assertions passed** across 663 test blocks; zero failures, errors, warnings or skips. |
| Final R3STEP edge-case follow-up | **56 assertions passed**, including six new one-group-class checks; zero failures or warnings. This reruns the affected controls file after the full suite. |
| External testthat equivalence suite | **503 assertions passed**; zero failures, warnings or skips. |
| Broader report harness | **1,018/1,018 comparisons agree**, all eight suites ran. The final table combines the complete run with the successful 421-comparison poLCA rerun after refining its reference. |
| Latent GOLD, freshly rebuilt package targets | **372/372 comparable quantities agree**; seven explicitly non-comparable, none missing or unparsed. Original input fingerprint verified. |
| Houle published-table arithmetic | **24/24 rounded table rows**, **96/96 package/formula values** agree. This is arithmetic validation, not replication of fitted analyses. |
| Harness self-tests, roxygen, lint, pkgdown | Passed; documentation generation clean, no lints, no pkgdown problems. |
| Source build and full `R CMD check --as-cran --no-manual` | Build and all vignette rendering/rebuilding passed. **Zero errors, zero warnings, two notes**: network/URL checks unavailable and current time unverifiable. |
| Final package follow-up after the R3STEP fix | **Zero errors, zero warnings, the same two environmental notes**. Installed tests: 3,000 passing assertions and 60 intentional CRAN skips. Vignette rebuilding was omitted here because it passed in the full check above. |

The final package was verified byte-for-byte against every current R source
and test file. The full source suite above runs the long tests omitted by
CRAN-mode checks; the affected controls file was rerun after the final fix.

Evidence is retained in [revision-review-2026-09-28/](revision-review-2026-09-28/):
[test results](revision-review-2026-09-28/tests.csv),
[follow-up tests](revision-review-2026-09-28/controls-followup.csv),
[initial harness findings](revision-review-2026-09-28/harness-before.csv),
[final comparisons](revision-review-2026-09-28/harness-after.csv),
[suite status](revision-review-2026-09-28/harness-suites.csv),
[Latent GOLD comparisons](revision-review-2026-09-28/latentgold.csv),
[full package check](revision-review-2026-09-28/R-CMD-check.log),
[final package follow-up](revision-review-2026-09-28/R-CMD-check-followup.log), and the final
[source fingerprint](revision-review-2026-09-28/SOURCE.md5).

Reproduction commands, from the repository root:

```sh
NOT_CRAN=true Rscript -e 'devtools::test(stop_on_failure = TRUE)'
Rscript equivalence/run.R
Rscript equivalence/harness/test-harness.R
Rscript equivalence/harness/run.R
Rscript -e 'devtools::document(); lintr::lint_package(); pkgdown::check_pkgdown()'
R CMD build --no-manual .
R CMD check --as-cran --no-manual latents_0.8.5.tar.gz
```

The final package follow-up additionally used `--no-vignettes`, retaining the
vignettes rendered by the preceding full build.

The report harness writes generated reports. For Latent GOLD, run
`equivalence/latentgold/make-kit.R` and `compare.R` in a scratch copy to preserve
the retained original kit/targets. The external executable is not required to
compare against the retained output, provided its input fingerprint matches.

## Limits that remain relevant to serious review

- Existing Monte Carlo coverage CSVs were inspected with their generating
  scripts; the large simulation studies were not all rerun. They are evidence
  for specified designs, not uniform guarantees. The LRT size study has only
  100 data sets per arm and 39 bootstrap samples per data set; its reported
  rejection rates condition on the 94/99 usable fits. Report failures and
  Monte Carlo uncertainty alongside rates.
- The known-parameter imputation arm is an oracle experiment, not proper MI
  with estimated imputation parameters. Its 0.995 coverage is conservative,
  not evidence of nominal 0.95 coverage. Rubin's scalar algebra is separately
  verified. See the original [mice paper](https://www.jstatsoft.org/article/view/v045i03)
  and the [study of MI under uncongeniality](https://arxiv.org/abs/1911.09980).
- No new Mplus or Latent GOLD executable runs were performed. Agreement is
  with retained external outputs or installed R implementations. Printed
  references and full-precision computations have different tolerances and
  must not be collapsed into a single numerical precision claim.
- Local maxima, near-identification and active constraints remain intrinsic
  to mixture models. Passing these checks does not guarantee the global
  maximum or reliable Wald inference for arbitrary inputs.
- Validation here uses macOS and R 4.5.2. Windows, Linux and other R releases
  still require the configured CI matrix. This review makes no claim that
  those remote jobs ran.
- The transition model remains homogeneous, first-order and measurement
  invariant. Transition-aware three-step methods and general MAR bootstrap
  calibration are not established by the tests above.

No release, commit or push was performed by this review.
