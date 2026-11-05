# Regression growth and article repairs

The regression and growth counterexamples are repaired, with 65 focused expectations passing. The affected existing regression and correlated-growth files pass another 88 expectations. Cross-level warning fixtures now explicitly assert the intended warnings and pass 98 expectations. Production changes are complete; no R process from this work remains active.

## Regression and growth behavior

| Change | Result |
|---|---|
| Drop unused membership, group-membership and random-factor levels | Zero design columns no longer increase the parameter count or invalidate otherwise available standard errors. Prediction retains the fitted factor coding. |
| Reject structurally unidentified covariance designs | Common two-wave unrestricted random intercept/slope fits with class-specific residual variances raise `latents_not_identified`. Diagonal covariance, staggered two-wave designs, and identified three-wave controls remain accepted. Multiple proportional class covariances with a shared residual variance retain their valid identification path. |
| Check repeated categorical capacity | Binary and ordinal mixtures without predictors are refused when their exchangeable response-count patterns cannot identify the requested parameters. The check uses positive-weight support. |
| Preserve identifiers that collide with output names | A source identifier named `class` is returned as `id_class`, beside the latent `class`. Cluster identifiers use `cluster_class`; an occupied prefix receives a numeric suffix. Assignment, individual, random-effect, prediction and trajectory readers use the same mapping. Ordinary identifier names keep their existing output names. |
| Validate recovery data | Reordered or modified fitted model columns raise `latents_bad_inference_data` rather than silently changing the recovery cross-tabulation. |
| Correct class-size labels | Two-level regression class counts print as rows; group-level recovery and grouped-class counts retain group labels. |
| Remove intercepts from ratio columns | `exp_estimate` and membership `odds_ratio` are `NA` on intercept rows. A printed ratio column containing no slope ratio is omitted. Raw intercept estimates and their tests remain available. |
| Report correlated-residual starts honestly | `starts` records the final correlated fit and its actual likelihood. `n_best_replicated` is `NA` because the earlier uncorrelated search does not establish replication of this model. Earlier starts are retained as `initialization_starts`, on the original likelihood scale. |
| Handle fixed one-class shares | A share of one with zero uncertainty receives the point interval `[1, 1]`, avoiding a logit-of-one `NaN`. |

The covariance check examines symmetric random-design products and residual identity in standardized coordinates. A row-space basis plus one extra row preserves the relevant rank while bounding each person's covariance block by the number of random effects, rather than allocating a full covariance matrix for a long series. Structural capacity checks identify the demonstrated deficiencies; passing them does not certify every mixture specification as identified.

Files changed for these repairs are `R/mixture-regression.R`, `R/mixture-regression-methods.R`, `R/growth-mixture-engine.R`, `R/growth-mixture-methods.R`, `R/growth-mixture-predict.R`, `R/growth-mixture-plot.R`, `R/growth-mixture-correlated.R`, and `R/growth-mixture.R`. Public documentation was updated in the source roxygen. Root handles generated help, NEWS, the conditions catalogue, README regeneration, and the shared session logs.

## Test evidence

| Validation | Result | Evidence |
|---|---|---|
| Pre-fix counterexamples against the isolated original 1.1.16 source | 42 failed and 22 passed expectations in the earlier 64-expectation version of the new file | [before.log](before.log), `before-results.rds` |
| Final `test-fixes-1105-regression.R` | 65 passed; zero failures, errors, warnings or skips | [focused.log](focused.log), [focused summary](focused-summary.log) |
| New regression tests plus existing `test-growth-correlated.R` and `test-mixture-regression.R` | 153 passed; zero failures, errors, warnings or skips | [affected.log](affected.log), `affected-results.rds` |
| Broader regression, growth, ordinal and related review files | 565 passed; one old snapshot expectation failed and was subsequently corrected and verified in the affected run | [broad.log](broad.log), `broad-results.rds` |
| Cross-level and Gaussian diagnostic fixtures | 98 passed; zero failures, errors, warnings or skips | [cross-warning.log](cross-warning.log), `cross-warning-results.rds` |
| Final installed-source review warning fixtures | 152 passed across `review-1026-multilevel-r2`, `review-families`, and `review-followup`; zero failures, errors, warnings or skips | [final-warning-fixtures-all.log](final-warning-fixtures-all.log), `final-warning-fixtures-all-results.rds` |
| Identifier trajectory preview | Rendered successfully; assertions preserve original IDs and select requested persons | [HTML preview](../../../tmp/regression-identifiers-1105.html), [render log](plot-render.log) |

The snapshot change removes the former membership-intercept “odds ratio” of 0.77. Its displayed raw logit, interval and p-value remain the same. The new numerical tests check parameter counts, likelihood equality, BIC equality, prediction coding, zero-weight support, covariance restrictions, interval endpoints and original identifier values. The expected cross-level criterion warning and the boundary warning in the homogeneous-transition Gaussian fixture are asserted by class.

The focused and affected runs overlap; their totals should not be added as independent coverage counts. All runs used R 4.5.2, `NOT_CRAN=true`, and single-threaded BLAS settings. Root performs the complete final-source suite and package check.

The final warning-fixture run used the already checked installed package at
`/private/tmp/latents-fixes-1105-vtcym3fq/final-check/latents.Rcheck/latents`,
with `_R_CHECK_LIMIT_CORES_=true`. The multilevel round-two helper now allows
the exact `latents_cross_level_criteria` condition alongside its existing
fixture conditions. The families and follow-up tests explicitly expect that
condition at their two restricted-cross-level fits each. These changes
preserve unrelated warnings, touch tests only, and do not reinstall or change
the checked production package.

## Article corrections

Nine source documents were corrected: `README.Rmd`, `vignettes/lca.Rmd`, and the articles `case-esm-lca`, `model-selection`, `class-enumeration`, `trajectory-classes`, `workflow-lpa`, `case-engagement-evaluation`, and `enumeration`.

The ESM article now assigns the 41.6% media-heavy group to student class 1 and the 58.4% low-media group to class 2. A fresh exact fit gives class-1 share `0.415974`, its media-heavy probability `0.565426`, class-2 share `0.584026`, and its low-media probability `0.973949`. Hidden guards now check those model quantities, so a future relabelling cannot silently invalidate the prose. See [ESM check](esm-doc-check.log).

Bootstrap explanations now use the number of valid replicates in the p-value denominator and state that the p-value is withheld below half the requested replicates. A direct check with 15 valid of 19 requested replicates gives `p=0.125`, matching the new explanation. AIC is no longer described as universally having the lightest penalty or decreasing at every grid increment. The AvePP explanation compares probability vectors at the same class count and distinguishes concentrated versus diffuse uncertainty. The LCGA article states that a bootstrap from a misspecified model can also favor extra classes.

Other corrections remove a universal claim of better group-class separation, replace an invalid yes/no profile-score example with the response-probability guidance, identify the residual-test z statistic correctly, remove a hardcoded thin-profile warning count, soften an inaccurate decimal-agreement claim, and distinguish requested candidates from successfully fitted candidates. Root handles the coordinated `drop_empty` explanation.

All nine edited Rmd sources pass extracted-R syntax checks. Five documents completed full execution and HTML rendering:

- [ESM latent class article](../../../tmp/nov5-docs/case-esm-lca.html)
- [Latent class vignette](../../../tmp/nov5-docs/lca.html)
- [Model selection article](../../../tmp/nov5-docs/model-selection.html)
- [Class enumeration article](../../../tmp/nov5-docs/class-enumeration.html)
- [Trajectory classes article](../../../tmp/nov5-docs/trajectory-classes.html)

At root's request, the long workflow render was stopped during its existing 199-replicate bootstrap to free resources for the complete test suite and package check. The remaining three corrections have explicitly labelled prose-only previews, with code, figures and computed outputs omitted:

- [Workflow LPA prose preview](../../../tmp/nov5-docs/workflow-lpa-prose-preview.html)
- [Engagement evaluation prose preview](../../../tmp/nov5-docs/case-engagement-evaluation-prose-preview.html)
- [Enumeration prose preview](../../../tmp/nov5-docs/enumeration-prose-preview.html)

No canonical article iteration counts or computation settings were changed for those previews. Full execution of these three articles remains outside this slice's validation. README execution belongs to root's regeneration step. [docs-check.csv](docs-check.csv) records the exact statuses, [docs-syntax-check.log](docs-syntax-check.log) records clean syntax checks, and [docs-check.log](docs-check.log) records full renders. The first purl pre-pass printed `has_ggplot` lookup messages because it did not execute setup; the subsequent clean syntax pass and the complete trajectory HTML confirmed that this was harness noise rather than an article runtime error.

## Independent audit of root changes

Read-only checks covered weighted grand moments, singleton refusal, standardized references, boundary diagnostics, weighted classification gaps, cross-level prevalences and counts, and categorical simulation names. A freshly shuffled, unbalanced cross-level dataset agrees with raw-data weighted individual prevalences to `1.1e-16`, with the correct positive-weight sample counts. The audit found that weighted population versus unweighted sample standard deviations still broke standardized-reference replication parity; root corrected both paths to population moments. The subsequent independent check confirms integer expansion and all-one/no-weight equality in [standard-reference-after.log](standard-reference-after.log).

No further confirmed scientific or API defect was found in those reviewed changes. No production code, existing tests or canonical article files are still being edited by this agent.
