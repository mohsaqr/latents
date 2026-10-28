# Pair preparation and weighted trees: fixes and repeated review

Date: 2026-10-20.

## Corrections

- F02: Generated categorical joint names cannot overwrite existing raw input columns, including other indicators, identifiers or covariates. A generated factor carries its source-item provenance; re-preparing it regenerates from the original columns, preserving idempotence while recognizing genuine name conflicts. Serialization retains provenance. Prepared joint-only model frames remain supported for refitting and prediction. An unmarked joint column alongside only one source item cannot masquerade as a prepared pair.
- F07: Extracted the existing native categorical validation into `.multilpa_categorical_levels()`. The encoder and pair preparer share the same dimension, type, finiteness, distinct printed label and category-order checks. Fitting retains the ordinary two-observed-category requirement for each source item. Prediction preparation permits constant, single-row and entirely missing subsets. Unsupported column types receive the categorical condition class explicitly.
- F03: Regression-tree adapters now use row weights for row posteriors, group weights for grouped class posteriors and group classes, and person weights for growth class posteriors. Clustered growth currently refuses sampling weights at its public API; its unweighted trees correctly count persons and clusters separately.

Changed production files: `R/categorical.R`, `R/local-dependence.R`, `R/tree.R`, `R/accessors.R`, `R/predict.R`. Dedicated regressions: `tests/testthat/test-review-pairs-trees.R`.

## Reproduction and verification

- `pre-fix-final-regressions.log`: Against the installed reviewed revision, 17 assertions failed and 17 passed. Includes the original 60 percentage point weighted tree reversal and the same weight-loss defect for grouped, two-level and growth posteriors. This uses the then-final regression suite before four additional identity regressions were added.
- `targeted.log`, `targeted-results.rds`: Final tests of the dedicated regressions plus existing local-dependence, tree and categorical files: **324 passing expectations; no failures, errors, warnings or skips** (58 + 60 + 131 + 75).
- `lint.log`: All three changed production files pass configured lint.
- `adversarial-review.R`, `adversarial-review.log`: Additional checks cover zero sampling weights, unequal group lengths, complete-case omission, row/group/two-level posteriors, six supported native categorical representations, single-row prediction and unweighted clustered growth (18 persons versus 6 clusters). Final checks pass.
- `tmp/weighted-tree-review.html`, source `tmp/weighted-tree-review.Rmd`: Rendered with `rmarkdown::render()`, showing corrected 80%/20% weighted bars and bands alongside numerical fit/tree agreement within 1e-12. User visual review remains for the coordinating agent to request.

## Additional issue found during repeated analysis

The original pair preparer also hid duplicate input indicator names, duplicate categorical declarations and an identifier reused as a paired indicator. Removing paired names through `setdiff()` let those inputs bypass the ordinary identity checks. Four new regressions first failed (`adversarial-pre-fix.log`). The local-dependence entry now performs these checks before joining; all four pass. The independent probe now reports `latents_bad_data` for duplicate original indicators.

## Practical limits

This work verifies the repaired pair-preparation and tree-adapter seams. It does not rerun external Latent GOLD, Windows/Linux checks or the entire package suite. The coordinating agent is running consolidated validation and owns the other review findings.

## Independent cross-review of other repairs

Read the changes to `R/inference.R`, `R/inference-service.R`, `R/lmr.R` and `R/three-step-ml.R`. No reproduced coordinate or score-selection regression was identified in those repairs. `cross-review.R` independently checks a mixed categorical-pair/direct-effect/ordinal/count fit after adding 2026 to its predictor: likelihood difference is zero and direct-effect SEs agree for observed, OPG and robust inference. Its conditioned score and information dimensions agree (60 x 17 and 17 x 17). A nonsingular shear of numerical Vuong coordinates changes its eigenvalues by only 3.04e-12.

This review did reproduce a regression in the pair provenance repair: base `rbind()` strips custom factor attributes. `provenance-pre-fix.log` preserves the failure. Added data-frame provenance as a fallback, accepted only for a factor whose labels still equal the joined source responses. Raw columns with coincidentally correct labels remain collisions. Stale frame metadata does not authorize replacement of an unrelated reassigned column. Six new assertions cover rbind, column selection, serialization and stale metadata; all pass. Final `cross-review.log` confirms rbind/subset/id-renaming provenance. No further production edits are pending in this scope.

The same independent integration probe also fits categorical and continuous local-dependence pairs together with ordinal/count blocks. Raw and prepared posterior predictions agree within 1e-12; refitting the prepared frame retains the covariance pattern and reproduces the likelihood exactly.

## Final preparation-context correction

Independent review by the numerical agent found that ordinary `subset()`, `merge()` and rbind followed by column selection can remove both provenance attributes. Retained the reproducer in `context-pre-fix.log`, then made prediction and known-fit data resolution explicitly request rebuilding existing joint columns from their native source items. The fitted model establishes those pair semantics, so public prediction no longer depends on preservation of arbitrary R attributes. Initial fitting remains strict even when an unmarked caller column happens to contain identical joint labels. Initial fitting also now refuses absent original indicators with `latents_bad_data`, before category-table construction. Joint-only prediction remains supported. Eleven further public assertions pass; the final four-file total is 324 expectations. `context-lint.log` is clean for all three files modified in this final correction. Source is frozen for consolidated package validation.
