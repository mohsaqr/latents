# Latents review fixes and repeated analysis

**Release verification complete:** Version 0.96.1 is committed and pushed at `d4c9f69`. All six full platform checks, equivalence and lint pass; the final unit suite passes 8,144 expectations. See the [final HTML release report](../release-0.96.1-2026-10-20/RELEASE.html) for platform follow-ups, final archive and exact-commit evidence. The results below preserve the earlier repair-phase evidence.

This working tree repairs the ten findings in the [review of e29a310](../review-2026-10-20/REVIEW.md), with regression tests and retained before/after evidence. The repeated analysis also identified and repaired related failures in contributing-unit counts, direct-effect centering, Vuong eigenvalues and paired-data preparation. These repairs and the subsequent release follow-ups are now committed and pushed; the final release record above supersedes the repair-phase status.

## Closure of the original findings

| Finding | Correction | Regression evidence |
| --- | --- | --- |
| F01 Empty clusters permit invalid robust inference | Count observed positive-weight groups for hard validity and score-rank checks. Use contributing units for CR1, with the explicit legacy empty-count convention confined to CR1. Apply the same principle to covariate, transition, additive, regression and growth inference. | `test-review-1020-inference.R`, `test-review-1020-zero-weights.R`; [inference evidence](inference/REPORT.txt) |
| F02 Pair names overwrite caller indicators | Refuse initial-fitting collisions before replacing a column. Preserve prepared-data identity and permit known-model prediction/data adapters to rebuild pairs from original items. | `test-review-pairs-trees.R`; [pair evidence](pairs-trees/REPORT.md) |
| F03 Weighted regression trees use sample shares | Select weights in the same row, group or person units as the posterior matrix. | `test-review-pairs-trees.R`; [tree evidence](pairs-trees/REPORT.md) |
| F04 Direct-effect uncertainty depends on predictor origin | Evaluate information in centered predictor coordinates, transform complete covariance matrices, information and scores back to the reported chart, and expose conditioned information to VLMR. | `test-review-1020-inference.R`; independent Gaussian `lm()` references |
| F05 ML distal uncertainty depends on outcome origin | Fit and differentiate on centered, scaled outcomes, retain the accurate internal coordinates, then translate reported means and rescale influences. Check positive information before inversion. | `test-three-step-ml.R`; [numerical evidence](numerics/REPORT.md) |
| F06 VLMR raw inversion loses scale protection | Condition information and scores consistently, then obtain reference eigenvalues through an equivalent symmetric calculation. | `test-lmr.R`; public model comparisons across indicator units and direct-effect origins |
| F07 Pair joining merges native numeric categories | Share native categorical validation between the ordinary encoder and pair preparation, before character conversion. | `test-review-pairs-trees.R`; ambiguous-label parity and supported-type cases |
| F08 Independent mover-stayer reference overflows | Subtract the largest logit before exponentiating, keeping the independent exhaustive path reference. Require finite reference probabilities. | [original failing run](reference-before.log), [repaired run](reference-after.log); `tests/equivalence/test-lta-general.R` |
| N01 Negative Imhof statistic follows the wrong regime | Use the sign/complement identity to handle negative statistics, with exact single-sign cases and degenerate-zero distributions. | `test-lmr.R`; exact chi-square and convolution references |
| N02 Relative cutoff removes real small weights | Retain nonzero distribution weights; handle numerical null directions at the symmetric eigenvalue producer using its floating-point error scale. | `test-lmr.R`; exact F probabilities across the former cutoff |

The new tests exercise the original public paths and transformations. The initial reviewed package is preserved in `validation/review-2026-10-20/rlib`; several before-fix runs load that installed copy explicitly. Existing fixtures and numerical tolerances were preserved. The F08 reference formula was stabilized, and its likelihood assertion remains independent of the production implementation.

## Additional defects found during the repeated analysis

### Contributing units across model families

Zero-weight groups were incorrectly counted in the finite-sample correction of additive, regression and growth models. The root regression file demonstrates that appending such groups changed the whole covariance by the old CR1 ratio, even though the estimated model was unchanged. The repaired covariance is invariant. Equivalent empty-group checks cover membership-covariate, homogeneous-transition and general LTA inference.

The first rank safeguard still allowed a legacy `drop_empty = FALSE` fit with one informative group and one free parameter to return a tiny nonzero variance from floating-point score noise. Validity now always requires genuinely observed positive-weight groups, independently of the legacy reporting convention. This case has a separate fixed-means regression.

### Direct effects on missing targets

Centering a predictor using rows whose target indicator was missing—or whose sampling weight was zero—could move the optimization chart far from the data that actually identify its effect. Extreme values on those rows caused a failed fit and likelihood drift. Internal predictor columns now receive target-specific centers and scales computed from observed positive-weight target rows. The public fit keeps the original covariate names and data. Shared-predictor, different-target missingness cases are checked against independent Gaussian regressions.

### Vuong numerical null directions and adjusted sample size

The old unsymmetric eigenvalue calculation could generate spurious positive/negative weights around `1e-7` for an exactly nested known-answer model. Its probability at zero was about `0.99981` instead of one. An information-whitened signed Gram matrix gives an equivalent symmetric problem and permits numerical null cleanup at the matrix-construction error scale, while retaining genuine small weights.

The adjusted LMR sample size also counted all score rows, including empty groups. In a reproduced 300-row example, adding 95 entirely missing rows changed the adjusted LR from 317.5581 to 318.3669 although the LR and reference moments were unchanged. The adjustment now uses the fit's contributing-unit convention; full score rows remain aligned for cross-model products.

### Paired input validation and prepared frames

Joining names could remove duplicate original `vars` or `categorical` declarations and conceal an identifier reused as a paired item. These identities are checked before transformation. Native types, dimensions, nonfinite values and ambiguous numeric labels receive consistent categorical validation.

The initial collision repair also exposed a compatibility regression: ordinary data-frame operations can discard custom factor or frame attributes. The repeated review exercised serialization, subsetting, row binding, merging and column selection. Initial fitting retains strict collision checks; fitted-model prediction and data reconstruction can rebuild a pair from its original source columns without relying on attributes surviving the caller's transformations. Stale metadata cannot authorize overwriting an unrelated raw column during initial fitting.

## Validation

| Check | Result | Retained evidence |
| --- | --- | --- |
| Complete development unit suite, with corrected catalogue file rerun | **1,209 cases; 8,104 passing expectations; zero failures, errors, warnings or skips** | [summary](validation-summary.log), [final cases](unit-final-cases.csv) |
| Complete development equivalence suite | **64 cases; 405 passing expectations; zero failures, errors, warnings or skips** | [log](equivalence-tests.log), [cases](equivalence-tests.csv) |
| Fresh package fits against retained Latent GOLD results | **376 agreements, zero disagreements; 11 not comparable** across 12 cases with 50 starts each | [comparison](latentgold/comparison.csv), [log](latentgold.log) |
| Independent core references | **13 passing expectations** | [log](independent-core-final.log) |
| Independent numerical probability stress study | **70 references; maximum absolute error 5.48e-12** | [results](numerics/adversarial.csv) |
| Configured package lint | **Clean** | [log](lint-final.log) |
| Final source build and native package check | **Zero errors, zero warnings, one NOTE: unable to verify current time** | [build](build-catalogue-final.log), [check](check-catalogue-final.log) |

The complete unit run initially had one failing documentation-contract assertion: the new `latents_invalid_lmr` condition was absent from the condition catalogue. Its catalogue entry and generated help were corrected, then the entire affected test file passed all 34 expectations. The final unit total replaces that file's original results, without double counting. The [original unit log](unit-tests.log) retains the failure and its source-change flag. The source change during that run was the catalogue comment; [parsed-code verification](catalogue-documentation.log) and [the aggregation script](summarize-final.R) establish that all 91 production files have identical executable expressions before and after that correction. No numerical implementation changed during the full run.

The package check rebuilds vignettes and runs examples, including `donttest` examples. Its test stage and PDF manual are omitted; unit and development equivalence tests run separately as shown above. The final check log records the environment-related time-verification NOTE explicitly.

The package was installed into an isolated library for socket workers. Full development runs use `NOT_CRAN=true`. The full native package check runs tests separately from its build/example/vignette checks. Source fingerprints are retained so that a result can be tied to the implementation actually tested.

Retained external Latent GOLD output is compared with newly computed package targets, without regenerating or overwriting the external fixtures. This verifies continued agreement with those reference cases; it is not a new run of the external program.

The old fixtures passed because their reference datasets did not exercise the triggering transformations and combinations: empty or zero-weight groups, large origins or changed measurement units, weighted tree adapters, and colliding or ambiguous categorical inputs. The new regressions test these properties directly and retain independent known-answer references. Passing the fixtures alone would still be insufficient evidence for those properties.

## Independent repeated review

The repeated analysis includes cross-review by reviewers who did not author the respective repairs, known-answer numerical distributions, full public transformation paths and combined-feature probes. The new integration test combines paired categorical data, ordinal and count indicators, direct effects, sampling weights, a large predictor origin, prediction, refitting and empty zero-weight groups. Its Gaussian coefficient block is checked against an independent weighted regression.

The numerical stress study covers 70 independent probability references; its largest measured absolute error is approximately `5.48e-12`. Tests distinguish representable data precision from algorithm errors by comparing large-origin outcomes with a centered version of those same floating-point values. Rescaling Gaussian indicators also requires scaling an explicit variance floor quadratically to preserve the same statistical model.

All original findings and additional reproduced issues listed above are closed. Cross-review and the repeated checks found no further reproducible defect within this review's scope.

Review limits remain concrete: local execution is on R 4.5.2/macOS, remote Windows/Linux/minimum-R CI was not run here, and no formal correctness proof or fresh branch-coverage claim is made. Unsupported model combinations remain explicitly refused. Ratios beyond representable floating-point distribution coordinates receive a typed unavailable-probability condition.

## Files and reproduction

- [run-suite.R](run-suite.R): complete unit or development equivalence execution, results and source fingerprints.
- [run-retained-latentgold.R](run-retained-latentgold.R): recompute package targets and compare retained external output without changing reference artifacts.
- [final delivery manifest](final-delivery-manifest.csv) and [integrity verification](delivery-integrity.log): fingerprints for 319 source/test/documentation files, with archive matching and the existing development-equivalence build exclusion recorded.
- [checked source archive](latents_0.96.0.tar.gz): the repaired package, built and checked locally; SHA256 `fa7d2d8f4eaaf8f964e52afee3ac473e1d3edf9ec123c2182b729dc4bae20ddd`.
- [inference/REPORT.txt](inference/REPORT.txt): contributing-unit, covariance-transformation and target-centering evidence.
- [pairs-trees/REPORT.md](pairs-trees/REPORT.md): pair preparation, weighted adapters and independent cross-review.
- [numerics/REPORT.md](numerics/REPORT.md): ML, VLMR, Imhof and numerical known-answer evidence.
- [weighted tree preview](../../tmp/weighted-tree-review.html): rendered corrected weighted shares and bands, with a numerical agreement check.

Run from the repository root. Set `R_LIBS` to `validation/fixes-2026-10-20/rlib` after installing the final source there, and `NOT_CRAN=true` for development suites. Local socket permissions are needed by `test-parallel.R` on this macOS configuration.
