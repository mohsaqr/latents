# Latents latest code review

**Historical review of `e29a310`.** The findings below describe the original reviewed code. All ten findings have since been repaired in the working tree; see the [fixes and repeated analysis](../fixes-2026-10-20/FIXES.md) for changes, regressions and final validation. Original reproduction evidence is preserved here.

The current `features` branch needs corrections before release. Seven defects are reproduced through public APIs, including one high-priority inference defect: entirely missing clusters can bypass the minimum-cluster guard and produce effectively zero robust standard errors. Two further defects silently change categorical measurement data. The ordinary estimation checks pass, but those successes do not establish the validity of these affected combinations.

The findings below identify concrete triggers, numerical consequences, source locations, repair directions and regression tests. Production source, existing tests and documentation were preserved. Reproduction scripts and outputs are retained beside this report.

## Reviewed state and scope

| Item | Reviewed state |
| --- | --- |
| Date | 20 October 2026 |
| Branch and revision | `features`, `e29a310` |
| Package version | `0.96.0` |
| Production source | 91 R files, 44,986 lines |
| Parsed top-level functions | 1,212 |
| Namespace | 38 exports and 172 S3 registrations |
| Existing tests | 141 unit-test files and 29 equivalence-test files |
| Latest feature comparison | `dfddf58..e29a310`, affecting 25 production files |
| Broader recent comparison | Changes after the 17 October review fixes at `e70ecb5` |
| Local platform | R 4.5.2, Apple Silicon, macOS Tahoe 26.3.1 |

Every production R file was parsed and fingerprinted. Detailed semantic review concentrated on the latest direct effects, categorical and continuous local dependence, ML three-step outcomes, VLMR/Imhof calculations, and their inference, prediction, enumeration, bootstrap, comparison and results integrations. The broader pass examined tree adapters, parallel starts and condition handling, informative-unit accounting, weights, inference caches, package metadata and CI configuration. Independent Gaussian and multilevel likelihood calculations supplement the existing suites.

This is a bounded code and numerical review, not a proof that every branch of all 1,212 functions is correct. Fresh line or branch coverage was not measured. Cross-platform CI was inspected but was not executed on Windows, Linux or R 4.1. Historical Mplus and Latent GOLD fixtures were exercised where included in tests; those external programs were not rerun. Production vignettes were built and checked; the separate website articles were not all rendered independently.

The source inventory is [function-inventory.csv](function-inventory.csv), with file fingerprints in [source-manifest.csv](source-manifest.csv). [inventory.R](inventory.R) reproduces the inventory.

## Findings by priority

P1 means a correction should precede release because the demonstrated result can invalidate an analysis. P2 means a reproducible defect in a supported case or a failing development check. The numerical-helper findings have narrower demonstrated exposure and are listed separately.

| ID | Priority | Finding | Demonstrated consequence |
| --- | --- | --- | --- |
| F01 | P1 | Empty clusters bypass robust-inference validity checks | One informative cluster produces SE `1.32e-17` and `p = 0` after adding empty clusters |
| F02 | P2 | Generated pair names overwrite existing indicators | Four parameters fitted where five are required; original data disappear |
| F03 | P2 | Regression trees discard sampling weights | Tree class shares differ from fitted weighted shares by 60 percentage points |
| F04 | P2 | Direct-effect inference depends on predictor origin | Adding 1,000 or 2,026 to a predictor causes false singular-information errors |
| F05 | P2 | ML distal uncertainty depends on outcome origin | A location shift changes SEs from about `0.09/0.08` to `2178/4165` |
| F06 | P2 | VLMR reinverts information without coordinate scaling | Rescaling indicators by `1e8` aborts a previously valid comparison |
| F07 | P2 | Pair joining bypasses categorical label validation | Four joint categories silently become two |
| F08 | P2 | Mover-stayer equivalence reference overflows | Development equivalence suite fails despite an exactly matching production likelihood |

F01 concerns older inference code interacting with empty-group handling. F03 was introduced with regression-tree support in `dfddf58`. F02 and F04–F07 concern the latest feature work. F08 is in an older equivalence test. These are findings against the current revision; they are not all regressions introduced by its most recent commit.

## Confirmed public API defects

### F01 Empty clusters can produce false certainty in robust inference

**P1.** The validity guard in [R/inference.R](../../R/inference.R), lines 1155–1159, counts `object$n_groups`. The sandwich calculation at line 942 also counts every score row. Neither excludes groups with no observed indicators, although `drop_empty = TRUE` records that they contribute no observed information.

A synthetic one-profile fit with 50 observed values in **one informative group** correctly rejects robust inference as `latents_too_few_groups`: it has one independent group and two free parameters. Appending 95 groups whose indicators are entirely missing preserves its likelihood and estimates. The padded fit nevertheless returns:

| Quantity | Padded result |
| --- | --- |
| Informative groups | 1 |
| Mean estimate | −0.1974855 |
| Robust mean SE | `1.320771e-17` |
| Mean z statistic | `−1.495229e16` |
| Mean p-value | 0 |
| Mean confidence interval | Collapses to the estimate at printed precision |

The lone informative group's score is essentially zero at the estimate. Treating its squared score as an adequate sandwich variance creates numerical false certainty. Two informative groups also change from correctly refused to accepted after adding empty groups. Padded OPG inference detects singularity later; that does not protect the robust path.

Even with enough informative groups, the cluster correction changes incorrectly. With five informative groups, adding 95 empty groups changes the mean SE from **0.1083709 to 0.09741823**, while likelihood and point estimates remain unchanged. The ratio is exactly the change from the CR1 factors for 100 versus five groups: `sqrt((100/99)/(5/4)) = 0.8989331`.

**Repair:** use the informative independent-unit count consistently in validity checks and the CR1 correction; filter score rows to those units where appropriate. Audit the equivalent covariate and transition paths. The public reproduction establishes the ordinary `multilpa` defect; the other paths require their own regression cases.

**Regression requirement:** empty-padding invariance for estimates, scores and robust uncertainty when inference is valid; unchanged refusal with one or two informative groups; coverage of positive and zero sampling weights and the legacy `drop_empty = FALSE` convention.

Evidence: [minimum-group probe](api/empty-group-minimum-guard.R), [output](api/empty-group-minimum-guard.log), [five-group probe](api/empty-group-robust.R), [output](api/empty-group-robust.log).

### F02 A generated joint indicator can overwrite a separate input indicator

**P2.** [R/local-dependence.R](../../R/local-dependence.R), line 214, assigns `frame[[pair$name]]` without reserving that name. The name rewrite at lines 155–157 removes duplicates, so the collision can silently remove a requested indicator.

Supply four independent binary columns named `a`, `b`, `a:b` and `c`, and request `local_dependence = "a:b"`. The fit replaces the existing third column with the joint responses of the first two and retains only `c("a:b", "c")` as its measurement variables. The caller's original `a:b` observations have disappeared.

| Fit | Parameters | Log likelihood |
| --- | --- | --- |
| Current paired implementation | 4 | −828.8379 |
| Collision-free joint indicator retaining the original third variable | 5 | −1105.5966 |

These likelihoods describe different measurement data; the difference is evidence of data loss, not a model-fit improvement. Parameter counts and information criteria consequently refer to the wrong requested analysis.

**Repair:** validate all generated names against existing column roles before mutation, or use reserved unique internal names with an explicit public-label map. Preserve intentional reuse of already prepared joint columns without mistaking an unrelated raw column for one.

**Regression requirement:** collisions with another indicator, identifier and covariate; multiple disjoint pairs; ordinary refit and stored-data prediction. Either retain all requested data or reject the collision clearly.

Evidence: [confirmation script](effects/confirm.R), [output](effects/confirm.log).

### F03 Weighted regression trees report unweighted class shares and flows

**P2.** [R/tree.R](../../R/tree.R), line 358, returns unit weights for regression and growth fits. Regression grids publicly advertise trees as supported, and `get_results(grid, "tree")` reaches this branch.

The reproduction has 160 observations near −4 with weight 1, and 40 near +4 with weight 16. A one-to-two-class regression grid converges. In the fit's class order, weighted posterior shares are **80% and 20%**. The tree reports **20% and 80%** for its destination shares and masses: a **60 percentage point discrepancy**. Its bands describe the unweighted sample composition, despite the model having been fitted with sampling weights.

**Repair:** match weights to the posterior units. Row-level regression posteriors need row weights; grouped regression posteriors and group classes need independent-group weights; growth posteriors need person-level weights. The shared regression/growth branch is implicated on inspection; the numerical public reproduction is for regression.

**Regression requirement:** independently calculated weighted node shares and edge masses through the public results API, plus grouped regression and growth cases. Confirm that every node's incoming and outgoing mass agrees with its weighted share. Rendering changes should receive visual review separately.

Evidence: [probe](api/weighted-regression-tree.R), [output](api/weighted-regression-tree.log).

### F04 Direct-effect Wald inference fails after a predictor translation

**P2.** [R/inference.R](../../R/inference.R), lines 864–881, centers the outcomes but evaluates the likelihood using raw direct-effect covariates. Scaling the effect coordinates at lines 903–914 does not remove the intercept/slope collinearity caused by a large predictor origin.

A 400-row, one-profile Gaussian fit with `direct_effects = "y1 ~ z + w"` has likelihood **−1116.26508119** and a slope SE of **0.04878294499**. Replacing `z` by `z + 1000` or `z + 2026` preserves the likelihood to about `5e-13` and slopes to about `3e-15`, but `parameter_inference()` raises `latents_singular_information`.

An independent Gaussian MLE reference based on `lm()` retains the same finite slope SE. The intercept SE at the shifted origin is also finite, about **98.8372**. The design is mathematically identified; the failure is in the information coordinates. Predictors such as calendar years can have this origin-to-spread relationship.

**Repair:** construct the information and group scores in centered predictor coordinates, as fitting already does, and transform the full covariance back, including intercept/effect cross-covariances.

**Regression requirement:** predictor-translation invariance for slopes, their SEs and tests, with independently verified intercept transformations; standard and robust inference; several direct covariates.

Evidence: [confirmation script](effects/confirm.R), [output and independent reference](effects/confirm.log).

### F05 ML three-step standard errors depend on the distal outcome origin

**P2.** [R/three-step-ml.R](../../R/three-step-ml.R), lines 227–233, adds finite-difference steps based on outcome spread directly to absolute class-mean coordinates. At a sufficiently large origin, floating-point rounding distorts or eliminates those steps.

On the public synthetic example, `three_step(..., method = "ml", vcov_type = "independent")` initially returns SEs of about **0.090 and 0.084**. Adding `1e12` to the outcome preserves its spread and the fitted means after reversing the translation, within representable data precision, but returns SEs **2177.97 and 4165.42**. At an offset of `1e10`, errors already approach **10%**. The incorrect uncertainty is returned without a diagnostic.

A separate hard-classification example supplies a known normal-mean reference and confirms the same failure. This is an extreme origin-to-spread case, rather than evidence that ordinary-scale ML outcomes are generally wrong. It extends the numerical issue addressed by the recent coordinate-scaling fixes.

**Repair:** estimate and differentiate in centered outcome coordinates, preferably scaled as well, then translate means and intervals back. Validate the conditioned information before using its inverse.

**Regression requirement:** exact hard-classification uncertainty and public outcome-translation invariance across large origin-to-spread ratios, both variance options, and both independent and clustered contrasts.

Evidence: [public and numerical probes](ml-lmr/probes.R), [output](ml-lmr/probes-output.txt), [additional probes](ml-lmr/additional-probes.R), [output](ml-lmr/additional-output.txt).

### F06 VLMR loses the information engine coordinate scaling

**P2.** [R/lmr.R](../../R/lmr.R), lines 48–49, calls `solve()` on raw information matrices returned by the inference engine. Those raw coordinates can be badly scaled even when the engine's scaled information is regular.

A converged EEI grid with one and two profiles gives LR **336.1164**, reference mean **4.418811**, and reference SD **3.286323**. Multiplying both indicators by `1e8` leaves both candidates converged, and their robust `parameter_inference()` calls remain valid. Nevertheless, `get_results(grid, "lmr")` aborts with:

```text
system is computationally singular: reciprocal condition number = 7.87581e-17
```

Applying consistent diagonal coordinate scaling to each information matrix and its score columns recovers the original reference mean, SD and LR from the rescaled public fits. This confirms a coordinate-conditioning failure rather than a newly unidentified model.

**Repair:** form the Vuong matrix from consistently scaled scores and information, reusing the inference service's scales where possible. Also preserve the table's intended per-comparison failure behavior if a comparison genuinely cannot be calculated.

**Regression requirement:** public LMR invariance under small and large indicator units, with full-covariance and direct-effect coordinates where supported.

Evidence: [probes](ml-lmr/probes.R), [output](ml-lmr/probes-output.txt), [whitened comparison](ml-lmr/additional-probes.R), [output](ml-lmr/additional-output.txt).

### F07 Joining pairs can silently merge distinct numeric categories

**P2.** [R/local-dependence.R](../../R/local-dependence.R), lines 205–214, converts pair members to character before the ordinary categorical encoder can validate their labels. [R/categorical.R](../../R/categorical.R), lines 32–46, already rejects distinct values that print alike.

The values `.3` and `.1 + .2` are distinct doubles but both print as `"0.3"`. Ordinary `lca()` correctly rejects a categorical column containing both as `latents_bad_categorical`, asking the caller to resolve the labels. Adding `local_dependence = "a:b"` bypasses the check, merges four joint categories into two, and succeeds.

In the deterministic example, the paired fit reports **two parameters and likelihood −277.2589**; explicit distinct category coding reports **four parameters and likelihood −415.8883**. A caller may reasonably choose to round these values, but that choice should be explicit and consistent across paired and unpaired fitting.

**Repair:** validate the original categorical columns before joining, or join validated native category codes while keeping a one-to-one label map.

**Regression requirement:** paired/unpaired parity for ambiguous numeric labels, invalid categorical types and nonfinite values, while retaining valid factors and numeric category labels.

Evidence: [probe](effects/labels.R), [output](effects/labels.log).

## Development check defect

### F08 The mover stayer equivalence reference overflows

**P2, test infrastructure.** The full development equivalence suite has one failure at [tests/equivalence/test-lta-general.R](../../tests/equivalence/test-lta-general.R), line 293. Its independent reference calculates initial probabilities using `exp(logits) / sum(exp(logits))` at lines 277–279.

The reproduced fitted stayer logit is **1313.061**. Exponentiating it produces infinity, so the reference probability becomes `Inf / Inf = NaN`. The assertion compares a finite production likelihood with `NaN`.

Replacing only the reference normalization with `exp(logits - max(logits))` gives the exact path-sum likelihood **−899.06784279732847**, equal to the production result with a difference of **zero** in this run. This establishes the observed failure as a reference-calculation defect, not evidence that the production forward likelihood is wrong.

**Repair:** stabilize the independent reference normalization without calling the production likelihood routine. Keep the large-logit case as a regression so that the reference remains independent and numerically usable.

Evidence: [equivalence suite output](equivalence-full.log), [isolated diagnosis](equivalence-mover-stayer.R), [original and stabilized reference values](equivalence-mover-stayer.log).

## Numerical helper defects with narrower demonstrated exposure

### N01 Negative statistics use the wrong Imhof integration regime

[R/lmr.R](../../R/lmr.R), lines 125 and 133, treats all `q <= 0` as nonoscillatory. Only `q == 0` has that property. `.latents_imhof_upper(-1, 1)` and the versions with weights `-1` or `c(1, -1)` raise `latents_inaccurate_probability` while integrating a well-defined distribution. For weight 1, the exact upper-tail probability at −1 is 1.

The helper defect is reproduced. A public neighboring-model grid that reaches it was **not** reproduced in the attempted synthetic search. A globally optimized nested likelihood ratio is nonnegative; an EM-converged local optimum does not by itself guarantee that ordering. Either handle negative arguments by an appropriate distributional symmetry or explicitly return an invalid-comparison result for a reversed likelihood. Add exact negative-statistic tests and a reversed-fit policy test. Evidence: [probe output](ml-lmr/probes-output.txt).

### N02 Relative eigenvalue filtering can erase a nonnegligible tail

[R/lmr.R](../../R/lmr.R), line 88, discards weights at or below `1e-10 * max(abs(lambda))`. Small relative magnitude does not establish small tail-probability error near zero.

For `q = 0` and weights `c(1, -1e10)`, the helper returns approximately **2.25e-13**. The exact F-distribution identity gives **6.3661977e-6**. Just below the filtering threshold, the corresponding test agrees to about `1e-16`. The error therefore comes from discarding a genuine small weight.

This is a confirmed helper-accuracy limitation; no public fitted grid producing that particular imbalance and near-zero statistic was reproduced. Use a justified numerical-rank or probability-error treatment, or signal the unsupported case. Add exact F cases on both sides of the cutoff. Evidence: [probe](ml-lmr/extreme-weights.R), [output](ml-lmr/extreme-weights-output.txt).

## Validation results

| Check | Result | Evidence |
| --- | --- | --- |
| Combined unit suite | 968 cases; **6,263 passed expectations**, zero failures/errors/warnings/skips | [summary](validation-summary.log), [case results](unit-combined.csv) |
| Full development equivalence suite | 64 cases; **403 passed, one failed expectation**, zero errors/warnings/skips | [log](equivalence-full.log), [case results](equivalence-full.csv) |
| Added independent core references | Four cases, **13 passing expectations** | [log](independent-core.log), [script](test-independent-core.R) |
| Parallel/socket tests | **28 passing expectations**, included in combined unit count | [log](parallel-tests.log) |
| Source installation | Passed | [log](install.log) |
| Configured lint | No lints found | [log](lint.log) |
| Full source build | Passed, including vignette creation | [log](build-full.log) |
| Native `R CMD check --as-cran --no-manual --no-tests` | **Zero errors, zero warnings, one NOTE**: unable to verify current time | [log](check-full.log), [native check log](full-check/latents.Rcheck/00check.log) |
| Examples, including `donttest` examples | Passed | [native check log](full-check/latents.Rcheck/00check.log) |
| Installed vignettes and rebuilding their outputs | Passed | [native check log](full-check/latents.Rcheck/00check.log) |
| Production source fingerprints | Unchanged after review and checks | [verification](validation-summary.log) |

The equivalence failure is F08. Its stabilized independent reference matches the production likelihood exactly; the original failing test was preserved. All seven public defects remain unfixed and are demonstrated by the separate retained probes, despite the existing unit suite passing.

The original unit run used the default CRAN setting. Every file containing a skipped case was subsequently selected for a development run with `NOT_CRAN=true`; the combined summary replaces those files' initial results rather than double-counting their passing cases. The separate development equivalence run also uses `NOT_CRAN=true`.

Six development test errors initially came from sandbox restrictions on local worker sockets. The complete parallel test file was rerun with the required socket access and passed. The combined summary replaces that file's blocked run with the successful run. [summarize-results.R](summarize-results.R) produces the aggregate and rechecks source fingerprints.

The current source was installed into an isolated library before tests, and that library was supplied through `R_LIBS`, so socket workers used this revision. The native package check skips its own tests because source suites ran separately. Remote incoming checks and forced optional dependencies were disabled in the restricted environment. Network repository-index messages are retained in the logs; this is not an online CRAN submission certification.

An initial build deliberately omitted vignette compilation and consequently produced expected missing-vignette warnings. The full build/check supersedes that diagnostic run; its results are the release-check results used above. Both logs are retained to make the distinction explicit.

## What the passing checks establish

The package has substantial numerical and integration testing. Independent checks here verified the following without using the production likelihood as the reference:

- One-profile full-covariance Gaussian MLE means, covariance, likelihood and mean standard errors against closed forms.
- Weighted Gaussian means and covariance against the weighted MLE formulas.
- A tiny multilevel mixture likelihood against exact marginalization over group classes and row mixtures.
- Missing-data Gaussian likelihood against independently evaluated observed margins.
- Nonadjacent continuous covariance pairs against the masked sample covariance and expected parameter count.
- Direct effects with missing responses and sampling weights against Gaussian regression references.
- A mixed model combining a categorical dependence pair, direct effects, an ordinal indicator and a count indicator through fit, prediction and inference.
- Public VLMR reference moments after consistent coordinate scaling.

The shared measurement-density hook, shared bootstrap measurement arguments and explicit unsupported-combination checks are useful architectural strengths. The observed defects repeatedly arise when an adapter reconstructs information already handled elsewhere: unit weights in trees, raw coordinates in VLMR, raw predictor origins in inference, or labels generated before categorical validation.

Existing feature tests pass on ordinary cases while the new adversarial probes fail. The missing test dimensions are transformations and integration invariants: adding empty units, translating predictors/outcomes, changing units, colliding names and crossing numerical-rank thresholds. Adding more fixed output snapshots alone would leave those dimensions exposed.

## Recommended repair order

1. Correct F01's informative-unit guard and sandwich correction. Verify refusal and padding invariance before trusting robust uncertainty for incomplete clustered data.
2. Correct F02 and F07 before paired categorical measurement can silently alter the input. Retain a clear distinction between raw and already prepared joint columns.
3. Correct F03's weighted tree adapters, with explicit node/edge mass references for each posterior unit.
4. Correct F04–F06 using centered and scaled information coordinates and complete covariance transformations. Keep fitting, scores, uncertainty and reporting in a coherent chart.
5. Correct F08's independent test reference, then require a clean full development equivalence run. Harden N01/N02 with known-answer distribution tests.
6. Run the regression cases, the complete unit and equivalence suites, a vignette-producing build and native package check; confirm supported platforms through CI before release.

The corrections should be separate, reviewable changes with tests demonstrating the requested statistical invariants. No fixes were applied as part of this review.

## Reproduction guide

Run from the repository root with the dependencies in `DESCRIPTION` and the development equivalence dependencies installed. The isolated library must contain the reviewed revision; recreate it with:

```sh
R CMD INSTALL --library=validation/review-2026-10-20/rlib .
```

Set `R_LIBS` to that directory and `NOT_CRAN=true` for development runs. The key evidence scripts are:

| Area | Script |
| --- | --- |
| Core closed-form and marginal-likelihood references | [test-independent-core.R](test-independent-core.R), run with `testthat::test_file()` |
| Empty-cluster validity and CR1 | [minimum guard](api/empty-group-minimum-guard.R), [five groups](api/empty-group-robust.R) |
| Weighted regression tree | [weighted-regression-tree.R](api/weighted-regression-tree.R) |
| Direct-effect translation and pair-name collision | [confirm.R](effects/confirm.R) |
| Pair-category validation | [labels.R](effects/labels.R) |
| Covariance and regression references | [invariants.R](effects/invariants.R) |
| Public ML and VLMR edges | [probes.R](ml-lmr/probes.R), [additional-probes.R](ml-lmr/additional-probes.R) |
| Imhof relative cutoff | [extreme-weights.R](ml-lmr/extreme-weights.R) |
| Equivalence reference diagnosis | [equivalence-mover-stayer.R](equivalence-mover-stayer.R) |

The specialist notes in [effects/review-report.txt](effects/review-report.txt) and [ml-lmr/REPORT.md](ml-lmr/REPORT.md) contain additional inspected paths, positive controls and numerical detail. All substantive conclusions in this report are supported by local source or retained executable evidence.
