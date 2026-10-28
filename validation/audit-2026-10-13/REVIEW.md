# Latents function and documentation audit

> Historical pre-repair audit. All 27 confirmed findings were subsequently addressed at the user's request; see the [repair closure and final checks](../fixes-2026-10-13/FIXES.md). Original evidence and descriptions are retained below.

The current working tree has **27 verified findings and contract gaps**, including **two high-priority paths that produce invalid uncertainty or classification results**. Several other findings change the requested model, misreport weighted quantities, or break combinations supported by the new API. The existing suites pass, so a passing package check alone is insufficient to approve this working tree for release.

The latest enumeration, named plotting, comparison and help consolidation work is the main focus. Earlier model and inference defects are identified separately. This is an audit: production code, existing tests, documentation sources and pending user edits were preserved. The proposed corrections below have not been applied.

## Scope and evidence

| Item | Audited state |
| --- | --- |
| Date | 13 October 2026 |
| Committed baseline | `618f8de`, release 0.9.23, 7 October 2026 |
| Compatibility baseline | `a10a919`, release 0.9.21, 6 October 2026 |
| Production source | 80 R files, 37,181 lines |
| Function definitions | 1,023, plus four public aliases |
| Public definitions and aliases | 197 exports or registered S3 methods |
| New or changed function bodies since HEAD | 92; comparison ignores source formatting/comments |
| Public documentation sources and metadata | 62 files, including every Rd and vignette/article Rmd |
| Supporting documentation | 46 validation reports, reference-study documents, plans and tutorial sources inventoried separately |
| Main local suite | 902 cases, 5,973 passed expectations, zero failures/errors/warnings/skips |
| Local equivalence suite | 56 cases, 376 passed expectations, zero failures/errors/warnings/skips |
| Installed unit-suite line coverage | 95.2468% |
| Functions with hits in that coverage run | 994 of 1,023; 29 without hits listed explicitly |
| Static usage checks | All 1,023 function definitions checked; raw diagnostics retained in the inventory |
| Configured lint | One style lint for numeric-category coercion using `suppressWarnings()` |
| Native package build | Passed, including vignette compilation |
| Native check | `--as-cran --no-manual --no-tests`: zero errors, zero warnings, two environment-related notes |
| Executable documentation | Regular and `donttest` examples passed; installed vignettes and rebuild checks passed; both new enumeration documents rendered with hidden numerical guards |

The two check notes concern network link/ORCID validation and inability to verify the current time. They are not evidence that the referenced links are invalid. The check deliberately skipped its own test execution because the complete source suites ran separately. An additional installed unit-suite execution supplied coverage. R was 4.5.2 on this machine; this audit does not establish current Windows, Linux or minimum-R behavior.

Every production function has an inventory row with source span, purpose where available, arguments, called functions, validation/condition calls, RNG/state effects, direct test-symbol references, semantic recency and measured runtime coverage. Direct test references are a navigation aid, not a coverage claim. All definitions received structural screening. Specialist inventories document semantic review of 402 core functions, 163 API/accessor functions and 101 plotting/prediction functions. Shared-service review and targeted probes supplement those tracks. No claim of a formal correctness proof for every input or every mathematical branch is made.

The raw per-function usage checker reports diagnostics for 41 functions, principally ggplot/formula evaluation names and some unused local bindings. These are retained for inspection rather than promoted to confirmed runtime defects; native package code checks pass. The single configured lint is `R/plot-gg.R:313`, where numeric category conversion suppresses coercion warnings before validating the result. It is a style-policy discrepancy, not evidence that the conversion's checked output is wrong.

The 29 functions without hits refer to the installed unit-suite coverage run only. Some are exercised by source-only equivalence tests, help examples, aliases or audit probes. Line coverage establishes execution, not complete branch coverage or correct scientific interpretation. Historical supporting documents were inventoried and consulted where relevant; they were not all rerendered or independently recalibrated.

## Highest priority findings

### SV01 Bootstrap inference accepts data from a different model

**P1, pre-existing.** In [R/inference.R](../../R/inference.R), lines 791–808, the bootstrap branch returns before the data-identity and likelihood checks used by Wald inference. [R/bootstrap-inference.R](../../R/bootstrap-inference.R) then resamples the supplied frame using the original fit's grouping metadata without checking that its indicators reproduce the fitted data.

A converged one-profile synthetic model has variance 0.8511. Supplying the same identifiers with indicators multiplied by ten is correctly rejected by Wald inference as `latents_bad_inference_data`. Bootstrap inference accepts it, retains the original point estimate, multiplies the mean SEs by ten and variance SEs by 100, and returns a variance interval of **66.10–107.20 around the reported original variance 0.85**. That table combines estimates and uncertainty from different datasets.

Validate the supplied frame before either inference branch, including original rows, identifiers, indicators, preprocessing and relevant weights. Bootstrap tests should reject modified indicators and changed row/group layouts, and establish equivalence between stored data and an explicit identical frame. Evidence: [service-probes.log](service-probes.log), `BOOTSTRAP DATA OVERRIDE`; [test-service-findings.R](test-service-findings.R).

### PD01 Named diagnostic plots bypass the noise restriction

**P1, new uncommitted work.** [R/plot-named.R](../../R/plot-named.R), lines 92–96, calls the diagnostic drawing helpers directly. It omits the noise guard used by `plot(fit, what = "entropy"/"avepp"/"posteriors")`.

The synthetic noise fit's Gaussian posterior row sums range from 0.0547 to 0.9885 because the noise column is stored separately. Named entropy plots treat those incomplete rows as complete probability distributions and produce uncertainty up to **1.0498 on a declared zero-to-one scale**. AvePP similarly assigns noise cases to Gaussian profiles. The generic plotting method correctly refuses the same requests.

Apply the generic's restrictions to named views, including the default diagnostics collection, while retaining any separately supported size view. A regression should compare named and generic behavior for every diagnostic on a noise fit. Evidence and exact call sites: [plot audit](plots-docs/review.txt), PD-01; [plot evidence tests](plots-docs/findings-tests.R).

## Recent API and comparison findings

The IDs below map to [the detailed API review](recent-api/findings.txt), which supplies exact reproductions, causes, repair directions and suggested regressions.

| ID | Priority and origin | Finding and consequence | Main location |
| --- | --- | --- | --- |
| RA01 | P2, new | Subsetting or reordering a comparison preserves its original models attribute. Classification output includes removed models or stale ordering. | `R/compare-models.R:434` |
| RA02 | P2, new | Comparing a general LTA enumeration, demonstrated with occasion-varying transitions, crashes in cross-sectional information-criteria extraction. | `R/compare-models.R:336`; `R/enumeration.R:1533` |
| RA03 | P2, new | Direct homogeneous transition fits and membership-covariate fits are rejected even where comparisons through a grid work or the new help accepts fits from those fitting verbs. | `R/compare-models.R:110` |
| RA04 | P2, new | The smallest-class plot drops covariance structure. EEI and VVI candidates become one curve with duplicated profile-count coordinates and no structure labels. | `R/selection.R:344` |
| RA05 | P2, new | `get_results(regression_comparison, "all")` always requests unsupported classification and errors instead of returning available tables. | `R/get-data.R:712`; `R/compare-models.R:425` |
| RA06 | P3, new documentation | Shared help promises `"all"` for family/transition enumerations, but their methods still accept only `"candidates"` and `"best"`. | `R/get-data.R:211`; `R/additive-compare.R:254`; `R/lta-compare.R:444` |
| RA07 | P3, new contract gap | The same-data signature compares centered indicators and categorical codes without offsets or category labels. Different raw frames can pass the stated identity guarantee. The verified shift/relabel examples have equivalent likelihoods; incorrect BIC is not demonstrated. | `R/compare-models.R:302` |
| RA08 | P2, new scientific naming | The stated Kneedle method is a simplified maximum perpendicular chord-distance rule, without the cited algorithm's smoothing, local maxima and sensitivity/crossing steps. | `R/selection.R:110`; `R/get-data.R:191` |
| RA09 | P2, new comparison and older enumeration | Noise-model entropy uses Gaussian posterior columns without noise mass. A one-Gaussian-profile plus noise candidate reports NA despite finite full entropy; multi-profile values also disagree. | `R/compare-models.R:337`; `R/enumeration.R:218` |
| RA10 | P2, pre-existing | Enumeration with membership covariates fits successfully but crashes when assembling candidate rows because `n_best_replicated` is NULL for that family. | `R/enumeration.R:224` |

For RA04, the sample EEI/VVI smallest shares at two profiles are 0.00660 and 0.09650, drawn as a single series. Preserve covariance structure in the plot data and its series/facets. For RA09, the complete fit entropy is 0.877411 for one Gaussian profile plus noise, while enumeration and comparison report NA. Use the existing complete-posterior helper consistently.

For RA08, retain the owner's intentionally calibrated statistic but describe it accurately, or implement the full intended algorithm deliberately. The [original Satopaa paper](https://lass.cs.umass.edu/papers/pdf/simplex11.pdf), Section III, describes additional smoothing and threshold steps. A near-linear synthetic curve currently reports a knee with distance 2.36e-7. Agreement with one published knee value validates that particular distance calculation, not equivalence to the full algorithm.

## Model and inference findings

### MC01 Requested centering is ignored with membership covariates

**P2, pre-existing.** `multilpa()` accepts `centering`, but its covariate dispatch at [R/multilpa.R](../../R/multilpa.R), lines 1280–1292, does not forward it. The covariate engine fits raw indicators without a notice or refusal.

On 30 synthetic persons with six observations each, requesting none, person or grand centering gives exactly the same means **96.42578 and 103.06085**, and likelihood **−476.2394**. Explicitly person-centering the data gives means **−2.391015 and 2.555912**, and likelihood **−288.3874**. The accepted argument changes neither the fit nor its recorded metadata, despite requesting a different analysis.

Carry centering and offsets through fitting, data reconstruction, inference and prediction, or refuse the unsupported combination before fitting. [Core review](model-core/README.md), MC-1, contains the reproduction and source provenance.

### MC02 Variance floors break constrained covariance models

**P2, pre-existing.** [R/covariance-structure.R](../../R/covariance-structure.R), lines 175 and 403, clips variances/eigenvalues separately after constrained updates. This can violate equal-volume or shared-shape restrictions across profiles.

With a valid floor of 0.1, an EVI example returns determinants **0.10 and 0.01 despite requiring equal determinants**. A VEI example returns profile variance ratios 0.4805 and 0.1519 across dimensions despite requiring a common shape. An EVV example likewise violates equal volume. Public synthetic fits with `min_variance = 0.05` fail every start for EVI, VEI, EVV and VEV after an EM likelihood decrease.

Solve the lower-bound constraints within the structured optimization and verify both model constraints and ascent. Componentwise clipping after an unconstrained solution is insufficient. Unconstrained VVI/EEE/VVV floors are not implicated by these reproductions. [Core review](model-core/README.md), MC-2, and [failing invariant tests](model-core/test-findings.R) retain the evidence.

### MC03 Growth likelihood still loses information at extreme scales

**P2, numerical limitation originating in 0.9.12.** [R/growth-mixture-engine.R](../../R/growth-mixture-engine.R), line 148, subtracts large quadratic terms in its Woodbury calculation. The later fix for unstable residual cross-products does not stabilize this second subtraction.

For a random-intercept example with `y = c(1e8-1, 1e8, 1e8+1)`, random-effect variance 1e16 and residual variance one, an independent scalar Gaussian calculation gives log density **−23.22680**; the current function gives **−21.72680**. Its implied quadratic is zero rather than three. Ordinary-scale controls agree; the discrepancy at amplitude 1e7 is 0.03125.

This is a verified internal likelihood edge, not a demonstration of routine public fits failing at default settings. A stable projection or posterior-residual fallback should preserve residual information without sacrificing the low-dimensional calculation. [Core review](model-core/README.md), MC-3, explains practical exposure and limits.

### SV02 Weighted residual diagnostics describe the wrong distribution

**P2, pre-existing.** [R/bivariate-residuals.R](../../R/bivariate-residuals.R), especially lines 110–129 and the categorical path, uses unweighted posteriors/counts. Sampling weights used to fit the model do not enter those diagnostic moments.

A converged one-class weighted synthetic fit reports residual correlation **+0.1657**. The independent cross-product with its actual sampling weights gives **−0.7504**. This reverses the sign of association. The current approximate Fisher/chi-square inference also does not implement a sampling-weight design correction.

Use consistent weighted descriptive moments and define the appropriate inference, or explicitly refuse weighted residual tests until supported. Simply multiplying counts by weights does not establish valid survey/cluster inference. Test weighted moments against independent calculations, plus the unweighted special case. Evidence: [service-probes.log](service-probes.log), `WEIGHTED RESIDUALS`.

### SV03 Pooling constant probabilities produces undefined missing information

**P2, pre-existing; recent edits touch other parts of the file.** [R/pool-imputations.R](../../R/pool-imputations.R), lines 247–249, calculates relative variance as between/within. A structural probability of one has both variances zero.

Pooling two identical one-profile fits gives correct point estimates and zero-width intervals for the structural probabilities, but **`riv` and `fmi` become NaN**. `summary()` then reports `fmi_max = NaN` and `fmi_median = NA`, hiding useful information for the estimable parameters. Existing identical-imputation tests use a covariate family whose table avoids this particular constant-row case.

Define missing-information output for structural constants and the zero-within/positive-between case explicitly, and prevent constant rows from contaminating aggregate summaries. Preserve valid estimates and intervals. Evidence: [service-probes.log](service-probes.log), `ONE CLASS IDENTICAL IMPUTATIONS` and `RUBIN ZERO WITHIN VARIANCE`.

### SV04 Matrix columns corrupt descriptive counts

**P2, pre-existing.** [R/descriptives.R](../../R/descriptives.R), lines 225–243, accepts a numeric matrix stored as one data-frame column. Logical extraction flattens it, and matrix missingness recycles the group selector.

A six-row frame with three actual groups and a two-column numeric matrix indicator reports **12 observations and four groups**, silently inventing a group from recycled/missing identifier indexing. Reject matrix-valued descriptive variables or deliberately expand them into named scalar variables with aligned rows. Evidence: [service-probes.log](service-probes.log), `DESCRIPTIVES MATRIX COLUMN`.

### SV05 R3STEP depends on predictor units

**P2, pre-existing numerical robustness issue.** The ordinary `r3step()` branch optimizes unscaled predictors and supplies an all-one Hessian scale at [R/three-step.R](../../R/three-step.R), lines 669–680. The separate group-class branch already uses a predictor-scale helper.

In the 400-row synthetic example, an ordinary predictor gives slope **−1.5063**, SE **0.1638**. Multiplying only that predictor by 1e-9, which preserves the statistical model after coefficient rescaling, causes `latents_singular_information`. The predictor is finite and the design remains mathematically full rank.

Optimize and evaluate curvature in appropriately scaled coordinates, then transform coefficients and covariance back. Verify fitted probabilities and inference under predictor-unit changes. This finding concerns extreme small units; ordinary-scale controls pass. Evidence: [service-probes.log](service-probes.log), `R3STEP PREDICTOR UNITS`.

## Plotting and documentation findings

PD01 is described above. [The plot and documentation review](plots-docs/review.txt) contains the remaining seven numbered findings and their detailed regressions.

| ID | Priority and origin | Finding | Location |
| --- | --- | --- | --- |
| PD02 | P2, pre-existing | Weighted tree flows/shares are unweighted. A fit with 9.75%/90.25% effective shares is plotted as 47.21%/52.79%. | `R/plot-gg.R:1393` |
| PD03 | P2, new compatibility regression | Released abbreviations such as `what = "enum"` and `"tr"` now become unknown criteria, while 0.9.21 resolved them with `match.arg`. | `R/plot.R:210` |
| PD04 | P2, new documentation regression | Shared plot help loses the released method aliases/usages, including `?plot.multilpa`; method-specific defaults are no longer inspectable there. | `man/latents-plots.Rd:4`; plotting methods marked `@noRd` |
| PD05 | P2, pre-existing | Aggregated binomial `predict(fit)` returns expected successes, but `predict(fit, training_data)` returns per-trial probabilities. Ratios equal the row trial counts, 10 or 20. | `R/mixture-regression-predict.R:46` and `:127` |
| PD06 | P2, new documentation | The new enumeration article describes the BLRT p-value as b/B, while implementation uses (b+1)/(B+1); its stated positive minimum conflicts with its prose formula. | `vignettes/articles/class-enumeration.Rmd:455` |
| PD07 | P3, carried-forward and new documentation | Shared prediction help unconditionally promises RNG restoration even though `seed = NULL` advances it, and broadly promises classed errors where regression prediction can return ordinary model-frame errors. | `R/predict-docs.R:31` and `:38` |
| PD08 | P2, intentional recent migration with unresolved compatibility | Ported base plots return ggplot objects and no longer draw when called inside a loop without printing. This contradicts unrestricted compatibility wording. | `R/mixture-regression-predict.R:224`; pooled/bootstrap ports |

For PD05, decide whether the public result is expected counts or per-trial probabilities, make training and new-data behavior consistent, and provide/document trial counts if new-data counts are supported. For PD08, decide the drawing contract explicitly. The probe confirmed an unchanged grid display list after two bare loop calls, and a populated display list after `print(plot(fit))`. Keeping the new ggplot design requires documenting that migration instead of claiming unchanged behavior.

### MC04 Public prior help omits intentional exceptions

**P3, pre-existing documentation.** The internal source in [R/prior.R](../../R/prior.R), lines 8–24, intentionally follows reference-program conventions: EII ignores prior scale, EEI does not shrink means and VEV ignores degrees of freedom. The public help presents these controls as generic conjugate-prior inputs without the exceptions.

The synthetic knob probes verify each exception. They do not establish an accidental estimator regression: the implementation intentionally preserves reference parity. Document the structure-specific behavior and qualify generic MAP wording. [Core review](model-core/README.md), MC-4, records the evidence.

### Additional documentation maintenance

The README's rendered enumeration output predates candidate names and AvePP columns, although its Rmd source calls the current API. Regenerate the rendered README after settling the API findings. The source data-generation comment that shipped engagement indicators are `log1p(count)` is too strong: the script directly simulates continuous values representing a log-count scale. Use the explicit README clarification consistently in case-study introductions, without rewriting historical release notes indiscriminately.

The shared plotting page should put argument-scope exceptions beside the general promise: some legacy methods warn about ignored arguments, others refuse them; the historical `intervals` scope excludes bars. Named wrappers also inspect `class(x)[1]` rather than the full S3 chain, so downstream subclasses can work with `plot()` and fail with a named wrapper. No package-built subclass failure was demonstrated.

The help-alias scan finds 102 registered/public definitions without their own Rd alias. Many are intentionally consolidated S3 methods, not 102 independently undocumented public APIs. Restore the old help entry points on shared topics where compatibility is required, while retaining the smaller manual.

## Assessment of the latest work

The candidate naming, compact recorded calls, lazy classification columns, selection tables and shared help are coherent improvements with substantial existing test support. The newly rendered enumeration vignette/article execute their numerical guards successfully. The remaining recent failures concentrate at boundaries between families: direct versus grid comparisons, general versus homogeneous transitions, regression versus profile classification, complete versus noise-truncated posteriors, and covariance structure as an extra enumeration dimension.

The recent test reorganization preserves a separate CI equivalence workflow and local execution of 376 independent expectations. Its separation from shipped tests is not itself a finding. The gap is that existing combinations and compatibility checks do not cover the reproductions above. Export presence, argument names and successful examples cannot establish unchanged drawing, abbreviation or help-navigation behavior.

The recent growth introduction and later stabilization have retained meaningful corrections, including stable count densities, better starts, outcome-location residual safeguards and restricted inference paths. MC03 identifies a remaining extreme-scale calculation, with its exposure stated narrowly. Earlier dispatch/floor/weight/pooling problems should not be attributed to the newest changes merely because this audit discovered them now.

## Repair order and acceptance checks

1. Close SV01 and PD01 with desired-behavior regressions before relying on their scientific output.
2. Correct accepted model arguments and structured fitting: MC01 centering, MC02 floor constraints, RA10 covariate enumeration.
3. Make comparison behavior consistent across families and subsets: RA01–RA03, RA05 and RA09. Verify that direct and grid comparisons give the same numbers and reject incompatible data.
4. Align weighted quantities: SV02 residual moments, PD02 tree shares/flows. State the supported inference design explicitly.
5. Correct smaller numerical/reporting edges: SV03, SV04, SV05 and contextual MC03.
6. Preserve structure in RA04 plots, settle PD05 prediction units and PD08 drawing contracts, then restore PD03 abbreviations and PD04 help aliases/usages.
7. Correct RA06/RA08/PD06/PD07/MC04 help claims, clarify RA07's data-equivalence contract, and regenerate affected executable documentation.
8. Rerun focused regressions, the two full suites and the native build/check after fixes. Verify numerical invariants independently rather than only asserting current implementation output.

The core's four invariant tests intentionally **fail on the current implementation**. The service/plot evidence tests intentionally **pass by asserting the observed defects**. Neither group is added to the shipped suite. Convert them into correct-behavior regression tests when implementing repairs; do not interpret passing characterization evidence as validation of the behavior.

## Audit artifacts

- [Function inventory](function-inventory.csv): every definition/alias with static analysis, semantic recency, runtime expression hits and review destination.
- [Source identity](source-files.csv): exact current production source hashes and counts.
- [Public documentation inventory](documentation-inventory.csv) and [supporting documentation inventory](supporting-documentation-inventory.csv).
- [Coverage summary](coverage-summary.txt), [file coverage](file-coverage.csv), [raw expression coverage](line-coverage.csv) and [functions without unit hits](functions-without-unit-hits.csv).
- [Per-function usage checks](usage-check.csv) and [configured lint output](lint.log).
- [Main/equivalence suite log](test-run.log), result CSV/RDS files, [build log](build.log), [native check log](package-check.log) and [session details](SESSION.txt).
- [API findings](recent-api/findings.txt), [plot/documentation findings](plots-docs/review.txt) and [model-core findings](model-core/README.md), each with independent scripts and outputs.
- [Shared-service probes](service-probes.R), their [outputs](service-probes.log), [characterization tests](test-service-findings.R), and [inventory tests](test-inventory.R).

`run-tests.R`, `inventory.R`, `run-coverage.R` and `summarize-coverage.R` reproduce the repository-level evidence from the package root. Coverage installs an isolated copy under `tmp/`; build/check artifacts are also there. The first coverage wrapper saved successful coverage but ended on an unsupported summary option; the option was corrected, and summary extraction reused the saved result without rerunning tests. The source and documentation hash checks pass.
