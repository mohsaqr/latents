# Latents code review for October 30 to November 1 2026

The review identified **ten issues: four high priority and six medium priority**. The high priority findings concern misleading numerical results: likelihood-ratio testing accepts incompatible predictor values or ordinal category orders as nested models, uses an unjustified chi-square reference for a singular covariance comparison, and a reordered bootstrap table can attach replicates to the wrong parameter labels. These should be repaired before releasing the new functionality.

The independent likelihood checks and retained external comparisons provide substantial support for the central estimation and uncertainty calculations. Those checks do not cover the data identity, table manipulation and boundary cases identified below.

Only review material was added. Package implementation, existing tests, version numbers and Git history were left unchanged.

## Scope

**Revision:** `4c8f5aa70873b54b648f6c7f184333d3ac4d94ae`, branch `features`, version `1.1.11`. **Comparison base:** `23a3de23`, the last commit before October 30. The primary period is October 31 and November 1; October 30 is included in full to avoid excluding work at the two-day boundary. The working tree was clean when the review began.

The comparison contains nine commits and 590 changed paths: 22 runtime R files, 25 equivalence R or shell scripts, seven validation R scripts, 107 test files, 26 documentation or metadata paths and 403 reference or validation evidence artifacts. Runtime changes total 1,638 added and 180 deleted lines. The complete inventory is [changed-files.csv](changed-files.csv).

All changed runtime code and substantive changes to tests, validation, comparison harnesses and documentation were read. Surrounding implementations were traced where they determine the new behavior. Recorded data tables and binary golden files were treated as evidence artifacts, rather than counted as independently reviewed executable code. Ignored scratch work and local planning files are outside the reviewed implementation.

| Commit | Date | Change |
| --- | --- | --- |
| `ab5bcef8` | October 30 | Preserve zero-weight classifications |
| `243503af` | October 30 | Allow platform-dependent reference differences |
| `0d0775bd` | October 30 | Sensitivity refits with supplied single-level data |
| `311f041e` | October 30 | Skip slower test files on CRAN |
| `0f79c7ed` | October 30 | Local-independence wording and harness evidence |
| `c979dce9` | October 31 | Likelihood-ratio tests in model comparison |
| `8518dffe` | October 31 | Consolidated release notes |
| `1b17adde` | November 1 | Alignment, categorical outcomes, covariate structures and BH correction |
| `4c8f5aa7` | November 1 | Step-one uncertainty, bootstrap accessors and plots |

## Findings by priority

| ID | Priority | Defect | Consequence |
| --- | --- | --- | --- |
| F1 | High | Shared predictor values are absent from LRT nesting checks | Invalid comparisons return highly significant p-values |
| F2 | High | Bootstrap replicate identity relies on row names | Sorting silently mislabels replicate values |
| F3 | High | Ordinal category order is absent from LRT nesting checks | Different ordinal models receive a nested-model p-value |
| F4 | High | Chi-square reference is used at a singular covariance restriction | Nominal 5% tests reject 8% in the targeted null study |
| F5 | Medium | Fixed inference and confidence intervals disagree about rows | `confint()` crashes |
| F6 | Medium | Slope nesting uses the wrong option name | Valid group-specific slope comparisons crash |
| F7 | Medium | Invalid BCH probability policy is bypassed for pairs | Intervals and tests contradict the documented policy |
| F8 | Medium | Alignment omits ordinal category definitions | Incompatible codings are accepted |
| F9 | Medium | Plot spread detection depends on measurement units | Nonconstant parameters are classified as constant |
| F10 | Medium | Categories are selected before excluded units are removed | An excluded row eliminates all categorical ML errors |

High priority means a supported operation silently returns materially misleading numerical output. Medium priority means a reproducible failure of supported behavior or an inconsistency affecting inference or interpretation.

The nine deterministic findings have runnable fixtures in [reproduce.R](reproduce.R), with observed output in [reproduce.log](reproduce.log). F4 has a separate 500-dataset calibration script and retained results. The script catches errors so subsequent probes still run; its zero process exit status does not mean the defects are absent.

From the repository root, reproduce them with:

```sh
OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 Rscript validation/review-2026-11-01/reproduce.R
```

### F1 LRT accepts different values for a shared predictor

**Location:** [R/compare-lrt.R line 184](/Users/mohammedsaqr/Documents/Github/latents/R/compare-lrt.R:184), especially the membership checks at line 190. Introduced in `c979dce9`.

The nesting check requires the null's covariate names to be included in the alternative's covariate names. It does not require shared covariates to contain the same values. The comparison signature checks indicators, grouping, weights and transformations, but omits membership predictor values.

A model using shuffled `z` can therefore be compared with a model using the original `z` plus `w`. Setting the coefficient of `w` to zero cannot recover the first model's design. Nevertheless, both converged fits in the 600-row reproduction pass the public checks:

```text
compare_models(small, large, test = "lrt")
LR = 222.0624
df = 1
p = 3.210245e-50
```

This is an invalid nested-model test of adding `w`: the models also replace the values of `z`. The extremely small p-value is the material failure.

**Repair:** Check retained values or design columns of every shared predictor before accepting nesting. Additional columns must retain the null's design. Audit common group covariates, direct-effect predictors and transition covariates for the same omission; the reproduced failure here is specifically profile membership.

**Regression checks:** Reject a shuffled shared predictor with a classed error. Accept adding a predictor when shared columns are unchanged. Cover group-level predictors and equivalent coding cases.

### F2 Sorting a bootstrap table silently corrupts replicate labels

**Location:** [R/bootstrap-inference.R line 300](/Users/mohammedsaqr/Documents/Github/latents/R/bootstrap-inference.R:300); incorrect pairing occurs at line 323. Introduced in `4c8f5aa7`.

The replicate matrix retains its original parameter order. Its guard checks only the table's row count and saved row names. A table can be reordered while retaining or regenerating those ordinary row names. The matrix then remains in its old order while the accessor reads labels from the reordered table.

Both operations reproduce the problem:

```r
changed <- boot[rev(seq_len(nrow(boot))), , drop = FALSE]
rownames(changed) <- NULL
as.data.frame(changed, what = "replicates")

changed <- dplyr::arrange(boot, dplyr::desc(estimate))
as.data.frame(changed, what = "replicates")
```

The guard accepts both tables. Standard deviations calculated from the exported, supposedly labelled replicates no longer reproduce the matching table standard errors. In the six-resample Faithful fixture, the maximum discrepancies are **3.693022** for the base R operation and **4.768007** for `dplyr::arrange()`.

This corrupts downstream analyses of exported replicates. Plotting also trusts the same guard when deciding which parameters have spread. The helper is shared by the profile, transition and growth bootstrap tables.

**Repair:** Store the original identities of the parameters, including every label column, and validate their actual ordered values before using the matrix. Alternatively, map replicate columns explicitly by stable parameter identities. Row names cannot establish correspondence. Preserve the intended refusal of modified tables unless reordering and subsetting are fully supported.

**Regression checks:** Include `arrange()`, reversing rows followed by resetting row names, repeated rows, changed parameter labels and column removal. A permitted table must reproduce its errors and quantiles from its labelled replicates; unsupported modifications must raise a classed error.

### F3 LRT accepts incompatible ordinal category orders

**Location:** [R/compare-lrt.R line 153](/Users/mohammedsaqr/Documents/Github/latents/R/compare-lrt.R:153), particularly the `ordinal` entry in the settings at line 157. Introduced in `c979dce9`.

The settings check records the names of ordinal indicators, but not their category order. The comparison signature converts ordinal observations to character labels and therefore also discards their factor level order. Two fits with the same response labels but different ordinal codings can pass both checks.

For an ordinal model, a nonmonotonic reordering changes the class-specific response constraints. The [ordinal likelihood](/Users/mohammedsaqr/Documents/Github/latents/R/extra-indicators.R:196) assigns successive categories scores 0, 1 and 2. Swapping the middle and last categories changes their scores to 0, 2 and 1, which cannot be absorbed by a linear change of the class location. Adding a membership coefficient does not recover the original measurement restriction.

The reproduction fits a two-profile model with ordinal levels `low`, `middle`, `high`, then a model with an added membership covariate and levels `low`, `high`, `middle`. The response values themselves are unchanged. The private nesting checker returns `NULL`, meaning accepted, and the public comparison reports:

```text
LR = 225.0226
df = 1
p = 7.25905e-51
```

This p-value combines a changed measurement model with the added covariate while treating the comparison as a one-parameter nested restriction. The missing coding check is separate from F1's missing shared-predictor values and from F8's public alignment compatibility check.

**Repair:** Include `extra_data$ordinal_levels`, with category order preserved, in profile-model nesting compatibility. The transition nesting helper already checks those definitions and provides a useful consistency reference.

**Regression checks:** Reject a changed ordinal order, including a nonmonotonic permutation with the same observed labels. Continue to accept added membership covariates when ordinal coding is unchanged. Confirm the same requirement for mixed and all-ordinal measurement models.

### F4 A singular covariance comparison receives a chi-square reference

**Location:** [R/compare-lrt.R line 46](/Users/mohammedsaqr/Documents/Github/latents/R/compare-lrt.R:46), together with the covariance inclusion rule at [line 222](/Users/mohammedsaqr/Documents/Github/latents/R/compare-lrt.R:222). Introduced in `c979dce9`.

The code accepts covariance restrictions by their nesting and computes a chi-square p-value using the parameter-count difference. That is insufficient for every accepted structure pair. In particular, `EII` within `EVE` is singular at the spherical null: the alternative's common orientation has no identified value when both profiles are spherical. The alternatives near that point do not form the regular linear parameter extension assumed by Wilks' reference.

This is distinct from changing the number of profiles, and from variances reaching `min_variance`. Neither the class-count refusal nor the variance-bound warning detects it. General likelihood-ratio theory explicitly distinguishes singularities from regular restrictions; see [Drton, Likelihood ratio tests and singularities](https://arxiv.org/abs/math/0703360). This application's singularity diagnosis is an inference from the covariance model's parameterization, checked by the targeted experiment below.

The new [lrt-orientation.R](lrt-orientation.R) study generated 500 datasets under the EII null, each with 1,000 observations, two Gaussian indicators and two profiles separated by 12 standard deviations. This separation keeps class overlap negligible. Each pair was fitted and tested through the public API. All pairs returned tests; no fit was on a variance bound. Results were:

```text
Reported degrees of freedom: 3
Mean LR statistic: 3.687609       (chi-square(3) has mean 3)
Rejection at nominal 5%: 0.080    (40 of 500 datasets)
P-value uniformity KS D: 0.1582
KS p-value: 2.706e-11
```

The study does not establish calibration for other pairs or quantify every singular case. It supplies a concrete counterexample to the current blanket interior-restriction claim, consistent with the nonidentified orientation at the null. Related work on Gaussian covariance-model tests also cautions against assuming chi-square calibration throughout the family: [Punzo, Browne and McNicholas, Hypothesis Testing for Parsimonious Gaussian Mixture Models](https://arxiv.org/abs/1405.0377).

**Repair:** Distinguish regular covariance comparisons from singular ones. Refuse the ordinary chi-square test for unsupported singular pairs, or calibrate those comparisons with a suitable parametric bootstrap. Revise the help's statement that all listed restrictions have an ordinary chi-square reference. Do not infer validity solely from equal class counts, an interior variance estimate and a positive parameter-count difference.

**Regression checks:** Include the spherical-to-common-orientation case in regularity refusals or bootstrap tests. Add retained calibration evidence for the supported comparison categories. Agreement with another program on the LR statistic and the arithmetic chi-square p-value does not validate the reference distribution.

Run the targeted study from the repository root with `Rscript validation/review-2026-11-01/lrt-orientation.R 500`. Its full results are in [lrt-orientation.csv](lrt-orientation.csv), and the summary is in [lrt-orientation.log](lrt-orientation.log).

### F5 Confidence intervals crash for fixed covariate models

**Location:** [R/covariate-inference.R line 749](/Users/mohammedsaqr/Documents/Github/latents/R/covariate-inference.R:749); the incompatible consumer is `confint.multilpa_covariates()` at [line 1107](/Users/mohammedsaqr/Documents/Github/latents/R/covariate-inference.R:1107). Introduced by fixed-model support in `1b17adde`.

The inference frame now removes held measurement parameters. `confint.multilpa_covariates()` still builds names for every natural parameter and assigns them to the shorter interval matrix.

The reproduced three-profile VVI two-step model has four free membership coefficients and twelve held measurement parameters. `parameter_inference()` correctly returns four rows. `confint()` gives those rows sixteen names and fails:

```text
length of 'dimnames' [1] not equal to array extent
```

The failure occurs for the default call and before `parm` selection. The numerical membership standard errors themselves are not implicated by this reproduction.

**Repair:** Derive interval names from the rows actually returned by inference. Default to estimated parameters and explicitly refuse requested held parameters with `latents_held_parameter`, consistently with the covariate-free method. Keep coefficient and covariance naming aligned with that policy.

**Regression checks:** Default, named and numeric `parm` requests with `fixed = "measurement"` and partially held blocks. Estimated coefficient intervals must match `parameter_inference()`; held-parameter requests must give the intended classed refusal.

### F6 Group-specific slope LRTs use the wrong option vocabulary

**Location:** [R/compare-lrt.R line 194](/Users/mohammedsaqr/Documents/Github/latents/R/compare-lrt.R:194). Introduced in `c979dce9`.

The fitter records `profile_slopes = "shared"` or `"group_class"`. The nesting code ranks values against `c("shared", "varying")`. Matching `"group_class"` yields `NA`, which enters an `if` condition.

Two converged two-level fits with a profile membership covariate reproduce the failure when shared slopes are tested against slopes by group class:

```text
simpleError: missing value where TRUE/FALSE needed
```

**Repair:** Use `c("shared", "group_class")` consistently with the public argument and stored field. Reject unknown values explicitly rather than allowing `NA` into control flow.

**Regression checks:** Shared-to-group-class comparisons must succeed when otherwise nested. Reverse comparisons must raise `latents_bad_nesting`. Comparing two group-class slope fits under another valid restriction must also work.

### F7 BCH pairwise tests ignore the invalid-probability policy

**Location:** [R/three-step-categorical.R line 126](/Users/mohammedsaqr/Documents/Github/latents/R/three-step-categorical.R:126), after invalid probabilities are identified at line 91. Introduced in `1b17adde`.

The help says BCH probabilities outside `[0, 1]` are reported without an interval or test. The individual table suppresses intervals, and the overall test suppresses its statistic. The pairwise path passes the invalid cells to its helper without a mask.

In the existing out-of-range fixture from the new tests, class 1 probabilities are **1.01069914** and **-0.01069914**. Individual intervals are `NA`, but comparisons involving those probabilities return:

```text
Difference for the rare category = 0.2364026
95% CI = [0.07550661, 0.39729855]
p = 0.003979877
```

Even the warning says the probabilities are reported without an interval. This finding concerns the package's explicit policy; it does not assert that every unbounded BCH contrast lacks an asymptotic interpretation.

**Repair:** Pass the invalid-cell mask into pairwise inference and withhold intervals, statistics and p-values for affected contrasts, preserving unaffected contrasts where intended. If pairwise inference is deliberately allowed, document that exception and revise the warning instead.

**Regression checks:** Add `contrast = "pairs"` to the existing out-of-range fixture. Include partially affected class-category cells and both `step_one` choices where supported.

### F8 Ordinal alignment accepts incompatible category definitions

**Location:** [R/align-labels.R line 56](/Users/mohammedsaqr/Documents/Github/latents/R/align-labels.R:56). Introduced with the public wrapper in `1b17adde`.

The compatibility layout checks ordinal indicator names but omits the ordered category definitions in `extra_data$ordinal_levels`. The signature then compares category probabilities by position.

The reproduction fits the same named ordinal indicator with levels `low`, `medium`, `high` and with levels `red`, `green`, `blue`. `align_labels(second, first)` accepts the pair, despite the documented requirement that indicators be coded the same way. Equal category counts do not establish equal coding.

**Repair:** Include the complete ordered category definitions in compatibility checks. Name incompatible indicators in a `latents_bad_argument` refusal.

**Regression checks:** Refuse different labels, reversed order and different category counts; accept exact coding matches. Exercise all-ordinal and mixed-indicator models.

### F9 Bootstrap plots treat small units as zero variation

**Location:** [R/bootstrap-inference.R line 360](/Users/mohammedsaqr/Documents/Github/latents/R/bootstrap-inference.R:360). Introduced in `4c8f5aa7`.

Spread is defined by `diff(range(values)) > 1e-10 * max(1, abs(values))`. The floor of one makes this an absolute `1e-10` cutoff for small-valued parameters. It can suppress nonconstant parameters solely because of their measurement units. Dependence on absolute magnitude also makes the rule vulnerable to large additive origins.

Rescaling both Faithful indicators by `1e-12` and requesting six bootstrap resamples gives nonzero mean replicate ranges:

```text
9.898963e-14, 7.471107e-13, 6.718518e-14, 2.554754e-12
```

Nevertheless, `plot(boot, what = "mean")` raises `latents_nothing_to_plot` and claims no chosen parameter varies. Estimates and percentile intervals are present.

**Repair:** Identify structurally constant parameters from model information, or use a numerical equality rule tied to representational precision. Do not decide whether uncertainty exists using a fixed threshold in measurement units.

**Regression checks:** Equivalent rescalings must retain the same parameter identities and panel structure, with scaled estimates and intervals. Include a large location shift and a genuinely constant probability.

### F10 An excluded row removes every categorical ML standard error

**Location:** [R/three-step.R line 528](/Users/mohammedsaqr/Documents/Github/latents/R/three-step.R:528), before kept units are selected at line 564. Introduced in `1b17adde`.

Categories are built from the full outcome column before `pieces$kept` removes unclassified units. A category seen only in excluded rows therefore enters step three as an all-zero indicator column. Its estimated probability is zero in every class, and the ML boundary guard suppresses uncertainty for the entire outcome table.

A binary-outcome fixture initially has four finite ML standard errors. Appending one row with all measurement indicators missing and a new outcome category called `excluded` changes that to:

```text
Finite standard errors before: 4
Finite standard errors after:  0
```

The row is excluded by the default classification policy. Its outcome label still changes the analyzed category space and disables the retained categories' intervals and tests.

**Repair:** Determine observed categories from the units actually entering step three, after alignment, level aggregation and kept-unit selection. Preserve factor order within those categories, recode consistently and perform the one-category check on the same analysis population.

**Regression checks:** An excluded observation or wholly excluded group carrying an otherwise unseen outcome category must leave retained estimates and inference unchanged. Cover weighting and ML paths.

## Verification

The current package installed successfully into an isolated library under `tmp/review-2026-11-01/library`. The user's R library was not modified. Runs used R 4.5.2 on Apple Silicon macOS with Accelerate BLAS and single-threaded numerical libraries. The full unit run sets `NOT_CRAN=true` and the two-core check configuration.

The five selected external comparison and claim suites completed with **250 passing expectations, zero failures, zero errors, zero warnings and zero skips**:

| Suite | Passing expectations | Coverage |
| --- | --- | --- |
| Documentation claims | 179 | Registered wording and referenced tests |
| Categorical distal outcomes | 6 | Input identity and retained Latent GOLD comparisons |
| Step-one uncertainty | 8 | Input identity and retained Latent GOLD correction comparisons |
| Likelihood-ratio tests | 2 | 27 recorded numerical comparisons with mclust, glca and Mplus |
| Covariate structures | 55 | Fits in all fourteen structures against MoEClust |

Expectations are not counts of independent numerical quantities: several check a complete comparison table. Commercial-program checks read saved outputs; Latent GOLD and Mplus were not rerun. See [equivalence-results.csv](equivalence-results.csv) and [equivalence-tests.log](equivalence-tests.log).

Verification safeguards also passed: **30 R tests with 152 expectations**, including the new optional-difference register behavior, and **four Python tests**. See [verification-safeguards.log](verification-safeguards.log) and [verification-python.log](verification-python.log).

The full unit suite completed across **233 files and 1,914 test cases** in about 31 minutes:

| Result | Count |
| --- | --- |
| Passing expectations | 13,184 |
| Failed assertions | 0 |
| Errors | 6 |
| Warnings | 2 |
| Skipped test cases | 23 |

All six errors are worker startup failures because a local server socket could not be opened: five in `test-parallel.R` and one in `test-search.R`. The isolated [parallel recheck](parallel-recheck.log) reproduced the five parallel errors on another port. They limit verification of live parallel execution here; they are not attributed to these commits. The skipped cases concern Accelerate fork safety, socket availability, the two-core check limit or unavailable machine-core detection.

The two warnings both come from the new alignment test fixture: its likelihood retains a score worth 1.6% of a standard error when result tables request Wald inference. Tighten that fixture's fit tolerance if a warning-free run is required. The alignment assertions still pass; the warning is retained rather than hidden.

The changed feature files include 172 passing expectations for step-one uncertainty, 50 for bootstrap tables, 34 for LRTs, 102 for alignment, 16 for fixed covariate models, 65 for constrained covariate inference and 35 for categorical three-step outcomes. The findings above exercise cases missing from those passing tests. Full outcomes are in [unit-results.csv](unit-results.csv) and [unit-tests.log](unit-tests.log).

The archived step-one coverage study meets its stated 0.92–0.98 criterion in **47 of 64 corrected settings**, versus **40 of 64 fixed-step-one settings**. Some poorly separated, three-class, small-sample settings still have coverage near 0.83–0.89. The new help already acknowledges this; it is not counted as a newly discovered defect. The scripts calculate coverage over runs with an available interval, so the figures are not an unconditional coverage guarantee. See [environment-and-coverage.txt](environment-and-coverage.txt).

These are local regression results on one platform. The archived large simulation studies were not repeated. A separate targeted 500-dataset null study was run for F4. The review did not regenerate the full golden catalogue, rebuild the website or perform a new release-wide `R CMD check --as-cran`.

## Review coverage

| Runtime area | Review focus | Assessment |
| --- | --- | --- |
| `compare-lrt.R`, `compare-models.R` | Data checks, nesting, degrees of freedom, reversed likelihoods, estimator restrictions and ordering | F1, F3, F4 and F6; external numerical comparisons pass |
| `align-labels.R`, `alignment-service.R` | Compatibility, signatures, permutations, rebased logits and Gaussian diagnostic rows | F8; substantial permuted-fit regression coverage |
| `covariates.R`, `covariate-inference.R`, `multilpa.R` | Structure forwarding, held blocks, counts, chart scores and delta maps | F5; independent likelihood and MoEClust checks support new structure paths |
| `three-step.R`, `three-step-categorical.R`, `three-step-ml.R` | Coding, weighting, EM, information, contrasts, clustering and correction | F7 and F10; retained Latent GOLD comparisons pass |
| `inference.R` | Shared coordinate scales and inference integration | No additional reproduced defect |
| `bootstrap-inference.R`, `lta-bootstrap.R`, `growth-mixture-parameter-inference.R` | Table construction, replicate export, identity and plot layout | F2 and F9 |
| `noise.R`, `accessors.R`, `predict.R`, `lta-methods.R`, `get-data.R` | Reporting versus counting masks and supplied single-level refit identifiers | No additional reproduced defect in changed paths |
| `report.R`, `diagnostics-checks.R`, `conditions.R` | BH adjustment, check wording and condition documentation | BH is explicitly passed; wording preserves check computation |
| Tests and development scripts | Independent references, optional differences, provenance, simulations and sabotage additions | Useful numerical evidence; CI registration gap below |

## Additional observations

**The new external suites need CI registration.** The new tests `lrt`, `latentgold-categorical-distal`, `latentgold-step-one` and `moeclust-covariate-structures` are in `equivalence/`. Current workflow commands run `tests/equivalence/`, a fixed retained-reference expression and ten named package kits. Those expressions and the kit allowlist omit these four suites. The full-output harness discovers its own suite directory, not these test files. Add explicit workflow coverage, dependencies and failure-on-skip checks. See [.github/workflows/equivalence.yaml line 89](/Users/mohammedsaqr/Documents/Github/latents/.github/workflows/equivalence.yaml:89) and [equivalence/run-package-kit-tests.R line 29](/Users/mohammedsaqr/Documents/Github/latents/equivalence/run-package-kit-tests.R:29).

**CRAN skips do not by themselves remove GitHub unit coverage.** Of the current 233 unit test files, 193 have a top-level `skip_on_cran()` call; 90 changed files add only that skip and a blank line. However, `setup-r@v2` exports `NOT_CRAN=true` when unset, so the ordinary GitHub check retains those tests. This suspected issue was investigated and excluded from the findings. [Upstream setup implementation](https://raw.githubusercontent.com/r-lib/actions/v2/setup-r/src/installer.ts#L777).

**Some descriptions lag the expanded API.** The `three_step()` interpretation still says maximum likelihood additionally assumes a normal outcome distribution, although categorical ML now uses a multinomial distribution. Restrict that sentence to numeric outcomes. The opening commentary in `bootstrap-inference.R` describes constrained-structure Wald inference as unavailable, although current charts support those structures. These are low-priority clarity issues.

**The categorical local-independence calibration remains a limitation.** The newly stored [bivariate residual simulation](../bivariate-residual-size.csv) marks all four LCA settings as failing its declared calibration criterion. Their per-test rejection rates range from 0.00067 to 0.00675 at a nominal 0.05, and none of the simulated datasets is flagged after BH adjustment. The four LPA settings meet that study's criterion. Passing BH-adjusted p-values through `diagnostics()` fixes the reported multiplicity error, but does not establish calibration for categorical residual tests. Treat the LCA results as an unresolved validation limitation when interpreting an `ok` check; investigate calibration or present those residuals descriptively. This is archived simulation evidence about an existing calculation, not an additional reproduced regression from these commits.

## Repair order

1. Repair shared-predictor and ordinal-coding identity in LRT nesting, the reference distribution for singular covariance comparisons, and parameter identity in bootstrap tables. These silently return misleading numerical results.
2. Repair fixed-model `confint()` and the slope option mismatch. Exercise complete public workflows.
3. Select categorical outcome support from analyzed units and apply the invalid-probability policy consistently across contrasts.
4. Add ordinal coding checks and make plotting invariant to units.
5. Register external suites in CI, add focused regressions and rerun affected checks followed by the full unit suite.

Retain the independent numerical checks. The uncovered cases mostly arise where tested components meet: model designs, labels and stored matrices, held rows and generic methods, or the full input and effective analysis population.
