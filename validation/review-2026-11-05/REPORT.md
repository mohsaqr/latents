# Latents changes since CRAN review

**Repair update:** the findings below describe the pre-fix 1.1.16 snapshot. The authorized repairs and final verification for local 1.1.17 are documented in the [completed repair report](../fixes-2026-11-05/REPORT.html).

Reviewed 5 November 2026. The working tree reports version **1.1.16**, on branch `features`, at committed HEAD `a555108111cfb8394dcd2ff47f02b45d7f0db671` plus the existing uncommitted changes.

The recent work adds useful inference, comparison and reporting features, and the latest numerical repairs pass their focused checks. **The current candidate still needs repairs before release:** the documentation evidence check fails, and important combinations continue to return misleading uncertainty, model-selection results or assignments. One still-open confirmed inference defect was introduced after 1.1.9; most remaining numerical defects are inherited.

## Scope and release baselines

There are two release references in the repository. NEWS and `cran-comments.md` describe changes from **0.9.21**, submitted in commit `a10a919b` on 6 October and recorded as reaching CRAN on 16 October. The latest handoff and November review records identify **1.1.9**, commit `8518dffe` on 31 October, as the newer CRAN version. Live CRAN metadata could not be retrieved during this review, so these are explicit repository baselines rather than a claim about the current online listing.

The review covers the complete change inventory from 0.9.21, with particular attention to the five commits after 1.1.9 and the uncommitted 1.1.15–1.1.16 work. It includes untracked help pages and regression tests already present at the start. Existing package changes were preserved; no fixes were applied.

Parsing function definitions, excluding comments and source formatting, gives:

| Comparison to current tree | Added function definitions | Modified existing definitions | Removed definitions | Added exports | Removed exports or S3 registrations |
|---|---:|---:|---:|---:|---:|
| 0.9.21 | 622 | 482 | 6 | 40 | 0 |
| 1.1.9 | 86 | 95 | 0 | 30 | 0 |
| Committed 1.1.14 HEAD | 1 | 4 | 0 | 0 | 0 |

These count top-level named function assignments, including internal helpers; they are not a claim that every definition received a fresh mathematical audit. Since 1.1.9, the only changed formal list among existing exported functions is `three_step()`: it adds `"overall"` as an accepted contrast value and appends `step_one`. The large uncommitted source diff is principally documentation: executable changes are the new comparison-count helper and four existing functions handling comparison counts, held response logits, cross-level warnings and residual adjustment.

Evidence: [all 86 commits since 0.9.21](commit-history.txt), [source audit](source-audit.log), [function inventory](source-function-inventory.txt), and [source manifest](source-manifest.json).

## What changed

| Period | Main changes | Practical consequence |
|---|---|---|
| After 0.9.21 through 20 October | Growth and transition extensions; scaling and empty-row handling; parallel enumeration; direct effects and local dependence; LMR and ML three-step; numerical and inference repairs | Substantially broader model coverage. Likelihoods, criteria, classifications and uncertainty can intentionally differ from 0.9.21. |
| 22–28 October | emEM becomes the default for EII, EEI and EEE; unit-aware stopping and acceleration; shared worker pools and worker recovery; zero-weight and covariance repairs; expanded independent comparisons | Fits can reach better optima. Reproducibility depends on retaining search settings and comparing fitted likelihoods, rather than expecting identical class labels. |
| 29–31 October through 1.1.9 | CRAN preparation and shorter checks; empty-row and zero-weight corrections; supplied-data sensitivity fix; local-independence adjustment; nested analytic LRT in `compare_models()` | Selection and diagnostic behavior changes. Many development tests move out of CRAN execution and rely on CI. |
| 1 November, 1.1.10–1.1.11 | `align_labels()`; categorical distal outcomes and overall tests; covariates across covariance structures; step-one uncertainty; bootstrap interval plots and replicate extraction | Useful new inference paths, including the staged uncertainty defect below. |
| 2 November, 1.1.12–1.1.14 | Stricter nesting, ordinal coding and bootstrap identity checks; count-argument repair; weighted plot repairs; processor-count fallback; named plot/table verbs and compact grid summary | Broader API coverage and stronger validation. The new verbs retain the earlier exported API. |
| Uncommitted 1.1.15–1.1.16 | Individual help pages and rewritten long pages; held categorical SE repair; staged-comparison parameter counts; BH residual default; cross-level criterion warning; dependency and character-encoding cleanup | Three numerical repairs are verified. Help changes break the documentation claim registry, and older tests need reconciliation with intentional behavior changes. |

Users revisiting 0.9.21 results should especially check small Gaussian units, occasion ordering, zero sampling weights, emEM defaults, class relabelling and entropy-based criteria. The current NEWS also correctly calls out staged/joint comparisons, held categorical inference and the new residual adjustment default.

## New regressions and validation failures

### Estimated step one uncertainty is incomplete for staged fits

**High priority.** `R/three-step.R:593` and `605` use only the second-stage free coordinates of a staged fit. `three_step(step_one = "estimated")` therefore omits uncertainty in the held measurement while labelling the result as estimated-step-one inference. This option arrived in **1.1.11**, after 1.1.9.

On a fresh synthetic example, the joint fit supplies a 15 by 15 covariance, while the staged fit supplies only 3 by 3. The staged outcome SE rises from 0.056501 to 0.056599, but measurement coordinates never enter the propagation. The joint fit is a different estimator, so its SE is not an exact reference; the omitted coordinates establish the defect. Refuse this combination until the required first-stage covariance and propagation are implemented. [Evidence and reproduction](inference/REPORT.md).

### The help rewrite breaks the documentation evidence registry

**Release validation failure.** `equivalence/test-claims.R:43` fails **23 expectations**: 23 of 44 registered phrases in `equivalence/CLAIMS.csv` no longer occur in their source. The same current registry has zero missing phrases at committed HEAD, isolating this to the uncommitted help rewrite. The equivalence workflow explicitly runs this check at `.github/workflows/equivalence.yaml:90`.

This is a documentation-to-evidence tracking failure; it does not establish disagreement with external estimates. Reconcile the rewritten prose and registry while retaining the underlying substantive comparisons. [Actual test output](claims-test.log), [missing phrases](claims-text-check.txt).

### Existing tests still encode superseded output

`tests/testthat/test-audit-1101-covariates-diagnostics.R:207` expects default residual `p_adjusted` to equal the raw p-value, contradicting the intentional BH change. `tests/testthat/test-plot-polish.R:146` expects the removed text “returns all of them” in the compact grid summary. These expectations should be updated to assert the intended numerical adjustment and summary behavior. Full-suite outcomes are recorded below.

### Parallel tests assume that processor detection succeeds

This sandbox returns `NA` from `parallel::detectCores()`. The new 1.1.13 fallback correctly uses one core and warns, but existing tests of mocked parallel dispatch still expect multiple cores. A targeted natural run of `test-lta-parallel.R` gives 79 passed expectations, 16 failures, two warnings, two skips and no errors. Eleven failures compare condition lists that now include the fallback warning; five assert task phases or costs that disappear on the serial path. Fitted-object equality checks pass.

With only the reported processor count replaced in memory, the same file gives **99 passed, zero failures/errors/warnings and two skips**. Worker paths exercised in that control remain mocked serial; Accelerate and core-limit skips remain active. This establishes test portability failures for those sixteen expectations, rather than evidence of incorrect transition estimates. Mock processor availability explicitly in dispatch unit tests, and skip live multicore tests where availability is unknown. [Diagnosis and controlled results](engines/lta-parallel-diagnosis.md).

## Remaining numerical and reporting defects

All fifteen major findings left open from the 4 November review were reproduced on current code with fresh synthetic probes. The table gives representative measurements; the linked area reports contain scripts, exact source locations, history and limitations. The staged step-one item is new after 1.1.9; the other fourteen already exist at that baseline.

| Finding and location | Fresh evidence and impact | Recommended response |
|---|---|---|
| Separated `r3step()` logits, `R/three-step.R:1125` | Slope 302.38, default robust SE 10.03 and p about 1e-199, without warning; observed-information SE is about two million. Neither Wald result is trustworthy at separation. | Detect separation or boundary behavior and withhold affected tests. |
| External-measurement staged bootstrap, `R/bootstrap-inference.R:137`, `154` | Held estimate 4.4851 receives interval [1.2876, 1.7242]; 16/16 replicates accepted. Replicates refit measurement on the current population rather than the external one. | Refuse unsupported provenance or implement the intended conditional/external-sample resampling. |
| OPG under grand transforms, `R/inference.R:1165`, `R/covariate-inference.R:595` | Invariant probability SEs change by roughly 8–12% after grand centering; observed and robust controls agree. | Use observed information for transform sensitivity while retaining the selected covariance estimator. |
| Person centering with singletons, `R/centering.R:203`, `315` | Ten singletons become identical zero vectors. BIC favors three profiles, 2919.77 versus 3229.00; removing them favors two, 3155.76 versus 3157.45. | Warn or refuse unsupported groups; evaluate treating singly observed deviations as missing. |
| Weights with grand transforms, `R/centering.R:116`, `192` | Weighted EII raw-scale means differ by 0.29145 from integer-row duplication; the untransformed control differs by about 2e-8. | Use weighted transform moments and verify constrained covariance structures. |
| Unused factor levels, `R/mixture-regression.R:690`, `891` | An unused membership level adds a parameter and BIC penalty log(70) = 4.2485 at equal likelihood; membership SEs become unavailable. | Remove unused design levels or reject rank-deficient designs before counting parameters. |
| Identifier named `class`, `R/mixture-regression-methods.R:283` and growth counterpart | IDs 1001–1070 are replaced with latent labels 1/2 in the executed group-table probe; assignment output has the same collision by source inspection. | Preserve identifiers using collision-safe output columns and joins. |
| Unidentified growth covariance, `R/mixture-regression.R:315` | With two occasions, random intercept, slope and residual variance lie on an exact likelihood ridge; fit reports convergence without warning. | Validate identifiability before reporting successful estimation. |
| Unidentified stayer model, `R/lta-engine.R:1472` and transition setup | With two occasions, reported stayer share 0.4299 and alternatives 0.01, 0.10, 0.20 have the same likelihood; 3/3 starts replicate. | Warn or refuse structurally unidentified specifications. Replication is insufficient. |
| Effectively empty transition origins, `R/lta-engine.R:1370`, `R/transitions.R:972` | Outgoing total 2.92e-7, or 1e-17 and smaller in the general engine, still produces rows marked estimated and charged in parameter counts. | Define an evidence-based effective-support rule; do not adopt an arbitrary 1e-6 threshold without checking it. |
| Categorical boundary diagnostics and LRT, `R/diagnostics-checks.R:68`, `R/compare-lrt.R:48` | Response probability 1e-10; check says no bound/ok, Wald inference refuses, analytic LRT gives a p-value without a boundary warning. | Share a boundary predicate for diagnostics and cautionary inference, accounting for bootstrap refusal semantics. |
| Categorical joint-pair bootstrap LRT, `R/enumeration.R:1102` | All 3/3 replicates fail with “undefined columns selected”; name repair converts `a:b` to `a.b`. | Preserve exact names and test non-syntactic names and partially observed pairs. |
| Ragged general-transition bootstrap, `R/lta-bootstrap.R:184`, `197` | A single longer person can disappear in a resample, causing a signature/shape error that aborts the bootstrap. | Preserve the fitted occasion schema or treat incompatible replicates as counted failures. |
| Cross-level prevalence interpretation, `R/cross-level.R:452` | Full-form table reports 0.5800/0.4200, versus individual posterior shares 0.3256/0.6744. Restricted form reports individual shares. | Name the population represented or report a consistent individual-based quantity. |
| Staged estimated-step-one correction | The new regression described above is also one of the fifteen outstanding findings. | Refuse or implement complete propagation. |

Area evidence: [engines and transformations](engines/REPORT.md), [inference and three step](inference/REPORT.md), [selection and reporting](surface/REPORT.md).

One additional older reporting defect was freshly verified: a cross-level fit with six zero-weight groups reports **60 groups and 1170 rows**, while `nobs()` and the used sample give **54 and 1053**. Likelihood and BIC match the positive-weight subset to about 2e-12. This affects displayed counts, not the fitted likelihood. Location: `R/cross-level.R:588`.

The r3step, external-measurement bootstrap, singleton, unused-level, identifier and reviewed growth/stayer ridge mechanisms exist in 0.9.21 source. Grand scaling, OPG transform handling, cross-level prevalence and several diagnostic paths were added after 0.9.21 but before 1.1.9. The categorical simulation name-repair bug was already present in 0.9.21; later joint-item support exposes it through a newer feature. All eight engine counterexamples were rerun against isolated 1.1.9 source with byte-identical logs.

## Repairs that are verified

The uncommitted numerical fixes are supported by fresh checks:

- **Held categorical covariate inference:** SEs 0.1275603 and 0.1401586 agree with an independently reconstructed likelihood and numerical Hessian. `vcov()` is correctly 2 by 2; `confint()` has two rows.
- **Staged versus joint comparison:** criteria count the measurement parameters when required. Tests also protect the conditional comparison when both fits hold the same measurement and leave joint-only comparisons unchanged.
- **Residual adjustment:** `get_results(fit, "residuals")` and `get_residuals()` now use BH by default; explicit `adjust = "none"` returns raw values.
- **Cross-level selection warning:** `latents_cross_level_criteria` now warns against selecting class count or full/restricted form from its criteria. This mitigates misuse; it does not make the duplicated-data likelihood suitable for that selection.
- **Bootstrap identity protection:** reordered parameter tables correctly lose access to their original replicate matrix through a classed refusal.

Focused package tests pass **237 expectations** with no failures, errors, warnings or skips: 23 latest-fix checks, 207 wildcard checks and seven help-page checks. Independent controls also reproduce a Gaussian likelihood by hand, the small-unit likelihood Jacobian to 5.46e-12 and one-class regression ML likelihood against `lm()` to 1.14e-13.

The previous summary listed 19 major findings. Three are numerically repaired (held categorical SEs, staged comparison counts, residual BH default), one is mitigated by a warning (cross-level criteria), and **15 remain open**. The old handoff's “14” understates that remainder; its “five fixes” also counted check cleanup outside those nineteen.

## Validation status

The full development suite completed across **254 files and 2064 test cases**, using `NOT_CRAN=true`:

| Run | Passed expectations | Failed expectations | Errors | Warnings | Skipped cases |
|---|---:|---:|---:|---:|---:|
| Full package suite | 14692 | 38 | 1 | 36 | 28 |
| Focused latest fixes, wildcards and help contracts | 237 | 0 | 0 | 0 | 0 |
| Documentation claim registry | 156 | 23 | 0 | 0 | 0 |
| Controlled LTA dispatch test | 99 | 0 | 0 | 0 | 2 |
| Supplemental development-source checks | 153 | 0 | 0 | 0 | 0 |

These runs overlap; their totals should not be added. The full suite is **not green**. Two failures are the stale BH and summary-text expectations described above. The other 36 failures and the one error are in parallel-dispatch expectations affected by unknown processor availability: 16 LTA, 14 profile-pool, five general-parallel and one sensitivity expectation. The error is `calls[[1L]]` at `test-parallel-profile-pool.R:247`, after serial fallback leaves the mocked parallel-call list empty. The controlled rerun isolates the LTA subset; the remaining classification follows the recorded traces and source paths rather than independent controlled reruns of every affected file.

All 36 warnings belong to the two intentional conditions: 26 cross-level criterion warnings and ten unknown-core warnings. The 28 skips comprise 15 Accelerate fork exclusions, four unavailable socket-port checks, three core-limit exclusions, one unknown-core exclusion and five checks whose development-study sources were initially absent from the isolated package copy. With those sources restored, the supplemental run passes 153 expectations, including **97 expectations in the five previously skipped cases**, with no failures, errors, warnings or skips. The other 23 exclusions remain limitations of this environment. [Supplemental output](supplemental-validation.log).

Evidence: [full results](full-suite-results.csv), [results by file](full-suite-by-file.csv), [failure details](full-suite-failures.log), [error trace](full-suite-error-trace.log), [warning and skip summary](full-suite-conditions-summary.log), and [complete console log](full-suite.log). The console abbreviates its failure display; the structured results retain every file. These test failures are separate from the fifteen numerical findings demonstrated by the synthetic probes.

Package installation into a private temporary library passed. All R files and 77 Rd files parsed and passed the Rd checks. The source/test/documentation manifest showed no drift during the review. Tests run from an isolated source copy, with an installation of that same copy placed first in `R_LIBS` for worker loading.

No new full `R CMD check --as-cran`, external-software comparison matrix, coverage simulation, website build or vignette rendering was completed. The earlier independent-comparison evidence remains useful but does not certify the changed source. Prior minor and edge findings in [the 4 November review](../review-2026-11-04/SUMMARY.md) remain a separate backlog unless reproduced here.

## Recommended order of work

1. Refuse incomplete staged step-one inference and prevent unguarded separation tests and mismatched external-measurement bootstrap intervals.
2. Reconcile the documentation registry and stale tests, then obtain a clean current-source suite and package check.
3. Repair weighted transforms, singleton handling, unused levels and identifier collisions; add meaningful counterexample regressions.
4. Add identification and boundary checks, then repair categorical and ragged transition bootstrap paths.
5. Settle cross-level prevalence semantics and effective transition-support rules. Retain the new cross-level criterion warning.
6. Rerun affected independent comparison kits, followed by the broader equivalence matrix and rendered articles, before publishing a release.

Reproduction scripts, structured results and logs are retained beside this report. Production code and existing tests remain unchanged by this review.
