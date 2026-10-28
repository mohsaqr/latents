# latents 0.96.1 release verification

**Ready for release: all configured verification gates passed.** Version 0.96.1 is committed and pushed to `features` at `d4c9f699525d9e07959463010e004b6362dcdaaf`. All six full platform checks, equivalence and lint passed on this exact commit. No tag, main-branch merge, GitHub release or CRAN submission has been performed.

## Changes and review findings

Version 0.96.1 repairs all ten findings in the [original review](../review-2026-10-20/REVIEW.md). The [repair report](../fixes-2026-10-20/FIXES.md) preserves each reproduction, correction and regression result. The changes cover observed positive-weight unit counts in robust inference; generated categorical pair collisions and native-category validation; weighted classification trees; direct-effect predictor origins; ML distal-outcome origins; conditioned VLMR information; stable independent mover–stayer references; and negative or genuinely small-weight Imhof distributions.

Further release testing found and corrected numerical and platform defects. No dependency or exported function was added. Retained external fixtures, reference tolerances and the EM likelihood-decrease guard were preserved. One unit-test fixture was replaced with a controlled interior fit after its old data proved to favor a boundary optimum; the existing assertions and tolerance remain intact.

| Commit | Release change |
| --- | --- |
| [`c4c1df2`](https://github.com/mohsaqr/latents/commit/c4c1df23a0083dd6b084714b3cae7b749587414f) | Reviewed repairs, four new regression files, version/NEWS and retained-reference CI coverage |
| [`a5f211f`](https://github.com/mohsaqr/latents/commit/a5f211f46fb93f4290032775e3d47af48f8f1717) | Reconciled documentation-claim registrations and checked them in CI |
| [`a40d5c7`](https://github.com/mohsaqr/latents/commit/a40d5c7810e616b3d162972506c2e0f43d82b845) | Numerically certified covariance-floor update, stable empirical-Gram Vuong calculation, portable parallel regressions and early platform gates |
| [`c065dec`](https://github.com/mohsaqr/latents/commit/c065deca5bddda3ee4d504965a903d082207cc48) | Explicit Bash installation step on Windows, avoiding the PowerShell `R` alias |
| [`cf1127f`](https://github.com/mohsaqr/latents/commit/cf1127fa66283826ed421d028999b5ac23e1901d) | Independently reviewed per-node socket cleanup, failure preservation and leak regressions |
| [`d4c9f69`](https://github.com/mohsaqr/latents/commit/d4c9f699525d9e07959463010e004b6362dcdaaf) | Controlled interior transition fixture, stronger eligibility assertions and early engine-regression coverage; no production-code change |

## Final-commit GitHub gates

All three workflow `headSha` values match `d4c9f69`. The six package-check jobs cover macOS and Windows R-release and Linux R-release, R-devel, R-oldrel-1 and R 4.1. Every matrix job installs the current package for socket workers and runs the sensitive LMR, parallel and covariance regressions before the full package check.

| Gate | Run | Result |
| --- | --- | --- |
| Six-platform package matrix | [R-CMD-check](https://github.com/mohsaqr/latents/actions/runs/39456155958) | Passed: all six full checks; zero errors and warnings |
| Development equivalence and retained inference/claim references | [equivalence](https://github.com/mohsaqr/latents/actions/runs/39456160271) | Passed: 405 + 390 assertions; zero failures/warnings/skips |
| Package lint | [lint](https://github.com/mohsaqr/latents/actions/runs/39456165092) | Passed; no lints found |

Machine-readable job records and complete logs are retained under [ci-final-d4c9f69/](ci-final-d4c9f69/).

| Platform | R | Errors | Warnings | NOTEs | Passed assertions | Skipped cases |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| [macOS release](https://github.com/mohsaqr/latents/actions/runs/39456155958/job/118743461680) | 4.6.1 | 0 | 0 | 0 | 8,136 | 4 |
| [Windows release](https://github.com/mohsaqr/latents/actions/runs/39456155958/job/118743461986) | 4.6.1 | 0 | 0 | 1 | 8,124 | 9 |
| [Linux devel](https://github.com/mohsaqr/latents/actions/runs/39456155958/job/118743462041) | devel r90719 | 0 | 0 | 0 | 8,133 | 4 |
| [Linux release](https://github.com/mohsaqr/latents/actions/runs/39456155958/job/118743462229) | 4.6.1 | 0 | 0 | 0 | 8,133 | 4 |
| [Linux R 4.1](https://github.com/mohsaqr/latents/actions/runs/39456155958/job/118743462311) | 4.1.3 | 0 | 0 | 2 | 8,090 | 6 |
| [Linux oldrel-1](https://github.com/mohsaqr/latents/actions/runs/39456155958/job/118743462399) | 4.5.3 | 0 | 0 | 0 | 8,133 | 4 |

All platform test runs report zero failures and unexpected warnings. Windows reports one temporary-directory NOTE for `Rscript123023e62f`.

R 4.1 reports zero errors and warnings with two NOTEs: optional `mice` is intentionally unavailable in the minimum-version job, and its installed package size is 8.1 MB (R 3.7 MB, documentation 3.0 MB, help 1.0 MB). The completed macOS and other Linux checks report zero NOTEs.

The complete local source-unit run has zero skips. Installed-package CI omits source-inspection cases and the Latent GOLD kit excluded from the built package; R 4.1 additionally omits optional `mice` tests, and Windows omits fork-specific tests. The exact per-platform passing counts and skip reasons are retained in the [platform summary](ci-final-d4c9f69/platform-summary.csv). The table reports each platform's actual counts. Feature-branch pushes do not automatically run these workflows; all three were explicitly dispatched against the final commit.

## Fresh local validation

| Check | Result | Evidence |
| --- | --- | --- |
| Complete development unit suite | 8,144 assertions, 1,224 cases; zero failures/errors/warnings/skips | [Log](unit-final.log) |
| Development equivalence | 405 assertions, 64 cases; zero failures/errors/warnings/skips | [Log](equivalence-final.log) |
| Affected retained and independent references | 566 assertions, 33 cases; zero failures/errors/warnings/skips | [Log](references-final.log) |
| Complete external-equivalence suite | 5,419 assertions, 198 cases; zero failures/errors/warnings/skips | [Log](outer-final/full-suite.log) |
| Complete parallel regression file | 37 assertions; zero failures/errors/warnings/skips | [Log](parallel-cleanup-followup/targeted-final.log) |
| Cleanup implementation lint | No lints found | [Log](parallel-cleanup-followup/lint.log) |
| Source build and vignette generation | Passed | [Log](build-final.log) |
| Native package check, examples and vignette rebuilds | Zero errors/warnings; one current-time verification NOTE | [Log](check-final.log) |

The native check uses `--as-cran --no-manual --no-tests`; regular examples, `donttest` examples and vignette rebuilds run there. Full tests run separately and within each CI package check. PDF manual generation is outside this validation.

The final full-unit run was executed at `d4c9f69`. The local development-equivalence, affected-reference and complete outer runs were executed at `cf1127f`; all 91 production files and all external-suite runtime inputs are byte-identical at the final revision. The only two changed files are the unit fixture and its CI early-gate filter. [Independent provenance proof](outer-final/source-equivalence-d4c9f69.json). Earlier full-unit results are preserved in the candidate directories and are not substituted for the new full-unit run.

The affected references include the claim registry, retained Latent GOLD and Mplus inference, weighted/robust inference, LMR adjustment, mclust covariance/prior/noise comparisons, missing-pattern conventions and independent numerical references. They use retained external outputs rather than claiming new external-program executions.

### Complete final external suite

The complete top-level `equivalence/` suite passed **5,419 expectations in 198 cases across 39 files**, with zero failures, errors, warnings or skips. It completed in **30 minutes 27 seconds** on `cf1127f`, in one continuous clean run without replacement or aggregation of earlier results. This suite is distinct from `tests/equivalence/`. The 566 targeted reference assertions are included in this 5,419-assertion complete-suite result.

The run includes all 160 Monte Carlo datasets, 40 JStats/Mplus 50-start refits, retained Latent GOLD feature and family comparisons, mclust and growth/transition references, and the 199-draw TECH14 bootstrap comparison. It used a fresh isolated copy of all tracked sources, six required sibling JStats inputs, and its own correctly installed worker package. The worker path was asserted and all 91 production-source fingerprints remained unchanged. Original repositories and their retained reports were preserved.

[Complete-suite report](outer-final/REPORT.txt), [measured summary](outer-final/summary.log), [complete log](outer-final/full-suite.log), [case-level CSV](outer-final/cases.csv).

The earlier `a5f211f` broad run and its explicit complete-file replacements remain under [outer/](outer/) as historical diagnostic evidence. They are not substituted for the final run above.

Earlier repair-phase studies separately recorded 376 agreements with retained Latent GOLD quantities, zero disagreements and 11 not-comparable quantities, plus 70 independent probability references with maximum absolute error approximately `5.48e-12`. These are supporting historical studies, not additional fresh unit counts.

## Installed-archive public reproductions

The final `d4c9f69` archive was installed into an isolated library after the final full unit run. Its package version is explicitly checked as 0.96.1. The original public probes show:

- Five observed groups and the same data padded with 95 empty groups give identical mean SE `0.1083709` and variance SE `0.1514179`.
- Weighted tree shares and masses both equal the fitted weighted proportions, `0.8` and `0.2`.
- Raw/generated pair-name collisions raise `latents_bad_data`; native numeric categories that print alike raise `latents_bad_categorical`.
- A real terminated socket worker raises `latents_parallel_failure`, emits no cleanup warnings and restores the original open-connection count.

[Installed-archive log](installed-public-proofs.log), [reproduction script](installed-final-smoke.R).

An independent comparison also verifies that all **1,219 installed function formals/bodies** match the committed runtime source and all **six installed datasets** match their source objects. [Function comparison](final-independent/installed-runtime-comparison.log), [dataset comparison](final-independent/installed-data-comparison.log).

## Why the earlier fixtures passed

The original equivalence coverage concentrated on ordinary fitted cases. It did not establish invariance to translated/scaled coordinates, correct contributing-unit counts after padding, preservation of native category identity, or adapter weights. The corrected tests exercise those public paths directly and compare against independent known answers. Existing fixture tolerances remain unchanged.

The release audit also found 26 stale exact phrases in `equivalence/CLAIMS.csv`. All 67 claims, source references, backing tests and notes remain. Reconciliation with actual documentation text makes all 271 registry assertions pass. [Claim evidence](claims/REPORT.txt). CI now runs this registry and additional retained inference references alongside the existing development-equivalence suite.

## Failures found during release validation

The first complete matrix at `a5f211f` failed on all six platforms; equivalence and lint passed. [Original matrix summary](ci/SUMMARY-a5f211f.txt). Those failures are preserved and are not counted as release sign-off.

- Exact-nested Vuong problems retained numerical null eigenvalues on different BLAS implementations. The correction identifies exact empirical-Gram information before scaling and uses an orthogonal score basis. For general metrics, propagated whitening uncertainty produces a typed unavailable result when genuine small weights cannot be resolved. An overbroad initial cutoff was independently rejected before commit. The revised implementation retains genuine small weights and passes independent rotated and ill-conditioned counterexamples. [Numerical review](numerical-followup/cross-review-v2/REPORT.txt).
- The covariance-floor barrier optimizer could fail with a nonfinite objective or an EM decrease near active bounds. A feasible active-set Newton update now certifies feasibility, projected stationarity, multiplier signs and complementarity. Independent analytical, scalar-optimization, scaling and KKT comparisons passed. A nonbinding-constraint defect in the initial proposal was rejected and corrected before commit. [Covariance review](covariance-followup/cross-review-v2/REPORT.md).
- Parallel regressions now respect the two-worker package-check limit and use the Windows-supported termination signal. Their original equality assertions are preserved.
- At `a40d5c7`, all five non-Windows early regression gates passed. Windows stopped before executing R because PowerShell interpreted `R` as its history alias. Commit `c065dec` explicitly selects Bash for installation; all workflows were redispatched.

- At `c065dec`, the numerical and covariance regressions also passed on Windows. The newly active worker-death test exposed a cleanup defect: a dead socket could make `stopCluster()` mask the original typed worker failure and leave other connections open. Per-node cleanup now closes every owned connection independently; deterministic failure injection and mixed dead/live socket regressions protect this path. Independent review accepted the repair, and the final full matrix passed.

- At `cf1127f`, the full R 4.1 check exposed a brittle occasion-LTA fixture. The source fitter and boundary guard are unchanged from the earlier candidate; a transition probability is near the existing boundary cutoff. The tests intend to check finite-difference step propagation and unsupported network conversion, both of which require a suitable fitted fixture. Independent profiling shows that reducing the old rare-move probability from `6.426e-6` to `8.697e-7` crosses the guard and improves likelihood; the data favor a boundary. The accepted replacement uses all 16 binary four-occasion paths with symmetric residuals, initial probabilities near `0.5` and transitions near `0.496`. Coarse/fine step SEs differ by approximately `5.65e-7`, exceeding the additional `1e-8` margin. All original assertions and their `1e-12` tolerance remain; an ignored-step mutation still fails three assertions. No production fitter or boundary guard changed. [Fixture review](engine-fixture-followup/REPORT.md).

These are finite-precision numerical algorithms. An accurately typed refusal on unresolved information is intentional; the evidence does not establish arbitrary-precision behavior for every ill-conditioned input.

## Source archive and provenance

[Download latents_0.96.1.tar.gz](latents_0.96.1.tar.gz)

SHA256: `033aedba1e4b2a88764dbe9fc6615a010093c29d5e48fab479cea071afbafe1f`.

The [integrity check](release-integrity.log) confirms a clean tracked working tree and verifies all **324 manifest files** against their committed Git blobs. It also checks **286 packaged files** against the archive. A separate inventory also verifies all 306 tracked files present in the archive against Git. The other 38 manifest files are development equivalence/CI files excluded by the existing build rules; they remain covered by the committed-source manifest. `DESCRIPTION` is compared by its original DCF values because package building wraps fields and adds standard metadata. [Verification script](verify-release.py), [manifest](release-source-manifest.csv).

Local `cran-comments.md` follows the repository's existing ignore convention. It records the final measured results and NOTE explanations, replacing obsolete 0.9.14 resubmission and historical win-builder claims.


## Final verification record

The [final gate verifier](verify-final-gates.py) completed successfully. It checks the current commit, clean tracked tree, the exact two-file test/CI delta from the external-run snapshot, all unchanged runtime inputs, local results and source fingerprints, all three exact-commit workflow conclusions, all six successful platform jobs, native check results, installed-archive smoke evidence and archive/Git identity. The remote `features` ref also matches the checked commit.

[Machine-readable verdict](final-validation.json), [gate summary](final-validation.csv), [complete verification log](final-gates.log), [independent final audit](final-independent/REPORT.txt).
