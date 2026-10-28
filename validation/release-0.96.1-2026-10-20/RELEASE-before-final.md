# latents 0.96.1 release preparation

**Release is held pending platform corrections and a fresh CI run.** The first complete matrix on commit `a5f211f46fb93f4290032775e3d47af48f8f1717`, pushed to `features`, exposed numerical and platform-specific test failures. The numerical repairs, four new regression files and finalized package metadata are committed. The tracked working tree is clean.

## Release contents

Version 0.96.1 closes all ten findings in the preceding review and the related defects found during repeated analysis. The [repair report](../fixes-2026-10-20/FIXES.md) records each numerical or API correction, independent references and regression evidence. No new dependency or exported function was introduced.

Two commits comprise this release preparation:

- [`c4c1df2`](https://github.com/mohsaqr/latents/commit/c4c1df23a0083dd6b084714b3cae7b749587414f): reviewed fixes, four regression files, version/NEWS finalization and additional retained-reference CI coverage.
- [`a5f211f`](https://github.com/mohsaqr/latents/commit/a5f211f46fb93f4290032775e3d47af48f8f1717): corrected documentation-claim registrations and their automatic CI check.

The independent release audit found 26 stale exact phrases in `equivalence/CLAIMS.csv`. Documentation wording and line wrapping had changed. All 67 claims remain, with their sources, reference programs, backing tests and explanatory notes preserved. The complete registry test now passes all 271 assertions. [Reconciliation evidence](claims/REPORT.txt) records the original failures and the precise source lines supporting the corrected registrations.

The existing equivalence workflow now runs both the 405-assertion development suite and an additional 390 assertions covering the claim registry, retained Mplus standard/robust/weighted inference, LMR adjustment and missing-data conventions. The additional step passed locally without failures, warnings or skips. [Local step log](ci-retained-step.log).

## Exact-commit GitHub validation

All three manually dispatched workflows have their `headSha` verified against `a5f211f`. Feature-branch pushes alone do not trigger these workflows, whose automatic push filters name `main` and `master`.

| Gate | Run | Result |
| --- | --- | --- |
| Package matrix: macOS and Windows release; Linux release, devel, oldrel and R 4.1 | [R-CMD-check](https://github.com/mohsaqr/latents/actions/runs/39439978015) | Failed; follow-up corrections and rerun required |
| Development equivalence plus retained inference/claim references | [equivalence](https://github.com/mohsaqr/latents/actions/runs/39439985705) | Passed: 405 + 390 assertions; no failures, warnings or skips |
| Configured package lint | [lint](https://github.com/mohsaqr/latents/actions/runs/39439992164) | Passed |

The superseded check/equivalence runs on `c4c1df2` were canceled so that the final gates correspond to the registry-inclusive revision. The earlier lint also passed. The final run/job records are retained under [ci/](ci/). The completed [equivalence log](ci/39439985705-complete.log) confirms both assertion counts, and the [lint log](ci/39439992164-complete.log) reports no lints.

## Local package and reference validation

The release source builds successfully. Its native package check has **zero errors, zero warnings and one environment NOTE: unable to verify current time**. Regular examples, `donttest` examples and vignette rebuilds pass. Tests and the PDF manual were omitted from this local check; the earlier complete unit/development runs and the exact-commit CI supply the test evidence. [Build log](build.log), [check log](check.log).

The underlying repaired source passed 8,104 unit expectations and 405 development equivalence expectations before the release metadata was finalized. One missing condition-catalogue entry was corrected and its complete test file rerun; all executable production expressions were unchanged. The repair report retains the unaggregated failing log and explicit replacement aggregation. Fresh package fits against retained Latent GOLD output produced 376 agreements, no disagreements and 11 not-comparable quantities; the external program was not rerun. Seventy independent probability references had maximum absolute error approximately `5.48e-12`.

The installed **0.96.1 archive** also passes the original public reproduction probes: appending empty groups preserves robust uncertainty, the weighted tree reports the correct 80%/20% shares, initial pair-name collisions are refused, and ambiguous native numeric categories receive categorical validation. [Installed-archive evidence](installed-public-proofs.log).

### Separate complete external suite

The top-level `equivalence/` harness is distinct from `tests/equivalence/`. Its complete suite finished in an isolated scratch copy: **198 cases, 5,419 passing expectations, zero failures, errors, warnings or skips** across 39 files, in 28 minutes 14 seconds. This is broad evidence for the `a5f211f` implementation before the platform follow-ups below; relevant references will be rerun after those corrections.

Six required sibling JStats inputs were copied into an isolated sibling tree to exercise external scripts without changing the original repositories or their retained reports. All 160 Monte Carlo datasets, growth/transition references and the 199-draw TECH14 comparison passed.

The initial registry failures are replaced by the complete corrected registry file (271 passing expectations). An initial enumeration error came from a socket worker loading the ambient installed 0.95.1 package; the complete file passes all 49 expectations after installing 0.96.1 into an isolated inherited `R_LIBS`. The aggregation replaces complete files without double counting. [Outer suite report](outer/REPORT.txt), [measured summary](outer/summary.log).

## Platform follow-ups in progress

The matrix completed with failures on all six platforms. The full [CI summary](ci/SUMMARY-a5f211f.txt) and per-job logs preserve the exact results. No release sign-off is based on that failed matrix.

- A small numerical null eigenvalue survives the original construction-only cleanup under the platform BLAS implementations. Revised inference recognizes exact empirical-Gram information before scaling and uses an orthogonal score basis; unresolved genuine weights receive a typed unavailable result. Independent review found and rejected an initial overly broad conditioning cutoff before it was committed. The revised numerical follow-up passes 134 assertions and its independent counterexamples. [Numerical review](numerical-followup/cross-review-v2/REPORT.txt).
- One parallel regression asked for three workers under the two-worker package-check limit. It now uses two workers and preserves its uneven five-start split and equality assertions.
- The Windows worker-death regression used `SIGKILL`, which is undefined there. It now uses the supported `SIGTERM`. The complete parallel file passes all 28 assertions locally under the check core-limit setting; native Windows rerun remains required.
- The covariance-floor barrier optimizer fails at active bounds: macOS reports a nonfinite initial barrier objective, and Windows reports an EM decrease. A feasible active-set correction is under independent KKT and objective review. The likelihood monotonicity guard and existing tolerances remain unchanged.

The matrix workflow now installs the current package for socket workers and runs the affected regression files before the full package check, providing earlier platform diagnostics on the next revision.

## Checked source artifact

[latents_0.96.1.tar.gz](latents_0.96.1.tar.gz)

SHA256: `5b5df7d726ec80789a7ba1a4217dc975fdbc96025ab2741b1f67cbe239207f61`.

The [integrity check](release-integrity.log) verifies 321 source/test/documentation/CI manifest files against the committed Git contents and verifies 286 packaged files against the source archive. The existing build exclusions omit development equivalence files and GitHub workflows; their committed contents are still included in the release manifest. `DESCRIPTION` is compared by its original DCF values because `R CMD build` wraps fields and adds normal packaging metadata. [Verification script](verify-release.py), [manifest](release-source-manifest.csv).

`cran-comments.md` has been rewritten with measured current results and the release's CI targets, removing the stale 0.9.14 resubmission and historical win-builder claims. It remains local under the repository's existing ignore convention.
