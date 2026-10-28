# Release audit — latents 0.98.0 — 2026-10-28

**Subsequent repair:** all six REL-1 through REL-6 defects below have been repaired and validated in the working tree. See the [repair report](../release-fixes-2026-10-28/REPORT.md) for current results and remaining release work. This audit preserves the original findings and release recommendation for the source reviewed before those repairs.

**Recommendation: hold this release.** Four concurrent reviewers (the coordinator and three specialist agents, the available maximum) confirmed two high-priority defects, four further correctness defects, and one lower-priority diagnostic issue. The two high-priority defects affect reported statistical results. The latest external-equivalence CI failure is a real regression, not merely a tolerance problem.

Scope: current working tree at `e907bfffe8f89217ec992188ce7c95f52b4843bf`, recent changes since `d5d4ae8` (68 changed R source files), and relevant surrounding implementation. The pre-existing uncommitted change to `tests/testthat/test-release-engine.R` was included in local testing and left intact. Production source, tests, version, and Git history were not changed by this audit. The retained scripts are synthetic reproductions and independent checks.

## Findings, ordered by priority

### REL-1 — P1: shared regression dispersions disappear from results

Locations: [mixture-regression-starts.R:395](../../R/mixture-regression-starts.R#L395), [mixture-regression-inference.R:31](../../R/mixture-regression-inference.R#L31) and line 34; downstream lookup in `R/mixture-regression-methods.R:108-125`.

The latest commit names `params$sigma2`. Packing an equal variance with `c(log_sigma.all = named_scalar)` now produces `log_sigma.all.class_1`; the result table still looks up `log_sigma.all`. The negative-binomial equal-dispersion branch has the same problem. Converged models report missing SD/dispersion and standard errors, including ordinary one-class fits. The Gaussian reproduction should report SD 0.844399; the negative-binomial reproduction should report dispersion 0.2821447 and SE 0.06338804.

This reproduces the failed `flexmix` CI suite. Removing incidental scalar names **in memory only** makes its entire 433-quantity / 29,933-value suite pass. That diagnostic run is not a pass of the unmodified release. Preserve explicit packed parameter names and add public-table/name regression tests for both families, one and multiple classes.

Evidence and scripts: [equivalence report](equivalence/REPORT.md), `equivalence/probe-flexmix.R`, `probe-gaussian-one.R`, `probe-equal-dispersion.R`, and `probe-flexmix-patched.R`.

### REL-2 — P1: general LTA covariance is not invariant to outcome units

Location: [lta-inference.R:523](../../R/lta-inference.R#L523), through line 525.

The general-engine LTA information calculation differences continuous means with a minimum absolute step of `1e-4` and tests conditioning on the unstandardized information. For an otherwise identical well-identified model, multiplying outcomes by `1e-4` changes reported SEs by up to **2.77%**; multipliers `3e-5`, `1e-5`, and `1e6` falsely raise `latents_singular_information`. The effect is present in both `vcov()` and `parameter_inference()`.

Independent refits have equivalent likelihoods within 3.64e-12 after the density Jacobian. Exact rescaling of the same fitted object isolates inference from optimization. At scale `1e-5`, appropriate derivative steps and correlation-scaled inversion restore SEs within 1.18e-8 relative. Standardize outcome coordinates for differentiation and inversion, including robust/OPG paths. This predates the latest commit; the new LTA `coef()` location correction itself passed review.

Evidence: [inference report](inference/REPORT.md), `inference/units-probe.R`, `units-inference.R`, and `units-inference.csv`.

### REL-3 — P2: zero-weight rows can erase regression inference

Location: [mixture-regression-inference.R:183](../../R/mixture-regression-inference.R#L183), through line 187; callers at 202-205.

The inference coordinate transformation computes predictor centers/scales from excluded rows. Appending one row with `x = y = 1e8, w = 0` to a 100-observation regression leaves coefficients and likelihood exactly unchanged but replaces both coefficient SEs, p-values, and confidence limits with `NA`. Values `1e7` also fail in this fixture; smaller tested values through `1e6` cause negligible errors, not substantive discrepancies. Recomputing inference on positive-weight units restores the entire baseline covariance exactly. Determine transformation moments using contributing rows/groups.

Evidence: [inference report](inference/REPORT.md), `inference/zero-weight-inference.R` and its log.

### REL-4 — P2: excluded sentinel groups can abort an additive fit

Location: [additive.R:274](../../R/additive.R#L274), through line 275; assertion at 300-302.

Posterior normalization subtracts two very large negative log values and loses the normalizing term. Two groups of weight zero with finite values around `100000` or `999999` cause `posterior rows must sum to one`, although the 30 contributing groups fit successfully with log likelihood -178.4893059150279. These are fixture-specific thresholds. The profile and cross-level engines preserve the baseline result through excluded values of `1e9` in the same probe. Use stable row-max-shifted probability normalization and preserve zero-weight exclusion throughout.

Evidence: [core report](core/report.md), `core/reproduce-additive-zero-weight.R` and its log.

### REL-5 — P2: comparison wrappers can silently omit a model

Location: [compare-models.R:115](../../R/compare-models.R#L115), through line 120.

The compatibility logic for fits labelled `n_profiles`, `n_group_classes`, or `cutoff` derives ordering from `sys.call()`, which counts forwarded `...` as one expression. The actual fit list has already expanded those dots. Subsetting by the shorter ordering drops a model without warning and changes BIC weights. With three distinct fitted candidates, a direct call keeps all three, but `function(...) compare_models(..., n_profiles = third)` keeps only the first two; a generic forwarding wrapper also fails. Expand dots consistently or restore models without a mismatched position vector.

Evidence: [public API report](equivalence/PUBLIC-API.md), `equivalence/probe-public-api.R` and its log.

### REL-6 — P2: noise classification is incorrectly declared trivial

Location: [diagnostics-checks.R:220](../../R/diagnostics-checks.R#L220) and the one-class shortcut at line 264.

A model with one Gaussian profile **and a noise component** has two classification categories, but the checks count only `n_profiles`. They report "one class: nothing to classify" and omit class-size, entropy, and average-posterior checks. The reproduction has 252 profile and 48 noise assignments and nontrivial reported entropy, so this is an incorrect diagnostic conclusion. Count the noise category consistently with the actual posterior/entropy tables.

Evidence: [diagnostics follow-up](core/DIAGNOSTICS.md) and its synthetic probes.

## Additional observations

- **P3, new: weighted modal/estimated gap mixes populations.** The informational diagnostic at `R/diagnostics-checks.R:287-291` compares raw modal proportions with weighted estimated proportions. A perfectly separated 150/150 sample with weights 1/9 reports a gap of 0.400, whereas weighted modal and posterior shares both equal 0.1/0.9. The underlying classification table deliberately retains raw modal counts and is not itself wrong. This check is `info`, with no warning threshold; align the diagnostic quantities or make its interpretation explicit. See `core/DIAGNOSTICS.md`.
- **Known, reproduced: sensitivity labels unfinished refits as best.** At `R/sensitivity.R:160-168`, `optimum` and `best` ignore convergence. Three `max_iter = 0` refits all return `best = TRUE`, `optimum = 1`, `converged = FALSE`. The plot uses a different, convergence-aware rule. This was already open in the handoff and is not counted as newly discovered. See `release/probe-sensitivity.R` and `sensitivity.csv`.
- **CI coverage gap:** the new `equivalence/test-latentgold-wider.R` is not invoked by the current equivalence workflow. Its provenance test also omits output hashes and the correspondence of `cases.rds` to hashed `.dat` inputs. This audit independently verified that all retained artifacts currently agree; this is a future regression-detection gap, not evidence of corrupt references.
- The tracked harness report identifies itself as version 0.97.2 with 15 suites. It is historical evidence; the current registry contains 24 suites. Current evidence is identified separately below.

## Validation evidence

- Full source build **with vignettes** passes (`release/build-vignettes.log`). The built tarball passes `R CMD check --as-cran --no-tests --no-manual`: **0 errors, 0 warnings, 1 NOTE** from network-dependent incoming/URL checks unavailable in this environment (`release/check-vignettes.log`). Examples, `--run-donttest` examples, vignette code, documentation, and package/code checks pass. Tests run separately below; the PDF manual is not checked.
- Seven parallel-related files rerun with local socket access: **819 passed, 0 failures/errors/warnings, 15 expected skips** for fork operations unsafe with Apple Accelerate (`release/parallel.log`). The remaining socket-dependent search file separately passes all 234 assertions (`release/search.log`).
- Full unit suite, consolidated after replacing nine environment-dependent files with their unrestricted reruns: **220 files, 1,953 test blocks, 13,340 passed, 0 failures, 0 errors, 0 warnings, 16 skips**. The initial sandbox run had six worker-socket errors (five in `test-parallel.R`, one in `test-search.R`); all pass with socket access. The CPU-count refusal file also passes all 17 assertions outside the sandbox; it had skipped one test because CPU detection was unavailable. The nine rerun files total 1,070 passed and 15 expected skips. The final 16 skips all concern fork execution unsafe with Apple Accelerate. Counts are not double-counted: `release/summarize-units.R` produces `unit-summary.log` and `consolidated-unit-results.csv` from the retained result objects. The existing uncommitted release-engine correction passes locally; cross-platform validation of that correction remains outstanding.

Fresh checks on unchanged production source:

- Nine current equivalence suites: **3,715/3,715 quantities, 298,945 individual values**, no skipped suites, warnings, or uncovered outputs. Detailed suite/coverage/source records: [harness report](equivalence/harness/REPORT.md).
- Wider Latent GOLD: **172/172 comparisons** against retained genuine listings. All 16 input and 16 output hashes match; all eight RDS fitting inputs correspond to the hashed data files. An injected discrepancy is rejected. This is agreement within the declared printed-reference tolerance, not full-double identity. No new proprietary-software run was performed.
- Harness self-tests: **201 expectations in 42 successful blocks**.
- Gaussian unit invariance: all **14 covariance structures × 3 scales = 42 fits plus inference** pass. Worst estimate discrepancy 7.78e-7 of its own SE; worst relative SE discrepancy 1.18e-7.
- Weighted Gaussian FIML: analytical ML means/variances agree within 1e-10 and likelihood within 1e-9.
- Six multilevel families × two scales: all **12 fits converge**, with adjusted likelihood discrepancies at most 7.28e-12; boundary and unsupported inference refusals remain consistent.
- Independent one-class growth-model reference: fixed effects agree within 8e-15, variances within 6.5e-10, fixed-effect covariance within 7.7e-12.
- Independent modal three-step reference: means and CR1 sandwich SEs agree exactly; contrast SE within 1e-12.
- LTA public coefficient means, coefficient/covariance names, pooled coefficients/covariances/intervals, and inspected result-table consistency checks pass.

Current remote evidence:

- [Equivalence CI on e907bff](https://github.com/mohsaqr/latents/actions/runs/40725430305) **fails** at `flexmix` for REL-1. The completed comparisons agree on 12,035 quantities / 3,302,066 values; this is not a passing complete run, and later CI test steps were skipped after failure.
- [R package CI on e907bff](https://github.com/mohsaqr/latents/actions/runs/40725420226): macOS passes; Windows and all four Linux jobs fail. The completed job logs show the same five failures in the grand-scaled SQUAREM release-engine test, the area of the already-present uncommitted correction. Oldrel-1 additionally reports four warnings from a GLM reference and a boundary enumeration fixture; they are recorded in `release/check-all-failed.log`. Local success of that correction does not establish cross-platform success. Final status snapshot is in `release/check-ci-jobs.json`.
- [Lint CI on e907bff](https://github.com/mohsaqr/latents/actions/runs/40725425361) passes. Production source is unchanged from that commit.

## Release follow-through

Repair REL-1 and REL-2 before release; add regression cases that fail on this source. Resolve the other confirmed correctness defects and rerun their reproductions, the relevant equivalence suites, and the full test suite on the repaired source. Add the wider Latent GOLD test and provenance checks to CI. Include the existing release-engine test correction in the candidate to be checked, then obtain green cross-platform checks on that exact candidate. This audit did not commit, publish, or submit a release.

## Scope and reproducibility limits

The audit combines source review, complete local unit testing, targeted independent numerical references, nine fresh equivalence suites, and current CI evidence. It does not prove every family/option combination correct. The other fifteen equivalence suites were not all rerun locally; remote evidence is distinguished from fresh local evidence. Existing documented gaps (including group-class Latent GOLD parameter/SE equivalence and robust Step-3 distal SE coverage) remain open. PDF-manual validation is excluded.

All R computations use one BLAS thread. The coordinator installed current source into `/private/tmp/latents-release-audit-20261028/lib` so socket workers can load this candidate. Reproduction scripts use that installation or `pkgload::load_all('.')`; run from the repository root. `release/source-manifest.json` records the reviewed HEAD and SHA-256 fingerprints of 332 source/test/metadata files. Specialist reports list their inspected areas, exact commands, and limits. Logs, CSV summaries, and reproductions are retained under this directory.
