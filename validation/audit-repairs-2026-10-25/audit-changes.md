# Audit change summary — 2026-10-26

### 2026-10-26 — Git delivery of the completed audit

- Publish the audit source, regression tests, corrected reports and portable audit handoff together on features. Root session notes remain local according to .gitignore.
- Preserve the ten pre-existing tracked result/report changes and unrelated untracked artifacts outside the audit commit.
- Reuse the completed numerical verification; check staged file coverage, evidence hashes, portable links and whitespace before publishing. Remove trailing spaces from one new regression fixture with identical parsed R expressions and 27 focused assertions passing; preserve archived logs verbatim.

### 2026-10-26 — What was wrong and what was fixed report

- validation/audit-repairs-2026-10-25/defects-and-fixes.html: readable problem/fix/verification report, organized into code repairs, validation/tooling repairs, and documentation/limits; 41 entries cover all 53 ledger groups.
- defects-and-fixes.json and defects-and-fixes-checks.json: structured entries and checked coverage, links, hashes and interactive controls. Independent wording review distinguishes fitted-prior checks, error units, fork timing and the fitted drop_empty policy.
- Checks: all 86 local links resolve; all groups appear exactly once; search, category filters, empty results, reset and print actions pass the JavaScript smoke check. Production code and audited numerical results are unchanged.

### 2026-10-26 — Audit repairs and completed verification

- R/covariates.R, R/kernel-blocks.R, R/kernel-structures.R, R/extra-indicators.R, R/centering.R and R/weights.R: condition on contributing rows, preserve native category distinctions, stabilize extreme count/ordinal and weight calculations, and retain ordinary categorical EM arithmetic.
- R/alignment-service.R, R/accessors.R, R/inference.R and R/covariate-inference.R: use covariance/dispersion information for class alignment, encode reserved category-name characters, and distinguish assignment-file roundtrips from strict objective-data identity.
- R/additive-inference.R, R/lta-three-step.R, R/lta-compare.R and related methods: condition inference coordinates, exclude empty units from robust corrections and informative sample counts, and compare actual data, normalized weights and fitting bounds consistently.
- R/three-step-ml.R: supplement existing equal-variance starts with guarded BFGS candidates. The retained regular counterexample now reaches the independent higher likelihood and covariance; finite search does not guarantee a global optimum.
- R/mixture-regression-*.R and R/growth-mixture-*.R: stabilize posterior normalization and predictor coordinates, preserve fitted random-design bases and ordered contrasts, reject negative Brownian prediction times, and enforce prediction/plot data contracts.
- R/descriptives.R, R/tna.R, R/tree.R and R/selection.R: retain distinct native numeric IDs, count mover/stayer group classes correctly, and show unavailable replication evidence explicitly.
- equivalence/harness/: reject nonfinite comparisons, validate tolerances and seeds, check convergence in stability results, and report mixed precision without treating every quantity as printed-rounding units.
- equivalence/features/: declare uniform direct-effects dual-search and joint-item strict-refinement comparison recipes, preserving original default diagnostics and unchanged reference tolerances. The independent joint-item EM reference retains coherent best states and requires nonnegative convergence gain; poLCA comparisons uniformly retain the original 20,000-iteration fits and verify with a 200,000-iteration cap using the same seeds, starts and tolerance. Seeded ML3/VLMR validation dispatch avoids unsafe forks. Supporting study runners receive explicit fork preflight guards without changing their scientific protocols.
- equivalence/*/run.sh and validation/audit-repairs-2026-10-25/test-external-runners.py: fix shell entry/export behavior, require normal/fresh external output, and verify failure and success paths with fake executables.
- future/R/: adapt deferred random-intercept reporting/plots to current services and validate deferred LMR data identity. The deferred LMR p-value remains unavailable.
- vignettes/articles/, Tutorial_docs/two-level-lpa-engagement.Rmd and man/: correct model/inference explanations, migrate retired tutorial calls and document finite-search, data-identity and RNG limits. DESCRIPTION declares numDeriv as an optional test dependency for the independent numerical covariance oracle.
- validation/audit-repairs-2026-10-25/: correct the supplied HTML report's counts, model descriptions and unsupported certifications; retain per-function source reviews, before/after evidence and hash manifests. The supplied Downloads HTML and the ten pre-existing tracked report/result artifacts remain unchanged.
- Verification: all 1,527 unit cases executed across matched-BLAS/safe-fork runs, 10,948 distinct case assertions, zero failures/errors/warnings; native check including installed tests Status OK; all 20 guides render; deferred tests 78, harness self-tests 56 and external-runner contracts four pass. Development equivalence has 405 passing assertions and one known thin-support fixture warning. The eight-suite harness has 1,021/1,021 quantities agreed. The full affected joint-item rerun has 164 passing assertions with unchanged numerical oracles. The consolidated external ledger passes 21,159 assertions across 52 files and 400 cases after complete affected-file reruns. Four in-case and 25 VLMR setup warnings are retained thin-support diagnostics; see HANDOFF.md for scope and provenance.
