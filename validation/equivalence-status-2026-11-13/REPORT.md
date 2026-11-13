# Equivalence status of latents, 13 November 2026

**Verdict: yes, we are good.** On the current code (1.1.20) nothing disagrees
with a reference anywhere. That covers 12,580 harness quantities (3.34 million
individual values against 37 sources), 8,806 argument-sweep checks and 21,448
equivalence-test expectations. The package's own tests pass (15,410
expectations) and `R CMD check --as-cran` is clean.

What does not count as agreement is listed below. It is 25 registered known
differences, 38 sweep checks with no reference able to judge them, and a
handful of skips. Every known difference is on the reference's side: a
reference stops short of the maximum, bounds its parameters, or settles at a
lower local maximum. None of them records latents computing a wrong value.

Every run is from today. Two things are not finished: the verification
matrix was last built on 1 November, and the harness report committed in the
repository is still the 1.1.6 run.

## Which code the evidence is for

Every run below was made today on the current code, 1.1.20 plus the
working tree, which adds only the CVS reporter under `equivalence/`. The
runtime code (`R/`, `NAMESPACE`) is also unchanged since the Nov 12 runs:
1.1.20 changed only `DESCRIPTION` (`git diff a26107dc 5138a96d -- R NAMESPACE`
is empty), and the harness, the sweep and shard A give the same results today
as on Nov 12. All runs were made on the ci VM
(Ubuntu 24.04 arm64, R 4.5.2).

## Evidence

| Run | When | Result | Failed |
|---|---|---|---|
| Comparison harness, 24 suites | **today** | 12,580 of 12,580 quantities agree (3,340,488 values, 37 sources); every produced output compared or exempted; no open defects registered | 0 |
| Argument sweep, 12 modules, 1,044 cases | **today** | 8,767 agree, 1 known difference, 38 no reference | 0 |
| Equivalence tests, shard C (`m`–`z`, 38 files) | **today** | 4,819 expectations pass, 0 skipped | 0 |
| Equivalence tests, shard B (`g`–`l`, 20 files) | **today** | 706 pass, 4 skipped (see below) | 0 |
| Equivalence tests, shard A (`a`–`f`, 13 files) | **today** | 15,923 expectations pass, 0 skipped (15,883 as on Nov 12, plus 40 in the new `test-cvs-report.R`) | 0 |
| Package tests (`tests/testthat`) | **today** | 15,410 pass, 0 skipped, 0 warnings | 0 |
| `R CMD check --as-cran --no-manual`, 1.1.20 | **today** | 0 errors, 0 warnings, 0 notes | — |

Largest difference as a fraction of its stated tolerance in the harness: 0.668
against references computed at full double precision, and 1 against references
read from published digits (a fraction above 1 fails).

## What each function is compared with

From the verification matrix (built 1 November; see "Open items"). A source is
an independent program or package with complete, agreeing comparisons at two
or more classes.

| Function | Sources | Independent sources | Published papers reproduced |
|---|---|---|---|
| `lpa()` | 7 | Mplus, Latent GOLD, Stata (gsem), mclust, Rmixmod, mixtools, EMCluster | 3 |
| `lca()` | 9 | Mplus, Latent GOLD, Stata (gsem), poLCA, depmixS4, glca, BayesLCA, randomLCA, multilevLCA | 7 |
| `multilpa()` | 4 | Mplus, Latent GOLD, depmixS4, mclust | 0 |
| `multilca()` | 4 | Mplus, Latent GOLD, glca, multilevLCA | 1 |
| `lta()` | 5 | Latent GOLD, Mplus, depmixS4, LMest, seqHMM | 3 |
| `mixture_regression()` | 5 | Latent GOLD, Mplus, flexmix, mixtools, MoEClust | 3 |
| `mixture_regression(random =)` | 3 | Mplus, Latent GOLD, lcmm | 1 |
| `three_step()` / `r3step()` | 4 | Latent GOLD, Mplus, depmixS4, sandwich | 0 |
| `bootstrap_lrt()` / VLMR | 3 | Latent GOLD, Mplus, mclust | 0 |
| `enumerate_*()` / `compare_models()` | 3 | Mplus, mclust, poLCA | 0 |
| `pool_imputations()` | 3 | mice, mitml, mitools | 0 |

The argument sweep adds a check of a different kind. Every argument of the
fitting and post-fit functions is run at each of its options on several
datasets, against references that do not use the package's likelihood code.

## What does not count as agreement

**Known differences (25), all on the reference's side.** Each is registered
with its reason in a `known-differences.csv` beside its kit. A registered
difference that stops occurring fails the run, so the register cannot outlive
what it describes.

| Kit | Count | What happens |
|---|---|---|
| `lca-packages` | 7 | multilevLCA bounds response probabilities to [1e-5, 1 − 1e-5], so at boundary maxima its likelihood is lower by 4e-5 to 0.006; one k-means start stops lower; one Linux run reached another solution of a flat model |
| `enumeration-packages` | 4 | mclust's vignette prints statistics at its default EM tolerance, short of the maximum latents reproduces; different bootstrap replicate refits |
| `mplus-growth` | 3 | Over-extracted models: Mplus stops at a lower maximum, or at a one-person class or a negative variance (Heywood case) |
| `stata-published` | 3 | Stata's printed digits are not converged (scaled gradient 1e-5); one unbounded diagonal model |
| `growth-published` | 2 | lcmm's grid search stops at a lower local maximum, or at a quasi-separated membership |
| `regression-packages` | 2 | MoEClust stops below the maximum on tonedata; one Linux run that did not reproduce itself |
| `latentgold-coverage`, `lpa-packages`, `lta-packages` | 1 each | Latent GOLD's errors for one zero-variance slope differ from the observed information by up to 0.8%; Rmixmod returns no estimates for one model; the paper's values are LMest's default stopping point |
| argument sweep | 1 | On five binary items with missing values, ten k-means starts (the default search for categorical indicators) stop 1.82 below the maximum poLCA reaches; `search = "emEM"` reaches it. A search limit, not a likelihood error. The default was kept deliberately (owner, 12 Nov) because changing it would move seeded results of the CRAN release |

**Sweep checks with no reference (38).** A reference that cannot judge a case
is recorded as such, with its reason, rather than passed. Thirty-four fits
(32 in families, 2 in gaussian) lie on a parameter boundary. One covariance
structure is not closed under rescaling, so its rescaling check does not
apply. Three transition cases set a variance option (equal variance, full
covariance) on data with no continuous indicator or only one, where the
option has nothing to act on.

**Skips.**

- Shard B, 4: lcmm 2.2.2 no longer returns subject-specific predictions, so
  three conditional-prediction expectations of `test-growth-predict-lcmm.R`
  cannot run. One check in `test-joint-items-reference-optimizer.R` depends on
  macOS rounding.
- Shard C: none. Its first run today skipped two checks of
  `test-review-1026-oracles.R` because `Rmpfr` was missing from that ci
  workspace; after installing it there, the rerun ran both, with 0 skips.

**Warnings.** Shard A's 482 warnings, the same as on Nov 12, are almost all (480) the
package's own boundary warning in `test-enumeration-packages.R`: a covariance
eigenvalue reaches its floor during enumeration grids. They are expected and
asserted behaviour, not failures. Two more come from
`test-features-direct-effects.R`. Shards B and C have 5 and 8.

## Open items

1. **Verification matrix not rebuilt since 1 November.** Its per-function
   sources above are from that build. Under its provenance rules, results made
   on older source count as historical until rerun and re-recorded with
   sidecars.
2. **Committed harness report is from 1.1.6.**
   `equivalence/harness/REPORT.md` and `report.csv` are the 30 October run. The
   current run's report exists only in the ci results.
3. **Not in scope, by decision (owner, 12 Nov).** `predict()` for transition
   and covariate fits, `three_step()` for transitions, `r3step()` with group
   classes, and homogeneous-transition standard errors under constrained
   structures are not built and therefore not compared.

## Reporting to CVS

Every run now reports to CVS, the validation service, under the target
`latents` (registered in validation PR #22, live since 13 Nov 17:30). Each
file declares its tier, `equiv` or `unit` (validation PR #23); CVS otherwise
classifies by file name, and an R file name carries no tier. The
reporters cover the testthat equivalence files, the package's own tests, the
harness and the argument sweep. Each sends one assertion per test block,
compared quantity or check, and failing runs report too. `.mini-ci` has
`report cvs`, so every named ci job delivers its report automatically. The
reporter and its tests are in `equivalence/cvs/` and
`equivalence/test-cvs-report.R`; see "Reporting to CVS" in
`equivalence/README.md`.

Delivered, with tiers:

- today: the harness, the argument sweep, the package tests, shards A and B,
  and both shard C runs;
- the last full runs before today: the Nov 12 harness and sweep, and the three
  Nov 5 package-test runs behind the November review and fixes reports, the
  last of them labelled `unit-reconciled` because it merges reruns.

The first deliveries, before tiers existed, were filed as `unit`; those runs
were deleted from the server and delivered again.
