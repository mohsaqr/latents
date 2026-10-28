# Equivalence re-run after the spike/fitting corrections — 2026-10-22

Current = working tree (R sources match ../source-sha256.txt); baseline =
`git archive 12e0145`, both installed into isolated scratch libraries.

## Harness (`equivalence/harness/run.R`)
- First run refused to start: `check_validation_api()` scanned archived
  `*.Rcheck` trees under validation/ (copies of the package's own tests that
  pass removed arguments on purpose). Fixed in `equivalence/harness/registry.R`
  (skip `.Rcheck/`), with a regression test in `test-harness.R` (15 tests pass;
  fails with the fix removed).
- Re-run on 0.96.2 working tree: 1018/1018 agree, 8 suites, largest
  difference 5.276e-05 (unchanged from 0.95.1).

## Golden fingerprints (260 cases, recorded under both versions)
- 121 identical; 139 differ only by: the new print line, the new
  `gaussian_diagnostics` table (and the summary's table count 19 -> 20),
  two new enumeration columns (35 -> 37), reworded thin-profile warnings and
  the reworded boundary refusal ("variance or covariance eigenvalue"). No
  estimate, likelihood, table value or refusal/acceptance status changed.
- The stored 2026-10-03 golden baseline is stale for both versions (218/260
  DIFFERENT at 12e0145 and now); not re-recorded.

## equivalence/run.R (current code, max_fails = Inf)
- claims 271, depmixs4-SE 15, enumeration-mplus 49, external-scripts 20,
  features-pair-covariance 1242: pass. enumeration-mplus needs R_LIBS to point
  at the matching install: socket workers otherwise load the global 0.95.1.
- features-direct-effects 2 fail (heavy_02), features-joint-items 9 fail
  (twopair_mod_05): failure details identical under 12e0145 (fail-*.txt).
- features-ml-three-step 2 fail (sweep_108), as documented 2026-10-21.
- features-vlmr: stopped by request (about 50 min; 13 documented failures).
- Then 11 more files passed (growth x3, jstats-lta, jstats-mplus 450,
  latentgold x3, mclust x3). The remaining 24 (Monte Carlo 2.0, Mplus) were
  stopped by the user's decision: the golden fingerprint diff already shows
  no numerical change, and these suites were run many times before.
- equivalence/helper-montecarlo.R: montecarlo_fit_all() capped at 2 workers
  (was up to 8, which froze the machine before).

## heavy_02 with search = "emEM" (re-confirms LEARNINGS.md:52)
See equivalence/features/direct-effects/FINDINGS.md finding 1: emEM reaches
Latent GOLD's and Mplus's -3098.487 in 8/8 runs; the default search misses in 8/8.
