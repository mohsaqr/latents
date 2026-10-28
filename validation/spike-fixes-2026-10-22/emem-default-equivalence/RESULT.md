# Full equivalence run: emEM default for shared-covariance models

Run 2026-10-22 23:47 to 2026-10-23 02:36, R_LIBS pointed at an install of the
working tree (emEM default for EII/EEI/EEE, kmeans otherwise), two workers,
single-threaded BLAS. Scripts: run-full.sh, eqrun.R. Logs: 30-harness.log,
31-equivalence-full.log. Per-test table: eq-full.rds.

## Harness (equivalence/harness/run.R)

1018 quantities compared, 1018 agreed. Every obtained value equals the
committed 0.96.2 report (largest change 0, compared to 1e-6); only timings
differ. The worst difference, 0.938, is a published-rounding quantity with
tolerance 1 and is the same as before.

## Equivalence suite (44 files)

375 tests, 20348 expectations, 13 failed, 0 errors, 0 skipped, 4 warnings
(the thin-profile diagnostic on small synthetic fixtures).

| File | Before (12e0145) | Now |
|---|---|---|
| features-vlmr | 13 failures (EEI/EEE search misses) | 0 failures, 5462 expectations |
| features-direct-effects | 2 (heavy_02, VVI) | 2, same case (VVI still kmeans) |
| features-joint-items | 9 (twopair_mod_05, LCA) | 9, same case |
| features-ml-three-step | 2 (sweep_108, step-three optimizer) | 2, same case |
| all other 40 files | pass | pass |

No new failure. The only change is that the VLMR file passes.
