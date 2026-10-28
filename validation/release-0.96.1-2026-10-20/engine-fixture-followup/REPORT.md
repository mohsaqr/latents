# Occasion-transition inference fixture follow-up

The R 4.1 full-check failure is a valid production refusal applied to an unsuitable fixture. Two API regression tests fitted the entire `course_engagement` panel with one start, then assumed that the occasion-specific transition fit was interior. The new CI fit correctly reports a transition/initial probability boundary; Wald inference must refuse it. The network-API test also emits that unrelated boundary warning before testing its intended unsupported-network condition.

Independent local inspection of the old fixture finds 1,422 observations, 106 persons and 15 occasions. Its smallest fitted transition probability is about 6.426e-6, only 6.4 times the existing 1e-6 boundary threshold; the associated logit is about -11.955. A numerical optimization trajectory can reach farther along that nearly boundary solution. The exact reason the old and new R 4.1 runs stop at different points has not been established. The general-LTA/quasi-Newton/boundary-guard production code is unchanged between the compared release revisions, and the newly repaired structured covariance-floor solver is off this default varying-diagonal path. No production guard or tolerance has been changed here.

Changed only `tests/testthat/test-engine-fixes.R`. A deterministic fixture now contains every four-occasion binary state path equally often, with the same symmetric residual quantiles in each path. Means are approximately +/-0.1 and within-state SD is approximately 0.05. This keeps the fitted transition probabilities near 0.496 and initial probabilities near 0.5, far from the probability boundary and variance floor. Explicit convergence/interior assertions establish that the fixture is eligible for the test it exercises.

All original step sizes, 1e-12 consistency/change assertions and network-refusal assertions remain. Added a stronger step-sensitivity margin: the coarse/fine standard errors must differ by more than 1e-8. Actual differences are about 5.64747e-7 in default, internal and BLAS matrix-product modes. This is numerical curvature sensitivity, with a margin beyond floating-point rounding, rather than dependence on an almost-boundary fitted sample.

Verification:

- Full engine-fixes file: **26 passing expectations**, zero failures, errors, warnings or skips (`targeted-final.log`, `results.rds`).
- Three matrix-product modes: converged/interior fits and matching step-sensitivity magnitude (`overlap-probe.log`, spread 0.5 rows). The probe also intentionally examined less-separated candidates; the spread=1 candidate was unsuitable and refused, so it was not selected.
- Memory-only mutation reproducing the original defect, in which the coefficient table discards the caller's information and recomputes at its default step: **three regression assertions fail** (`ignored-step-mutation.R/.log`). Thus the new fixture still detects the original bug.

Independent review is coordinated with `release_audit`; cross-platform confirmation remains the responsibility of the fresh CI matrix.

Independent review accepted the final fixture. It also directly verified that moving the old fixture's small transition probability from 6.426e-6 to 8.697e-7 improves the likelihood by 3.811e-5, and moving it to about 7.23e-13 improves the likelihood by 4.408e-5. The boundary is therefore supported by an improving likelihood direction, rather than caused by the boundary guard. Mocking the changed VLMR/floor/cleanup helpers to throw did not affect the old fixture fit, independently confirming that those repairs are off its execution path. Full independent evidence: `../r41-fixture-followup/cross-review/`.
