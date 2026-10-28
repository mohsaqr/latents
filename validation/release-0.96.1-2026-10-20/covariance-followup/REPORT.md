# Active covariance floors: release follow-up

The Windows release log reports an EM decrease inside the seven-code floor-active profile regression. It does not identify the individual covariance code. The exact synthetic fixture passes locally under default, internal and BLAS matrix-product modes.

Two local counterexamples expose a failure in the same bounded-eigenvalue solver without changing the statistical problem: rounding the fixture to nine significant digits fails EVE with `initial value in 'vmmin' is not finite`, and permutation 21 of the identical rows fails EVE with the same message. Reproductions and saved data are retained in `search.R`, `permutation-search.R`, their baseline logs and RDS files. `new-regressions-baseline.log` preserves failures against the original bounded solver.

The bounded eigenvalue problem remains convex in its existing log-volume/log-shape coordinates, with the same linear floor inequalities. Replaced the numerical interior barrier with a feasible active-set Newton solve. It reaches active floors directly, searches along their null space, and releases constraints whose Lagrange multipliers have the wrong sign. Singular-curvature directions use projected gradient descent. Backtracking evaluates the true objective difference through its gradient and an `expm1` remainder, avoiding subtraction of almost equal totals at convergence. Rounding-sized velocities in constraints dependent on the active set are treated as algebraic zeros.

No changes to the likelihood, variance floor, covariance-family constraints, EM decrease guard or existing test tolerances were made. No new dependency or for-loop was introduced. Changed source: `R/covariance-structure.R`. Changed tests: `tests/testthat/test-review-core-audit.R`. The existing seven-code test now labels each covariance code separately, so a future CI failure identifies the responsible model. Added regressions for the rounded-data failure, an independent analytic EVI active-floor optimum and the all-zero scatter oracle with dependent active constraints.

Verification:

- Six covariance/prior/missing-data/inference test files: **475 passing expectations**, zero failures, errors, warnings or skips. Results and full output: `targeted-results.rds`, `targeted-v2.log`.
- All 56 combinations of seven codes and eight rounding precisions pass (`search-v2.log`).
- All 210 fits from 30 row permutations and seven codes pass (`permutation-search-v2.log`).
- Sixty independently generated constrained problems satisfy their floor/family constraints and attain costs no worse than the original solver where it returns an answer. Largest normalized objective improvement: 1.54e-4. An analytic active EVI solution agrees within 1.12e-16 (`reference-v2.log`).
- Configured lint passes (`lint-v2.log`).

A new Windows CI run is still required to establish that its exact unlabeled EM-decrease failure is resolved. The local counterexamples establish a genuine platform-sensitive bounded-optimization defect and directly verify its correction; they do not recreate that exact Windows error message.


## Independent review and v2 correction

Independent review blocked v1 on a further active-set error, despite its passing targeted suite. With scatter `rbind(c(1e-5, 1), c(1e-4, 100))`, weights `c(1,1)` and floor 1, v1 falsely froze a variance above its floor and differed from the exact optimum by 0.00149794. The absolute tolerance comparing step lengths could treat a nonbinding constraint as active when the allowable step was very small. The retained failure and review are under `cross-review/`; `kkt-v1-regressions.log` records the new permanent tests failing before correction.

V2 activates only the first genuinely limiting constraint. Simultaneously binding independent constraints are added on subsequent zero-length steps, and retained active constraints are checked by their actual slack. Convergence explicitly certifies primal feasibility, projected stationarity, multiplier signs and complementary slackness. A tiny step alone never establishes convergence; exhausting the iteration budget without certification raises `latents_no_converge`.

The independent counterexample now matches its exact solution within 7.11e-15, with objective excess zero and active log slacks about 2.22e-16 (`kkt-v2-independent.log`). Five nearby positive scatter scales form a permanent oracle regression. The final six-file run passes 475 expectations; no existing tolerance or EM decrease rule was relaxed. A fresh macOS CI log also identifies the original `constrOptim` nonfinite failure in the primitive floor test, directly confirming the failure path reproduced locally.
