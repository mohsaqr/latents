# Model-core audit fixes — 2026-10-13

All four confirmed model-core findings are addressed. The implementation changes are confined to the assigned R sources and two new regression files. No generated help pages or namespace files were edited in this track.

| Finding | Closure | Validation |
|---|---|---|
| MC-1: ignored covariate centering | The covariate dispatch, base start, measurement design, stored indicators, reconstructed data and Wald-inference inputs now use the requested centering. The fit records the scheme and offsets. | Person-centered means and likelihood reproduce a manual transform; grand-centered means shift back to the ordinary fit; raw-data reconstruction and inference agree; incorrect indicator values are refused. Mixed missing indicators, full covariance and sampling weights pass together. |
| MC-2: covariance floors broke model constraints | Active floors are solved within a linear chart of log eigenvalues, common/varying volume and determinant-one shape. Full-covariance models retain their allowed orientations; common orientations use objective-improving pairwise rotations. A previous covariance is retained only if it satisfies the floor and the requested family, including after SQUAREM extrapolation. The same repair applies to floor-active EVI/VEI/VEV prior updates. | Floor, volume, shape and common-orientation properties pass across all 14 structures; the covariance objective does not deteriorate from a feasible previous covariance. Seven floor-active public models converge without an EM decrease. Six malformed cross-structure starts cannot be retained as invalid fallbacks. |
| MC-3: growth quadratic cancellation | At large random/residual variance ratios, the likelihood uses the equivalent sum of a conditional residual quadratic and standardized random-effect penalty. Direct conditional residuals also stabilize the variance M-step and scores. Random-space precision and projected fixed-effect scores avoid the corresponding Woodbury subtraction. | Independent scalar Gaussian densities agree across amplitudes 1e3–1e8. Extreme fixed-effect, residual-SD and random-SD score references pass. Existing growth and multilevel-growth equivalence suites pass. |
| MC-4: unexplained prior conventions | Public prior documentation now states the EII scale, EEI mean/shrinkage, VEV dof and mean-prior denominator conventions explicitly. Intentional mclust behavior is retained. | Existing prior unit and equivalence suites pass; new prior floor tests retain equal-volume/shared-shape constraints. Root regenerates and checks the public help. |

The bounded solver changes only cases where the raw or proposed covariance approaches the configured floor. The ordinary, floor-inactive update arithmetic remains on its established paths. Common-orientation searches remain generalized EM updates when their configured inner budget stops before the orientation settles. Active-floor updates can take longer than the old componentwise clipping, because they now solve the requested constrained problem.

Growth fallbacks retain the established arithmetic when the relevant subtraction is well conditioned. The independent fixed-effect score test allows the package's pre-existing relative Cholesky jitter; it still distinguishes the formerly incorrect score by several-fold.

Validation completed in this track:

- New regression suite: 96 successful unit expectations; independent growth suite: 7 successful expectations.
- Existing unit suites: covariance structures, covariate inference, covariate measurement, covariate missingness, growth mixtures, weights and priors all pass.
- Existing equivalence suites: growth mixtures, multilevel growth, missing covariance and priors all pass.
- Namespace usage checks of all five new helpers report no diagnostics; `git diff --check` passes for the changed sources and tests.

Artifacts: `new-tests.log`, `targeted-tests.log`, `equivalence-targeted.log`, `dispatch-fixed.log`, `covariance-fixed.log`, `growth-fixed.log`, and `fix-validation.log`. The original audit failure logs remain as the historical evidence; rerunning the original probe scripts refreshes their RDS snapshots with corrected results.

Root owns the final full-suite run, documentation generation, package check, and consolidated learning/handoff files. Individual targeted suites were not expanded into a second full package run in this track.
