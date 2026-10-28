# Independent covariance v2 cross-review — accepted

**Disposition: the reproduced blocker is repaired.** No further blocker was found in this review. No tracked source or tests were edited.

The revised solver activates only the first limiting inequality, rather than treating small absolute differences between step limits as a shared boundary. It rechecks retained active rows against their actual log slack. Convergence now requires primal feasibility, projected gradient stationarity, multiplier signs and complementary slackness. A small accepted step alone no longer establishes convergence; exhausting 500 iterations raises a typed nonconvergence condition.

The original counterexample now returns the exact oracle to maximum variance error `7.11e-15`, with zero displayed objective excess. Its retained active slacks are both `-2.22e-16`, its multipliers are positive, and it settles after eight iterations. The original 23-scale sweep from `1e-2` through `1e-24` also agrees with the exact oracle to at most `4.98e-14`.

Independent checks retained in this folder:

- **107 core-audit expectations**, no failures, warnings or skips (`core-tests.log`).
- **100 scalar-reduced EVI/VEI reference problems** in two dimensions. The independent reference conditions on shared log volume or shared log shape and performs a separate one-dimensional optimization. Production objective excess is at most `7.56e-16` after normalization. Some production answers improve on the scalar optimizer's finite kink/stopping precision; no claim of bitwise equality to that reference is made.
- **100 independently recomputed KKT certificates** across two to four profiles, two to five dimensions, mixed zero/positive scatters, unequal weights and active/inactive floors. The external check reconstructs the volume/shape design and verifies the final primal feasibility, stationarity, dual sign and complementarity from returned variances. Maximum observed stationarity residual is `1.74e-11`; maximum complementarity residual is `1.60e-15`. All use at most 21 iterations. External stationarity/dual checks use `1e-10 * (1 + max(abs(gradient)))`.
- **40 all-zero-scatter problems** with two to five profiles and two to six dimensions return the floor, exercising singular curvature and redundant active constraints.
- **15 unit/count transformations** of the original exact oracle, spanning variance units and effective counts from `1e-100` to `1e100`, preserve the answer; maximum relative error is below `2.26e-14` (`scaling.log`).

Scripts and CSVs: `verify.R`, `scalar-references.csv`, `kkt-certificates.csv`, `original-oracle-sweep.csv`, `scaling.R`. The `expm1` objective-difference expression and the null-space/multiplier signs remain consistent with the convex log-eigenvalue objective.

This acceptance covers the numerical implementation and the reproduced local failure. Exact-commit Windows and other supported-platform CI remain separate release gates; local cross-review does not establish their completion.
