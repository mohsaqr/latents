# Conditioning-aware Vuong numerical-zero follow-up

The R 4.1 CI job retained a null eigenvalue of `-5.54537857482984e-13` in the exact nested reference. Windows retained `+2.5350041423043e-12`. Their nonzero probabilities near a zero statistic exposed a missing error term in the producer's numerical-rank treatment.

The previous numerical-zero bound covered only construction/eigendecomposition of the symmetric signed-Gram matrix. It omitted the information factorization and score whitening that produced that matrix. A near-collinear information chart amplifies a componentwise rounding error before its apparently well-scaled score matrix reaches the eigendecomposition.

`R/lmr.R` now:

1. Solves directly with the triangular information factor rather than multiplying scores by a generic explicit inverse.
2. Bounds the Cholesky information perturbation by `gamma[p+1] |R'| |R|`, including the two-product diagonal scaling bound. Here `gamma[k] = k*eps / (1-k*eps)`.
3. Propagates that componentwise bound through `R^-1`. Its spectral norm `eta` bounds the whitened information perturbation; information whose bound cannot resolve the chart (`eta >= 1`) is refused.
4. Bounds the score perturbation from the triangular solve and original column scaling. Its effect on the score Gram matrix is bounded by `2 ||S|| ||dS|| + ||dS||^2`; the factorization perturbation contributes `||S||^2 eta/(1-eta)`.
5. Adds the two models' bounds to the **unchanged** symmetric-matrix construction bound. Weyl's eigenvalue perturbation bound supplies the numerical-zero criterion.

This accounts for measured conditioning through the actual triangular factors. It does not enlarge a fixed constant or relax any test tolerance. The Imhof helper remains unchanged and continues retaining genuine nonzero weights.

## Deterministic local reproduction

`probe.R` and `probe-before.log` use a null score chart with columns `s1`, `s1 + separation*s2`, `s3`, and `s4`, nested in a five-coordinate alternative. Exactly one chi-square(1) weight must remain. Separations `1e-2`, `1e-3`, and `1e-4` give information condition numbers about `4.6e4`, `4.6e6`, and `4.6e8`.

The former producer retained spurious null weights of `1.6397e-12`, `6.7185e-10`, and `7.1994e-8`. The added permanent regression failed **12 expectations** before the repair (`regression-before.log`). It now returns one unit weight and numerical zeros in all three charts (`probe-after.log`), including correct exact chi-square probabilities at zero and `1e-12`.

## Targeted validation

`regression-after.log`: complete development run of `tests/testthat/test-lmr.R`, **132 passing expectations, zero failures/warnings/skips**. This includes:

- The original R 4.1/Windows-sensitive fixture, with unchanged assertions.
- The three deterministic near-collinear charts and unchanged exact-reference tolerances.
- The independent genuine `1e-10` producer eigenvalue.
- Direct Imhof repeated weights spanning `1e200` in relative magnitude.
- Existing Latent GOLD probability fixtures, public indicator-unit transformations, direct-effect translations, empty-unit/zero-weight adjustment, and invalid-comparison policy.

Only `R/lmr.R` and `tests/testthat/test-lmr.R` were modified. Platform CI reruns remain the parent's release gate; these local results do not claim a completed R 4.1/Windows execution.
