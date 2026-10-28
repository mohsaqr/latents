# Vuong follow-up v2 — awaiting independent acceptance

This supersedes `REPORT.md`. The first conditioning-bound repair was rejected by independent review: enlarging the numerical-zero cutoff deleted genuine small negative weights, and an overflowing bound could delete every finite weight. `lmr-before-v2.R` preserves that rejected implementation; `cross-review/` retains the independent counterexamples and evidence.

## Revised implementation

- Recognize **exact supplied empirical-Gram information** (`information` equals a fresh `crossprod(scores)`, with dimensions preserved and names ignored) before chart scaling. For that explicit identity, obtain whitened scores directly as an orthogonal SVD basis. This avoids converting score information to ill-conditioned normal equations and back. It does **not** presume the alternative/null score spaces nest: the subsequent signed-Gram calculation still uses their actual relative orientations.
- Other information matrices continue through triangular whitening and its propagated conditioning/error bounds.
- Restore the original `64 * eps * dimension * largest_singular_value^2` matrix-construction cutoff. The additional whitening uncertainty is **not** evidence of a zero eigenvalue. A nonzero weight between that construction cutoff and the additional uncertainty causes `latents_inaccurate_probability`, preserving an unavailable comparison with its reason.
- Nonfinite signed matrices/error bounds and bounds overwhelming the computed reference scale cause the same typed refusal. They cannot zero all finite weights.
- The Imhof helper is unchanged and retains genuine nonzero weights.

## Tests and independent references

`regression-after-v2.log`: full development run of `tests/testthat/test-lmr.R`, **134 passing expectations**, no failures, warnings, or skips. Original seed-8 exact-nested assertions are unchanged. Added ill-conditioned synthetic nesting cases explicitly permit a typed inability to resolve their chart; two of the three are refused on this platform.

New permanent regressions cover both independent genuine-negative-weight counterexamples, their accurate well-conditioned counterparts, and both finite-overwhelming/nonfinite uncertainty examples. The unresolved examples now produce a typed refusal rather than probability one or an all-zero reference. Existing genuine `1e-10` producer weights, near-zero chi-square/F references, Imhof ratios through `1e200`, public scale/translation checks, and Latent GOLD fixtures pass.

`gram-span-verification.R/.log` independently verifies that empirical-Gram recognition does not force a nested chi-square answer. Rotating four null-score directions out of a five-dimensional alternative gives known weights `1` and four pairs `+/- sin(angle)`. Angles `0.2`, `1e-6`, and `1e-10` retain all nine weights, with maximum absolute weight error `2e-15`.

The same artifact verifies the blocking ill-conditioned references with true probabilities `0.99985000000125` and `0.99998500000000123` are honestly refused with `latents_inaccurate_probability`.

## Remaining gates

Independent release cross-review and exact-commit supported-platform CI must accept v2 before release. This report does not claim those checks are complete. The numerical author froze `R/lmr.R` and `tests/testthat/test-lmr.R` after the successful targeted run and made no staging, commit, push, or CI changes.
