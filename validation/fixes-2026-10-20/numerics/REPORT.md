# Numerical repairs and repeat review — 2026-10-20

## Changes

- **F05, ML distal outcome origin:** `R/three-step-ml.R` estimates in centered, standardized outcome coordinates, retains those means before reported-coordinate rounding, calculates influence contributions in that chart, and restores the original units. Mean contrasts cancel the origin before subtraction. Conditioned information must be positive definite before inversion. Both variance options, independent/clustered uncertainty, and mean/pair output are covered.
- **F06, VLMR coordinate units:** `R/lmr.R` uses the inference engine's conditioned Hessian and score attributes (added by the inference repair), with consistent diagonal equilibration for standalone score/information inputs. Cholesky whitening and a symmetric signed-Gram eigendecomposition preserve the Vuong distribution under unit changes and translated direct-effect covariates.
- **N01, negative Imhof arguments:** continuous-distribution reflection handles negative statistics; one-sided support and the all-zero degenerate distribution are explicit. A genuinely worse alternative is an unavailable comparison with an explanatory `note`; near-equal likelihoods within numerical tolerance use LR zero. Expected per-comparison numerical failures preserve the table row.
- **N02, genuine tiny Imhof weights:** the distribution helper removes only exactly zero weights. It retains nonzero weights across the old `1e-10` relative cutoff. Log-domain factors/truncation bounds prevent products and squares from overflowing at extreme ratios. Unrepresentable relative weight/integration ranges raise `latents_inaccurate_probability`.
- **Additional repeat-review defect, null Vuong eigenvalues:** the former unsymmetric matrix is defective at exact regular nesting and generates spurious paired eigenvalues of order `1e-7`. In the saved independent 5-versus-4-dimensional known-answer case, its probability at zero was `0.9998101` instead of 1. The equivalent symmetric signed-Gram calculation, with a machine-backward-error zero criterion at its producer, gives one unit weight and exact numerical zeros. A genuine `1e-10` weight survives. The helper itself applies no relative magnitude cutoff.
- **Additional repeat-review defect, aLMR empty units:** `n_units` now counts the inference service's informative positive-weight independent-unit mask. Score matrices keep their original aligned rows. Empty padding leaves the adjustment unchanged under the default convention; `drop_empty = FALSE` preserves its explicit legacy count. Observed zero-weight units are excluded under either convention.

## Regression evidence

- `pre-fix.log`: original implementation, new LMR regressions demonstrate raw-information inversion and negative-argument defects. Some initial tiny-unit public inputs hit the declared fixed variance floor; the final unit-invariance tests transform `min_variance` quadratically with the indicators, preserving the same statistical model.
- `pre-fix-ml.log`: original implementation, exact hard-classification mean uncertainty and public large-origin ML cases fail (including false singular information).
- `nested-eigen-reference.R/.log`: independent reconstruction of the previous matrix formula alongside repaired weights and the exact chi-square reference.
- Root-level `../lmr-empty-probe.log`: pre-repair empty-unit adjustment defect; `lmr-empty-fixed.log` reruns that public probe after the repair.

## Verification

- `final-tests.log`: full development runs of `tests/testthat/test-lmr.R` (**120 passing expectations**) and `tests/testthat/test-three-step-ml.R` (**69 passing expectations**), with zero failures, warnings, or skips. Existing Latent GOLD fixture comparisons remain passing.
- Public LMR tests cover EEI and EEE structures, indicators at `1e-8` and `1e8` units, direct-effect predictor offset `1e6`, and simultaneous indicator scaling. For changed units, the variance lower bound changes quadratically so it defines the same model.
- ML known-answer tests use hard assignments and independently calculated per-class normal-mean influence variance. Large-origin public tests compare against the centered version of the **actual representable shifted dataset**, across both variance specifications and both variance estimators, including pair contrasts.
- `adversarial.R/.csv/.log`: 70 independent F/chi-square/convolution probabilities, with zero-statistic repeated mixed-sign weights at ratios `1e-12` to `1e24`, near-zero positive and negative statistics, and one-sided single weights. **Zero errors; maximum absolute probability error 5.475176e-12.** Permanent tests additionally cover repeated weights spanning `1e200` and a typed refusal beyond representable relative precision.
- `direct-lmr-probe.R/.log`: public translated-direct-effect and scaled-indicator comparison preserves LR `452.1687`, reference mean `4.220817`, and reference SD `3.793443`.
- `original-probes.R` copies the original review's public/numerical probe without its artifact-overwrite statement. `original-probes-fixed.log` retains the repeated results; original review evidence remains unchanged.

## Limits

- Extremely imbalanced nonzero weights beyond representable relative precision are refused with a numerical-probability condition; they are not silently dropped.
- A reversed-likelihood public table policy is tested by injecting the reversed value at the scoring adapter after genuine fitted models pass data validation. No naturally occurring reversed, regular neighboring public fit was claimed or required to establish the helper/policy repair.
- Reported absolute means and confidence endpoints near `1e12` necessarily have the precision of doubles at that origin; uncertainty and pair contrasts retain their conditioned-coordinate precision.
- This scope report covers the two numerical source files and their tests. The parent review owns full-suite, full-package, integration, and tracking-log verification.
