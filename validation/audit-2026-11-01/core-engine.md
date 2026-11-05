# Audit 2026-11-01 — slice `core-engine`

Files: R/kernel-em.R, R/kernel-blocks.R, R/kernel-structures.R,
R/covariance-structure.R, R/gaussian-moments.R, R/start-pool.R,
R/multilpa-starts.R, R/multilpa-pool.R, R/acceleration.R, R/prior.R,
R/scaling.R, R/centering.R (121 top-level functions).

Checking scripts (synthetic data, run with `OMP_NUM_THREADS=1
VECLIB_MAXIMUM_THREADS=1 NOT_CRAN=true`, `devtools::load_all()`) are in the
session scratchpad
/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/f561d1fc-3262-4e98-94dd-fd3bab536edf/scratchpad/core/
(not kept in the repository): `mstep.R` (M-step vs `mclust::mstep`, objective
lattice), `eve.R` (EVE/VVE full solve), `lattice-fit.R` (whole fits, 23
nesting pairs, warm-started), `likelihood.R` (two-level likelihood and
posteriors by brute force, missing values, full covariances), `structures.R`,
`structures2.R` (mixture structures, scores vs numerical gradient, mixing
M-step vs `optim`), `fb.R` (forward-backward vs path enumeration),
`missing.R` (FIML M-step vs `optim` of the observed-data likelihood),
`bounded.R` (bounded M-steps vs `constrOptim` and a stricter floor),
`extra.R` (Gaussian + ordinal + count likelihood by hand), `logit.R`,
`misc1.R`–`misc3.R` (validation paths, centring/scaling by hand, helpers,
EM driver on a toy fixed point). The mclust equivalence files were also run
(see the end).

## 1. Function table

| function | file:line | checked how | verdict |
|---|---|---|---|
| `.latents_em_settings` | R/kernel-em.R:48 | read; covered by test-kernel.R | ok |
| `.latents_observed_values` | R/kernel-em.R:88 | ran synthetic (`misc2.R`: NA cells and a zero-weight row → 1) | ok |
| `.latents_replication_tolerance` | R/kernel-em.R:111 | read | ok |
| `.latents_em_threshold` | R/kernel-em.R:122 | read; covered by test-review-1026-kernel-r2.R (K2-6) | ok |
| `.latents_em_decrease` | R/kernel-em.R:131 | read; covered by test-kernel.R | ok |
| `.latents_em_step` | R/kernel-em.R:151 | ran synthetic (`misc3.R`, NaN step) | logged |
| `.latents_em` | R/kernel-em.R:197 | ran synthetic (`misc3.R`: plain EM 24 steps to x = 2; SQUAREM 5; pause at 3 + resume gives identical history and count) | ok |
| `.latents_squarem_cycle` | R/kernel-em.R:289 | read (S3 step length, backtracking, fallback); ran via `misc3.R` | ok |
| `.latents_run_starts` | R/kernel-em.R:367 | read; covered by test-start-helpers.R, test-kernel.R | ok |
| `.latents_prepare_starts` | R/kernel-em.R:397 | read; covered by test-start-helpers.R | ok |
| `.latents_start_attempt` | R/kernel-em.R:415 | read; covered by test-start-helpers.R | ok |
| `.latents_start_failure` | R/kernel-em.R:440 | read; covered by test-start-helpers.R | ok |
| `.latents_select_starts` | R/kernel-em.R:455 | read; covered by test-kernel.R | ok |
| `.latents_run_screened` | R/kernel-em.R:511 | read; covered by test-start-helpers.R (mocked mapper) | ok |
| `.latents_screen_attempt` | R/kernel-em.R:552 | read; covered by test-start-helpers.R | ok |
| `.latents_rank_screened` | R/kernel-em.R:576 | read; covered by test-start-helpers.R | logged |
| `.latents_select_screened` | R/kernel-em.R:609 | read; covered by test-start-helpers.R | ok |
| `.latents_quasi_newton` | R/kernel-em.R:685 | ran synthetic (`misc3.R`: quadratic with a lower bound → (1, 0.5), converged) | ok |
| `.latents_gaussian_log_density` | R/kernel-blocks.R:28 | ran synthetic (`likelihood.R`: equals the moments path to 1.8e-15) | ok |
| `.latents_measurement_log_density` | R/kernel-blocks.R:66 | ran synthetic (`likelihood.R`, `extra.R`: Gaussian + ordinal + count likelihood equals a hand sum to 2e-13; zero-column case gives zeros) | ok |
| `.mixture_row_lse` | R/kernel-structures.R:34 | ran synthetic (`fb.R`: a -Inf row and a 1000 row) | ok |
| `.mixture_log_softmax` | R/kernel-structures.R:54 | ran via `structures*.R` | ok |
| `.mixture_normalize_log_weights` | R/kernel-structures.R:64 | ran via `structures*.R` | ok |
| `.mixture_weigh` | R/kernel-structures.R:81 | ran synthetic (`structures2.R`: integer weights = duplication, diff 0) | ok |
| `.mixture_structure_expectation` | R/kernel-structures.R:114 | ran synthetic (`structures.R`: two-level and group likelihood vs brute force, diff 0) | ok |
| `.mixture_two_level_design` | R/kernel-structures.R:189 | ran via `structures.R` | ok |
| `.mixture_solve` | R/kernel-structures.R:201 | ran synthetic (`misc2.R`: singular system) | logged |
| `.mixture_multinomial` | R/kernel-structures.R:235 | ran synthetic (`structures.R`: covariate at 2000, objective -54.874 vs BFGS -54.880) | ok |
| `.mixture_multinomial_newton` | R/kernel-structures.R:267 | read (gradient, information blocks, ordering); ran via `.mixture_multinomial` | ok |
| `.mixture_update_mixing` | R/kernel-structures.R:332 | ran synthetic (`structures*.R`: Q at the M-step equals `optim` Q to 6e-14, all three nestings) | ok |
| `.mixture_structure_scores` | R/kernel-structures.R:389 | ran synthetic (score sums vs central differences ≤ 1.3e-9, all three nestings) | ok |
| `.latents_forward_backward` | R/kernel-structures.R:454 | ran synthetic (`fb.R`: log likelihood vs path enumeration 2e-16, unequal lengths; posteriors sum to 1) | ok |
| `.latents_scaled_forward_backward` | R/kernel-structures.R:503 | ran synthetic (`fb.R`: 4e-16 vs paths; messages equal the log-space ones to 9e-16) | ok |
| `.latents_log_sum_exp` | R/kernel-structures.R:561 | ran synthetic (`fb.R`) | ok |
| `.multilpa_nested_structure` | R/kernel-structures.R:593 | ran synthetic (`likelihood.R`: likelihood and posteriors vs brute force; renumbered groups) | fixed |
| `.multilpa_logit_structure` | R/kernel-structures.R:660 | ran synthetic (`logit.R`: likelihood vs brute force; renumbered groups) | fixed |
| `.multilpa_weigh` | R/kernel-structures.R:711 | ran synthetic (test-audit-1101-core-engine.R) | ok |
| `.mixture_structure_pack` | R/kernel-structures.R:734 | ran synthetic (`structures.R`: names, round trip exact) | ok |
| `.mixture_structure_unpack` | R/kernel-structures.R:769 | ran synthetic (round trip exact) | ok |
| `.multilpa_structures` | R/covariance-structure.R:21 | read | ok |
| `.multilpa_is_ellipsoidal` | R/covariance-structure.R:35 | read | ok |
| `.multilpa_inferable_structures` | R/covariance-structure.R:47 | read | ok |
| `.multilpa_resolve_structure` | R/covariance-structure.R:62 | ran synthetic (`misc2.R`: round trip with `.multilpa_structure_arguments` for all 14) | ok |
| `.multilpa_structure_parameters` | R/covariance-structure.R:124 | ran synthetic (= `mclust::nVarParams` for 14 structures x d 1..4 x k 1..4, 224 cells) | ok |
| `.multilpa_structure_variances` | R/covariance-structure.R:158 | ran synthetic (`mstep.R`: = `mclust::mstep` ≤ 2e-12; `bounded.R`: EVI floor vs `constrOptim`) | ok |
| `.multilpa_evi_variances` | R/covariance-structure.R:202 | read (Celeux & Govaert [λB_k]); ran via `mstep.R` | ok |
| `.multilpa_vei_variances` | R/covariance-structure.R:227 | read ([λ_k B]); ran via `mstep.R` (2e-12 from mclust) | ok |
| `.multilpa_scatter_diagonals` | R/covariance-structure.R:267 | read | ok |
| `.multilpa_check_structure_inference` | R/covariance-structure.R:284 | read | ok |
| `.multilpa_project_start` | R/covariance-structure.R:313 | ran synthetic (`lattice-fit.R`, `misc1.R`: looser start into EEI/EEE refused by class `latents_bad_start`) | ok |
| `.multilpa_scaled_determinant` | R/covariance-structure.R:364 | read; ran via EVV/VEE/EVE | ok |
| `.multilpa_structure_covariances` | R/covariance-structure.R:384 | ran synthetic (`mstep.R`, `eve.R`, `bounded.R`) | ok |
| `.multilpa_eev_covariances` | R/covariance-structure.R:449 | ran synthetic (= mclust 6e-15) | ok |
| `.multilpa_vev_covariances` | R/covariance-structure.R:472 | ran synthetic (= mclust 4e-12) | ok |
| `.multilpa_evv_covariances` | R/covariance-structure.R:512 | ran synthetic (= mclust 4e-15) | ok |
| `.multilpa_vee_covariances` | R/covariance-structure.R:534 | ran synthetic (= mclust 9e-9) | ok |
| `.multilpa_eve_covariances` | R/covariance-structure.R:578 | ran synthetic (`eve.R`: full solve below mclust's objective by 7.7e-6 (EVE) and 5.8e-5 (VVE); the 5-step GEM above it by 1.5e-4 / 3.6e-4, by design) | ok |
| `.multilpa_reset_structure_log` | R/covariance-structure.R:662 | read | ok |
| `.multilpa_warn_structure` | R/covariance-structure.R:674 | read | ok |
| `.multilpa_report_structure_log` | R/covariance-structure.R:690 | read; covered by test-covariance-structure.R | ok |
| `.multilpa_structure_arguments` | R/covariance-structure.R:714 | ran synthetic (round trip, 14 structures) | ok |
| `.multilpa_seed_block` | R/covariance-structure.R:741 | read | ok |
| `.multilpa_orientation_steps` | R/covariance-structure.R:763 | read | ok |
| `.multilpa_eigenvalue_design` | R/covariance-structure.R:770 | ran synthetic (dims VVI 9x9, EVI 9x7, VEI 9x5 for k = d = 3) | ok |
| `.multilpa_bounded_variances` | R/covariance-structure.R:800 | ran synthetic (`bounded.R`: EVI with a binding floor 78.05757 vs `constrOptim` 78.05761) | ok |
| `.multilpa_bounded_orientation` | R/covariance-structure.R:939 | read (Jacobi angle derived by hand: exact pairwise minimiser); ran via `bounded.R` | ok |
| `.multilpa_bounded_covariances` | R/covariance-structure.R:971 | ran synthetic (`bounded.R`: VEV/EVV/VEE/EVE/VVE feasible at the floor, family kept, objective below the 1.5x-floor solution) | ok |
| `.multilpa_gaussian_moments` | R/gaussian-moments.R:6 | ran synthetic (`likelihood.R`: marginal densities by hand with missing patterns and an all-missing row) | ok |
| `.multilpa_bound_covariance` | R/gaussian-moments.R:58 | ran synthetic (eigenvalues 2, 0.1 from a singular matrix) | ok |
| `.multilpa_maximize_moments` | R/gaussian-moments.R:86 | ran synthetic (`mstep.R` all 14 structures vs mclust; `missing.R`: FIML fixed point, `optim` improves by 0) | ok |
| `.multilpa_covariance_boundary` | R/gaussian-moments.R:175 | ran synthetic | ok |
| `.latents_resolve_min_variance` | R/gaussian-moments.R:204 | ran synthetic (`missing.R`: 1e-6 x smallest variance, zero-weight rows ignored, invariant to a 1e9 shift) | ok |
| `.latents_execute_start` | R/start-pool.R:5 | read; covered by test-parallel-candidates / profile-pool (not run: they start workers) | ok |
| `.latents_start_tasks` | R/start-pool.R:20 | read | ok |
| `.latents_pool_phase` | R/start-pool.R:53 | read | ok |
| `.latents_pool_setup` | R/start-pool.R:79 | read | ok |
| `.latents_pool_progress` | R/start-pool.R:103 | read | ok |
| `.latents_count` | R/start-pool.R:117 | ran synthetic | ok |
| `.latents_format_duration` | R/start-pool.R:126 | ran synthetic (0, 59.6, 3599.9, 7300 s) | ok |
| `.multilpa_run_start` | R/multilpa-starts.R:8 | ran via every `multilpa()` fit above | ok |
| `.multilpa_start_interior` | R/multilpa-starts.R:30 | read | ok |
| `.multilpa_run_plan` | R/multilpa-starts.R:39 | ran via every `multilpa()` fit above | ok |
| `.multilpa_report_failed_plan` | R/multilpa-starts.R:57 | read | ok |
| `.multilpa_finish_plan` | R/multilpa-starts.R:66 | read | ok |
| `.multilpa_finalise` | R/multilpa-starts.R:75 | read (criteria, counts, weights, noise); ran via fits; `missing.R`/`extra.R` log likelihoods recomputed by hand | ok |
| `.multilpa_pool_fallback` | R/multilpa-pool.R:5 | read | ok |
| `.multilpa_pool_width` | R/multilpa-pool.R:20 | ran synthetic (continuous 2 + categorical 1 + ordinal 2 = 5) | ok |
| `.multilpa_execute_start` | R/multilpa-pool.R:36 | read | ok |
| `.multilpa_pool_loss` | R/multilpa-pool.R:50 | read | ok |
| `.multilpa_setup_task` | R/multilpa-pool.R:60 | read | ok |
| `.multilpa_start_tasks` | R/multilpa-pool.R:79 | read | ok |
| `.multilpa_enumerate_fits` | R/multilpa-pool.R:98 | read; covered by golden enumeration cases | ok |
| `.multilpa_accelerated_blocks` | R/acceleration.R:13 | read | ok |
| `.multilpa_flatten_parameters` | R/acceleration.R:25 | ran synthetic (round trip 9e-16) | ok |
| `.multilpa_accelerator` | R/acceleration.R:43 | ran synthetic (round trip; reference scales = column sd) | ok |
| `.latents_reference_scales` | R/acceleration.R:65 | ran synthetic (zero-weight rows at 1e6 ignored) | ok |
| `.latents_scale_chart` | R/acceleration.R:87 | ran via `.multilpa_accelerator` | ok |
| `.multilpa_restore_parameters` | R/acceleration.R:113 | ran synthetic (round trip) | ok |
| `.multilpa_valid_parameters` | R/acceleration.R:139 | ran synthetic (negative variance and non-normalised probabilities refused) | ok |
| `.lta_accelerator` | R/acceleration.R:199 | ran synthetic (`misc2.R`) | fixed |
| `prior_control` | R/prior.R:96 | ran synthetic (defaults print; asymmetric and vector `scale` refused with `latents_bad_argument`) | ok |
| `print.latents_prior` | R/prior.R:127 | ran synthetic | ok |
| `.multilpa_positive_scale` | R/prior.R:149 | ran synthetic | ok |
| `.multilpa_prior_structures` | R/prior.R:160 | read | ok |
| `.multilpa_prior_parameters` | R/prior.R:177 | read vs `mclust::defaultPrior()`; covered by equivalence/test-mclust-prior.R | ok |
| `.multilpa_prior_maximize` | R/prior.R:228 | covered by equivalence/test-mclust-prior.R (one M-step = `mclust::me` to 1e-8, 10 structures); bounded EVI/VEI pseudo-weights derived by hand | ok |
| `.multilpa_prior_vei` | R/prior.R:342 | covered by equivalence/test-mclust-prior.R | ok |
| `.multilpa_prior_vev` | R/prior.R:373 | covered by equivalence/test-mclust-prior.R | ok |
| `.multilpa_check_prior` | R/prior.R:416 | ran synthetic (`misc1.R`: VEE, bad class, FIML, ordinal-only, weights each refused by class) | ok |
| `.multilpa_design_scale` | R/scaling.R:22 | ran synthetic (RMS by hand; zero column → 1) | ok |
| `.multilpa_covariance_coordinate_scale` | R/scaling.R:62 | ran synthetic (1 1 1 / 1 4 1 by hand) | ok |
| `.multilpa_centering_schemes` | R/centering.R:17 | read | ok |
| `.multilpa_scaling_schemes` | R/centering.R:27 | read | ok |
| `.multilpa_center_indicators` | R/centering.R:56 | ran synthetic (`misc2.R`: person offsets/scales by hand ≤ 9e-16 with an unordered group index and NA; grand = `scale()`; zero-weight group excluded) | ok (comment corrected) |
| `.multilpa_indicator_scales` | R/centering.R:107 | ran synthetic (as above) | ok |
| `.multilpa_center_only` | R/centering.R:185 | ran synthetic (as above) | ok |
| `.multilpa_uncentered_indicators` | R/centering.R:223 | ran synthetic (round trip exact) | ok |
| `.multilpa_center_like` | R/centering.R:248 | ran synthetic (equals the fitted transform exactly) | ok |
| `.multilpa_is_transformed` | R/centering.R:270 | ran synthetic | ok |
| `.multilpa_check_centered` | R/centering.R:301 | ran synthetic (constant-within-unit → `latents_bad_data`) | ok |
| `.multilpa_check_atomic_groups` | R/centering.R:356 | ran synthetic (2 of 3 groups, 57.1% of rows) | ok |

Totals: 121 functions; 3 fixed, 3 logged, 115 ok (one with a corrected
comment), 0 needs-other-slice.

## 2. Fixed

**F1. Group-level rows of the nested and logit structures followed first
appearance, not the group index** (R/kernel-structures.R:610, :624, :676,
:685). `.multilpa_nested_structure()` and `.multilpa_logit_structure()` summed
rows to groups with `rowsum(..., reorder = FALSE)`, which orders the result by
first appearance of each index value, while every other group-level input
(`group_design`, `group_log_density`, `weights`) and every consumer
(`group_posteriors[group_index, ]`) reads row `j` as group `j`. With groups
numbered 1..5 but appearing as 3, 1, 5, 2, 4 (`likelihood.R`):
`group posteriors aligned to index (row j = group j)? max diff: 0.7258485`
(nested; likelihood unchanged, posteriors and subject posteriors wrong), and
for the logit structure (`logit.R`): `permuted index: loglik diff -0.03282125
group posterior alignment: 0.1079615 subject posteriors: 0.02783102` (wrong
likelihood). Every current caller builds its index with
`match(x, unique(x))`, i.e. numbered in order of appearance, so no fitted
number changes today; the functions documented "Group of each row" without
that requirement, and the rest of the slice (`.mixture_structure_expectation`,
`.multilpa_center_only`) already sums with `reorder = TRUE` for this reason.
Change: `reorder = TRUE` at the four `rowsum()` calls (bit-identical for
indices numbered in order of appearance) and the contract stated in both
roxygen blocks. Guarded by test-audit-1101-core-engine.R: "the nested
structure puts group j in row j for any group numbering", "renumbering the
groups only permutes the nested structure's rows", "the logit structure's
likelihood does not depend on the group numbering", "the full-covariance
likelihood with missing values is the marginal one" (brute-force references).
With the old lines restored, these 4 tests fail (8 expectations).

**F2. `.lta_accelerator(general = FALSE)` dropped `row_weights`**
(R/acceleration.R:204). The documented `row_weights` argument reached only
the general branch; the homogeneous branch called
`.multilpa_accelerator(x, parameters)`, so excluded rows set the SQUAREM
reference spread (`misc2.R`: `lta non-general honours row_weights: FALSE`
with ten zero-weight rows at 1e6). The only homogeneous caller
(`.multilpa_transition_em`, R/transitions.R:385) passes no weights, so no
fitted number changes today. Change: pass `row_weights` through. Guarded by
"the homogeneous transition accelerator measures rows of positive weight"
(expected chart units 1/sd of the 40 rows of positive weight, computed
independently); with the old line restored it fails (2 expectations).

**F3 (comment only). `.multilpa_center_indicators()`** (R/centering.R:72)
said "A standard deviation does not move with a shift, so the scales are read
off the raw values whatever the centring". Under `centering = "person",
scaling = "grand"` the raw grand SD includes the between-group spread
(`misc2.R`: scale = raw SD exactly; differs from the SD of the
person-centred values by 0.907). The behaviour matches the documented
`scaling` ("its standard deviation over all observations") and is frozen; the
comment now says what happens. No test (no behaviour change).

## 3. Logged, not fixed

- `.latents_em_step` (R/kernel-em.R:151): a NaN log likelihood reaches
  `if (NA)` and fails with the unclassed "missing value where TRUE/FALSE
  needed" (with or without the decrease guard). Every engine's E-step stops
  earlier on a non-finite likelihood (`.multilpa_nested_structure`, the LTA
  engines, `latents_vanished_density`) or uses `stop_on_nonfinite`, so this
  is unreachable from realistic data; a start failing this way is still
  recorded in the starts table.
- `.latents_rank_screened` (R/kernel-em.R:589): when every screened start is
  valid but scores non-finite, the "All n starts failed:" message lists no
  reason. Unreachable while the E-steps above refuse non-finite likelihoods.
- `.mixture_solve` (R/kernel-structures.R:201): a singular, unridged system
  ends in an unclassed `qr.solve` error. The Newton system carries a ridge and
  collinear designs are refused upstream.
- `.multilpa_nested_structure` (R/kernel-structures.R:633): "Non-finite
  likelihood or posterior probabilities." is unclassed; it is caught per start
  and surfaces in `latents_all_starts_failed`.
- `.multilpa_structure_covariances` (R/covariance-structure.R:422): the
  axis-parallel fallback recomputes `.multilpa_structure_variances()` once per
  profile (and would count a VEI cap k times). No caller reaches it: diagonal
  structures go to `.multilpa_structure_variances()` directly.
- EVE/VVE inside EM take 5 warm-started orientation steps
  (`.multilpa_orientation_steps()`), so one M-step from a cold start is above
  mclust's objective by 1.5e-4 / 3.6e-4; documented as a generalized EM step
  (frozen engine). Solved to convergence the same M-step is below mclust's.
- A default search does not always reach the best maximum of a looser
  structure: VEE with 3 starts ended at -1575.992, while VEE warm-started at
  the EEE fit reached -1575.013 (`lattice-fit.R`). A search-budget matter,
  not an M-step defect; every warm-started nesting pair gained (23/23).

## 4. Needs another slice

None found.

## 5. Test runs

New file and existing slice files, each through `testthat::test_file()`
after `devtools::load_all()`, `NOT_CRAN=true`, one process:

```
test-audit-1101-core-engine.R: PASS 74 FAIL 0 SKIP 0 ERR 0
test-engine-fixes.R: PASS 26 FAIL 0 SKIP 0 ERR 0
test-covariance-structure.R: PASS 167 FAIL 0 SKIP 0 ERR 0
test-centering.R: PASS 36 FAIL 0 SKIP 0 ERR 0
test-acceleration.R: PASS 29 FAIL 0 SKIP 0 ERR 0
test-lta-acceleration.R: PASS 39 FAIL 0 SKIP 0 ERR 0
test-kernel.R: PASS 50 FAIL 0 SKIP 0 ERR 0
test-review-1026-kernel.R: PASS 18 FAIL 0 SKIP 0 ERR 0
test-review-1026-kernel-r2.R: PASS 42 FAIL 0 SKIP 0 ERR 0
test-prior.R: PASS 25 FAIL 0 SKIP 0 ERR 0
test-scaling.R: PASS 31 FAIL 0 SKIP 0 ERR 0
test-scaling-invariance.R: PASS 37 FAIL 0 SKIP 0 ERR 0
test-start-helpers.R: PASS 23 FAIL 0 SKIP 0 ERR 0
test-default-search.R: PASS 28 FAIL 0 SKIP 0 ERR 0
test-api-compat-0921.R: PASS 44 FAIL 0 SKIP 0 ERR 0
test-compat-0921.R: PASS 10 FAIL 0 SKIP 0 ERR 0
test-covariates.R: PASS 23 FAIL 0 SKIP 0 ERR 0
test-weights.R: PASS 76 FAIL 0 SKIP 0 ERR 0
test-review-1020-zero-weights.R: PASS 18 FAIL 0 SKIP 0 ERR 0
```

Break checks (old lines restored in both files at once, then restored from a
saved copy and confirmed byte-identical with `cmp`):
`BROKEN: PASS 64 FAIL 10 ERR 0` — the four F1 tests (8 expectations) and
the F2 test (2) failed. The `sed` used for the break also flipped the two
`reorder = TRUE` rowsums of `.mixture_structure_expectation()` (lines 159,
163), which none of these tests calls; the file was restored from the saved
copy.

Not run: test-search.R and test-parallel-profile-pool.R, which start socket
workers (`cores = 2`), forbidden in this audit.

mclust equivalence files (not shipped), after `load_all()`, after the fixes:

```
equivalence/test-mclust.R PASS 41 FAIL 0 SKIP 0 ERR 0
equivalence/test-mclust-prior.R PASS 51 FAIL 0 SKIP 0 ERR 0
equivalence/test-mclust-noise.R PASS 14 FAIL 0 SKIP 0 ERR 0
```

(Run before the fixes and again after them, with the same counts.)

Engine fingerprint (`Rscript equivalence/engine-golden/compare.R` from the
repository root, whole working tree including other slices' edits):

```
identical: 247; within_tolerance: 0; DIFFERENT: 13; ERROR: 0
```

All 13 differences are one attribute of the bootstrap result renamed
`replicate_rows` -> `replicate_identity`, which only R/bootstrap-inference.R
(another slice, modified in the working tree) contains: profiles_two_level_VEI,
profiles_lca_single, profiles_count_two_level, lta_homogeneous_bootstrap,
lta_general_bootstrap, lta_general_weights_bootstrap and seven
svc_inference_* cases, each at `.../inference_bootstrap/...` or `probes`. To
separate them from this slice, those cases (pattern
`^(profiles_two_level_VEI|profiles_lca_single|profiles_count_two_level|lta_homogeneous_bootstrap|lta_general_bootstrap|lta_general_weights_bootstrap|svc_inference_)`,
15 cases) were rerun on a scratch copy of `git archive HEAD` with only this
slice's three edited files laid over it:

```
identical: 15; within_tolerance: 0; DIFFERENT: 0; ERROR: 0
```

So no fingerprint case changes because of this slice. The run rewrote
equivalence/engine-golden/last-compare.csv (a tracked file; another agent's
golden run writes it too).

roxygen: `roxygen2::roxygenise()` on two scratch copies of the current tree,
one with this slice's three files at HEAD, gave identical man/ and NAMESPACE
(the edited blocks are `@noRd`), and both equal the repository's man/ and
NAMESPACE, so `devtools::document()` was not run in place: it would change
nothing.
