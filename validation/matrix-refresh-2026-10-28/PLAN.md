# Equivalence extension run (started 2026-10-28 evening, ~9 h, no commits)

User goal: extend the equivalence to more packages, all Latent GOLD tests,
several datasets; test, fix bugs, re-test. No questions, no commits.
Constraints: at most 2 cores, single-threaded BLAS, one heavy job at a time.

## Phases
1. [running] Rerun the 40 stale/unsaved equivalence/ test files on current code
   -> validation/matrix-refresh-2026-10-28/equivalence-results.csv
   (excluded: test-vlmr-runtime.R (3 cores), 3 runtime/infrastructure files).
   After it ends: fix run-equivalence.R line 27 bracket style (hook).
2. [done 402 q: 396 agree, 6 known (g03 LG s.e.)] Latent GOLD coverage kit: equivalence/latentgold-coverage/
   (make-kit.R, run-kit.sh done; 20 cases run; t03 fixed: each transition term
   needs its own `| State[-1]`). TODO: compare-functions.R (robust listing
   reader: tables have 1/3/5 term columns; skip the duplicate
   "Variances / Covariances continuous latent" table; K = 1 has no
   "Paired Comparisons"), compare.R, test-latentgold-coverage.R, README.
   Then add more datasets per family (second seed / real data).
3. [done 186 q / 60,017 values: 179 agree, 7 known (multilevLCA bounds/start)] lca() R-package sources: glca, multilevLCA, randomLCA, BayesLCA (EM),
   depmixS4 (mix) on poLCA datasets (carcinoma, cheating, gss82, values,
   election) -> equivalence/lca-packages/.
4. More packages for other functions where they fit (mixtools, tidySEM,
   seqHMM, nlme identities, Rmixmod) and published data (Masyn/Stata,
   lcmm paquid, LMest datasets, randomLCA myocardial, flexmix NPreg).
5. Fix any latents bug found (with regression tests), re-test.
6. Matrix: read new evidence (per-quantity comparison CSVs + rerun results),
   rebuild PDF; LEARNINGS/CHANGES/HANDOFF.

## Notes
- latents coef names: lca measurement.response.profile_k.v:c (probabilities),
  profile.probability...; lca covariates profile.coefficient.profile_k.*;
  lpa measurement.mean/variance(.shared)/covariance.profile_k.v:w;
  lta (plain) like multilpa plus profile.initial_probability /
  transition_probability.profile_i.group_class_1:profile_j;
  lta with covariates: mean.all / log_variance.all / initial / transition
  logits; growth coefficient.class_k.*, log_sigma.class_k,
  random.log_cholesky.(class_k|shared).a:a, random.cholesky.class_k.b:a,
  membership.class_2.(Intercept).

## Progress log
- LG main kit and wider kit rerun on current code (376 + 172 agree).
- Matrix reads kit CSVs (kit-evidence.csv, kit-functions.csv); kit wrapper
  testthat rows dropped to avoid double counting.
- Strength after phases 2-3: lca very strong+++++++ PPPPPP; growth strong;
  enumeration weak. Rerun results pending.
- Next: lpa-packages (Rmixmod, EMCluster, mixtools, flexmix FLXMCmvnorm),
  regression-packages (mixtools regmixEM, MoEClust, CO2data published),
  growth paquid (lcmm JSS), lta seqHMM + LMest RLMSdat, pool mitml/mitools.

## Rerun triage (2026-10-29): 20,858 expectations, 96 failed, 1 error
- test-montecarlo.R, 90 (ml_lca/ml_lpa icl/awe/clc_individual): stale test.
  Since 0.97.2 ICL/AWE/CLC charge the complete-data entropy (CB classes +
  CW profiles); the test charged the CW posterior entropy. New helper
  montecarlo_complete_entropy() computes it from Mplus's joint CB x CW
  posteriors. Rerun: 3461/3461 pass, tolerance unchanged.
- test-enumeration-mplus.R, 1 error (`.latents_collect_conditions` not
  found): environment. Socket workers (Accelerate BLAS) load the installed
  latents (0.95.1), not load_all(). With the current code installed in a
  scratch library on R_LIBS: 49/49 pass. run-equivalence.R now installs the
  current code to a temporary library and sets R_LIBS for the workers.
- test-review-1026-oracles.R, 1: sigma2 is named by class since 0.98.0
  (e907bff, intentional); unname() in the oracle.
- test-mclust-prior.R, 1: stopping slack. At tol 1e-13 means up to 1.45e-8
  and posteriors 2.1e-8 from mclust::me (bound 1e-8); with both at 1e-16
  every structure within 2.8e-10 (EII 6e-15). Test now runs both at 1e-16:
  51/51 pass.
- test-features-pair-covariance.R, 1: vvi_one_k3_l_02 posterior gap to
  Mplus 2.0069e-5 (bound 2e-5; was 1.69e-5 at 0.96.1); LL 2.76e-8 below
  Latent GOLD. Stopping slack of a flat case: at tol 1e-12 (Latent GOLD's)
  5.4e-7 to LG, 1.0e-6 to Mplus, LL 3e-11. pc_fit() now fits at 1e-12.
- test-features-ml-three-step.R, 2: sweep_077 / sweep_095 (VVI, k 4, n 50,
  5 starts) refit to HIGHER maxima than the 2026-10-20 kit (+7.48, +3.13);
  the kit's step one was a local maximum. Test restated: refit never below,
  posteriors agree where the maximum is the same, higher cases named.
- test-features-vlmr.R, 1: vvi_small_01-03 k3/k4 — both programs at
  single-unit spikes (share 1/60, variance on the floor); latents' relative
  floor (2026-10-26) now puts its spike above LG's in some. Test restated:
  degenerate pairs are VVI, and where latents is higher its fit is
  bound-active (new boundary0/boundary1 columns).
- Retest: retest.R (merges into equivalence-results.csv, first run kept as
  equivalence-results-first-run.csv).

## Enumeration kit (2026-10-29)
- First run 127 q: 119 agree, 8 disagree. faithful EEE k2 BIC: tolerance
  slip (BIC tolerance must be 2x the LL's). diabetes VVV bootstrap: latents'
  maxima above mclust's hierarchical init (k2 -2343.012 vs -2355.638
  converged; k4 -2281.257 vs -2285.168; proper solutions). mclust's EM stops
  at relative 1e-5 (default), up to 0.04 short; its vignette's statistics
  are that stopping point. Kit now: mclust at tol 1e-12 everywhere, bootstrap
  compared at mclust's solutions (lpa(start =)), vignette-corrected diabetes
  in the grid too, bootstrap.csv (means of the bootstrap statistics).

## Retest (2026-10-29)
- Serial (retest-serial.log): test-features-ml-three-step, -vlmr, mclust-prior,
  montecarlo, review-1026-oracles: 11,188 expectations, 0 failed (first run:
  11,187, 95 failed). The 29 VLMR warnings are low-Gaussian-mass diagnostics
  raised outside test_that(), as in the first run (not in the results table).
- Two workers (retest-workers.log): test-features-pair-covariance,
  test-enumeration-mplus.
- Kits added after the triage: lta-packages + data_market_sim (Pennoni et al.
  2024): 52 q, 48 agree, 4 known; growth-published + lcmm data_hlme: 37 q, 33
  agree, 4 known (data_hlme k3: lcmm's quasi-separated solution 0.0017 higher).
