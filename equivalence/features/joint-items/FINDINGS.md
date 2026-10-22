# Findings: joint items (`local_dependence`, categorical)

## 1. EM reports convergence short of the maximum (not fixed: released engine)

- Tests: `equivalence/test-features-joint-items.R`, "latents stops at the
  maximum it reports as converged" and the per-quantity comparisons of case
  `twopair_mod_05` (8 failures), plus "the reference study" (1 failure).
- Repro: `Rscript equivalence/features/joint-items/findings/squarem-early-stop.R`
- On `twopair_mod_05` latents stops at -2916.757139, 1.28e-4 below the
  maximum (Latent GOLD prints -2916.7570); posteriors are off by up to
  7.8e-3. Seven of the 140 study fits stop 1.8e-6 to 1.0e-5 short.
- Cause: the stopping rule in `R/kernel-em.R` stops on a small gain per
  step; on a flat ridge EM (and SQUAREM) gains little per step long before
  the maximum. It is not specific to joint items: the same pair joined by
  hand into one item behaves the same.
- Not fixed, on purpose: the rule is shared by every released model, whose
  results were checked against Mplus, Latent GOLD and mclust. A scratch
  Aitken-extrapolated rule still stopped short on this case (-2916.757109
  with SQUAREM; plain EM reached -2916.757035 at 3.3 times the iterations).

## 2. Bootstrap inference stopped with an unclassed error (fixed)

- Test: `tests/testthat/test-features-joint-items.R`, "bootstrap inference
  survives a resample without a rare combination".
- Repro: `Rscript equivalence/features/joint-items/findings/bootstrap-rare-combination.R`
- A resample that missed a combination observed only twice (or a rare
  category of a single item) refitted that item with one category fewer,
  and aligning its parameters with the original's stopped with "the two
  signatures must have the same shape".
- Fix (2026-10-21): `.multilpa_bootstrap_replicates()`
  (`R/bootstrap-inference.R`) gives up such a replicate with the reason
  "the resample does not contain every category of ..."; it is dropped and
  counted by the existing `latents_bootstrap_dropped` warning. The test now
  expects that warning and finite standard errors (59 of 59 pass). A
  resample that keeps every category is extracted as before, so the
  existing bootstrap suites are unchanged (5419 equivalence expectations
  pass).
