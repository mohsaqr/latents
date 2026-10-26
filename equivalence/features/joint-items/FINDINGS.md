# Findings: joint items (`local_dependence`, categorical)

## 1. EM reports convergence short of the maximum (not fixed: released engine)

The production contract is a relative log-likelihood gain per step or
SQUAREM cycle. It does not bound the remaining objective gap, establish
stationarity in a flat direction, or guarantee a global maximum. That
contract and the production stopping rule remain unchanged.

- Original-recipe tests: `equivalence/test-features-joint-items.R`, formerly
  "latents stops at the maximum it reports as converged", with per-quantity
  comparisons of `twopair_mod_05` (8 historical failures), plus "the reference
  study" (1 historical failure). The current strict-recipe test retains the
  independent objective-gap oracle and separately records the original fits.
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
- Current exact-kit verification (2026-10-26): the original 15-start,
  two-seed recipe at tol = 1e-10 gives -2916.757139377, with independent
  remaining gap 0.0001276719795 and posterior gap 0.00782306434. Both original
  seed fits, their convergence flags and their controls remain separate
  diagnostics. Each seed's own selected parameters are then resumed with
  `starting_values()`, one supplied start, tol = min(original tol, 1e-13), and
  at least the original iteration budget (minimum 5000). No reference
  parameters or reference-dependent seeds enter this recipe.
- This strict verification recipe is uniform over every external case and
  seed, and over both joined and plain fits in the 140-dataset study. The
  independent 1e-9 relative objective-gap oracle and all external comparison
  tolerances are unchanged. The test names the strict recipe explicitly;
  it does not turn per-step convergence into a production accuracy promise.
  On the retained difficult case, a supplied-start resume takes about two
  seconds and gives -2916.757011942618, independent remaining relative gap
  8.69274e-11 and posterior gap 1.80174e-5; the retained external comparisons
  then pass without changing their tolerances.
- When the runner explicitly sets `latents.equivalence.output_dir`, fresh
  comparisons, strict-polish gaps, original-seed diagnostics and reference
  study results are written beneath that artifact directory. Normal tests
  do not overwrite retained reference reports. Diagnostics distinguish the
  independent categorical EM reference from the tighter production refit
  used for the two-level cases.

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

## 3. The poLCA reference stopped at its iteration cap (verification protocol)

The original poLCA recipe used 10 starts, the case seed, an absolute
`tol = 1e-12` and `maxiter = 20000`. On `twopair_mod_05`, all ten attempts
reached that cap; the selected attempt was still gaining about `5.7e-9`
per step. Its likelihood was `-2916.757049011803`, leaving five comparison
rows outside their unchanged tolerances. This was an iteration-limited
reference fit, rather than evidence that the two likelihoods differ.

Every one of the 49 eligible cases now uses both declared budgets: the
original 20000 iterations and a verification budget of 200000. Each run
resets the same case seed and retains the same 10 starts, absolute
`tol = 1e-12`, data, model and likelihood convention. Both use stock poLCA;
no production or independent-reference parameters initialize either run.
The verification fit supplies the comparisons. The complete original fits,
all attempt likelihoods, selected iteration counts, controls and original
comparison rows remain separate diagnostics. With an explicit artifact
output directory, the runner writes both comparison tables, cap diagnostics
and attempt tables without changing retained reports. poLCA solves a
one-class model directly; those diagnostics explicitly mark a closed-form
fit and do not invent random attempts.

On the difficult case, all verification attempts stopped after 54540--57705
iterations; the selected attempt stopped at 56767 with final gain
`9.09e-13`. Its likelihood was `-2916.757011694772`; all seven comparison
rows passed their original tolerances. An ordinary one-class fit returns
identical estimates and comparison rows under either cap. This protocol
still makes no claim that any fitted local maximum is globally optimal.
