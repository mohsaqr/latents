# Findings: direct effects (`direct_effects`)

Finding 1 keeps a failing assertion (a limitation of the released start
search, which is not changed: see below). Findings 2 and 3 were fixed on
2026-10-21 in `.latents_parse_direct_effects()` (`R/direct-effects.R`) and
their test now passes.

## 1. heavy_02: latents stops 8.50 log-likelihood units below both programs

- Test: `equivalence/test-features-direct-effects.R`, "latents reaches the best
  maximum either program found" (and, as a consequence, "every case is
  compared ... with at least one program").
- Repro: `Rscript equivalence/features/direct-effects/repro/heavy-02-lower-maximum.R`
- Data: 2 profiles, VVI, n = 500, t(3) noise (heavy tails), `y3 ~ z`.
- Latent GOLD (100 sets) and Mplus (400/80 starts, best value replicated) both
  report -3098.4868. latents reports -3106.9860 from 3 seeds x 200 starts;
  nearly every start converges there (198 of 200 for seed 1). The
  independent likelihood (reference.R) evaluated at Latent GOLD's estimates
  gives -3098.48679023, exactly Latent GOLD's value, and 500 further EM steps
  leave it there: it is a regular interior maximum (smallest variance 0.455),
  not a degenerate one.
- Diagnosis: a start-search limitation, not the direct-effect algebra. A plain
  VVI fit (no direct effect) on y3 - 0.8063 z misses the same maximum
  (-3106.9967 from 49 of 50 starts), and the independent EM with k-means and
  random-partition starts misses it too. latents' starts (converged k-means
  partitions) all fall into one basin on these heavy-tailed data; the
  programs' screened random starts reach the higher one. `n_best_replicated`
  then reports the lower solution as well replicated. The Monte Carlo kit
  recorded the same narrowness for over-extracted models
  (`equivalence/montecarlo/README.md`).
- Not fixed, on purpose: the start strategy is shared by every released model,
  whose results were checked against Mplus, Latent GOLD and mclust. A scratch
  copy with random-partition starts reached -3098.4868 from every seed; it
  also changes which start wins in other models, so it is recorded here
  instead.
- Opt-in emEM (2026-10-22, `repro/heavy-02-emem.R`, `heavy-02-emem.rds`):
  `search = "emEM"` reaches -3098.487 (Latent GOLD's and Mplus's value) in all
  8 runs, from 200 starts x seeds 1-3 and from 10 starts x seeds 1-5; the
  default search stays at -3106.986 in all 8. The winning finalists are
  random-centre and random-partition starts; the solution is regular
  (smallest standardized eigenvalue 0.153, no boundary, no warning). emEM
  took 5-6 times as long (11 s against 2 s at 200 starts). The equivalence
  test still runs the default search, so its assertion still fails.

## 2. A trailing `+` in a term is accepted silently (fixed)

- Test: `tests/testthat/test-features-direct-effects.R`, "every malformed
  `direct_effects` value is refused as a bad argument".
- Repro: `Rscript equivalence/features/direct-effects/repro/malformed-terms.R`
- `direct_effects = "y3 ~ z +"` fits `y3 ~ z` with no condition:
  `strsplit("z +", "+", fixed = TRUE)` drops the trailing empty field, so the
  empty-covariate check in `.latents_parse_direct_effects()` never sees it.
  The roxygen says a term that is not `"indicator ~ covariate"` raises
  `latents_bad_argument` (a leading `"y3 ~ + z"` does).

## 3. NA or non-character `direct_effects` stop with an unclassed error (fixed)

- Same test and repro as 2.
- `direct_effects = NA_character_` and `direct_effects = 1` stop in
  `stopifnot()` with a `simpleError`, not the documented
  `latents_bad_argument`, so they cannot be caught by class.
- Fix (2 and 3): the parser refuses a value that is not a character vector
  without NA as `latents_bad_argument`, and splits each term with a trailing
  space appended, so an empty field after a final `+` or `~` is seen and
  refused.

## Reference-side observations (not latents defects)

- **farmean_01, Latent GOLD**: with the covariate centred at 1000 Latent GOLD
  ends at a degenerate solution (4 activated constraints, class-2 variances
  0, 14 parameters, LL -3244.0504) although its own best start is reported;
  Mplus and latents agree at -3069.3959 with 18 parameters. The case is
  compared with Mplus only. farmean_02 is fine in both programs.
- **Mplus refuses the unit-scaled covariates** (variance above 1e6, or "a
  variance of zero" for units of 1e-4). The kit rescales them with
  `DEFINE: z = z * s` (s = 1e-4 and 1e4) and the reader converts the effect
  and its standard error back (b = b' s).
- **Mplus frees the residual covariance** of two indicators that are both
  regressed on covariates; the kit fixes it (`y1 WITH y3@0`) in every class to
  hold local independence.
- **Mplus's log likelihood** differs from Latent GOLD's and latents' at
  agreeing estimates by up to 4.3e-3 at n = 20000 (3.5e-8 relative), of
  either sign; Latent GOLD and latents agree to 8.1e-8 (2.4e-11 relative).
- **Latent GOLD writes -1** as the standard error of the class size when k = 1.
