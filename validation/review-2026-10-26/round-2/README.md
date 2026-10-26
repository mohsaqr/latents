# Round 2 of the 2026-10-26 review: every formula re-verified after the fixes

Nine reviewers (one per area, the partition of round 1) re-derived every
algorithm in `R/` against an independent check: brute-force enumeration,
numerical derivatives of verified analytic scores, closed forms, or an installed
reference package. Each then fixed what it found, with a regression test shown
to fail on the pre-fix code. The area reports in this folder hold the findings,
the checks that passed (with the largest discrepancy seen), and a "Fixes
applied" section.

## What round 2 found

No error in the model mathematics. Every finding was in the machinery around
it, and almost all of one of three kinds:

1. **A convention not carried everywhere.** The morning's "zero-weight units
   are absent" rule (ML-1) had not reached the cross-level fits (log
   likelihood +90, BIC −176; CL-1), grand centring and scaling, the regression
   and additive classification tables, `drop_empty = FALSE`, the LTA general
   path, or the additive identification checks and starting values.
2. **An absolute constant compared with a quantity in data units.** The VEI and
   EVE/VVE shape floors (`min_variance` against a unitless shape: 13x off,
   every start failing for a rate next to a population; K2-1), the default
   variance floor (iris in metres: VVV 1.0, VVI 414 and EEE 465 log-likelihood
   units low; K2-2), the EM stop `tol (1 + |LL|)`, the separation flag at
   |coefficient| > 30 (fires on a calendar year), the growth random-effect
   floors and degeneracy rule (time in days: slope variances pinned, SEs 2x
   off; GR2-1), the membership Newton step and the differencing chart of the
   regression information (N1, N6), the LTA covariate chart (LTA-F5).
3. **A quantity treated as known that is estimated.** Wald standard errors of
   centred fits treated the grand mean and SD as constants: 22-31% too small,
   and the robust path that would correct person centring was refused
   (INF2-1, INF2-2). Fixed by stacking the centring equations into the
   information; Monte Carlo calibration now 0.90-1.09 of the true spread.

Three were conventions that gave a wrong answer silently: ICL, AWE and CLC
under the "groups" convention charged only the group-class entropy, which is
zero for every single-level fit, so `icl_groups` equalled BIC and `clc_groups`
always chose the largest model (S-1); under sampling weights the classification
estimates mixed weighted and unweighted quantities (S-3, S-4); and a
half-observed categorical pair was dropped under FIML, a MAR bias of 50 SE at
n = 40,000 (C2-2). Two were copies of round-1 defects that the first fix had
missed in another engine (the shared negative-binomial dispersion rule, 23.5
log-likelihood units low through `lta()`, C2-1; the raw ridge in the membership
Newton step, N1).

| Area | Findings | High | Fixed | Open |
|---|---|---|---|---|
| kernel | 7 | 0 | 7 | SQUAREM steplength is unit-dependent (0.0003-0.018 SE); observation |
| multilevel | 12 | 1 (CL-1) | 12 | — |
| categorical | 6 | 1 (C2-1 through `lta()`) | 6 | parametric-bootstrap BVR p-value deferred; the chi-square is documented as conservative (0.25% / 0.08% at nominal 5%) |
| inference | 6 + round-1 INF-3 | 2 (INF2-1, INF2-2) | 7 | — |
| selection | 7 | 1 (S-1) | 7 | S-5 kept as a documented convention (KIC = −2L + 3(q + 1)); S-7 the soft-entropy convention, documented |
| mixture regression | 8 | 0 | 8 (+ the bootstrap rule) | — |
| growth | 8 | 1 (GR2-1) | 8 (+ a seed regression found and fixed) | — |
| additive | 8 | 0 | 8 | a weighted class collapse that stops 1.1e-4 apart is not flagged (its inference is refused); ADD2-6 boundary fits documented 0.0015 SE short |
| LTA | 6 | 0 (2 medium) | see lta.md | — |

Decisions made for the fixes (all by the project owner, 2026-10-26): zero-weight
units are absent from every count, estimate, start and refusal, in every family;
the bootstrap LRT p-value is taken over the valid replicates, with the invalid
count reported, and is NA only below half; ICL/AWE/CLC charge the complete-data
classification entropy under both conventions (transition fits keep one entropy
per level); every classification estimate is sampling-weighted; `min_variance`
defaults to 1e-6 times the smallest indicator variance in the fitted units;
growth mixtures are fitted on standardised outcome and predictors and mapped
back exactly; the EM stop is `tol × max(1, observed values)`.

## Verification

See the section of the same name in `../README.md`, written after the final
runs on the finished tree.
