# Protocol: diagnostic for rare, weakly separated group classes

Written 1 October 2026, before any diagnostic value was computed. Motivation:
in the `rare_weak` simulation cell ([ADDITIVE_SIMULATION.md](ADDITIVE_SIMULATION.md))
Wald intervals were available in 99.5% of datasets but covered as little as
0.544, with no refusal or warning.

## Statistic (chosen from theory, not from data)

For group class h with weight estimate w_h and observed-information variance
Var(w_h), the **effective number of groups** is

    effective_groups_h = w_h^2 (1 - w_h) / Var(w_h)

With every group's class known, Var(w_h) = w_h (1 - w_h) / J, and
effective_groups_h = J w_h, the expected number of groups in the class.
Classification uncertainty inflates Var(w_h) (the fraction of missing
information), so the statistic falls when a class is small, poorly separated,
or both. The fit-level statistic S is its minimum over classes. It needs no
extra model: only the fitted weights and the observed information already
used for standard errors. Not defined for one class.

Comparators, reported but not shipped unless they clearly dominate: relative
entropy over groups, and the smallest diagonal average posterior probability
(AvePP).

## Data

The existing registry cells, replicates 1-1000, refitted with identical seeds
(replicate 7 of `rare_weak` reproduced the stored estimates to 5.6e-16) and
joined to the stored coverage records in `additive-simulation-full-replicates.csv`.
Datasets whose intervals were withheld are excluded: there is nothing to warn
about.

- **Reference cells** (correctly specified, well separated): `regular_varying`,
  `regular_shared`, `unequal_sizes`, `singleton_admixture`, `three_classes`.
- **Target cell:** `rare_weak`.
- Other cells are reported descriptively.

## Threshold rule (training: replicates 1-500)

Grid for S: 5, 10, 15, 20, 25, 30, 40, 50, 75, 100. Choose the **largest**
grid value whose flag rate in the pooled reference cells is at most 1%.
Comparators use the same rule on their own grids (entropy: 0.5-0.9 by 0.05,
flag when below; AvePP: 0.6-0.95 by 0.05, flag when below).

## Evaluation (held out: replicates 501-1000)

Declared success criteria for shipping the warning:

1. flag rate in the pooled reference cells <= 2%;
2. flag rate in `rare_weak` >= 80%;
3. pooled over all cells, coverage of the class means and weights among
   **unflagged** datasets is closer to 0.95 than among flagged datasets.

Reported: per-cell flag rates, coverage of class means and weights by flag
status with Monte Carlo SEs. If a criterion fails, the result is reported and
the warning is not shipped as a validated diagnostic.
