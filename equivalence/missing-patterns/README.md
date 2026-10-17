# Missing indicators: five patterns against Mplus and Latent GOLD

Five simulated datasets, each with one pattern of missing indicators, fitted
by this package (`missing = "fiml"`), genuine **Mplus VERSION 9 DEMO (Mac)**
and **Latent GOLD 6.1** (under Wine), all on 2026-10-17. Every model has three
continuous indicators with class-specific variances and no covariances (VVI).
Profiles are only moderately separated, so entropy (0.55 to 0.65) and AvePP
are informative rather than near one.

| Case | Design | Missingness | Rows with nothing observed |
|---|---|---|---|
| m1 | single level, 3 profiles, 600 rows | 15% of cells completely at random | 0 |
| m2 | m1 | m1, plus 10% of rows with every indicator missing | 62 |
| m3 | single level, 3 profiles, 600 rows | `y2`, `y3` missing more often as `y1` rises (MAR) | 0 |
| m4 | single level, 3 profiles, 600 rows | `y3` missing for 70% of rows | 0 |
| m5 | two level, 2 profiles by 2 group classes, 60 groups of 15 | 10% of cells, 8% of rows, and groups 7 and 31 entirely | 110 (2 groups) |

## Files

- `make-kit.R` simulates the data and writes `kit/`: Mplus data (`-999` for
  missing) and `.inp`, Latent GOLD data (`.` for missing) and `.lgs`, and
  `cases.rds` with the data and design.
- `run.sh` runs both programs (`sh equivalence/missing-patterns/run.sh`, or
  with `mplus` or `latentgold`) and writes their output to `returned/`.
- `../test-missing-patterns.R` compares the package with the returned output;
  it runs offline.

## What agrees

| Quantity | Mplus | Latent GOLD |
|---|---|---|
| Log likelihood | every case, to the 3 decimals printed | every case, to the 4 decimals printed |
| Individual posteriors, row by row | largest gap 1.9e-4 (m4, 70% of `y3` missing) | largest gap 2.0e-4; rows with nothing observed 1e-5 |
| Sample size and individual BIC | every case | — (different convention, below) |
| Individual entropy and AvePP | every case, to the 3 decimals printed | — (Latent GOLD reports a different entropy R-squared) |

The m2 and m5 entropy and AvePP agree only because classification statistics
leave out rows with nothing observed (0.9.38); with them counted, the test
fails on both cases.

## Conventions the programs differ on

Rows and groups with nothing observed carry no information about their class:
their posterior is the prior their group implies. Latent GOLD, which keeps
them, gives them the same posterior as this package (to 1e-5); Mplus leaves
them out of its saved posteriors, so its output says nothing about them. What
differs is whether they are counted.

| | Mplus | Latent GOLD (`missing includeall`) | This package |
|---|---|---|---|
| Rows with nothing observed, in N and individual BIC | excluded ("These cases were not included in the analysis") | included (N = 600 for m2, 900 for m5) | excluded |
| The same rows, in entropy and AvePP | excluded | — | excluded |
| The same rows, in class counts | excluded | included | included |
| Groups with nothing observed (m5: 7 and 31) | dropped (58 clusters) | kept (60 groups) | kept, in `n_groups`, `bic_groups` and group-level entropy and AvePP |
| Group-level entropy and AvePP | averaged over observed rows, each carrying its group's posterior | — | averaged over groups |

The group-level statistics of m5 are therefore not compared. With the two
empty groups dropped and the posteriors averaged over observed rows, this
package's own group posteriors give Mplus's printed values exactly (entropy
0.5710, AvePP 0.8566 and 0.8923 against 0.571, 0.857 and 0.892); averaged over
all 60 groups they give 0.552, 0.855 and 0.870.
