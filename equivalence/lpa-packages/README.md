# lpa() against independent Gaussian mixture implementations

Every case is fitted by `lpa()` and by independent maximum-likelihood
implementations of the same covariance structure, and compared on the
likelihood, shares, means, (co)variances and every row's posteriors.

```
Rscript equivalence/lpa-packages/compare.R      # comparison.csv, search.csv, RUN.txt
```

`equivalence/test-lpa-packages.R` runs the same comparison under testthat.

**Result (2026-10-29): 156 quantities (16,048 values) compared, 155 agree,
1 known difference (Rmixmod returns no solution for one model), 0
unexplained.**

| structure | references (settings) |
|---|---|
| VVI | Rmixmod `Gaussian_pk_Lk_Bk` (20 tries, `epsilonInAlgo = 1e-12`) |
| EEI | Rmixmod `Gaussian_pk_L_B` |
| VVV | Rmixmod `Gaussian_pk_Lk_Ck`; EMCluster (20 `init.EM` starts, `em.eps = 1e-12`); mixtools `mvnormalmixEM(arbvar = TRUE)` (20 random starts, `epsilon = 1e-12`) |
| EEE | Rmixmod `Gaussian_pk_L_C`; mixtools `mvnormalmixEM(arbvar = FALSE)` |

Datasets: `faithful` (VVV and EEE at 2, VVI and EEI at 3), `iris` (VVV, EEE and
EEI at 3, VVI at 2), mclust's `diabetes` (VVV and VVI at 3), mclust's
`banknote` (VVV and EEE at 2), MASS's `crabs` (VVV and EEE at 2) and mclust's
`thyroid` (VVI and VVV at 3).

Tolerances, fixed before any comparison was run: the log likelihood to 1e-6
plus 1e-10 of its size; shares, means and (co)variances to 1e-4 of
max(1, |value|); posteriors to 1e-4.

## How the comparison is made

Gaussian mixtures have several maxima, and a VVV likelihood is unbounded at
covariances that collapse onto a few points. Each package's own search finds
its own maximum, so the estimates are compared **at each reference's
solution**: `lpa()` is started there (`start =`, one start) and must reach the
same stationary point, with the same likelihood, parameters and posteriors.
This tests the likelihood, the EM and the reported quantities without
confounding them with where two searches stopped. `lpa()`'s own search is
recorded beside it in `search.csv`.

- **EMCluster** keeps starts whose covariance collapses onto a few points
  (infinite or implausibly high likelihoods); only starts whose covariance
  eigenvalues all exceed 1e-8 of the largest variance compete.
- **flexmix** (`FLXMCmvnorm`) is not used: it estimates covariances with the
  unbiased divisor, so its fixed point is not the maximum likelihood estimate
  (faithful VVV: -1130.272 against every other program's -1130.264).

## What the search record shows (`search.csv`)

| case | lpa() default (kmeans, 50 starts) | lpa() emEM (50 starts) | best reference |
|---|---|---|---|
| banknote VVV, 2 | -729.952 (all 50 starts, and all of 500) | -718.396 | -718.396 (Rmixmod, EMCluster, mixtools) |
| crabs VVV, 2 | -1423.591 (all 50 starts, and all of 500) | -1354.157 | -1354.157 (EMCluster) |
| iris VVV, 3 | -180.185 | -75.151 (a covariance at the floor) | -180.185 |

On the two strongly correlated real datasets the default search returns the
same lower maximum from every start, so its replication count (all starts)
is no evidence of a global maximum; on crabs that maximum is even below the
nested EEE fit (-1413.511). `emEM` reaches the best maximum there, but on iris
it lands on a variance-floor spike. Started at the references' solutions,
`lpa()` reaches them exactly. The default search is frozen (2026-10-22); this
is recorded, not changed.

## Known differences (`known-differences.csv`)

- **thyroid VVV at three profiles:** Rmixmod returns no estimates (likelihood
  0, means NA) from 20 tries. EMCluster and mixtools reach latents' maximum.
