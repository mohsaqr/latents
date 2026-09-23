# Missing indicators and full covariance: genuine Mplus comparison

Both models were fitted independently using native R and the local genuine
**Mplus VERSION 9 DEMO (Mac)** executable on 2026-09-17. Every retained Mplus run
terminated normally and replicated its best likelihood. This extends the earlier
validation to full residual covariance matrices and observed-data maximum
likelihood for partially missing indicators together.

Data: 1,200 individuals nested in 60 observed groups, two continuous indicators,
two individual profiles, and two discrete group classes. Starting with
`../generated/synthetic.dat`, replace `y2` by `y2 + 0.35*y1`; make `y1` missing for
rows `seq(1, 1200, 5)`; make `y2` missing for rows `seq(2, 1200, 7)` only when `y1`
is observed. Thus 240 `y1` and 138 `y2` values are missing, and no row is entirely
missing. Saved data uses `-999` as the missing-value marker. The package separately
tests fully missing individuals; Mplus comparisons here avoid its differing
row-exclusion policy for such individuals.

`full.inp` fits profile-specific covariance matrices (13 free parameters), while
`equal.inp` labels both variances and their covariance equal across profiles
(10 free parameters). Both use ML, the same measurement invariance across group
classes, and `cw#1 ON cb` for group-specific profile proportions. R uses 20 starts,
and Mplus uses 100 initial / 20 final starts. R estimates are not supplied to
Mplus as starting values. The variance/eigenvalue bound is inactive.

| Maximum absolute difference | Varying covariance | Equal covariance |
|---|---:|---:|
| Log likelihood | 2.58e-5 | 2.16e-6 |
| Means | 5.25e-8 | 2.15e-8 |
| Covariance entries | 7.11e-8 | 3.61e-8 |
| Profile probabilities | 1.16e-8 | 1.20e-9 |
| Group probabilities | 8.87e-11 | 5.85e-12 |
| Individual posteriors | 1.60e-7 | 5.85e-8 |
| Group posteriors | 1.21e-10 | 1.90e-12 |

Mplus `SAVEDATA RESULTS` retains eight significant digits, explaining the larger
absolute likelihood differences. Posterior values are saved to twelve decimal
places. Parameter ordering was inspected in retained `TECH1` output. The four
saved probabilities are joint CB/CW patterns 11, 12, 21, 22; comparisons marginalize
these and align both latent-class label sets and explicit observation IDs.

Run from the project root:

```r
source("equivalence/mplus/extensions-missing-full/compare.R")
```

The script verifies retained Mplus outputs, fits the R models, prints differences,
and regenerates `tests/fixtures/mplus/twolevel-missing-full-{varying,equal}.rds`.
Fixtures include original data, Mplus targets, Mplus standard errors in TECH1
order, version, and source hashes. Offline regression coverage is in
`tests/testthat/test-mplus-missing-full.R`.

To rerun Mplus itself, execute `/Applications/MplusDemo/mpdemo full.inp full.out`
and the corresponding `equal` command from this directory, redirecting console
output to the retained `*-console.txt` files. R standalone comparisons require no
Mplus installation; rerunning Mplus requires a compatible licensed/demo binary.

These two matched models validate this implementation subset; they do not
establish general Mplus parity, missing-at-random assumptions, or identification
for arbitrary covariance/missingness patterns.
