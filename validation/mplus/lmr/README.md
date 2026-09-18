# Genuine Mplus TECH11 comparison

Two `TECH11` runs on the published Example 7.9 data (500 individuals, four
indicators, shared variances), providing external references for the
Lo-Mendell-Rubin adjusted likelihood-ratio statistic.

TECH11 is available only under `ESTIMATOR = MLR`; Mplus reports
`TECH11 option is available only for estimator MLR` and ignores the request
otherwise. The robust estimator therefore had to exist before this reference
could be generated.

| File | Contents |
|---|---|
| `lmr2.inp` / `lmr2.out` | One versus two classes; TECH11 block |
| `lmr3.inp` / `lmr3.out` | Two versus three classes; TECH11 block |
| `ex7.9.dat` | Input data, copied from `../public/` |
| `compare.R` | Reproduces the comparison from the project root |
| `comparison.csv` | Retained numerical result |

| Comparison | Mplus statistic | Mplus adjusted | Native adjusted | Difference |
|---|---:|---:|---:|---:|
| 1 versus 2 classes | 747.063 | 723.771 | 723.7704 | 5.53e-4 |
| 2 versus 3 classes | 8.065 | 7.814 | 7.8135 | 4.57e-4 |

The native fits reproduce the statistic end to end: fitting one- and two-profile
models in R on the same data gives `747.0633` against the `747.063` Mplus prints
to three decimals. The adjustment is `1 + 1 / (df log n)`, confirmed on both
runs rather than fitted to one.

## What is deliberately not reproduced

Mplus prints a Mean, a Standard Deviation and a p-value for the
Vuong-Lo-Mendell-Rubin reference distribution. None is implemented, and no
substitute is returned. A chi-square tail probability is not a usable stand-in:
on the second run it gives 0.153 where Mplus reports 0.582. Moment-matched
gamma and scaled chi-square approximations were also tested against that value
and did not reproduce it either. Use `bootstrap_lrt_multilpa()` for a calibrated
class-count p-value.

```sh
Rscript validation/mplus/lmr/compare.R
```
