# Genuine Mplus comparison, two-level latent class analysis

`CATEGORICAL = u1-u5` with `TYPE = TWOLEVEL MIXTURE`, on 1,200 simulated
individuals in 60 groups with five binary indicators.

| File | Contents |
|---|---|
| `categorical.inp` / `categorical.out` | Mplus specification and VERSION 9 DEMO (Mac) output |
| `categorical-results.dat` | Saved estimates, standard errors and fit statistics |
| `categorical-posteriors.dat` | Saved class probabilities |
| `categorical.dat` | Input data |
| `compare.R` | Reproduces the comparison from the project root |
| `comparison.csv` | Retained numerical result |

| Quantity | Largest absolute difference |
|---|---:|
| Thresholds (10) | 5.75e-6 |
| Profile probabilities within group class | 1.06e-6 |
| Group-class probabilities | 8.30e-7 |
| Log likelihood | 2.75e-6 |

## The `MODEL cw:` block is load-bearing

Without it Mplus estimates **23** free parameters, not 13: its default frees the
thresholds across *both* latent class variables, giving `2 x 2 x 5 = 20`
thresholds. That is a different model from this package's, which assumes
measurement invariance across the group classes.

The `MODEL cw:` block fixes the thresholds to vary with `cw` only, reducing the
count to `10 + 1 + 2 = 13` and matching the native parameterization. TECH1
confirms the layout: parameters 1-5 are the `cw#1` thresholds and 6-10 the
`cw#2` thresholds, shared across the `cb` patterns.

Anyone extending these comparisons should check the Mplus parameter count
against the native one before comparing estimates; a silent mismatch here would
compare two different models.

## Parameterization

Mplus reports `u$1` thresholds. For a binary indicator the threshold is
`qlogis(P(u = lowest category))`, which is exactly the first column of
`.ml_lpa_categorical_thresholds()`, so the two are directly comparable once
class labels are aligned.

```sh
Rscript validation/mplus/categorical/compare.R
```
