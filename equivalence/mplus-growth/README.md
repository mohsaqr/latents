# Growth mixtures against genuine Mplus 9 runs (JStats fixtures)

Ten Mplus 9 Demo runs kept in JStats (`tests/fixtures/mplus/`), copied to
`returned/` with their inputs and data (`PROVENANCE.txt`), are reproduced by
`mixture_regression(y ~ time, data, k, id = "id", class_level = "group", ...)`:

| series | model | estimator | classes |
|---|---|---|---|
| `lcga-ml` | latent class growth, five waves (growth-factor variances fixed at 0) | ML | 1, 2, 3 |
| `runner-lcga` | latent class growth, four waves | MLR (Mplus's default) | 1, 2, 3 |
| `gmm-ml` | growth mixture, class-specific random intercept and slope covariance (`random = ~time`, `random_covariance = "varying"`) | ML | 2 |
| `runner-gmm-growth` | the same, four waves | MLR | 1, 2, 3 |

Every class has its own residual variance, equal across waves (Mplus's labels
`(r1)`, `(r2)`, ...), which is latents' default. Compared: the log likelihood,
BIC, parameter count and every printed estimate and standard error (class
logits against Mplus's last class; growth-factor means, variances and
covariance; the residual variance).

```
Rscript equivalence/mplus-growth/compare.R      # comparison.csv, RUN.txt
```

`equivalence/test-mplus-growth.R` runs the same comparison under testthat.

**Result (2026-10-29): 212 quantities compared, 134 agree (every printed value
of the seven one- and two-class runs), 78 known differences (the three
three-class runs, below), 0 unexplained.**

Tolerances, fixed before the comparison: Mplus prints three decimals, so half a
unit (5e-4) plus 1e-6 of the value; counts exactly.

**MLR.** The runner series uses Mplus's default MLR: the person-level sandwich
without a small-sample factor. latents' robust covariance (`vcov_type =
"robust"`) carries the CR1 factor G / (G - 1); with it removed every MLR error
agrees, as `test-mplus-mixture-weights.R` documents for Mplus WEIGHT.

## Known differences (`known-differences.csv`)

The three-class runs over-extract these data, and the two programs report
different maxima:

- **lcga-ml, three classes.** Mplus stops at -2087.250 with a 1% class;
  latents reaches -2087.181 (0.069 higher, from 50 and from 300 starts) with a
  1.6% class.
- **runner-lcga, three classes.** Mplus's -2053.234 has a third class of about
  one person with residual variance 0.027 and standard errors of exactly zero:
  a degenerate spike with singular information.
- **runner-gmm-growth, three classes.** Mplus's -2715.247 has a negative
  intercept variance (-0.539). latents parameterizes every random-effect
  covariance by its Cholesky factor, so it cannot reach that point.

## Also seen

On runner-lcga with three classes latents' best maximum from 300 starts
(-2056.463) is below its best from 50 (-2055.762): starts are screened for 50
iterations and only the better half finished, so a start that would have won
can be screened out. Recorded, not changed (the start search is frozen).
