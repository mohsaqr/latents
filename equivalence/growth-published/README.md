# The paquid growth mixtures of Proust-Lima, Philipps and Liquet (2017)

`mixture_regression()` reproduces section 6 of Proust-Lima, Philipps and Liquet
(2017), *Estimation of extended mixed models using latent classes and latent
processes: the R package lcmm*, JSS 78(2), doi:10.18637/jss.v078.i02: normMMSE
on a quadratic in age65 with class-specific trajectories, CEP shared by the
classes, a correlated random intercept, slope and quadratic with one
covariance for every class, and one residual variance. The same three models
are fitted by the current lcmm as well.

```
Rscript equivalence/growth-published/compare.R      # comparison.csv, RUN.txt
```

`equivalence/test-growth-published.R` runs the same comparison under testthat.

**Result (2026-10-29): 37 quantities compared, 33 agree (all 19 published
values; lcmm at one and two classes on both datasets), 4 known differences
(lcmm's three-class grid search on both, below), 0 unexplained.**

In latents the model is

```
mixture_regression(normMMSE ~ age65 + age65sq + CEP, data, k, id = "ID",
                   class_level = "group", common = ~CEP,
                   random = ~ age65 + age65sq, random_covariance = "equal",
                   variance = "equal")
```

with `age65sq = age65^2` (the paper's `poly(age65, 2, raw = TRUE)`), 50 starts.

**Published values reproduced** (Table in section 6 and `postprob(m2)`): the
log likelihood, number of parameters and BIC at one, two and three classes
(-8920.623 / 11 / 17909.61; -8899.228 / 15 / 17891.67; -8890.148 / 19 /
17898.37, the grid search's best), the modal class percentages (12.4 / 87.6;
0.8 / 11.8 / 87.4), the 62 and 438 subjects of the two-class model and their
mean posteriors 0.8054 and 0.8730. All agree to their printed precision.

**The current lcmm** (2.2.1) fits the same models as the paper does (`B = m1`
for two classes; `gridsearch(rep = 30, maxiter = 15)` for three). Its maximized
likelihood, parameter count and every fixed effect are compared; its
likelihood to 1e-4 (its stopping criteria include 1e-4 of the likelihood) and
the fixed effects to 1e-3 of max(1, |value|).

**lcmm's `data_hlme`** (simulated, 100 subjects, 326 rows) adds a second
dataset with the model of `hlme()`'s help page: `Y ~ Time * X1` with the
intercept and Time class-specific, X1 and Time:X1 shared, a correlated random
intercept and slope with one covariance for every class, one residual
variance, and class membership on X2 + X3:

```
mixture_regression(Y ~ Time + X1 + Time:X1, data_hlme, k, id = "ID",
                   class_level = "group", common = ~ X1 + Time:X1,
                   membership = ~ X2 + X3, random = ~Time,
                   random_covariance = "equal", variance = "equal")
```

against `hlme(..., classmb = ~ X2 + X3)` by `gridsearch(rep = 30, maxiter =
15)` from the one-class fit. Compared at one to three classes: the likelihood,
parameter count and every estimate (class intercepts and slopes, shared
effects, membership logits against a common class, the random covariance and
the residual SD) to 1e-3 of max(1, |value|). At one and two classes the
likelihoods agree within 4.5e-6 and the estimates within 1.6e-5.

## Data

`paquid-normmmse.csv` is `lcmm::paquid` with `normMMSE = NormPsy::normMMSE(MMSE)`
and `age65 = (age - 65) / 10`, computed once with NormPsy 1.0.8 (installed in a
scratch library, not as a dependency): see `PROVENANCE.txt`. 2,214 of 2,250
rows have an MMSE.

## Known differences (`known-differences.csv`)

- **Three classes, lcmm, paquid.** lcmm's grid search (30 random vectors of
  15 iterations, seed 1) stops at the local maximum -8891.351, the paper's own
  default-start m3; `mixture_regression()` reaches -8890.148, the paper's best
  (m3b), against which all its published three-class values agree.
- **Three classes, lcmm, data_hlme.** lcmm's grid search reaches -769.6274 at
  a quasi-separated membership: one class's X2 logit is -21.8, so the class is
  all but excluded where X2 = 1 (classes 15/31/54). That is the supremum
  approached as the logit runs to minus infinity. `mixture_regression()`
  reaches the interior maximum -769.6292 (16/39/45) from 50 starts with two
  seeds and from 200 starts, 0.0017 lower.
