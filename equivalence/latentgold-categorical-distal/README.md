# Categorical distal outcomes against Latent GOLD 6.1

`three_step(x, data, outcome)` with a factor, character or logical outcome
reports the probability of each category in each class. Latent GOLD 6.1's
step-three model of a nominal dependent variable (`dependent distal nominal`,
`distal <- 1 + Cluster`) is the reference, for the four methods:

| `three_step(method =)` | Latent GOLD |
|---|---|
| `"bch"` | `step3 modal bch` |
| `"modal"` | `step3 modal none` |
| `"proportional"` | `step3 proportional none` |
| `"ml"` | `step3 modal ml` |

Latent GOLD reads latents' step-one posteriors from the kit
(`posterior = (post1 ...)`), assigns modally and builds the classification
errors itself. It was run in batch under Wine on this machine
(`~/.wine-latentgold`, `lg61.exe <case>.lgs /b /o <case>.lst`) on
31 October 2026 (the four ML files again on 1 November, without the size
equation); its listings are in `returned/`, so the comparison runs
offline. `returned/INPUTS.md5` fingerprints the data and syntax it was given,
and the test refuses a kit that has changed since.

```sh
Rscript equivalence/latentgold-categorical-distal/make-kit.R
sh equivalence/latentgold-categorical-distal/run-kit.sh
Rscript equivalence/latentgold-categorical-distal/compare.R   # -> comparison.csv
Rscript equivalence/run.R latentgold-categorical-distal       # the same, as a test
```

| case | first step | units | outcome |
|---|---|---|---|
| d01 | `lpa()`, 2 profiles, 2 Gaussian indicators | 600 | 3 categories |
| d02 | `lca()`, 3 classes, 5 binary items | 900 | 2 categories |

**Result: 144 quantities compared, 144 agree** (each within Latent GOLD's
printed rounding, four decimals): every class probability of every category
for all four methods, with model-based and with robust standard errors
requested (96); the robust standard errors of BCH, modal and ML (36); and the
robust Wald test of the class effect with its degrees of freedom (12).
latents' side uses `vcov_type = "independent"`, the variance Latent GOLD's
robust standard errors use; both carry the factor n / (n - 1).

## Not comparable (70 rows)

- **Model-based standard errors and Wald tests.** Latent GOLD's come from the
  model's information matrix; latents' are sandwich estimates, for every
  method.
- **Proportional standard errors and Wald test.** Latent GOLD treats the
  posterior weights as frequency weights: its robust standard error is
  `sqrt(p (1 - p) / T_k * n / (n - 1))`, with `T_k` the class's total
  posterior weight. That formula reproduces every printed value. latents'
  is the sandwich of the weighted mean over independent units, as for the
  other methods and for a numeric outcome.

## What the comparison established

- Latent GOLD's step-three model estimates the class sizes when the syntax
  has `Cluster <- 1;` and holds them at the average posterior probabilities
  when it does not. latents' ML model holds them, for a categorical outcome
  as for a numeric one (whose kit, `features/ml-three-step/`, omits the
  equation too, and agrees with Mplus's DU3STEP, which also holds them). The
  ML files here omit it; Latent GOLD then reports K (C - 1) parameters. With
  the sizes free the model would be saturated and reproduce BCH; an earlier
  run of this kit with them free showed the ML standard errors of the two
  conventions differ by up to 4%.
- The overall test is on the baseline-category logits, as Latent GOLD's
  Wald test of `distal <- Cluster` is. A test on the probability scale gives
  a different statistic (Wald tests depend on the parameterization).
