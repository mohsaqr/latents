# The Mplus User's Guide examples latents fits

The Mplus 8.8 User's Guide posts each example's input, data and Mplus's own
output (statmodel.com/usersguide, `chapter7.zip`, `chapter8.zip`; run
2022-04-20). The examples latents fits are kept unchanged in `listings/`
(`LISTINGS.md5`), and every compared value is read from the `.out` file.
**[CATALOGUE.md](CATALOGUE.md) maps every example of chapters 7, 8 and 10
to a latents function or to the feature latents would need: check it before
adding a mixture model.**

| Case | Example | latents |
|---|---|---|
| `ex7.1_regression` | mixture regression: common slope x1, class slope x2, class variances, x1 on membership | `mixture_regression(y ~ x1 + x2, common = ~ x1, membership = ~ x1)` |
| `ex7.3_lca`, `ex7.4_lca` | LCA, four binary indicators | `lca()` |
| `ex7.6_lca` | LCA, three-category indicators | `lca()` |
| `ex7.7_lca` | LCA, nominal indicators (Mplus prints logits against the last category) | `lca()` |
| `ex7.9_lpa`, `ex7.10_lpa`, `ex7.22_lpa` | LPA: EEI, VVI, EEE with three classes | `lpa(model =)` |
| `ex8.9_lcga` | LCGA, binary outcome, four occasions | `mixture_regression(family = "binomial", class_level = "group")` |
| `ex8.12_hmm` | hidden Markov model | `lta(transitions = "homogeneous")` |
| `ex8.14_covariate` | LTA, covariate on the states and each origin's transitions | `lta(initial_covariates =, transition_covariates =)` |
| `ex8.15_mover_stayer` | mover-stayer LTA | `lta(mover_stayer = TRUE, transitions = "occasion")` |

Compared: log likelihood, parameters, AIC, BIC; entropy, class proportions
and modal counts (chapter 7 and ex8.9); every estimate Mplus prints in MODEL
RESULTS, or RESULTS IN PROBABILITY SCALE for categorical indicators, with
its MLR standard error (delta method from latents' robust covariance); each
state's item probabilities and the transition tables of the latent
transition examples (see `compare-functions.R`).

```
Rscript equivalence/mplus-ug/compare.R      # comparison.csv, RUN.txt
```

`equivalence/test-mplus-ug.R` runs the same comparison under testthat.
Tolerances: half a unit of the last printed digit plus 1e-6 of the value;
counts exactly. Standard errors are latents' robust errors times
sqrt((n - 1) / n) (Mplus's MLR has no CR1 factor). Class sizes are compared
as the proportions Mplus prints.

## Known differences (`known-differences.csv`)

- **Rounding boundaries** (ex8.9's C#1; one ex8.14 transition on Linux).
- **Class proportions of ex7.4 and ex7.6** at the fifth decimal (Mplus's
  convergence from user-specified starts).
- **ex7.7's near-flat class**: one barely identified error (estimate 2.35,
  error 3.92; relative gap 2.5e-4), and on Linux only its class proportions,
  one logit at a rounding boundary and one error (the converged point moves
  by about 1e-5 between platforms).

## Not compared

The other examples of chapters 7, 8 and 10: [CATALOGUE.md](CATALOGUE.md)
gives each one's reason. Also Mplus's entropy of several latent class
variables, ex8.15's transition tables and class counts (Mplus mixes movers
and stayers), and the chi-square tests.
