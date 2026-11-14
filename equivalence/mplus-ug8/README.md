# Growth and latent transition models against the Mplus 8.8 User's Guide, chapter 8

The Mplus User's Guide posts each example's input, data and Mplus's own output
(statmodel.com/usersguide/chapter8.shtml, `chapter8.zip`, Mplus 8.8, run
2022-04-20). Four are models latents fits, kept unchanged in `listings/`
(`LISTINGS.md5`); every compared value is read from the `.out` file.

| Case | Example | latents |
|---|---|---|
| `ex8.9_lcga` | 8.9: latent class growth analysis, binary outcome, 4 occasions, 2 classes | `mixture_regression(u ~ time, family = "binomial", id = "id", class_level = "group")` |
| `ex8.12_hmm` | 8.12: hidden Markov model, one binary indicator, 4 occasions | `lta(transitions = "homogeneous")` |
| `ex8.14_covariate` | 8.14: two occasions, 5 binary indicators, 3 states, a covariate on the first state and on each origin's transitions | `lta(initial_covariates = "x", transition_covariates = "x")` |
| `ex8.15_mover_stayer` | 8.15: mover-stayer, three occasions, initial states differing between movers and stayers | `lta(mover_stayer = TRUE, transitions = "occasion")` |

Compared: log likelihood, parameters, AIC, BIC; ex8.9's entropy, class
proportions, modal counts and every estimate with its MLR standard error
(Mplus's threshold, I and S means and C#1 are linear in latents' class
intercepts, slopes and membership logit); each state's item probabilities,
with their errors in ex8.12; the transition tables of ex8.12 and ex8.14.

```
Rscript equivalence/mplus-ug8/compare.R      # comparison.csv, RUN.txt
```

`equivalence/test-mplus-ug8.R` runs the same comparison under testthat.
Tolerances: half a unit of the last printed digit plus 1e-6 of the value;
counts exactly. Standard errors are latents' robust errors times
sqrt((n - 1) / n) (Mplus's MLR has no CR1 factor).

## Known differences (`known-differences.csv`)

- **ex8.9's C#1** sits on a rounding boundary: Mplus's printed proportions
  put it at -0.15350 +- 2e-5; latents -0.153498, Mplus prints -0.154.

## Not compared

Example 8.10 (latents' ordinal mixture regression gives each class its own
thresholds: 7 parameters against Mplus's 6), 8.13 (Mplus models the known
class's proportion too), the growth mixtures with random effects (8.1 to
8.8: Mplus's residual variances differ by occasion), the zero-inflated
models (8.5, 8.11), Mplus's entropy of several latent class variables,
ex8.15's transition tables and class counts (Mplus mixes movers and
stayers), and the chi-square tests.
