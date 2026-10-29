# lta() against LMest, seqHMM and two published LMest analyses

LMest's `RLMSdat` (job satisfaction in the Russia Longitudinal Monitoring
Survey: one five-category item at seven occasions, 1,718 people) with a basic
latent Markov model: time-homogeneous transitions and the same response
probabilities at every occasion. `lta(data, "y", "id", k, time = "time",
categorical = "y")`, 50 starts, is compared with

- the current **LMest** (`lmest(modBasic = 1)`, 20 starts, `tol = 1e-12`) at
  two, three and four states;
- **seqHMM** (`build_hmm()` and `fit_model()` with 20 EM restarts) at two and
  three states;
- the published analysis of **Bartolucci, Pandolfi and Pennoni (2017)**, *LMest:
  An R package for latent Markov models for longitudinal categorical data*,
  JSS 81(4), section 3.1: the three-state log likelihood, parameter count, BIC
  and initial probabilities, and the four-state log likelihood and BIC.

And LMest's `data_market_sim` (200 customers, five occasions; `brand`, four
categories, and `price`, five ranges, as two categorical items) with the same
model: `lta(data, c("brand", "price"), "id", k, time = "time", categorical =
c("brand", "price"))` at two, three and four states against

- the current **LMest** (`lmest(responsesFormula = brand + price ~ NULL,
  modBasic = 1)`, 20 starts, `tol = 1e-12`);
- the hidden Markov models of **Pennoni, Pandolfi and Bartolucci (2024)**,
  *LMest: An R package for estimating generalized latent Markov models*, The R
  Journal, RJ-2024-036, section 3.3 (`lmestSearch()`, seed 12345): log
  likelihood, parameter count, AIC and BIC. The paper's one-state model is not
  a transition model (`lta()` refuses it). Its covariate model holds the
  response probabilities fixed at the three-state estimates and uses LMest's
  `difflogit` transition slopes, a parametrization `lta()` does not have, so
  it is not compared.

Compared: the maximized log likelihood and every initial, transition and
response probability (after matching state labels).

```
Rscript equivalence/lta-packages/compare.R      # comparison.csv, RUN.txt
```

`equivalence/test-lta-packages.R` runs the same comparison under testthat.

**Result (2026-10-29): 52 quantities compared, 48 agree, 4 known differences
(the 2017 paper's unconverged three-state values, below), 0 unexplained.** On
the market data all 24 agree: the log likelihoods within 6e-8 of LMest's and
every probability within 4e-5; the 12 published values to their printed
precision.

Tolerances, fixed before the comparison: software, the log likelihood to 1e-6
plus 1e-10 of its size and every probability to 1e-4; published values, half a
unit of the last printed digit.

## Known differences (`known-differences.csv`)

- **The paper's three-state values are LMest's default stopping point.** With
  `start = 0` and its default `tol = 1e-8` the current LMest reproduces them
  exactly (log likelihood -13557.20657, initial probabilities .1739 / .3638 /
  .4623, BIC 27263.391). Run to convergence (`tol` 1e-10 or tighter) LMest
  reaches -13557.20497 with .1737 / .3639 / .4624, the maximum `lta()` reaches.
  The published BIC and the four-state values agree.
