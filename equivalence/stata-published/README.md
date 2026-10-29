# lpa() and lca() against Stata's printed output and Masyn (2013)

Genuine Stata 19 output printed in its manual, and the published model-fit
table of Masyn (2013), are reproduced:

- **Stata [SEM] Example 52g** (latent profile models, `gsem ..., lclass()`),
  Reaven and Miller's diabetes data divided by 10, n = 145: the two-class model
  with class-invariant variances (EEI) and the three-class model with
  class-varying covariances (VVV): the log likelihood and every printed
  estimate and standard error (class logits, means, variances, covariances).
- **Stata [SEM] Example 50g** (a latent class model), Stouffer and Toby's four
  role-conflict items, n = 216, two classes: the log likelihood, the class
  logit, every item logit with its error, and the class probabilities with
  their error (`estat lcprob`).
- **Masyn (2013)**, *Latent class analysis and finite mixture modeling*, The
  Oxford Handbook of Quantitative Methods vol. 2, Table 25.12, the same data:
  log likelihood, parameter count and BIC for EEI and VVI at one to five
  profiles, EEE at one to four and VVV at one to three (Mplus).

```
Rscript equivalence/stata-published/compare.R      # comparison.csv, RUN.txt
```

`equivalence/test-stata-published.R` runs the same comparison under testthat.
The data are kept as CSV (`PROVENANCE.txt`): Stata's stored single-precision
values, which its fits used.

**Result (2026-10-29): 154 quantities compared, 145 agree, 9 known differences
(below), 0 unexplained.**

Tolerances, fixed before the comparison: half a unit of the last printed digit
plus 1e-6 of the value.

## Known differences (`known-differences.csv`)

- **Stata's convergence.** Stata stops by default when its scaled gradient
  g'H^-1 g falls below 1e-5, so the last of its seven printed digits are not all
  converged. Six values of the three-profile VVV model's second class and one
  error of the latent class model are outside the printed rounding, all by less
  than 1e-5 standard errors (largest 8.1e-6); the log likelihoods agree to 2e-6.
- **VVI at five profiles.** The likelihood is unbounded there: latents' search
  reaches a variance-floor spike (one member, flagged bound-active), 7.8 above
  Masyn's proper maximum. The parameter count agrees.

## Not compared

Masyn's adjusted LMR p-values come from Mplus TECH11, which fits the K-class
null model by deleting the first class of the K + 1 solution. latents' VLMR
(like Latent GOLD's, which `test-features-vlmr.R` compares) is Vuong's test
between the two separately fitted maxima, so the two p-values answer different
questions (EEI 4 against 5: 0.29 in Mplus, 0.007 in latents).
