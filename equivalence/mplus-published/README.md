# lca(), lpa(), lta() and VLMR against Mplus output printed by UCLA OARC

Genuine Mplus output printed on the UCLA OARC Statistical Methods and Data
Analytics pages, with the data files posted there (`data/`, `DATA.md5`,
`PROVENANCE.txt`), is reproduced:

| Case | Page | Compared |
|---|---|---|
| `mccutcheon_ego_k1`, `_k2` | Applied Latent Class Analysis (Hagenaars and McCutcheon 2002), ch. 2, McCutcheon: Tables 2 to 4, Ego's Dilemma (Stouffer and Toby; group 1 of `page59_a.dat`) | log likelihood, parameters, BIC, entropy, class counts, every probability with its error, the thresholds and class logit with errors |
| `mccutcheon_abortion_k2` | ch. 2, Table 8 model H1 (`page75.dat`) | log likelihood, parameters, BIC, entropy |
| `dayton_pleural_k2`, `dayton_cheating_k2`, `dayton_female_k2`, `dayton_male_k2` | Latent Class Scaling Analysis (Dayton 1998): Table 3.2 model I, the cheating example (p. 30), Table 6.2 models II and III | log likelihood, parameters, BIC, entropy, class counts, every probability with its error |
| `ucla_drinking_k3`, `ucla_drinking_vlmr_k2_k3` | Mplus data analysis example "Latent Class Analysis" (`lca1.dat`, 1000 persons, nine items) | each class's probability of "yes", modal and posterior class counts; TECH11's H0 log likelihood, LR, parameter difference and adjusted LRT |
| `vermunt_diabetes_*` | ALCA ch. 3, Vermunt and Magidson, Table 1 (`EXAM1.dat`) | log likelihood, parameters, BIC, entropy of ten models: VVV, EEE, VVI, EEI at one to three clusters |
| `langeheine_lm_*` | ALCA ch. 11, Langeheine and van de Pol, Table 3 models 7a and 7b (`chap11.dat`) | log likelihood, parameters, BIC of the latent Markov model with homogeneous and occasion-specific transitions |

```
Rscript equivalence/mplus-published/compare.R      # comparison.csv, RUN.txt
```

`equivalence/test-mplus-published.R` runs the same comparison under testthat.

Tolerances, fixed before the comparison: half a unit of the last printed
digit plus 1e-6 of the value; counts of parameters and modal counts exactly.

## Conventions

- **Standard errors.** Mplus fits mixtures by MLR: a sandwich without a
  small-sample factor. latents' robust covariance carries the CR1 factor
  n / (n - 1), so its errors are scaled by sqrt((n - 1) / n). Dayton's
  cheating page printed ML errors (observed information): its input names no
  estimator, but all sixteen agree with ML and four with MLR.
- **Frequency tables** are expanded to one row per person (Mplus's
  `weight is freq (freq)`).
- **Ego's Dilemma.** `page59_a.dat` holds two groups of 216 (for Table 10's
  two-group models); Tables 2 to 4 are of group 1, whose log likelihood is
  Mplus's -504.468.
- **Langeheine and van de Pol.** Mplus reads the file without its group
  column, so the two groups' counts are pooled (5147 persons).

## Known differences (`known-differences.csv`)

- **Diabetes VVV, two clusters.** Mplus stopped at a local maximum
  (-2355.905); 10 of latents' 50 starts reach -2354.647, no boundary fit.
- **Class counts to five decimals** (female and male cheating, drinking): not
  converged digits of Mplus; the gaps are 1.2e-4 to 0.017 persons.

## Not compared

TECH11's p-values, mean and SD (Mplus deletes the first class of the K-class
model for its H0 and builds the weights differently; `test-features-vlmr.R`),
the chi-square fit tests, Dayton's constrained scaling models and McCutcheon's
restricted models (latents fixes no parameter), the entropy of Mplus's
latent Markov models (over the joint pattern of five latent variables),
Vermunt and Magidson's models with one free covariance and their prostate
models.
