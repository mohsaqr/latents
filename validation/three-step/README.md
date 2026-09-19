# Three-step correction: independent checks

The BCH correction is validated three ways, because no seminal source publishes
a worked example that can be reproduced from public data.

## 1. Against an independent implementation

`compare-tidysem.R` hands this package's posterior probabilities to
`tidySEM:::classification_probs_mostlikely()` and rebuilds the weights and the
corrected means from them. tidySEM implements Bolck, Croon and Hagenaars (2004)
for OpenMx mixtures, and its weight construction is a pure function of the
posteriors, so the comparison is exact rather than approximate.

| quantity | max absolute difference |
|---|---|
| classification error matrix | 8.88e-16 |
| BCH weights | 1.11e-15 |
| corrected class means | 1.78e-15 |

## 2. Against the method's defining identities

Verified in `tests/testthat/test-three-step.R`, so they run on every check:

- the error matrix is row-stochastic;
- the weight matrix is exactly its inverse;
- the weights sum to one within a unit;
- **the weighted class sizes reproduce the model's own estimated class sizes**,
  which is the property that makes the third step unbiased;
- **the table formulation agrees with the per-unit weights**. Bolck et al.
  correct the assigned-by-external contingency table with the inverse error
  matrix; this package carries per-unit weights instead. The two routes are
  algebraically the same method and agree to 9e-13.

## 3. Against the paper's headline claim

Bolck, Croon and Hagenaars showed that naive three-step *underestimates* the
association between class membership and an external variable. Over 60
replications with 8% misclassification and a true class difference of 10:

| method | bias in the difference | 95% interval coverage |
|---|---|---|
| `bch` | -0.126 | 0.85 |
| `modal` | -1.285 | 0.00 |
| `proportional` | -1.821 | 0.00 |

The attenuation is reproduced, and the correction removes about nine tenths of
it. The 0.85 coverage is the known cost of treating the error matrix as known
rather than estimated, and is documented on `three_step()`.

## R3STEP: covariates predicting membership

`r3step()` implements the ML three-step correction of Vermunt (2010): the
assigned class is treated as an error-prone indicator of the true one, with the
error rates fixed at what step one found, and the multinomial logit is fitted
against that. Over 60 replications with a true log-odds slope of 1.2 and
intercept -0.3:

| parameter | naive bias | R3STEP bias |
|---|---|---|
| intercept | +0.012 | **+0.006** |
| slope | -0.163 | **-0.002** |

The attenuation the seminal paper describes is reproduced and removed. The
correction costs precision, as expected: the slope's standard deviation across
replications was 0.127 against 0.100 for the naive fit. The nominal 95%
intervals covered the true slope in 95% of the replications, measured where the
profiles separate well.

There is also an internal cross-check available that no external package
provides: `fit_covariates()` fits the same covariate model in one step. Three-
step should approximate it, and any large disagreement is a signal about one of
the two.

## What could not be done

No equivalence against a published table. Bolck et al. (2004), Vermunt (2010)
and Bakk & Vermunt (2016) are simulation studies rather than worked examples,
and the Latent GOLD Step-3 tutorial is image-based with no obtainable dataset.
An Mplus `R3STEP`/`BCH` run would be the natural addition whenever Mplus is at
hand; the existing `validation/mplus/` artifacts show the pattern to follow.

## References

Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent structure models
with categorical variables: One-step versus three-step estimators.
*Political Analysis*, 12, 3--27.

Vermunt, J. K. (2010). Latent class modeling with covariates: Two improved
three-step approaches. *Political Analysis*, 18, 450--469.

Bakk, Z., & Vermunt, J. K. (2016). Robustness of stepwise latent class modeling
with continuous distal outcomes. *Structural Equation Modeling*, 23, 20--31.
