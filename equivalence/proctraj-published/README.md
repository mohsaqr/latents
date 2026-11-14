# Latent class growth models against SAS PROC TRAJ (Andruff et al. 2009)

Andruff, H., Carraro, N., Thompson, A., Gaudreau, P., and Louvet, B. (2009).
Latent class growth modelling: a tutorial. *Tutorials in Quantitative Methods
for Psychology*, 5(1), 11-24. The tutorial prints its data (Appendix A:
disengagement coping of 107 persons at times 0, 1, 2; `disen.csv`) and SAS
PROC TRAJ's output (Jones, Nagin and Roeder 2001) for a censored normal model
with `MAX 5`. No score lies at a bound (1.00 to 4.63), so the likelihood is
the normal one: `mixture_regression(disengagement ~ time + I(time^2), id =
"id", class_level = "group", variance = "equal")`.

| Case | Compared |
|---|---|
| `proctraj_k1` | Figure 3, one quadratic trajectory: L, AIC, both BICs, every estimate and error, Sigma |
| `proctraj_k1_linear`, `proctraj_k2` | the BICs the text reports for one linear and two quadratic trajectories |
| `proctraj_k3` | Figure 4, three quadratic trajectories: L, AIC, both BICs, every estimate and error, Sigma, the group percentages and their errors |

```
Rscript equivalence/proctraj-published/compare.R      # comparison.csv, RUN.txt
```

`equivalence/test-proctraj-published.R` runs the same comparison under
testthat. Tolerances: half a unit of the last printed digit plus 1e-6 of the
value. PROC TRAJ's criteria are on the log-likelihood scale: AIC = L - p,
BIC = L - p log(N) / 2, N the 321 observations or the 107 persons.

**Standard errors.** Every error both figures print is latents' (observed
information) times sqrt(N / (N - q)), N = 321 and q the parameters other than
Sigma: 1.004706 (q = 3) and 1.017588 (q = 11), to the printed digits, a
degrees-of-freedom factor on the residual variance. It was measured on these
two outputs, not taken from PROC TRAJ's documentation.

## Known differences (`known-differences.csv`)

- **The text's linear BIC** (-307.79) contradicts the paper's Figure 3 (the
  quadratic term has p = .9159, so the linear BIC is -307.46).
- **Three-trajectory estimates past the fifth decimal**: gaps of 5.1e-6 to
  2.8e-5 (shares 2.4e-4 points), at most 6.2e-5 standard errors, mostly in the
  smallest group (3.4 persons); L, the criteria and every error agree.

## Not compared

Figure 5 and Table 2 (quadratic, linear, quadratic trajectories: one formula
serves every class in mixture_regression()), the four-trajectory BIC (its
orders are not given) and the posterior summaries of the final model.
