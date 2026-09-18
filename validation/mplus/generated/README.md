# Genuine Mplus Demo multilevel LPA comparisons

These are new Mplus VERSION 9 DEMO (Mac) runs, not published Mplus User's Guide results. Both use the existing project synthetic demonstration data: 1,200 individuals in 60 unbalanced groups (16, 20, or 24 people), with two continuous indicators, two individual profiles, and two discrete group classes. Original generation code is `validation/synthetic-demo.R`; `synthetic.dat` preserves the data supplied identically to both estimators, adds numeric group and subject IDs, and omits latent truth labels.

The command used in this directory was:

```sh
/Applications/MplusDemo/mpdemo varying.inp varying.out > varying-console.txt 2>&1
/Applications/MplusDemo/mpdemo equal.inp equal.out > equal-console.txt 2>&1
```

Both Mplus runs use independent 100 initial starts and 20 final optimizations, ML estimation, seed 20260917, and tight convergence settings. Both terminate normally with replicated best likelihood. `compare-generated.R` independently fits the native R model with 30 starts and seed 20260917. Run from the project root:

```sh
Rscript validation/mplus/generated/compare-generated.R
Rscript -e 'testthat::test_file("validation/mplus/generated/test-generated.R")'
```

Model equivalence is checked from the saved specification and TECH1: all indicators are within-only; means are profile-specific but invariant across group classes; covariance matrices are diagonal; group classes affect profile probabilities through `cw#1 ON cb`; no continuous random effects or covariates. The varying-variance model has 11 free parameters; the equal-variance model has 9. Mplus's warnings that indicators are uncorrelated within class describe the intended diagonal covariance model.

`*-results.dat` are genuine Mplus SAVEDATA RESULTS files: parameter estimates (TECH1 order), standard errors, parameter count, likelihood, AIC, BIC, adjusted BIC, entropy, and condition number. They use eight significant digits; the printed `.out` uses three decimals. `*-posteriors.dat` contains Mplus joint posteriors at 12 decimal places. The R comparison sums joint probabilities for each marginal level. It aligns subjects by saved ID (Mplus reorders rows), profiles by the first indicator mean, and group classes by the conditional probability of profile 1. The Mplus BIC matches the native R `bic_individual`, not its default group-count `bic`.

Both comparisons pass all script assertions and the accompanying testthat tests. Across the two scenarios, largest mean/variance discrepancy is 5.05e-8, mixing-probability discrepancy 5.12e-9, subject-posterior discrepancy 1.69e-10, and group-posterior discrepancy 5.10e-12. Likelihood differences are at most 4.45e-5, within saved-result rounding precision. Exact individual metrics are in `comparison.csv`.

These two well-separated synthetic scenarios validate this specific model calculation. They do not establish overall Mplus feature parity, difficult-case global optimization reliability, or equivalence of standard errors (the R package does not provide these).

The first experimental derivative tolerance of 1e-8 stalled at floating-point precision. That run was interrupted and replaced by final runs using 1e-6 derivative tolerance; log-likelihood criteria remain 1e-8 absolute and 1e-9 relative. Only final successful runs are retained. Mplus Demo limits printed in these outputs: 6 dependent variables, 2 independent variables, 2 between variables, 2 continuous latent variables for time-series analysis. These models fall within those limits.
