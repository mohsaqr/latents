# Membership covariates: genuine Mplus comparison

The retained Mplus 9 Demo input/output fits the same one-step ML model as
`fit_ml_lpa_covariates()`: two individual profiles, two group classes, two
conditionally independent Gaussian indicators, individual predictor `z` with
shared profile-logit slope, and group predictor `w` for group-class membership.
There are 720 individuals in 60 groups and 13 free parameters. Mplus used
100 initial/20 final starts, terminated normally and replicated its best
likelihood. The comparison R fit uses ten independent starts; neither engine
is initialized from the other's estimates. All constraints are inactive.

`covariates.dat` is the exact synthetic input; `r-fit.rds` preserves the initial
R probe and its data. Generation used seed 223: independent standard-normal
individual `z` and group `w`; group classes from logistic(-.2 + .7*w), individual
profiles from logistic(-1.5 + 3*group_class + .9*z), and indicators with means
(-3,-2)/(3,2) and residual SDs (.5,.6).

TECH1 order in `covariates-results.dat`: two means and two variances for each
profile (1:8), group intercept (9), profile intercept (10), group slope (11),
profile slope (12), group-class effect on profile logit (13). Conditional
profile intercepts are parameter 10+13 and parameter 10. Results then contain
13 standard errors and the summary values. Standard errors are retained but
are not compared because this covariate API does not yet expose inference.

Class labels are canonicalized by increasing first-indicator mean and then
decreasing first-profile intercept across group types. Profile label reversal
negates the profile logits; group label reversal negates the group logits.
Saved joint posterior patterns are marginalized and rows aligned by `id`.

| Quantity | Maximum absolute difference |
|---|---:|
| Log likelihood | 3.63e-5 |
| Means | 5.14e-8 |
| Variances | 1.76e-8 |
| Profile logit coefficients | 3.43e-6 |
| Group logit coefficients | 2.35e-6 |
| Individual posteriors | 3.34e-15 |
| Group posteriors | 5.90e-6 |

Eight-significant-digit Mplus RESULTS precision accounts for the likelihood
rounding allowance (5e-5). Parameter/posterior test tolerance is 1e-5, and
information-criterion tolerance is 5e-4. Input warnings about uncorrelated
variables are retained in the output.

Run `Rscript validation/mplus/covariates/compare.R` from the project root to
recreate `comparison.rds`, print the comparison, and refresh the offline
genuine-Mplus fixture. To rerun Mplus, invoke the installed executable on
`covariates.inp` with this directory as shell working directory.
