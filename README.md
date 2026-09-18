# multilpa

A dependency-light R package for **two-level latent profile and latent class
analysis**. All estimation runs in R; Carm, Mplus, and a server are not required.

Indicators may be continuous, categorical, or a mixture of the two. The model
identifies individual profiles from the indicators and latent group classes from
differences in profile proportions. For example, students
have achievement profiles, while schools differ in the prevalence of those
profiles. The implementation uses nested expectation-maximization (EM) with
multiple starts.

## Install and run

From this project directory:

```sh
R CMD INSTALL .
```

Then, with your own data frame:

```r
library(multilpa)

fit <- multilpa(
  data = students,
  indicators = c("reading", "maths", "engagement"),
  group = "school_id",
  n_profiles = 3,
  n_group_classes = 2,
  variance_model = "varying",
  n_starts = 20,
  seed = 42
)

summary(fit)
as.data.frame(fit)                                   # one row per profile and indicator
as.data.frame(fit, what = "profile_probabilities")   # one row per group class and profile
as.data.frame(fit, what = "posteriors")              # one row per input individual
as.data.frame(fit, what = "group_posteriors")        # one row per observed group
as.data.frame(fit, what = "starts")                  # one row per EM start
information_criteria(fit)                            # every information criterion
classification_table(fit, level = "both")            # classification quality
entropy_table(fit)                                   # entropy at each level

plot(fit)                                            # profile means by indicator
plot(fit, scale = "standardized")                    # comparable profile shapes
plot(fit, what = "probabilities")                    # prevalence by group class
```

Every result table is a base `data.frame`, so nothing needs to be pulled out of
the fitted object by hand.

`students` is a placeholder for your data. For a complete reproducible example
using simulated data, run:

```sh
Rscript validation/synthetic-demo.R
```

The script checks recovery, prints the data and fitted structure, and saves
the data, true classes, and fitted model to `tmp/synthetic-demo.rds`. It sources
the local implementation, so package installation is optional for this example.

## Model

For person `i` in group `j`, the vector of continuous indicators has a Gaussian
profile `k`. By default, the indicators are conditionally independent within a profile.
The latent group class `h` determines the prior profile probabilities.

```text
P(Y_j) = sum_h eta_h * product_i [sum_k pi_hk * product_d
                                Normal(y_ijd | mu_kd, variance_kd)]
```

- `means[k, d]` and `variances[k, d]` are shared across group classes.
- `profile_probabilities[h, k]` describes profile prevalence within group class h.
- `group_probabilities[h]` describes prevalence across groups, giving each
  observed group one contribution, regardless of its size.
- Individual posteriors average over uncertainty in the group's latent class.
- `variance_model = "equal"` shares each indicator's variance across profiles;
  different indicators can still have different variances.
- `n_group_classes = 1` reduces to ordinary pooled LPA.

The formulation extends the nonparametric multilevel mixture framework of
[Vermunt (2003)](https://jeroenvermunt.nl/sm2003.pdf) to Gaussian indicators.
It assumes measurement invariance across group classes. Full residual covariance,
membership covariates, and an alternative continuous random-intercept model are
available through the extensions below. Longitudinal transitions and random
slopes are not implemented.

## Categorical indicators

Name any indicators that should be treated as categorical, and they are modelled
by unrestricted, profile-specific response probabilities instead of Gaussian
means and variances. That is the latent class measurement model, so the same
verb fits two-level LPA, two-level LCA, and mixtures of the two.

```r
# Two-level latent class analysis: every indicator categorical.
lca <- multilpa(survey, c("u1", "u2", "u3", "u4", "u5"), "school_id",
                  n_profiles = 3, n_group_classes = 2,
                  categorical = c("u1", "u2", "u3", "u4", "u5"),
                  n_starts = 40, seed = 42)

as.data.frame(lca, what = "responses")   # probabilities and logit thresholds
plot(lca, what = "responses")            # response curves by profile

# Mixed measurement: name only the categorical ones.
mixed <- multilpa(survey, c("reading", "maths", "u1", "u2"), "school_id",
                    n_profiles = 3, n_group_classes = 2,
                    categorical = c("u1", "u2"), n_starts = 40, seed = 42)
```

Binary, ordinal and unordered indicators all use the same unrestricted
parameterization, which is what mixture software estimates by default when every
threshold is free. Numeric, integer, logical, character and factor columns are
accepted; factor levels keep their declared order, and everything else is
ordered by `sort()`, so the coding never depends on row order. Category counts
come from the observed values.

`missing = "fiml"` extends to categorical indicators: an unobserved cell simply
drops out of that individual's likelihood, which is the observed-data likelihood
under an ignorable mechanism, with no imputation.

Standard errors and the parametric bootstrap do **not** yet cover categorical
indicators; both refuse such a fit rather than returning a number from the
Gaussian score functions. `min_probability` bounds every response probability
away from zero, in the same way `min_variance` bounds variances, and the bound
is applied as the exact constrained solution rather than by rescaling.

## Additional methods

```r
# Observed-data ML for incomplete indicators; correlated residuals within profiles.
fit <- multilpa(students, c("reading", "maths", "engagement"), "school_id",
                  n_profiles = 3, n_group_classes = 2,
                  covariance_model = "full", missing = "fiml",
                  n_starts = 20, seed = 42)

# Observed-information standard errors and Wald intervals.
as.data.frame(parameter_inference(fit, data = students))
confint(fit, data = students)

# MLR robust sandwich errors, accumulating the score over independent groups.
robust <- parameter_inference(fit, data = students, vcov_type = "robust")
as.data.frame(robust)
vcov(fit, data = students, vcov_type = "robust")

# Likelihood-ratio statistic with the Lo-Mendell-Rubin adjustment.
# No p-value is returned; use the bootstrap below for a calibrated one.
lmr_lrt(smaller_fit, larger_fit)

# One-step membership regressions (numeric predictors, complete indicators).
with_predictors <- fit_covariates(
  students, c("reading", "maths", "engagement"), "school_id", 3, 2,
  profile_covariates = "age", group_covariates = "school_resources", seed = 42)
with_predictors$profile_coefficients   # covariate models have no tidy accessor yet
with_predictors$group_coefficients

# Alternative: one continuous group intercept, loading 1 on every indicator.
random_intercept <- fit_random_intercept(
  students, c("reading", "maths", "engagement"), "school_id", 3, seed = 42)
random_intercept$quadrature_log_likelihood_difference

# Compare profile/group-class counts; inspect diagnostics in every row.
candidates <- enumerate_classes(
  students, c("reading", "maths", "engagement"), "school_id",
  profiles = 2:4, group_classes = 1:3, n_starts = 20, seed = 42)
as.data.frame(candidates)
plot(candidates, criterion = "sabic_individual")

# Complete-data nested models differing by one class at one level.
# bootstrap_lrt(smaller, larger, students, n_boot = 199, seed = 42)
```

| Capability | Supported scope |
|---|---|
| Measurement model | Gaussian, categorical, or mixed, via `categorical`; categorical indicators use unrestricted profile-specific response probabilities |
| Missing indicators | `missing="fiml"`: observed Gaussian marginals and observed categorical responses, ignorable missingness assumption; no missing covariates |
| Residual covariance | `covariance_model="diagonal"` or `"full"`; shared or profile-specific via `variance_model` |
| SEs/CIs | Observed-Hessian ML inference for the Gaussian discrete model, including full covariance and FIML; not available for categorical indicators |
| Robust SEs | `vcov_type="robust"`: Huber-White sandwich over independent groups, plus the MLR scaling correction factor; matches Mplus `ESTIMATOR=MLR` |
| Membership covariates | Numeric predictors at both levels; shared individual-profile slopes across group classes; complete diagonal model |
| Continuous random effect | One shared Gaussian group intercept with unit indicator loadings; complete diagonal model, no discrete group classes |
| Class enumeration | Grid of discrete models, every information criterion under both sample-size conventions, entropy, failed-fit and convergence diagnostics |
| Information criteria | AIC, BIC, SABIC, CAIC, AWE and ICL, each reported under both the group-count and individual-count conventions |
| Classification quality | Modal and model-estimated class sizes, average posterior probabilities, odds of correct classification, relative entropy |
| LMR statistic | Likelihood-ratio statistic and the Lo-Mendell-Rubin adjustment; **no p-value**, because the VLMR reference distribution is not reproduced |
| Bootstrap LRT | Parametric bootstrap preserving group sizes; complete discrete models differing by one class; not Mplus TECH14 |
| Plots | Profile means (raw or standardized), prevalence by group class, and any enumeration criterion; base graphics only |

These are explicit model families, not every combination of Mplus options.
Covariate and random-intercept fits currently do not provide standard errors.
Random-intercept mixtures use Gaussian quadrature; increase node counts and
refit if the higher-order likelihood check fails. Inference rejects variance
boundaries and nonpositive information matrices. Natural-scale Wald intervals
can extend beyond parameter bounds; log-variance/logit coordinates are also
available. A bootstrap p-value is withheld if any replicate fails convergence
or has a reversed likelihood. Enumeration never automatically declares a winner.

## Plots

`plot()` methods are provided for fitted models and for enumeration grids. They
use base graphics only; no plotting package is added as a dependency.

Every series is distinguished by **colour, point symbol and line type together**
and labelled directly at the line end, so the plots stay readable in greyscale,
for colour-blind readers, and without a legend. Colours are Okabe-Ito.

```r
plot(fit)                                   # means in input units
plot(fit, scale = "standardized")           # shapes comparable across indicators
plot(fit, what = "probabilities")           # profile prevalence per group class
plot(candidates, criterion = "icl_groups")  # any column of the enumeration grid
```

Use `scale = "standardized"` whenever indicators are on different scales,
otherwise the largest-scale indicator dictates the apparent profile shape. The
divisor is each indicator's observed standard deviation, not its within-profile
residual standard deviation, so the values are comparable but are not effect
sizes.

Plots show point estimates only, carry no standard errors, and profile order is
arbitrary. In an enumeration plot the ringed point marks the lowest criterion
value, which identifies an extremum and does not select a model; candidates that
failed to converge appear as crosses on the baseline rather than being dropped.

No visual constant is hard-coded. Pass any of them inline:

```r
plot(fit, palette = c("#0072B2", "#D55E00", "#CC79A7"),
     panel_fill = "#FFFFFF", point_size = 1.8, line_width = 3,
     main = "My title", subtitle = "My subtitle", labels = FALSE)
```

## Diagnostics and interpretation

The result retains input row order. Group order follows first occurrence in
the data; `group_values` preserves the original identifiers, and `group_ids`
provides unique display labels. Profile and group-class numbers are arbitrary:
align labels before comparing separate fits.

Row order is part of the effective seed. Initialization runs k-means on the
rows as supplied, so the same data in a different row order can converge to a
different local optimum with the same `seed`. The likelihood itself is exactly
row-order invariant given identical starting values; it is the search that
differs. Fix the row order alongside the seed when a fit must be reproducible,
and raise `n_starts` until `n_best_replicated` is comfortably above one.

Inspect `converged`, `starts`, `n_best_replicated`, and
`log_likelihood_history`. The best finite likelihood is returned even when
that start did not converge, with a warning. Failed starts are recorded with
their error messages. More starts and replicated best likelihoods improve
confidence but do not establish a global maximum. Convergence is based on
relative likelihood change, not parameter stability or a Hessian test.

`min_variance` is an explicit lower bound in **squared input units**, defaulting
to `1e-6`. Gaussian mixtures with freely varying variances otherwise allow
singular solutions. The returned estimates solve a constrained likelihood
problem; `boundary = TRUE` and a warning identify a variance at the bound.
Choose a sensible bound for the indicator scales and examine sensitivity.
`small_classes` flags effective memberships below one; effective class counts
are also printed by `summary()`.

Missing indicators are rejected by default; use `missing="fiml"` to integrate
them out of the likelihood. Completely unobserved rows contribute no indicator
information, but are retained for posterior prediction. Entirely unobserved
indicators, infinite values, constant indicators, and missing group IDs are
rejected. Mean filling is used only to initialize incomplete-data fits.
Numeric, character,
and factor group identifiers are supported. Group-class identification requires
more than basic count checks: indistinguishable profile distributions or group
prevalence vectors can produce unidentified models. The package rejects some
obviously unidentified cases but does not certify general identification.

### Information criteria

Let `J` be the number of observed groups, `N` the number of individuals, and
`q` the number of free parameters:

| Field | Definition |
|---|---|
| `aic` / `AIC(fit)` | `-2 * logLik + 2 * q` |
| `bic` / `bic_groups` / `BIC(fit)` | `-2 * logLik + log(J) * q` |
| `bic_individual` | `-2 * logLik + log(N) * q` |

The BIC sample-size convention matters in multilevel mixtures.
[Lukočienė, Varriale, and Vermunt (2010)](https://doi.org/10.1111/j.1467-9531.2010.01231.x)
studied this issue. Both conventions are exposed; use one consistently when
comparing multilevel candidates fitted to the same observations and indicators.
For ordinary pooled-LPA enumeration (`n_group_classes = 1`), use
`bic_individual`; group-count BIC has zero penalty if there is only one observed
group. `nobs(fit)` and the `logLik` object's `nobs` attribute return `J`.

For full covariance, the variance constraint applies to every covariance
eigenvalue. New discrete fits retain `indicator_data` to verify the original
observations before inference and bootstrap comparisons.

## Verification

```sh
Rscript -e 'pkgload::load_all("."); testthat::test_dir("tests/testthat")'
R CMD build .
R CMD check --no-manual multilpa_0.4.0.tar.gz
```

`pkgload` is only a development convenience. Installed-package testing via
`R CMD check` does not require it. Estimation uses only base R and `stats`;
`testthat` and `mclust` are optional testing dependencies.

Verification includes:

| Check | Observed agreement on retained synthetic examples |
|---|---|
| Exhaustive latent-assignment likelihood and posteriors | posterior differences below `1e-15` |
| Independent first EM update, both variance constraints | within `1e-9` |
| Direct `stats::optim()` BFGS likelihood optimization | likelihood difference `7.82e-11`; parameter differences below `6e-7` |
| `mclust` single-level VVI and EEI limits | likelihood differences below `6e-14` |
| Simulated parameter/class recovery | tested after aligning arbitrary labels |
| Input, dimension, numerical, restart, and RNG edge cases | covered by regression tests |

The `mclust` limit checks initialize at the reference solution to test that it
is a fixed point of the new engine; they do not compare global searches.
Independent likelihood enumeration and direct optimization validate the
multilevel calculations on tested examples.

### Direct Mplus comparisons

The original six comparisons cover three datasets:

- Published Mplus examples **7.9 and 7.10**, using the same 500-person,
  four-indicator dataset: the two single-level variance models reproduce the
  published estimates within their displayed precision.
- New Mplus 9 Demo runs on **1,200 simulated individuals in 60 groups**:
  both two-level variance models agree on parameter estimates within `6e-8`
  and on marginal posterior probabilities within `2e-10`.
- New Mplus 9 Demo runs on the **public Example 10.4 dataset**, with 1,000
  individuals in 110 groups and five indicators: both matched two-level LPA
  models agree on parameter estimates within `1.1e-5` and posteriors within
  `3e-5`. These are new LPA fits; the original published CFA mixture is a
  different model and is not used as the numerical target.

A genuine `CATEGORICAL =` run validates the two-level latent class model on
1,200 individuals in 60 groups with five binary indicators: thresholds agree
within `5.8e-6`, profile probabilities within `1.1e-6`, group-class
probabilities within `8.3e-7`, and the likelihood within `2.8e-6`. The
single-level limit additionally agrees with `poLCA` to `2.7e-11` in likelihood
and `7.9e-8` in response probabilities, with identical free-parameter counts.

Additional genuine Mplus runs validate full covariance with missing indicators
(both equal and varying covariance), one-step membership covariates at both
levels, and the single-profile continuous random-intercept limit. Observed
information standard errors are also compared with Mplus.

A further genuine `ESTIMATOR = MLR` run validates the robust sandwich errors and
the new information criteria: robust standard errors agree within `1.2e-7`, the
MLR scaling correction factor within `7.1e-7`, and AIC, BIC and SABIC within
`4.1e-5`. Two genuine `TECH11` runs reproduce the Lo-Mendell-Rubin adjusted
statistic to the three decimals Mplus prints, end to end from native fits on the
published Example 7.9 data (`747.0633` against `747.063`, adjusted `723.7708`
against `723.771`). The Vuong-Lo-Mendell-Rubin reference distribution, that is
the mean, standard deviation and p-value Mplus reports beside it, is **not**
reproduced, and no p-value is returned in its place. See
[the extension report](validation/mplus/EXTENSIONS.md) for exact scope and
tolerances. The random-intercept mixture likelihood additionally agrees with
independent adaptive numerical integration; its Mplus comparison currently
covers one profile only.

All new two-level comparisons use independent starts in R and Mplus. Class labels and saved subject
IDs are aligned before comparison. Mplus's BIC is compared to
`bic_individual`, not this package's default group-count BIC. Saved Mplus
likelihoods and information criteria have finite output precision, accounted
for explicitly in the assertions. Tests run offline from retained fixtures;
neither Mplus nor internet access is needed to run them.

See [the comparison report](validation/mplus/COMPARISON.md) for numerical
tables, original URLs, raw Mplus artifacts, reproduction commands, and scope.
These checks establish agreement for the tested specifications, not full
Mplus feature parity or a guarantee for every dataset.

The independent reference implementation lives in
`tests/testthat/helper-independent-likelihood.R`; it uses explicit assignment
enumeration and direct probability arithmetic, independently of the fitting
engine's log-domain EM calculation.
