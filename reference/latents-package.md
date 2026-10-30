# Latent Profile, Class and Transition Models

latents is an R package for modelling unobserved heterogeneity in
response patterns, regression relationships and trajectories. It
supports latent profile and class analysis, multilevel models, latent
transition analysis, and mixture regression and growth models.

The package combines model estimation and comparison with tools for
assessing parameter uncertainty, classification quality and solution
stability. It also supports analyses relating latent membership to
external variables and provides tables and plots for interpreting and
reporting results.

## Details

Latent profiles describe heterogeneity in continuous indicators through
component-specific distributions. Latent classes describe heterogeneity
in discrete response patterns. In the two-level model, observation-level
profiles or classes are nested within groups, and latent group classes
differ in their conditional probabilities of observation-level
membership. The standard two-level formulation shares measurement
parameters across group classes. Transition models additionally estimate
dependence between latent states at successive occasions.

The number of latent classes and the measurement specification are
supplied by the analyst. Estimation uses multiple starting values
because different starts can converge to different local maxima.
Convergence, replication of the selected likelihood, class sizes,
parameter boundaries, and substantive interpretability should be
considered together when evaluating a solution. Classification summaries
describe uncertainty under the fitted model; they do not establish the
substantive validity of the classes.

## Model fitting

Use [`lpa`](https://pak.dynasite.org/latents/reference/lpa.md) or
[`lca`](https://pak.dynasite.org/latents/reference/lca.md) for
single-level models and
[`multilpa`](https://pak.dynasite.org/latents/reference/multilpa.md) or
[`multilca`](https://pak.dynasite.org/latents/reference/multilca.md) for
two-level models.
[`lta`](https://pak.dynasite.org/latents/reference/lta.md) fits latent
transition models, and
[`mixture_regression`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
fits finite mixtures of regression models. The `family` argument of
`multilpa` provides additional model families documented on that
function's help page.
[`fit_staged`](https://pak.dynasite.org/latents/reference/fit_staged.md)
estimates a two-level model with measurement parameters fixed at their
first-stage estimates.

## Model evaluation and inference

[`enumerate_classes`](https://pak.dynasite.org/latents/reference/enumerate_classes.md),
[`enumerate_lpa`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md),
[`enumerate_lca`](https://pak.dynasite.org/latents/reference/enumerate_lca.md),
and
[`enumerate_regressions`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md)
fit candidate model sets.
[`compare_models`](https://pak.dynasite.org/latents/reference/compare_models.md)
summarizes selected models, and
[`bootstrap_lrt`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
provides parametric bootstrap likelihood ratio tests for supported
comparisons.

[`get_results`](https://pak.dynasite.org/latents/reference/get_results.md)
extracts model-specific tables;
[`diagnostics`](https://pak.dynasite.org/latents/reference/diagnostics.md)
and
[`sensitivity`](https://pak.dynasite.org/latents/reference/sensitivity.md)
assess fitted solutions and sensitivity to starting values.
[`parameter_inference`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
provides parameter estimates, standard errors, and intervals for
supported model families.
[`three_step`](https://pak.dynasite.org/latents/reference/three_step.md)
and [`r3step`](https://pak.dynasite.org/latents/reference/r3step.md)
relate external variables to latent membership with correction for
classification error.
[`pool_imputations`](https://pak.dynasite.org/latents/reference/pool_imputations.md)
combines compatible analyses across imputed datasets.

## Interpretation and assumptions

The measurement distribution, conditional independence assumptions,
covariance constraints, and treatment of time depend on the selected
model. Consult the fitting function for these assumptions and for
restrictions on missing data, sampling weights, covariates, and
parameter inference. Posterior class assignments are estimates with
uncertainty. Class labels are arbitrary identifiers and should be
interpreted through the estimated response patterns. Numerical fit and
classification separation alone do not establish causal effects or
distinct populations.

## Getting started

The vignettes
[`vignette("lpa", package = "latents")`](https://pak.dynasite.org/latents/articles/lpa.md),
[`vignette("lca", package = "latents")`](https://pak.dynasite.org/latents/articles/lca.md),
[`vignette("enumeration", package = "latents")`](https://pak.dynasite.org/latents/articles/enumeration.md),
and
[`vignette("lta", package = "latents")`](https://pak.dynasite.org/latents/articles/lta.md)
introduce the main workflows. See
[`latents-plots`](https://pak.dynasite.org/latents/reference/latents-plots.md)
for graphical summaries,
[`latents-predict`](https://pak.dynasite.org/latents/reference/latents-predict.md)
for prediction, and
[`latents-conditions`](https://pak.dynasite.org/latents/reference/latents-conditions.md)
for diagnostic conditions.

## References

Vermunt, J. K. (2003). Multilevel latent class models. *Sociological
Methodology*, 33, 213–239.
[doi:10.1111/j.0081-1750.2003.t01-1-00131.x](https://doi.org/10.1111/j.0081-1750.2003.t01-1-00131.x)
.
