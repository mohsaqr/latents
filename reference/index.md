# Package index

## Fitting models

Two-level latent profile, latent class and latent transition models.

- [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  : Fit a two-level latent profile model
- [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
  : Fit a two-level latent class model
- [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) : Latent
  transition analysis
- [`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
  : Fit a two-level model in stages, holding the measurement solution
- [`starting_values()`](https://pak.dynasite.org/latents/reference/starting_values.md)
  : Build starting values for a multilevel latent profile fit

## Results

Every table of a fitted model, retrieved by name.

- [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
  : Every table a fitted model holds, from one verb
- [`descriptives()`](https://pak.dynasite.org/latents/reference/descriptives.md)
  : Describe the variables a model uses, or will use
- [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  [`as.data.frame(`*`<multilpa_diagnostics>`*`)`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  [`print(`*`<multilpa_diagnostics>`*`)`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  [`plot(`*`<multilpa_diagnostics>`*`)`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  : Every classification diagnostic, in one call
- [`report()`](https://pak.dynasite.org/latents/reference/report.md) :
  Everything about a fit, in one call

## Evaluating a model

- [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  : Enumerate numbers of individual profiles and group classes
- [`candidate_fit()`](https://pak.dynasite.org/latents/reference/candidate_fit.md)
  : Take one fitted model out of an enumeration grid
- [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
  : Parametric bootstrap likelihood-ratio comparison
- [`sensitivity()`](https://pak.dynasite.org/latents/reference/sensitivity.md)
  : Does the solution survive a different seed?
- [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  : Tidy inference for a covariate fit
- [`vcov(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/vcov.multilpa_transitions.md)
  [`parameter_inference(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/vcov.multilpa_transitions.md)
  [`confint(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/vcov.multilpa_transitions.md)
  : Standard errors are not available for a latent transition model

## External variables

Three-step analysis with a correction for classification error.

- [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
  : Relate a latent class to an outcome it did not help define
- [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md) :
  Covariates predicting class membership, corrected for
  misclassification

## Plots

- [`plot_views()`](https://pak.dynasite.org/latents/reference/plot_views.md)
  : The plots this package can draw
- [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  [`as.data.frame(`*`<multilpa_diagnostics>`*`)`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  [`print(`*`<multilpa_diagnostics>`*`)`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  [`plot(`*`<multilpa_diagnostics>`*`)`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  : Every classification diagnostic, in one call
- [`plot(`*`<multilpa>`*`)`](https://pak.dynasite.org/latents/reference/plot.multilpa.md)
  : Plot a fitted multilevel latent profile model
- [`plot(`*`<multilpa_bootstrap_lrt>`*`)`](https://pak.dynasite.org/latents/reference/plot.multilpa_bootstrap_lrt.md)
  : Plot a simulated bootstrap null distribution
- [`plot(`*`<multilpa_covariates>`*`)`](https://pak.dynasite.org/latents/reference/plot.multilpa_covariates.md)
  : Plot a covariate model
- [`plot(`*`<multilpa_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/plot.multilpa_enumeration.md)
  : Plot a class-enumeration grid
- [`plot(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/plot.multilpa_transitions.md)
  : Plot a fitted latent transition model

## Transition networks

Latent transitions handed to the tna package.

- [`get_tna()`](https://pak.dynasite.org/latents/reference/get_tna.md) :
  A transition network for the whole sample
- [`get_group_tna()`](https://pak.dynasite.org/latents/reference/get_group_tna.md)
  : A transition network for each latent group class

## Data

- [`course_engagement`](https://pak.dynasite.org/latents/reference/course_engagement.md)
  : Student engagement across a sequence of courses
- [`student_esm`](https://pak.dynasite.org/latents/reference/student_esm.md)
  : Leisure activities of university students in daily life

## Methods

Standard generics for fitted models, enumeration grids and summaries.

- [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  [`as.data.frame(`*`<multilpa_diagnostics>`*`)`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  [`print(`*`<multilpa_diagnostics>`*`)`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  [`plot(`*`<multilpa_diagnostics>`*`)`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  : Every classification diagnostic, in one call
- [`print(`*`<multilpa>`*`)`](https://pak.dynasite.org/latents/reference/print.multilpa.md)
  : Print a fitted multilevel latent profile model
- [`print(`*`<multilpa_bootstrap_lrt>`*`)`](https://pak.dynasite.org/latents/reference/print.multilpa_bootstrap_lrt.md)
  : Print a parametric bootstrap likelihood-ratio comparison
- [`print(`*`<multilpa_covariates>`*`)`](https://pak.dynasite.org/latents/reference/print.multilpa_covariates.md)
  : Print a covariate LPA fit
- [`print(`*`<multilpa_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/print.multilpa_enumeration.md)
  : Print a class-enumeration grid
- [`print(`*`<multilpa_start>`*`)`](https://pak.dynasite.org/latents/reference/print.multilpa_start.md)
  : Print a set of starting values
- [`print(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/print.multilpa_transitions.md)
  : Print a fitted latent transition model
- [`print(`*`<summary_multilpa>`*`)`](https://pak.dynasite.org/latents/reference/print.summary_multilpa.md)
  : Print a multilevel LPA summary
- [`print(`*`<summary_multilpa_bootstrap_lrt>`*`)`](https://pak.dynasite.org/latents/reference/print.summary_multilpa_bootstrap_lrt.md)
  : Print a bootstrap likelihood-ratio summary
- [`print(`*`<summary_multilpa_covariates>`*`)`](https://pak.dynasite.org/latents/reference/print.summary_multilpa_covariates.md)
  : Print a covariate LPA summary
- [`print(`*`<summary_multilpa_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/print.summary_multilpa_enumeration.md)
  : Print an enumeration summary
- [`print(`*`<summary_multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/print.summary_multilpa_transitions.md)
  : Print a latent transition summary
- [`summary(`*`<multilpa>`*`)`](https://pak.dynasite.org/latents/reference/summary.multilpa.md)
  : Summarize a fitted multilevel latent profile model
- [`summary(`*`<multilpa_bootstrap_lrt>`*`)`](https://pak.dynasite.org/latents/reference/summary.multilpa_bootstrap_lrt.md)
  : Summarise a parametric bootstrap likelihood-ratio comparison
- [`summary(`*`<multilpa_covariates>`*`)`](https://pak.dynasite.org/latents/reference/summary.multilpa_covariates.md)
  : Summarize a covariate LPA fit
- [`summary(`*`<multilpa_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/summary.multilpa_enumeration.md)
  : Summarise a class-enumeration grid
- [`summary(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/summary.multilpa_transitions.md)
  : Summarize a fitted latent transition model
- [`as.data.frame(`*`<multilpa>`*`)`](https://pak.dynasite.org/latents/reference/as.data.frame.multilpa.md)
  : Coerce a fitted multilevel latent profile model to its primary table
- [`as.data.frame(`*`<multilpa_bootstrap_lrt>`*`)`](https://pak.dynasite.org/latents/reference/as.data.frame.multilpa_bootstrap_lrt.md)
  : Coerce a bootstrap likelihood-ratio comparison to its primary table
- [`as.data.frame(`*`<multilpa_covariates>`*`)`](https://pak.dynasite.org/latents/reference/as.data.frame.multilpa_covariates.md)
  : Coerce a fitted covariate model to its primary table
- [`as.data.frame(`*`<multilpa_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/as.data.frame.multilpa_enumeration.md)
  : Coerce a class-enumeration grid to its primary table
- [`as.data.frame(`*`<multilpa_start>`*`)`](https://pak.dynasite.org/latents/reference/as.data.frame.multilpa_start.md)
  : Coerce a starting-value set to its primary table
- [`as.data.frame(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/as.data.frame.multilpa_transitions.md)
  : Coerce a fitted latent transition model to its primary table
- [`as.data.frame(`*`<summary_multilpa>`*`)`](https://pak.dynasite.org/latents/reference/as.data.frame.summary_multilpa.md)
  : Coerce a multilevel LPA summary to its primary table
- [`as.data.frame(`*`<summary_multilpa_bootstrap_lrt>`*`)`](https://pak.dynasite.org/latents/reference/as.data.frame.summary_multilpa_bootstrap_lrt.md)
  : Coerce a bootstrap comparison summary to its primary table
- [`as.data.frame(`*`<summary_multilpa_covariates>`*`)`](https://pak.dynasite.org/latents/reference/as.data.frame.summary_multilpa_covariates.md)
  : Coerce a covariate-model summary to its primary table
- [`as.data.frame(`*`<summary_multilpa_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/as.data.frame.summary_multilpa_enumeration.md)
  : Coerce a class-enumeration summary to its primary table
- [`as.data.frame(`*`<summary_multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/as.data.frame.summary_multilpa_transitions.md)
  : Coerce a latent transition model summary to its primary table
- [`coef(`*`<multilpa>`*`)`](https://pak.dynasite.org/latents/reference/coef.multilpa.md)
  : Extract multilevel LPA coefficients
- [`coef(`*`<multilpa_covariates>`*`)`](https://pak.dynasite.org/latents/reference/coef.multilpa_covariates.md)
  : Estimated parameters of a covariate fit
- [`coef(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/coef.multilpa_transitions.md)
  : Fitted parameters of a latent transition model
- [`confint(`*`<multilpa>`*`)`](https://pak.dynasite.org/latents/reference/confint.multilpa.md)
  : Wald confidence intervals for multilevel LPA coefficients
- [`confint(`*`<multilpa_covariates>`*`)`](https://pak.dynasite.org/latents/reference/confint.multilpa_covariates.md)
  : Wald confidence intervals for a covariate fit
- [`vcov(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/vcov.multilpa_transitions.md)
  [`parameter_inference(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/vcov.multilpa_transitions.md)
  [`confint(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/vcov.multilpa_transitions.md)
  : Standard errors are not available for a latent transition model
- [`vcov(`*`<multilpa>`*`)`](https://pak.dynasite.org/latents/reference/vcov.multilpa.md)
  : Extract multilevel LPA covariance estimates
- [`vcov(`*`<multilpa_covariates>`*`)`](https://pak.dynasite.org/latents/reference/vcov.multilpa_covariates.md)
  : Covariance matrix of a covariate fit
- [`logLik(`*`<multilpa>`*`)`](https://pak.dynasite.org/latents/reference/logLik.multilpa.md)
  : Extract the multilevel model log likelihood
- [`logLik(`*`<multilpa_covariates>`*`)`](https://pak.dynasite.org/latents/reference/logLik.multilpa_covariates.md)
  : Extract a covariate LPA log likelihood
- [`logLik(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/logLik.multilpa_transitions.md)
  : Log likelihood of a fitted latent transition model
- [`nobs(`*`<multilpa>`*`)`](https://pak.dynasite.org/latents/reference/nobs.multilpa.md)
  : Extract the number of independent groups
- [`nobs(`*`<multilpa_covariates>`*`)`](https://pak.dynasite.org/latents/reference/nobs.multilpa_covariates.md)
  : Count independent groups in a covariate LPA fit
- [`nobs(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/nobs.multilpa_transitions.md)
  : Number of independent units in a fitted latent transition model

## Conditions

- [`latents-conditions`](https://pak.dynasite.org/latents/reference/latents-conditions.md)
  : Conditions this package raises
