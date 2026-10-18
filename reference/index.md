# Package index

## Fitting models

Two-level latent profile, latent class and latent transition models.

- [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md) : Fit a
  latent profile model
- [`lca()`](https://pak.dynasite.org/latents/reference/lca.md) : Fit a
  latent class model
- [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  : Fit a two-level latent profile model
- [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
  : Fit a two-level latent class model
- [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) : Fit a
  latent transition model
- [`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
  : Fit a two-level model in stages, holding the measurement solution
- [`starting_values()`](https://pak.dynasite.org/latents/reference/starting_values.md)
  : Extract starting values for a latent profile fit
- [`prior_control()`](https://pak.dynasite.org/latents/reference/prior_control.md)
  [`print(`*`<latents_prior>`*`)`](https://pak.dynasite.org/latents/reference/prior_control.md)
  : Specify a conjugate prior for Gaussian mixture estimation

## Results

Every table of a fitted model, retrieved by name.

- [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
  : Extract the result tables of a fitted model
- [`descriptives()`](https://pak.dynasite.org/latents/reference/descriptives.md)
  : Descriptive statistics and intraclass correlations of model
  variables
- [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  [`as.data.frame(`*`<multilpa_diagnostics>`*`)`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  [`print(`*`<multilpa_diagnostics>`*`)`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  [`plot(`*`<multilpa_diagnostics>`*`)`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  : Classification diagnostics of a fitted model
- [`report()`](https://pak.dynasite.org/latents/reference/report.md) :
  Print a full report of a fitted model

## Evaluating a model

- [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md)
  : Compare latent profile models
- [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md)
  : Compare latent class models
- [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  : Enumerate numbers of individual profiles and group classes
- [`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md)
  [`candidate_fit()`](https://pak.dynasite.org/latents/reference/get_candidate.md)
  : Take one fitted model out of an enumeration grid or a comparison
- [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
  : Parametric bootstrap likelihood-ratio test of two nested models
- [`sensitivity()`](https://pak.dynasite.org/latents/reference/sensitivity.md)
  : Refit a model under several seeds
- [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  : Tidy inference for a covariate fit
- [`parameter_inference(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/parameter_inference.multilpa_transitions.md)
  [`vcov(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/parameter_inference.multilpa_transitions.md)
  [`confint(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/parameter_inference.multilpa_transitions.md)
  : Standard errors and intervals for a latent transition model
- [`pool_imputations()`](https://pak.dynasite.org/latents/reference/pool_imputations.md)
  : Fit a model to multiply imputed data and pool it with Rubin's rules

## Mixture regression and growth models

Regressions whose coefficients differ across latent classes, at one or
two levels, for continuous, binary and count outcomes; classes of
trajectories, with random effects within classes (growth mixture
models).

- [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
  : Fit a finite mixture of regressions
- [`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md)
  [`as.data.frame(`*`<latents_regression_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md)
  [`get_results(`*`<latents_regression_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md)
  [`print(`*`<latents_regression_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md)
  : Compare mixture regressions with different numbers of classes
- [`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)
  : Compare fitted models on the same data
- [`predict(`*`<latents_growth_mixture>`*`)`](https://pak.dynasite.org/latents/reference/latents-predict.md)
  [`simulate(`*`<latents_growth_mixture>`*`)`](https://pak.dynasite.org/latents/reference/latents-predict.md)
  [`predict(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/latents-predict.md)
  [`simulate(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/latents-predict.md)
  [`predict(`*`<multilpa>`*`)`](https://pak.dynasite.org/latents/reference/latents-predict.md)
  [`predict(`*`<multilpa_covariates>`*`)`](https://pak.dynasite.org/latents/reference/latents-predict.md)
  [`predict(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/latents-predict.md)
  : Predict from and simulate a fitted model

## External variables

Three-step analysis with a correction for classification error.

- [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
  : Estimate a distal outcome's mean in each latent class
- [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md) :
  Regress latent class membership on covariates by the three-step
  approach

## Plots

- [`plot(`*`<latents_family_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot(`*`<multilpa_additive>`*`)`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot(`*`<latents_comparison>`*`)`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot(`*`<multilpa_cross_level>`*`)`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot(`*`<multilpa_bootstrap_lrt>`*`)`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot(`*`<latents_growth_mixture>`*`)`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot(`*`<latents_transition_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot(`*`<multilpa_lta>`*`)`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot_profiles()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot_responses()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot_sequences()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot_prevalence()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot_diagnostics()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot_transitions()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot_tree()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot_selection()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot(`*`<multilpa>`*`)`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot(`*`<multilpa_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot_enumeration()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot(`*`<multilpa_covariates>`*`)`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot_views()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot(`*`<latents_pooled>`*`)`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  : Plot latents results

## Transition networks

Latent transitions handed to the tna package.

- [`get_tna()`](https://pak.dynasite.org/latents/reference/get_tna.md) :
  Build one transition network for the whole sample
- [`get_group_tna()`](https://pak.dynasite.org/latents/reference/get_group_tna.md)
  : Build one transition network per latent group class

## Data

- [`course_engagement`](https://pak.dynasite.org/latents/reference/course_engagement.md)
  : Student engagement across a sequence of courses
- [`student_esm`](https://pak.dynasite.org/latents/reference/student_esm.md)
  : Leisure activities of university students in daily life
- [`study_hours`](https://pak.dynasite.org/latents/reference/study_hours.md)
  : Study hours and quiz scores under two study strategies
- [`growth_scores`](https://pak.dynasite.org/latents/reference/growth_scores.md)
  : Achievement growth of students with three kinds of trajectory
- [`growth_schools`](https://pak.dynasite.org/latents/reference/growth_schools.md)
  : Reading growth of students nested in schools
- [`srl`](https://pak.dynasite.org/latents/reference/srl.md) :
  Self-regulated learning scores simulated by a large language model

## Methods

Coefficients, covariance, intervals and likelihood of a fit.

- [`vcov(`*`<multilpa_covariates>`*`)`](https://pak.dynasite.org/latents/reference/latents-model-methods.md)
  [`coef(`*`<multilpa_covariates>`*`)`](https://pak.dynasite.org/latents/reference/latents-model-methods.md)
  [`confint(`*`<multilpa_covariates>`*`)`](https://pak.dynasite.org/latents/reference/latents-model-methods.md)
  [`logLik(`*`<multilpa_covariates>`*`)`](https://pak.dynasite.org/latents/reference/latents-model-methods.md)
  [`nobs(`*`<multilpa_covariates>`*`)`](https://pak.dynasite.org/latents/reference/latents-model-methods.md)
  [`coef(`*`<multilpa>`*`)`](https://pak.dynasite.org/latents/reference/latents-model-methods.md)
  [`vcov(`*`<multilpa>`*`)`](https://pak.dynasite.org/latents/reference/latents-model-methods.md)
  [`confint(`*`<multilpa>`*`)`](https://pak.dynasite.org/latents/reference/latents-model-methods.md)
  [`logLik(`*`<multilpa>`*`)`](https://pak.dynasite.org/latents/reference/latents-model-methods.md)
  [`nobs(`*`<multilpa>`*`)`](https://pak.dynasite.org/latents/reference/latents-model-methods.md)
  [`vcov(`*`<latents_pooled>`*`)`](https://pak.dynasite.org/latents/reference/latents-model-methods.md)
  [`logLik(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/latents-model-methods.md)
  [`nobs(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/latents-model-methods.md)
  [`coef(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/latents-model-methods.md)
  : Coefficients, covariance, intervals and likelihood of latents fits

## Conditions

- [`latents-conditions`](https://pak.dynasite.org/latents/reference/latents-conditions.md)
  : Conditions raised by the package
