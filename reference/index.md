# Package index

## Fitting models

Two-level latent profile, latent class and latent transition models.

- [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md) : Latent
  profile analysis
- [`lca()`](https://pak.dynasite.org/latents/reference/lca.md) : Latent
  class analysis
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
- [`prior_control()`](https://pak.dynasite.org/latents/reference/prior_control.md)
  [`print(`*`<latents_prior>`*`)`](https://pak.dynasite.org/latents/reference/prior_control.md)
  : Conjugate prior for Gaussian mixture estimation

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

- [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md)
  : Compare latent profile models
- [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md)
  : Compare latent class models
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
- [`parameter_inference(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/parameter_inference.multilpa_transitions.md)
  [`vcov(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/parameter_inference.multilpa_transitions.md)
  [`confint(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/parameter_inference.multilpa_transitions.md)
  : Wald inference for a latent transition model
- [`pool_imputations()`](https://pak.dynasite.org/latents/reference/pool_imputations.md)
  : Fit a model to multiply imputed data and pool it with Rubin's rules
- [`get_results(`*`<latents_pooled>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_pooled.md)
  [`as.data.frame(`*`<latents_pooled>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_pooled.md)
  : Tables of a pooled multiply imputed fit
- [`get_results(`*`<multilpa_additive>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_additive.md)
  [`as.data.frame(`*`<multilpa_additive>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_additive.md)
  [`coef(`*`<multilpa_additive>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_additive.md)
  [`vcov(`*`<multilpa_additive>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_additive.md)
  [`confint(`*`<multilpa_additive>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_additive.md)
  [`logLik(`*`<multilpa_additive>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_additive.md)
  [`nobs(`*`<multilpa_additive>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_additive.md)
  : Tables of an additive group-class fit
- [`get_results(`*`<multilpa_cross_level>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_cross_level.md)
  [`as.data.frame(`*`<multilpa_cross_level>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_cross_level.md)
  [`print(`*`<multilpa_cross_level>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_cross_level.md)
  [`summary(`*`<multilpa_cross_level>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_cross_level.md)
  [`logLik(`*`<multilpa_cross_level>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_cross_level.md)
  [`nobs(`*`<multilpa_cross_level>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_cross_level.md)
  : Tables of a cross-level fit
- [`get_results(`*`<latents_family_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_family_enumeration.md)
  [`as.data.frame(`*`<latents_family_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_family_enumeration.md)
  [`print(`*`<latents_family_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_family_enumeration.md)
  [`summary(`*`<latents_family_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_family_enumeration.md)
  : Tables of a group-class family enumeration

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
- [`get_results(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_mixture_regression.md)
  [`as.data.frame(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_mixture_regression.md)
  [`coef(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_mixture_regression.md)
  [`vcov(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_mixture_regression.md)
  [`confint(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_mixture_regression.md)
  [`logLik(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_mixture_regression.md)
  [`nobs(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_mixture_regression.md)
  : Tables of a mixture-of-regressions fit
- [`get_results(`*`<latents_growth_mixture>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_growth_mixture.md)
  [`as.data.frame(`*`<latents_growth_mixture>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_growth_mixture.md)
  [`coef(`*`<latents_growth_mixture>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_growth_mixture.md)
  [`vcov(`*`<latents_growth_mixture>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_growth_mixture.md)
  [`confint(`*`<latents_growth_mixture>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_growth_mixture.md)
  [`logLik(`*`<latents_growth_mixture>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_growth_mixture.md)
  [`nobs(`*`<latents_growth_mixture>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_growth_mixture.md)
  : Tidy results of a growth mixture model
- [`predict(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/predict.latents_mixture_regression.md)
  : Predict from a mixture-of-regressions fit
- [`simulate(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/simulate.latents_mixture_regression.md)
  : Simulate outcomes from a mixture-of-regressions fit
- [`simulate(`*`<latents_growth_mixture>`*`)`](https://pak.dynasite.org/latents/reference/simulate.latents_growth_mixture.md)
  : Simulate outcomes from a growth mixture model

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
- [`plot(`*`<latents_comparison>`*`)`](https://pak.dynasite.org/latents/reference/plot.latents_comparison.md)
  : Plot a model comparison
- [`plot(`*`<latents_family_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/plot.latents_family_enumeration.md)
  : Plot information criteria across group-class families
- [`plot(`*`<latents_growth_mixture>`*`)`](https://pak.dynasite.org/latents/reference/plot.latents_growth_mixture.md)
  : Plot a growth mixture model
- [`plot(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/plot.latents_mixture_regression.md)
  : Plot a mixture-of-regressions fit
- [`plot(`*`<latents_pooled>`*`)`](https://pak.dynasite.org/latents/reference/plot.latents_pooled.md)
  : Plot a pooled multiply imputed fit
- [`plot(`*`<latents_transition_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/plot.latents_transition_enumeration.md)
  : Plot information criteria across transition models
- [`plot(`*`<multilpa>`*`)`](https://pak.dynasite.org/latents/reference/plot.multilpa.md)
  : Plot a fitted multilevel latent profile model
- [`plot(`*`<multilpa_additive>`*`)`](https://pak.dynasite.org/latents/reference/plot.multilpa_additive.md)
  : Plot an additive group-class fit
- [`plot(`*`<multilpa_bootstrap_lrt>`*`)`](https://pak.dynasite.org/latents/reference/plot.multilpa_bootstrap_lrt.md)
  : Plot a simulated bootstrap null distribution
- [`plot(`*`<multilpa_covariates>`*`)`](https://pak.dynasite.org/latents/reference/plot.multilpa_covariates.md)
  : Plot a covariate model
- [`plot(`*`<multilpa_cross_level>`*`)`](https://pak.dynasite.org/latents/reference/plot.multilpa_cross_level.md)
  : Plot a cross-level fit
- [`plot(`*`<multilpa_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/plot.multilpa_enumeration.md)
  : Plot a class-enumeration grid
- [`plot(`*`<multilpa_lta>`*`)`](https://pak.dynasite.org/latents/reference/plot.multilpa_lta.md)
  : Plot transition probabilities of a general transition fit
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
- [`study_hours`](https://pak.dynasite.org/latents/reference/study_hours.md)
  : Study hours and quiz scores under two study strategies
- [`growth_scores`](https://pak.dynasite.org/latents/reference/growth_scores.md)
  : Achievement growth of students with three kinds of trajectory
- [`growth_schools`](https://pak.dynasite.org/latents/reference/growth_schools.md)
  : Reading growth of students nested in schools
- [`srl`](https://pak.dynasite.org/latents/reference/srl.md) :
  Self-regulated learning scores simulated by a large language model

## Methods

Standard generics for fitted models, enumeration grids and summaries.

- [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  [`as.data.frame(`*`<multilpa_diagnostics>`*`)`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  [`print(`*`<multilpa_diagnostics>`*`)`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  [`plot(`*`<multilpa_diagnostics>`*`)`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  : Every classification diagnostic, in one call
- [`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md)
  [`as.data.frame(`*`<latents_regression_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md)
  [`get_results(`*`<latents_regression_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md)
  [`print(`*`<latents_regression_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md)
  : Compare mixture regressions with different numbers of classes
- [`get_results(`*`<latents_family_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_family_enumeration.md)
  [`as.data.frame(`*`<latents_family_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_family_enumeration.md)
  [`print(`*`<latents_family_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_family_enumeration.md)
  [`summary(`*`<latents_family_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_family_enumeration.md)
  : Tables of a group-class family enumeration
- [`get_results(`*`<latents_transition_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_transition_enumeration.md)
  [`as.data.frame(`*`<latents_transition_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_transition_enumeration.md)
  [`print(`*`<latents_transition_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_transition_enumeration.md)
  [`summary(`*`<latents_transition_enumeration>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_transition_enumeration.md)
  : Tables of a transition-model enumeration
- [`get_results(`*`<multilpa_cross_level>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_cross_level.md)
  [`as.data.frame(`*`<multilpa_cross_level>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_cross_level.md)
  [`print(`*`<multilpa_cross_level>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_cross_level.md)
  [`summary(`*`<multilpa_cross_level>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_cross_level.md)
  [`logLik(`*`<multilpa_cross_level>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_cross_level.md)
  [`nobs(`*`<multilpa_cross_level>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_cross_level.md)
  : Tables of a cross-level fit
- [`get_results(`*`<multilpa_lta>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_lta.md)
  [`as.data.frame(`*`<multilpa_lta>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_lta.md)
  [`print(`*`<multilpa_lta>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_lta.md)
  [`summary(`*`<multilpa_lta>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_lta.md)
  [`coef(`*`<multilpa_lta>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_lta.md)
  [`vcov(`*`<multilpa_lta>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_lta.md)
  [`logLik(`*`<multilpa_lta>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_lta.md)
  [`nobs(`*`<multilpa_lta>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_lta.md)
  : Tables of a latent transition fit with occasion- or
  covariate-dependent transitions or occasion-specific measurement
- [`print(`*`<latents_growth_mixture>`*`)`](https://pak.dynasite.org/latents/reference/print.latents_growth_mixture.md)
  : Print a growth mixture model
- [`print(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/print.latents_mixture_regression.md)
  : Print a mixture-of-regressions fit
- [`print(`*`<latents_plots>`*`)`](https://pak.dynasite.org/latents/reference/print.latents_plots.md)
  : Print several plots
- [`print(`*`<latents_pooled>`*`)`](https://pak.dynasite.org/latents/reference/print.latents_pooled.md)
  : Print a pooled multiply imputed fit
- [`print(`*`<latents_table>`*`)`](https://pak.dynasite.org/latents/reference/print.latents_table.md)
  [`as.data.frame(`*`<latents_table>`*`)`](https://pak.dynasite.org/latents/reference/print.latents_table.md)
  : Print a latents result table
- [`print(`*`<multilpa>`*`)`](https://pak.dynasite.org/latents/reference/print.multilpa.md)
  : Print a fitted multilevel latent profile model
- [`print(`*`<multilpa_additive>`*`)`](https://pak.dynasite.org/latents/reference/print.multilpa_additive.md)
  : Print a group-class fit
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
- [`prior_control()`](https://pak.dynasite.org/latents/reference/prior_control.md)
  [`print(`*`<latents_prior>`*`)`](https://pak.dynasite.org/latents/reference/prior_control.md)
  : Conjugate prior for Gaussian mixture estimation
- [`summary(`*`<latents_growth_mixture>`*`)`](https://pak.dynasite.org/latents/reference/summary.latents_growth_mixture.md)
  [`print(`*`<summary_latents_growth_mixture>`*`)`](https://pak.dynasite.org/latents/reference/summary.latents_growth_mixture.md)
  : Summarize a growth mixture model
- [`summary(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/summary.latents_mixture_regression.md)
  [`print(`*`<summary_latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/summary.latents_mixture_regression.md)
  : Summarize a mixture-of-regressions fit
- [`summary(`*`<multilpa_additive>`*`)`](https://pak.dynasite.org/latents/reference/summary.multilpa_additive.md)
  [`print(`*`<summary_multilpa_additive>`*`)`](https://pak.dynasite.org/latents/reference/summary.multilpa_additive.md)
  [`as.data.frame(`*`<summary_multilpa_additive>`*`)`](https://pak.dynasite.org/latents/reference/summary.multilpa_additive.md)
  : Summarize an additive group-class fit
- [`summary(`*`<latents_pooled>`*`)`](https://pak.dynasite.org/latents/reference/summary.latents_pooled.md)
  : Summarize a pooled multiply imputed fit
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
- [`get_results(`*`<latents_growth_mixture>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_growth_mixture.md)
  [`as.data.frame(`*`<latents_growth_mixture>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_growth_mixture.md)
  [`coef(`*`<latents_growth_mixture>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_growth_mixture.md)
  [`vcov(`*`<latents_growth_mixture>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_growth_mixture.md)
  [`confint(`*`<latents_growth_mixture>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_growth_mixture.md)
  [`logLik(`*`<latents_growth_mixture>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_growth_mixture.md)
  [`nobs(`*`<latents_growth_mixture>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_growth_mixture.md)
  : Tidy results of a growth mixture model
- [`get_results(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_mixture_regression.md)
  [`as.data.frame(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_mixture_regression.md)
  [`coef(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_mixture_regression.md)
  [`vcov(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_mixture_regression.md)
  [`confint(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_mixture_regression.md)
  [`logLik(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_mixture_regression.md)
  [`nobs(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_mixture_regression.md)
  : Tables of a mixture-of-regressions fit
- [`get_results(`*`<latents_pooled>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_pooled.md)
  [`as.data.frame(`*`<latents_pooled>`*`)`](https://pak.dynasite.org/latents/reference/get_results.latents_pooled.md)
  : Tables of a pooled multiply imputed fit
- [`get_results(`*`<multilpa_additive>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_additive.md)
  [`as.data.frame(`*`<multilpa_additive>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_additive.md)
  [`coef(`*`<multilpa_additive>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_additive.md)
  [`vcov(`*`<multilpa_additive>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_additive.md)
  [`confint(`*`<multilpa_additive>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_additive.md)
  [`logLik(`*`<multilpa_additive>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_additive.md)
  [`nobs(`*`<multilpa_additive>`*`)`](https://pak.dynasite.org/latents/reference/get_results.multilpa_additive.md)
  : Tables of an additive group-class fit
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
- [`parameter_inference(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/parameter_inference.multilpa_transitions.md)
  [`vcov(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/parameter_inference.multilpa_transitions.md)
  [`confint(`*`<multilpa_transitions>`*`)`](https://pak.dynasite.org/latents/reference/parameter_inference.multilpa_transitions.md)
  : Wald inference for a latent transition model
- [`vcov(`*`<latents_pooled>`*`)`](https://pak.dynasite.org/latents/reference/vcov.latents_pooled.md)
  : Pooled covariance of a multiply imputed fit
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
- [`predict(`*`<latents_growth_mixture>`*`)`](https://pak.dynasite.org/latents/reference/predict.latents_growth_mixture.md)
  : Prediction for a growth mixture model (not implemented)
- [`predict(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/predict.latents_mixture_regression.md)
  : Predict from a mixture-of-regressions fit
- [`predict(`*`<multilpa>`*`)`](https://pak.dynasite.org/latents/reference/predict.multilpa.md)
  : Classify new observations with a fitted model

## Conditions

- [`latents-conditions`](https://pak.dynasite.org/latents/reference/latents-conditions.md)
  : Conditions this package raises
