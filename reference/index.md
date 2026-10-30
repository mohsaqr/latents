# Package index

## Overview

- [`latents-package`](https://pak.dynasite.org/latents/reference/latents-package.md)
  [`latents`](https://pak.dynasite.org/latents/reference/latents-package.md)
  : Latent Profile, Class and Transition Models

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
  : Estimate group classes with fixed individual profile definitions
- [`starting_values()`](https://pak.dynasite.org/latents/reference/starting_values.md)
  : Extract starting values for a subsequent fit
- [`prior_control()`](https://pak.dynasite.org/latents/reference/prior_control.md)
  : Specify prior settings for Gaussian profile estimation

## Results

Every table of a fitted model, retrieved by name.

- [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
  : Extract model estimates and analysis results
- [`descriptives()`](https://pak.dynasite.org/latents/reference/descriptives.md)
  : Descriptive statistics and intraclass correlations of model
  variables
- [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  : Diagnostics of a fitted model
- [`report()`](https://pak.dynasite.org/latents/reference/report.md) :
  Display model estimates and diagnostic summaries

## Evaluating a model

- [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md)
  : Compare latent profile models
- [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md)
  : Compare latent class models
- [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  : Class enumeration and comparison of mixture models
- [`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md)
  [`candidate_fit()`](https://pak.dynasite.org/latents/reference/get_candidate.md)
  : Retrieve a fitted model from an enumeration or comparison
- [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
  : Parametric bootstrap likelihood-ratio test for nested mixture models
- [`sensitivity()`](https://pak.dynasite.org/latents/reference/sensitivity.md)
  : Assess sensitivity of a fitted mixture model to initialisation
- [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  : Standard errors, confidence intervals and tests for model parameters
- [`pool_imputations()`](https://pak.dynasite.org/latents/reference/pool_imputations.md)
  : Fit and pool mixture models across multiple imputations

## Mixture regression and growth models

Regressions whose coefficients differ across latent classes, at one or
two levels, for continuous, binary and count outcomes; classes of
trajectories, with random effects within classes (growth mixture
models).

- [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
  : Estimate regression relationships that vary across latent classes
- [`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md)
  : Assess the number of classes in mixture regression models
- [`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)
  : Compare the fit and classification of mixture models
- [`predict(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/latents-predict.md)
  [`simulate(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/latents-predict.md)
  [`predict(`*`<multilpa>`*`)`](https://pak.dynasite.org/latents/reference/latents-predict.md)
  : Predict from and simulate a fitted model

## External variables

Three-step analysis with a correction for classification error.

- [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
  : Compare latent classes on an external outcome
- [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md) :
  Relate latent class membership to external covariates

## Plots

- [`plot_profiles()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot_responses()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot_sequences()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot_prevalence()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot_diagnostics()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot_transitions()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot_tree()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot_selection()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot_enumeration()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  [`plot_views()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
  : Plot model estimates, classification summaries, and comparisons

## Transition networks

Latent transitions handed to the tna package.

- [`get_tna()`](https://pak.dynasite.org/latents/reference/get_tna.md) :
  Represent transitions between latent states as a network
- [`get_group_tna()`](https://pak.dynasite.org/latents/reference/get_group_tna.md)
  : Construct transition networks for latent group classes

## Data

- [`course_engagement`](https://pak.dynasite.org/latents/reference/course_engagement.md)
  : Student engagement across a sequence of courses
- [`student_esm`](https://pak.dynasite.org/latents/reference/student_esm.md)
  : Leisure activities of university students in daily life
- [`study_hours`](https://pak.dynasite.org/latents/reference/study_hours.md)
  : Study hours and quiz scores under two study strategies
- [`growth_scores`](https://pak.dynasite.org/latents/reference/growth_scores.md)
  : Simulated achievement scores with three growth trajectories
- [`growth_schools`](https://pak.dynasite.org/latents/reference/growth_schools.md)
  : Reading growth of students nested in schools
- [`srl`](https://pak.dynasite.org/latents/reference/srl.md) :
  Self-regulated learning scores generated by a large language model

## Methods

Coefficients, covariance, intervals and likelihood of a fit.

- [`vcov(`*`<multilpa_covariates>`*`)`](https://pak.dynasite.org/latents/reference/latents-model-methods.md)
  [`coef(`*`<multilpa>`*`)`](https://pak.dynasite.org/latents/reference/latents-model-methods.md)
  [`vcov(`*`<multilpa>`*`)`](https://pak.dynasite.org/latents/reference/latents-model-methods.md)
  [`confint(`*`<multilpa>`*`)`](https://pak.dynasite.org/latents/reference/latents-model-methods.md)
  : Extract model coefficients, covariance matrices and likelihood
  summaries

## Conditions

- [`latents-conditions`](https://pak.dynasite.org/latents/reference/latents-conditions.md)
  : Errors, warnings and diagnostic messages
