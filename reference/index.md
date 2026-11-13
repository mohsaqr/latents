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
- [`overview()`](https://pak.dynasite.org/latents/reference/overview.md)
  : Everything about any object, printed and drawn
- [`get_checks()`](https://pak.dynasite.org/latents/reference/get_checks.md)
  : Every diagnostic check of a fit
- [`get_entropy()`](https://pak.dynasite.org/latents/reference/get_entropy.md)
  : Relative entropy of the classification
- [`get_classification()`](https://pak.dynasite.org/latents/reference/get_classification.md)
  : Class sizes and classification certainty
- [`get_avepp()`](https://pak.dynasite.org/latents/reference/get_avepp.md)
  : Average posterior probability matrix
- [`get_residuals()`](https://pak.dynasite.org/latents/reference/get_residuals.md)
  : Bivariate residuals of the indicators
- [`get_replication()`](https://pak.dynasite.org/latents/reference/get_replication.md)
  : How many starts reached the best log likelihood

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
- [`get_candidates()`](https://pak.dynasite.org/latents/reference/get_candidates.md)
  : Every candidate of a grid
- [`get_best()`](https://pak.dynasite.org/latents/reference/get_best.md)
  : The best candidate under each criterion
- [`get_deltas()`](https://pak.dynasite.org/latents/reference/get_deltas.md)
  : Change in the criterion with each added class
- [`get_knees()`](https://pak.dynasite.org/latents/reference/get_knees.md)
  : Where each criterion line bends
- [`get_tree()`](https://pak.dynasite.org/latents/reference/get_tree.md)
  : How the classes of each candidate split into the next
- [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
  : Parametric bootstrap likelihood-ratio test for nested mixture models
- [`sensitivity()`](https://pak.dynasite.org/latents/reference/sensitivity.md)
  : Assess sensitivity of a fitted mixture model to initialisation
- [`align_labels()`](https://pak.dynasite.org/latents/reference/align_labels.md)
  : Number one fit's profiles like another's
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

## Plots of a fit

One function per picture; each takes any fit that has the view.

- [`plot_profiles()`](https://pak.dynasite.org/latents/reference/plot_profiles.md)
  : Profile means across the indicators
- [`plot_bars()`](https://pak.dynasite.org/latents/reference/plot_bars.md)
  : Profile means as bars
- [`plot_heatmap()`](https://pak.dynasite.org/latents/reference/plot_heatmap.md)
  : Profile means, or a criterion, as a heat map
- [`plot_raincloud()`](https://pak.dynasite.org/latents/reference/plot_raincloud.md)
  : Distribution of each indicator by profile
- [`plot_parallel()`](https://pak.dynasite.org/latents/reference/plot_parallel.md)
  : Every case as a line across the indicators
- [`plot_pairs()`](https://pak.dynasite.org/latents/reference/plot_pairs.md)
  : Scatter-plot matrix of the indicators by profile
- [`plot_responses()`](https://pak.dynasite.org/latents/reference/plot_responses.md)
  : Response probabilities of the categorical indicators
- [`plot_prevalence()`](https://pak.dynasite.org/latents/reference/plot_prevalence.md)
  : Profile prevalence within each group class
- [`plot_sequences()`](https://pak.dynasite.org/latents/reference/plot_sequences.md)
  : Assigned profiles over the occasions
- [`plot_transitions()`](https://pak.dynasite.org/latents/reference/plot_transitions.md)
  : Transition probabilities between profiles
- [`plot_trajectories()`](https://pak.dynasite.org/latents/reference/plot_trajectories.md)
  : Class trajectories of the outcome
- [`plot_individuals()`](https://pak.dynasite.org/latents/reference/plot_individuals.md)
  : Trajectories of individual persons
- [`plot_coefficients()`](https://pak.dynasite.org/latents/reference/plot_coefficients.md)
  : Class coefficients with intervals
- [`plot_classification()`](https://pak.dynasite.org/latents/reference/plot_classification.md)
  : Class probabilities of every person
- [`plot_all()`](https://pak.dynasite.org/latents/reference/plot_all.md)
  : Every plot of an object

## Classification plots

The four classification views, for a fit of any family or a grid.

- [`plot_sizes()`](https://pak.dynasite.org/latents/reference/plot_sizes.md)
  : Class sizes
- [`plot_posteriors()`](https://pak.dynasite.org/latents/reference/plot_posteriors.md)
  : Certainty of each unit's classification
- [`plot_entropy()`](https://pak.dynasite.org/latents/reference/plot_entropy.md)
  : Entropy of the classification
- [`plot_avepp()`](https://pak.dynasite.org/latents/reference/plot_avepp.md)
  : Average posterior probabilities
- [`plot_diagnostics()`](https://pak.dynasite.org/latents/reference/plot_diagnostics.md)
  : Classification summaries of a fit

## Plots of a candidate grid

- [`plot_criteria()`](https://pak.dynasite.org/latents/reference/plot_criteria.md)
  : Every information criterion
- [`plot_enumeration()`](https://pak.dynasite.org/latents/reference/plot_enumeration.md)
  : Information criteria across a candidate grid
- [`plot_selection()`](https://pak.dynasite.org/latents/reference/plot_selection.md)
  : Model-selection evidence across a candidate grid
- [`plot_deltas()`](https://pak.dynasite.org/latents/reference/plot_deltas.md)
  : Change in the criterion with each added class
- [`plot_trust()`](https://pak.dynasite.org/latents/reference/plot_trust.md)
  : Convergence, boundary and replication of every candidate
- [`plot_tree()`](https://pak.dynasite.org/latents/reference/plot_tree.md)
  : How classes split as more are added
- [`plot_views()`](https://pak.dynasite.org/latents/reference/plot_views.md)
  : The plot views of a fitted model, or of one grid

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
