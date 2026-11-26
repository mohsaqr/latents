# Package index

## Overview

- [`latents-package`](https://pak.dynasite.org/latents/reference/latents-package.md)
  [`latents`](https://pak.dynasite.org/latents/reference/latents-package.md)
  : Latent Profile, Class and Transition Models

## Fitting models

Two-level latent profile, latent class and latent transition models.

- [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md) : Fit a
  latent profile model that groups observations by their pattern on
  continuous indicators
- [`lca()`](https://pak.dynasite.org/latents/reference/lca.md) : Fit a
  latent class model that groups observations by their pattern of
  categorical responses
- [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  : Fit a two-level latent profile model: profiles of observations and
  classes of the groups they are nested in
- [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
  : Fit a two-level latent class model: classes of categorical responses
  and classes of the groups they are nested in
- [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) : Fit a
  latent transition model of how units move between latent states over
  repeated occasions
- [`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
  : Fit a two-level profile model in two stages, holding the profiles
  fixed while the group classes are estimated
- [`starting_values()`](https://pak.dynasite.org/latents/reference/starting_values.md)
  : Get the estimates of a fitted profile or class model as starting
  values for another fit
- [`prior_control()`](https://pak.dynasite.org/latents/reference/prior_control.md)
  : Set a prior that regularises the means and covariances of Gaussian
  profile models

## Results

Every table of a fitted model, retrieved by name.

- [`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
  : Get any table of a fitted model, grid, comparison, test or
  diagnostics result as a data frame
- [`descriptives()`](https://pak.dynasite.org/latents/reference/descriptives.md)
  : Describe each variable of a data frame or a fitted model, with its
  variation within and between groups
- [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
  : Diagnose a fitted model: convergence, replication, class sizes,
  separation and local independence
- [`report()`](https://pak.dynasite.org/latents/reference/report.md) :
  Print the estimates, descriptive statistics and diagnostic checks of a
  fitted model, and draw its plots
- [`overview()`](https://pak.dynasite.org/latents/reference/overview.md)
  : Print and plot everything about any object of the package: its
  summary, its checks and every plot
- [`get_checks()`](https://pak.dynasite.org/latents/reference/get_checks.md)
  : Get every diagnostic check of a fitted model, with its value, rule
  and status
- [`get_entropy()`](https://pak.dynasite.org/latents/reference/get_entropy.md)
  : Get the relative entropy, from 0 to 1, of how sharply a model
  separates its classes
- [`get_classification()`](https://pak.dynasite.org/latents/reference/get_classification.md)
  : Get how many units each class holds and how certain their assignment
  is
- [`get_avepp()`](https://pak.dynasite.org/latents/reference/get_avepp.md)
  : Get the average posterior probability of every class among the units
  assigned to each class
- [`get_residuals()`](https://pak.dynasite.org/latents/reference/get_residuals.md)
  : Get the bivariate residuals: the association between pairs of
  indicators that each class leaves unexplained
- [`get_replication()`](https://pak.dynasite.org/latents/reference/get_replication.md)
  : Get how many random starts converged, failed and reached the best
  log likelihood

## Evaluating a model

- [`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md)
  : Fit latent profile models with different numbers of profiles and
  covariance structures, and compare them
- [`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md)
  : Fit latent class models with different numbers of classes, and
  compare them
- [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
  : Fit models with different numbers of classes and collect their
  criteria and diagnostics for comparison
- [`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md)
  [`candidate_fit()`](https://pak.dynasite.org/latents/reference/get_candidate.md)
  : Get one fitted model out of an enumeration grid or a model
  comparison
- [`get_candidates()`](https://pak.dynasite.org/latents/reference/get_candidates.md)
  : Get the table of every candidate model in an enumeration grid, with
  its criteria and diagnostics
- [`get_best()`](https://pak.dynasite.org/latents/reference/get_best.md)
  : Get the converged candidate with the lowest value of each
  information criterion
- [`get_deltas()`](https://pak.dynasite.org/latents/reference/get_deltas.md)
  : Get how much each information criterion improves as a class or a
  group class is added
- [`get_knees()`](https://pak.dynasite.org/latents/reference/get_knees.md)
  : Get the candidate after which adding classes improves each
  information criterion much less
- [`get_tree()`](https://pak.dynasite.org/latents/reference/get_tree.md)
  : Get how the classes of each candidate split into the classes of the
  next, as a table
- [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
  : Test whether a model with one more class improves on a nested model,
  with a parametric bootstrap likelihood ratio
- [`sensitivity()`](https://pak.dynasite.org/latents/reference/sensitivity.md)
  : Refit a model from several random seeds and check whether it reaches
  the same solution
- [`align_labels()`](https://pak.dynasite.org/latents/reference/align_labels.md)
  : Renumber one fit's profiles and group classes to match the most
  similar ones of another fit
- [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  : Compute standard errors, confidence intervals and tests for the
  estimates of a fitted model
- [`pool_imputations()`](https://pak.dynasite.org/latents/reference/pool_imputations.md)
  : Fit a profile or class model to each imputed data set and pool the
  results across the imputations

## Mixture regression and growth models

Regressions whose coefficients differ across latent classes, at one or
two levels, for continuous, binary and count outcomes; classes of
trajectories, with random effects within classes (growth mixture
models).

- [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
  : Fit a mixture regression model in which the relationship between an
  outcome and its predictors differs across latent classes
- [`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md)
  : Fit mixture regression models with different numbers of classes, and
  compare them
- [`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)
  : Compare the fit, classification and convergence of several fitted
  models in one table
- [`predict(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/latents-predict.md)
  [`simulate(`*`<latents_mixture_regression>`*`)`](https://pak.dynasite.org/latents/reference/latents-predict.md)
  [`predict(`*`<multilpa>`*`)`](https://pak.dynasite.org/latents/reference/latents-predict.md)
  : Predict class membership or outcomes from a fitted model, and
  simulate data from it

## External variables

Three-step analysis with a correction for classification error.

- [`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
  : Compare latent classes on an outcome not used to define them,
  correcting for classification error
- [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md) :
  Regress latent class membership on covariates not used to define the
  classes, correcting for classification error

## Plots of a fit

One function per picture; each takes any fit that has the view.

- [`plot_profiles()`](https://pak.dynasite.org/latents/reference/plot_profiles.md)
  : Plot how the profiles differ on the indicators: each profile's mean,
  or another summary, on every indicator
- [`plot_bars()`](https://pak.dynasite.org/latents/reference/plot_bars.md)
  : Plot the mean of each continuous indicator in each profile as bars
- [`plot_heatmap()`](https://pak.dynasite.org/latents/reference/plot_heatmap.md)
  : Plot profile means as tiles of profiles by indicators, or a grid's
  criterion as tiles of candidates
- [`plot_raincloud()`](https://pak.dynasite.org/latents/reference/plot_raincloud.md)
  : Plot the distribution of each continuous indicator in each profile
  as a density, its quartiles and the points
- [`plot_parallel()`](https://pak.dynasite.org/latents/reference/plot_parallel.md)
  : Plot every case as a line across the continuous indicators, one
  panel per profile
- [`plot_pairs()`](https://pak.dynasite.org/latents/reference/plot_pairs.md)
  : Plot every pair of continuous indicators against each other, with
  the cases and an ellipse for each profile
- [`plot_responses()`](https://pak.dynasite.org/latents/reference/plot_responses.md)
  : Plot each profile's probability of one response category on every
  categorical or ordinal indicator
- [`plot_prevalence()`](https://pak.dynasite.org/latents/reference/plot_prevalence.md)
  : Plot how common each profile is within each group class of a
  two-level model
- [`plot_sequences()`](https://pak.dynasite.org/latents/reference/plot_sequences.md)
  : Plot each unit's sequence of assigned profiles across the occasions
  as a row of coloured cells
- [`plot_transitions()`](https://pak.dynasite.org/latents/reference/plot_transitions.md)
  : Plot the estimated probability of moving from each profile to each
  profile between occasions
- [`plot_trajectories()`](https://pak.dynasite.org/latents/reference/plot_trajectories.md)
  : Plot the curve of the outcome in each class of a mixture regression
  model
- [`plot_individuals()`](https://pak.dynasite.org/latents/reference/plot_individuals.md)
  : Plot individual persons' observations and predicted curves beside
  the trajectory of their class
- [`plot_coefficients()`](https://pak.dynasite.org/latents/reference/plot_coefficients.md)
  : Plot the regression coefficients of each class of a mixture
  regression, with 95% confidence intervals
- [`plot_classification()`](https://pak.dynasite.org/latents/reference/plot_classification.md)
  : Plot each person's posterior probability of every class in a growth
  model as a stacked bar
- [`plot_all()`](https://pak.dynasite.org/latents/reference/plot_all.md)
  : Plot every view an object has and return the plots as one list

## Classification plots

The four classification views, for a fit of any family or a grid.

- [`plot_sizes()`](https://pak.dynasite.org/latents/reference/plot_sizes.md)
  : Plot how many units each class holds, or the size of the smallest
  class of every candidate in a grid
- [`plot_posteriors()`](https://pak.dynasite.org/latents/reference/plot_posteriors.md)
  : Plot the posterior probability of each unit's assigned class, by
  class
- [`plot_entropy()`](https://pak.dynasite.org/latents/reference/plot_entropy.md)
  : Plot how uncertain the classification is, by unit for a fitted model
  or by candidate for a grid
- [`plot_avepp()`](https://pak.dynasite.org/latents/reference/plot_avepp.md)
  : Plot the average posterior probability matrix of a fit, or the
  lowest average posterior probability of every candidate
- [`plot_diagnostics()`](https://pak.dynasite.org/latents/reference/plot_diagnostics.md)
  : Plot how well a fitted model separates its classes: sizes,
  posteriors, entropy and average posterior probabilities

## Plots of a candidate grid

- [`plot_criteria()`](https://pak.dynasite.org/latents/reference/plot_criteria.md)
  : Plot every information criterion of a grid, or of a model
  comparison, against the number of classes
- [`plot_enumeration()`](https://pak.dynasite.org/latents/reference/plot_enumeration.md)
  : Plot the information criteria of an enumeration grid against the
  number of classes
- [`plot_selection()`](https://pak.dynasite.org/latents/reference/plot_selection.md)
  : Plot the evidence for choosing the number of classes from an
  enumeration grid
- [`plot_deltas()`](https://pak.dynasite.org/latents/reference/plot_deltas.md)
  : Plot how much the information criterion falls each time a class or a
  group class is added
- [`plot_trust()`](https://pak.dynasite.org/latents/reference/plot_trust.md)
  : Plot whether each candidate of a grid converged, lies on a boundary
  and replicated its best log likelihood
- [`plot_tree()`](https://pak.dynasite.org/latents/reference/plot_tree.md)
  : Plot how the classes of each solution split or rearrange as one more
  class is added
- [`plot_views()`](https://pak.dynasite.org/latents/reference/plot_views.md)
  : List the views a plot call can draw for a fitted model or a grid,
  with what each one shows

## Transition networks

Latent transitions handed to the tna package.

- [`get_tna()`](https://pak.dynasite.org/latents/reference/get_tna.md) :
  Get the transitions between a fitted model's profiles as one
  transition network for the whole sample
- [`get_group_tna()`](https://pak.dynasite.org/latents/reference/get_group_tna.md)
  : Get one transition network per latent group class, to compare how
  their members move between profiles

## Data

- [`course_engagement`](https://pak.dynasite.org/latents/reference/course_engagement.md)
  : Student engagement across a sequence of courses
- [`engagement_momentum`](https://pak.dynasite.org/latents/reference/engagement_momentum.md)
  : Simulated engagement states and trajectories over a study programme
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
  : Extract the coefficients, their covariance matrix and the log
  likelihood of a fitted model

## Conditions

- [`latents-conditions`](https://pak.dynasite.org/latents/reference/latents-conditions.md)
  : Errors, warnings and diagnostic messages
