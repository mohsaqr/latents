#' Get any table of a fitted model, grid, comparison, test or diagnostics result
#' as a data frame
#'
#' @description
#' `get_results()` returns the results of an object of this package as data
#' frames. Every fitted model, enumeration grid, model comparison, bootstrap
#' test, starting-value set and diagnostics result holds a set of named
#' tables, and `what` chooses among them. The sections below list every
#' table name for each kind of object, with what one row is, its columns and
#' the arguments it takes. Many tables also have a verb of their own that
#' takes any object, such as [get_entropy()], [get_best()] or [get_tree()].
#'
#' @param x A fitted model, an enumeration grid, a model comparison, a
#'   bootstrap test, a starting-value set, a diagnostics result or pooled
#'   imputations; or the `summary()` of a fitted model, of a profile or
#'   class grid or of a bootstrap test, which holds the same tables.
#' @param what Which table to return:
#'   * `NULL` (default): the object's primary table, the same table as
#'     `as.data.frame(x)` and the first one listed for the object below. A
#'     profile or class fit without continuous indicators returns the first
#'     of `"responses"`, `"ordinal"` and `"count_means"` that it has; a
#'     general transition fit (class `multilpa_lta`) returns `"transitions"`;
#'   * one table name, matched in full and not by its first letters;
#'   * `"all"`: every table the object builds by default, as a named list.
#' @param ... Arguments of the requested table, named in its entry below
#'   and described there or under *Table arguments*.
#'
#' @return For one table, a data frame: the table's entry below gives what
#'   one row is and the columns. For `what = "all"`, a named list of data
#'   frames in the order of the entries, each built with its default
#'   arguments. The list leaves out the tables this object cannot produce
#'   (the sequence tables of a fit made without `time`, for example) and
#'   those built only when asked for (`"lmr"`, `"recovery"`, and
#'   `"group_posteriors"` of a fit with fewer than two group classes).
#'
#'   Conditions raised:
#'   * `latents_bad_argument`: `what` is not a table of the object (the
#'     message lists the names it has), `...` holds an argument the table
#'     does not take, or `x` is not an object of this package;
#'   * `latents_unused_argument` (a warning): such an argument, given to a
#'     table of a cross-level fit or of a transition or family grid;
#'   * `latents_no_time`: a sequence table of a fit made without `time`;
#'   * `latents_no_local_dependence`, `latents_no_direct_effects`: the
#'     table of those terms, for a fit without them;
#'   * `latents_unsupported_noise`: the classification errors, BCH weights
#'     or residuals of a fit with a noise component.
#'
#' @section Tables of a profile, class or transition fit:
#' The tables of a fit from [lpa()], [lca()], [multilpa()] or [multilca()],
#' with or without membership covariates, and of an [lta()] fit of class
#' `multilpa_transitions` (the basic transition model). Any other [lta()]
#' fit has class `multilpa_lta` and the tables listed under *Tables of the
#' other model families*.
#' \describe{
#'   \item{`"profiles"`}{The Gaussian measurement model, one row per profile
#'     and continuous indicator: `profile`, `indicator`, `mean`, `variance`,
#'     `standard_deviation`, `mean_standard_error`,
#'     `variance_standard_error`. Takes `scale` and `data`.}
#'   \item{`"responses"`}{The categorical measurement model, one row per
#'     profile, indicator and category: `profile`, `indicator`, `category`,
#'     `probability`, `threshold` (`qlogis(P(y <= category))`, `NA` for the
#'     last category), `probability_standard_error`. Takes `data`.}
#'   \item{`"covariances"`}{One row per profile and ordered pair of continuous
#'     indicators: `profile`, `indicator`, `indicator_2`, `covariance` (exact
#'     zeros off the diagonal under a diagonal structure).}
#'   \item{`"ordinal"`, `"count_means"`}{Only for fits with ordinal or
#'     count indicators. `"ordinal"` has one row per profile, indicator and
#'     category: `profile`, `indicator`, `category`, `probability`,
#'     `location` (the profile's location on the indicator's latent scale)
#'     and `location_standard_error`. `"count_means"` has one row per
#'     profile and indicator: `profile`, `indicator`, `mean`,
#'     `standard_error`. Take `data`.}
#'   \item{`"local_dependence"`, `"direct_effects"`}{For fits with those
#'     arguments. Local dependence has `pair`, `kind`, `profile`, `item_1`,
#'     `category_1`, `item_2`, `category_2`, `probability`, `independent`,
#'     `ratio`, `covariance`, `correlation`: a categorical pair has one row per
#'     profile and category combination with the joint `probability`, the
#'     product of the marginals (`independent`) and their `ratio`; a continuous
#'     pair one row per profile with `covariance` and `correlation`. Direct
#'     effects have one row per effect: `indicator`, `covariate`, `estimate`,
#'     `standard_error`, `statistic`, `p_value`, `conf_low`, `conf_high`.}
#'   \item{`"counts"`}{Effective class sizes (summed posteriors), one row
#'     per level and class: `level`, `class`, `effective_count`,
#'     `effective_proportion`.}
#'   \item{`"model"`, `"information_criteria"`}{One row describing the
#'     fit: dimensions, likelihood, parameter counts, criteria, convergence
#'     and boundary status; `"information_criteria"` has the same columns
#'     for every family and takes `format` and `definitions` (see
#'     *Information criteria*).}
#'   \item{`"posteriors"`}{One row per observation and profile: `row`,
#'     `group`, `profile`, `posterior`, `modal`. Takes `format`. Excluded
#'     observations keep their rows with `NA`.}
#'   \item{`"group_posteriors"`}{One row per group and group class: `group`,
#'     `group_size`, `log_likelihood`, `group_class`, `posterior`, `modal`.
#'     Takes `format`.}
#'   \item{`"assignments"`}{One row per observation: the columns the fit
#'     was built from (or those of `data`), `profile`, `group_class` where
#'     there is one, `uncertainty` (one minus the posterior of the assigned
#'     profile), one `posterior_profile_*` column per profile, and
#'     `posterior_noise` for a noise fit (profile 0). Takes `data` and
#'     `truth`; with `truth` it returns a recovery table (*Table arguments*).
#'     Rows with no observed indicator have `NA` assignments under
#'     `drop_empty = TRUE`, but keep a group class supported by other rows.}
#'   \item{`"classification"`}{One row per level and class: `level`, `class`,
#'     `n_modal`, `proportion_modal`, `estimated_n`, `estimated_proportion`,
#'     `average_posterior`, `odds_correct_classification` (see *Classification
#'     quality*). Takes `level`.}
#'   \item{`"average_posteriors"`}{For the units assigned to each class,
#'     their mean posterior for every class: `level`, `assigned_class`,
#'     `class`, `n_assigned`, `average_posterior`; rows sum to one within
#'     an assigned class. Takes `level`.}
#'   \item{`"classification_errors"`, `"bch_weights"`}{The probability
#'     that a unit of one true class is assigned to another (`level`,
#'     `true_class`, `assigned_class`, `probability`), and the BCH weights
#'     that invert that matrix (`level`, `unit`, `assigned_class`, `class`,
#'     `weight`). Take `level`. Marginal over the units on a covariate fit.}
#'   \item{`"entropy"`}{One row per level: `level`, `n_classes`,
#'     `n_units`, `entropy_sum`, `relative_entropy`.}
#'   \item{`"replication"`}{One row: the numbers of starts run, screened,
#'     converged and failed, and `n_best_replicated`, the starts within
#'     `tolerance` of the best log likelihood. See [get_replication()].}
#'   \item{`"starts"`}{One row per start: `start`, `log_likelihood`,
#'     `converged`, `iterations`, `error`, `boundary` where recorded, and
#'     under `search = "emEM"` `kind`, `stage` (`"final"` or `"screened"`)
#'     and `screened_log_likelihood`.}
#'   \item{`"residuals"`}{Bivariate residuals of every indicator pair:
#'     `profile`, `indicator_1`, `indicator_2`, `kind`, `observed`, `expected`,
#'     `residual`, `effective_n`, `statistic`, `df`, `p_value`, `p_adjusted`.
#'     One row per pair and profile, or per pair under `by = "overall"`. Takes
#'     `data`, `by` and `adjust`. See *Bivariate residuals*.}
#'   \item{`"gaussian_diagnostics"`}{Only for fits with continuous
#'     indicators. One row per Gaussian component: `context`, `profile`,
#'     `gaussian_mass`, `positive_rows`, `weight_effective_n`,
#'     `min_variance_ratio`, `min_standardized_eigenvalue`, `boundary`: how
#'     much posterior mass a component carries, how concentrated it is, and
#'     how small its variances are relative to the data. Descriptive, with
#'     no cutoff and no effect on start selection.}
#'   \item{`"sequences"`, `"sequence_summary"`}{Fits with `time`: the
#'     assignments in occasion order (`group`, `group_class`, `time`,
#'     `profile`; takes `format`) and one row per group class
#'     (`group_class`, `groups`, `observations`, `mean_length`,
#'     `median_length`, `shortest`, `longest`, `complete`, `gaps`).}
#'   \item{`"data"`}{The data the model was fitted to, in input order: the
#'     identifier, the occasion, the rescaled weights and the indicators
#'     under their original names; not the covariates, so [three_step()]
#'     and [r3step()] take `data`.}
#'   \item{`"profile_probabilities"`, `"stages"`}{A fit without membership
#'     covariates or transitions: one row per group class and profile
#'     (`group_class`, `profile`, `probability`, `group_class_probability`),
#'     and one row per estimation stage (`stage`, `group_classes`, `fixed`,
#'     `log_likelihood`, `parameters`, `parameters_with_measurement`,
#'     `converged`; two rows for a [fit_staged()] result).}
#'   \item{`"coefficients"`}{A fit with membership covariates
#'     (`profile_covariates` or `group_covariates`): one row per membership
#'     coefficient (`level` `"profile"` or `"group"`, `outcome`, `term`,
#'     `parameter`, always `"logit"`, `estimate`); standard errors come from
#'     [parameter_inference()].}
#'   \item{`"transitions"`, `"initial"`, `"sequence_lengths"`}{An [lta()]
#'     fit: one row per group class and ordered pair of profiles
#'     (`group_class`, `from`, `to`, `probability`, `expected_count`,
#'     `origin_expected_count`, `supported`, `stable`, `estimated`,
#'     `group_class_probability`; takes `estimated` and `stable`); one row
#'     per group class and profile (`group_class`, `profile`, `probability`,
#'     `prevalence`, `group_class_probability`); one row per group (`group`,
#'     `group_class`, `occasions`, `observations`, `complete`).}
#' }
#' Under sampling weights every population quantity, here and in the grid
#' tables, is weighted: effective and estimated counts and proportions,
#' average posteriors, the odds of correct classification, entropy,
#' classification errors, BCH weights, the criteria and the tree's `mass`.
#' Modal counts stay unweighted counts of the units of positive weight; a
#' unit of weight zero keeps its row and posteriors but enters no statistic.
#'
#' @section Tables of an enumeration grid:
#' A profile or class grid from [enumerate_lpa()], [enumerate_lca()] or
#' [enumerate_classes()] has the tables below. [get_candidates()],
#' [get_best()], [get_deltas()], [get_knees()] and [get_tree()] return them
#' too, and their pages describe the columns one by one.
#' \describe{
#'   \item{`"candidates"`}{One row per candidate, failed ones included:
#'     `candidate`, its class counts and `model`, log likelihood, parameter
#'     count, every criterion under both sample-size conventions, both
#'     entropies, the AvePP summaries of each level, convergence, boundary,
#'     replication, the Gaussian diagnostics, `warnings` and `error`. Takes
#'     `cutoff` (default 0.8) for `avepp_below`; a class with no assigned
#'     unit counts as below it. [get_candidates()] drops the columns empty
#'     in every row.}
#'   \item{`"best"`}{One row per criterion and sample-size convention, with
#'     the converged candidate that has the lowest value. `"criteria"` is
#'     its older name, still in `"all"`.}
#'   \item{`"classification"`}{One row per fitted candidate, level and
#'     class: `candidate`, its class counts and `model`, the columns of a
#'     fit's `"classification"` table and the level's `relative_entropy`.}
#'   \item{`"deltas"`, `"knees"`}{Model-selection evidence, taking `using`, one
#'     criterion of the grid (default `"bic_groups"`); neither selects a model.
#'     `"deltas"` has one row per candidate, with the decrease in the criterion
#'     from the converged candidate with one profile fewer (`delta_profile`) and
#'     with one group class fewer (`delta_group_class`): positive when this
#'     candidate is lower, `NA` where the neighbour is missing or unconverged.
#'     `"knees"` has one row per line of the grid (along the profiles, the group
#'     classes or the parameters), with the candidate furthest below the chord
#'     of the line (`knee_candidate`), a simplified Kneedle rule; a line needs
#'     three candidates for a knee and four for `stable`.}
#'   \item{`"tree"`}{The flows of [plot_tree()], one row per pair of classes
#'     in consecutive rows of a tree, with `mass`, the share of all units
#'     the two classes share (summed over the products of their
#'     posteriors). Takes the arguments of [plot_tree()]: `n_group_classes`,
#'     `n_profiles` (`n_classes` on a regression grid), `model`,
#'     `candidates`, `from`, `to` and `level` (default `"individuals"`).}
#'   \item{`"lmr"`}{The Vuong-Lo-Mendell-Rubin and adjusted Lo-Mendell-Rubin
#'     tests, one row per pair of converged candidates differing by one
#'     profile or one group class at the same structure: `direction`,
#'     `null_candidate`, `alternative_candidate`, the alternative's class
#'     counts and `model`, `n_units`, `parameters_added`, `lr`, `vlmr_p`,
#'     `adjusted_lr` (`lr / (1 + 1 / (parameters_added * log(n_units)))`),
#'     `almr_p`, `reference_mean`, `reference_sd` and `note`. The reference
#'     is a weighted chi-square with the eigenvalues of Vuong's matrix as
#'     weights, evaluated by Imhof's method; a pair with a singular
#'     information matrix or a boundary candidate keeps its row with `NA`
#'     and the reason in `note`. Built only when asked for.}
#' }
#' A transition grid (`enumerate_classes(time = )`) or a family grid
#' (`enumerate_classes(family = )`) has three tables:
#' \describe{
#'   \item{`"candidates"`}{One row per candidate, with `model`, `warnings`,
#'     `error` and the columns [get_candidates()] lists for the grid.}
#'   \item{`"best"`}{One row: the row of `"candidates"` with the lowest
#'     value of one criterion among the converged candidates. Takes
#'     `criterion`: `"bic"` (default), `"aic"` or `"bic_individual"`.}
#'   \item{`"tree"`}{As for a profile or class grid. Not in `"all"`.}
#' }
#' A regression grid from [enumerate_regressions()] has:
#' \describe{
#'   \item{`"fit"`}{One row per candidate, with the columns
#'     [get_candidates()] lists for a regression grid and `note`; `best_bic`
#'     is `TRUE` for the converged candidate with the lowest BIC.}
#'   \item{`"tree"`}{As for a profile or class grid, with `n_classes` in
#'     place of `n_profiles`.}
#'   \item{`"model"`}{Not a table: the fitted model of one candidate, chosen
#'     with `n_classes` and `n_group_classes` (default 1). Not in `"all"`.}
#' }
#' The fitted models of the other grids and of a comparison are extracted
#' with [get_candidate()].
#'
#' @section Tables of the other objects:
#' \describe{
#'   \item{A [compare_models()] result}{`"comparison"`: one row per model,
#'     with the columns [compare_models()] lists. `"classification"`: one
#'     row per model, level and class, with the columns of a fit's
#'     `"classification"` table and `relative_entropy`; only when every
#'     model compared is a profile, class, covariate or transition model.}
#'   \item{A [bootstrap_lrt()] result}{`"test"`: one row, the test.
#'     `"replicates"`: one row per simulated data set.}
#'   \item{A starting-value set ([starting_values()])}{`"profiles"`,
#'     `"responses"`, `"profile_probabilities"` and `"covariances"`: the
#'     rows of a fit's tables of those names, with the estimates only (no
#'     standard errors or thresholds).}
#'   \item{A [diagnostics()] result}{`"entropy"`, `"classification"`,
#'     `"average_posteriors"`, `"residuals"` and `"replication"`, each as
#'     the table of that name of a fit, and `"checks"`, one row per check
#'     (columns in [get_checks()]). `"residuals"` exists only for a profile,
#'     class or `multilpa_transitions` fit.}
#' }
#'
#' @section Tables of the other model families:
#' Each family has its own tables; the first is the primary one. Tables
#' with intervals take `level` (the confidence level, default 0.95) and
#' `vcov_type` (`"observed"`, `"robust"` or `"opg"`; `NULL` reuses the type
#' stored with a regression fit).
#'
#' A mixture regression ([mixture_regression()]):
#' \describe{
#'   \item{`"coefficients"`}{One row per class and term: `class`, `term`,
#'     `estimate`, `std_error`, `statistic`, `p_value`, `conf_low`,
#'     `conf_high`, `p_adjusted` (Benjamini-Hochberg over the
#'     non-intercepts) and, where the link makes it an odds or rate ratio,
#'     `exp_estimate` (`NA` for intercepts).}
#'   \item{`"classes"`}{One row per class: `class`, `share`, `count`,
#'     `n_assigned`, `mean_posterior`, `prior`, `prior_std_error`, and for
#'     a Gaussian model the residual standard deviation `sigma` with
#'     `sigma_std_error`.}
#'   \item{`"membership"`}{One row per coefficient of the class
#'     membership model: the logit `estimate` against a `reference` class,
#'     with its standard error, test, interval and `odds_ratio`.}
#'   \item{`"group_classes"`}{Two-level fits (empty otherwise): one row
#'     per group class and class, with the `probability` of the class
#'     within the group class and `group_share`.}
#'   \item{`"fit"`}{One row: the family, nesting, class counts, numbers of
#'     observations and groups, log likelihood, number of parameters, `aic`,
#'     `bic` (on the rows under observation nesting, on the groups under
#'     group or two-level nesting), `bic_rows`, `sabic`, `icl`, `entropy`,
#'     `group_entropy`, `smallest_share`, convergence and replication.}
#'   \item{`"assignments"`}{One row per observation: `row`, the identifier,
#'     `class`, `posterior` and one `probability_class_*` column per class.}
#'   \item{`"groups"`}{Fits with `id` (empty otherwise): one row per group,
#'     with its number of rows, its class or group class and the
#'     posterior probabilities.}
#'   \item{`"fitted"`}{One row per observation: `row`, `observed`,
#'     `fitted`, `fitted_modal`, `residual`, `class`. An ordinal outcome
#'     is scored `1,...,C`.}
#'   \item{`"starts"`}{One row per start: `start`, `kind`,
#'     `log_likelihood`, `converged`, `iterations`, `stage` (`"screened"`
#'     or `"completed"`), `degenerate`, `selected`.}
#'   \item{`"classification"`}{For the units assigned to each class,
#'     their mean posterior of every class: `assigned`, `class`,
#'     `mean_posterior`, `n_assigned`.}
#'   \item{`"recovery"`}{The assigned classes against known ones: one row
#'     per assigned class and truth value, with `n` and `share`. Needs
#'     `data` and `truth` (one column of `data`), and takes
#'     `by = "class"` (default) or `"group_class"`. Not in `"all"`.}
#' }
#' Identifier columns that collide with output names use an `id_` or
#' `cluster_` prefix, with a unique suffix if needed.
#'
#' A growth mixture ([mixture_regression()] with `random`):
#' \describe{
#'   \item{`"coefficients"`}{One row per class and term of the growth
#'     curve, with standard error, test and interval.}
#'   \item{`"classes"`}{One row per class: `share` with its standard error
#'     and interval, `expected_persons`, `assigned_persons`,
#'     `average_posterior`, `residual_sd`.}
#'   \item{`"random"`}{Variances, standard deviations and correlations of
#'     the random effects, one row per class and parameter, with intervals
#'     on the log and Fisher scales.}
#'   \item{`"membership"`}{The coefficients of the class membership model.}
#'   \item{`"fit"`}{One row: class count, persons, observations, the
#'     random-effect structure, number of parameters, log likelihood, `aic`,
#'     `bic`, `bic_rows`, `sabic`, `icl`, `entropy`, `smallest_share`,
#'     convergence and replication. In a multilevel growth mixture `bic` is
#'     on the clusters and `bic_persons` on the persons.}
#'   \item{`"assignments"`}{One row per person: the identifier, `class`,
#'     `probability` and one `posterior_*` column per class.}
#'   \item{`"random_effects"`}{One row per person: the identifier,
#'     `class` and the predicted random effects.}
#'   \item{`"trajectories"`, `"individual"`}{Each class's mean curve with
#'     its band (`class`, the time variable, `estimate`, `std_error`,
#'     `conf_low`, `conf_high`), and each person's curve (the identifier,
#'     the time variable, `class`, `observed`, `class_mean`, `predicted`).
#'     Both take `time`, the name of the numeric time variable.}
#'   \item{`"starts"`}{One row per start, as for a mixture regression.}
#'   \item{`"group_classes"`, `"clusters"`}{Only in a multilevel growth
#'     mixture: its group classes and its clusters.}
#'   \item{`"recovery"`}{As for a mixture regression, without `by`. Not
#'     in `"all"`.}
#' }
#'
#' An additive fit (`multilpa()` with `family = "additive"`, `"dispersion"`
#' or `"additive_dispersion"`):
#' \describe{
#'   \item{`"parameters"`}{Class means, between- and within-group
#'     variances and class weights: `level`, `group_class`, `indicator`,
#'     `parameter`, `estimate`, `standard_error`, `statistic`, `p_value`,
#'     `p_adjusted`, `conf_low`, `conf_high`. Variance intervals are
#'     formed on the log scale and weight intervals on the logit scale.
#'     Takes `adjust` (`"none"` by default).}
#'   \item{`"group_classes"`}{One row per group class: `group_class`,
#'     `weight`, `count`, `n_assigned`, `mean_posterior` and
#'     `effective_groups`, the number of perfectly classified groups that
#'     would estimate the weight as precisely (below 50, intervals warn).}
#'   \item{`"groups"`}{One row per group: the identifier, `n`, `group_class`,
#'     `posterior` and one `probability_group_class_*` column per group class.}
#'   \item{`"intercepts"`, `"conditional_intercepts"`}{The predicted
#'     intercept of each group on each indicator (`intercept_mean`,
#'     `intercept_variance`): one row per group and indicator, and one row
#'     per group, indicator and group class.}
#'   \item{`"assignments"`}{One row per observation, with the group class
#'     of its group.}
#'   \item{`"classification"`}{For the groups assigned to each group
#'     class, their mean posterior of every group class.}
#'   \item{`"fit"`, `"starts"`}{One row describing the fit; one row per start.}
#'   \item{`"recovery"`}{The assigned group classes against known ones.
#'     Needs `data` and `truth`. Not in `"all"`.}
#' }
#'
#' A cross-level fit (`multilpa(family = "restricted_cross_level")` or
#' `"full_cross_level"`):
#' \describe{
#'   \item{`"profiles"`}{The profiles, then the group classes' means of
#'     the group means: `level`, `profile`, `indicator`, `mean`,
#'     `variance`, `prevalence`.}
#'   \item{`"group_classes"`}{One row per group class: `group_class`,
#'     `weight`, `count`, `n_assigned`.}
#'   \item{`"composition"`}{The share of each profile within each group class
#'     (`group_class`, `profile`, `share`): model-implied for the full family, a
#'     posterior-weighted description for the restricted one.}
#'   \item{`"groups"`, `"assignments"`}{One row per group, and one row per
#'     observation, with their classes and posterior probabilities.}
#'   \item{`"fit"`, `"starts"`, `"gaussian_diagnostics"`}{One row describing the
#'     fit; the record of the starts (for the restricted family, one row per
#'     part of the search); one row per Gaussian component.}
#' }
#' The criteria of a cross-level fit come from a working likelihood and do
#' not choose the number of group classes or the form; see [multilpa()].
#'
#' A general transition fit, the [lta()] fit of class `multilpa_lta` that
#' the extensions of the basic transition model listed in [lta()] return:
#' \describe{
#'   \item{`"transitions"`, `"initial"`}{The transition and initial
#'     probabilities the model implies, averaged over the groups and the
#'     occasions they apply to. `"transitions"` has one row per group class,
#'     occasion (`NA` when the transitions do not differ by occasion) and
#'     ordered pair of profiles: `group_class`, `occasion`, `from`, `to`,
#'     `probability`, `origin_expected_count`, `supported`, `estimated`.
#'     `"initial"` has `group_class`, `profile`, `probability`.}
#'   \item{`"transition_coefficients"`, `"initial_coefficients"`}{The same
#'     parameters as log odds, one row per coefficient, with standard error,
#'     test, interval and `odds_ratio`.} \item{`"second_order_transitions"`,
#'     `"second_order_coefficients"`}{Only for `order = 2`.} \item{`"profiles"`,
#'     `"covariances"`, `"responses"`, `"ordinal"`, `"count_means"`}{The
#'     measurement model, each with an `occasion` column: profile means and
#'     variances, covariances (under `covariance_model = "full"`), and the
#'     categorical, ordinal and count models of fits with such indicators.}
#'   \item{`"group_classes"`}{One row per group class: `group_class`, `weight`,
#'     `count`.} \item{`"assignments"`, `"posteriors"`, `"entropy"`,
#'     `"information_criteria"`}{As for the other transition fits, with the same
#'     arguments.}
#'   \item{`"fit"`, `"starts"`, `"gaussian_diagnostics"`}{One row
#'     describing the fit; one row per start; and, for continuous
#'     indicators, one row per Gaussian component.}
#' }
#' `profile` is `"profile_1"`, ... in `"profiles"`, `"responses"` and the
#' transition tables of a general transition fit, and the profile's number
#' in its other tables.
#'
#' Pooled imputations ([pool_imputations()]):
#' \describe{
#'   \item{`"estimates"`}{One row per pooled parameter: `level`, `outcome`,
#'     `term`, `parameter`, `estimate`, `standard_error`, `statistic`,
#'     `p_value`, `p_adjusted`, `conf_low`, `conf_high`, `df`, and the
#'     `within` and `between` variances with `riv` and `fmi`.}
#'   \item{`"imputations"`}{Each imputation's aligned estimates, one row
#'     per imputation and parameter.}
#'   \item{`"fits"`}{One row per imputation: `imputation`, `log_likelihood`,
#'     `n_parameters`, `converged` and the label orders used to align it.}
#' }
#'
#' @section Table arguments:
#' \describe{
#'   \item{`data`}{One row per observation of the fit, in the fitted order;
#'     `NULL` uses the columns the fit carries. The standard errors of
#'     `"profiles"` and `"responses"` come from [parameter_inference()] on
#'     these data; a fit that cannot supply them gives `NA` with a message
#'     when `data` is omitted, and raises the condition when it is given.}
#'   \item{`scale`}{For `"profiles"`: `"raw"` (default) or
#'     `"standardized"`, each indicator's deviation from its grand mean in
#'     units of the observed population standard deviation (a weighted mean
#'     and spread, without zero-weight rows, for a weighted fit).}
#'   \item{`format`}{For the posterior tables and `"sequences"`, `"long"`
#'     (default, one row per unit and class) or `"wide"` (one row per unit).
#'     For `"information_criteria"`, `"wide"` (default, one row) or `"long"`
#'     (one row per criterion and convention).}
#'   \item{`level`}{For the classification tables of a profile, class or
#'     transition fit: `"individuals"`, `"groups"` or `"both"`; `NULL`
#'     (default) is `"both"` for a fit with two or more group classes and
#'     `"individuals"` otherwise. For the other model families `level` is
#'     the confidence level of the intervals (default 0.95).}
#'   \item{`truth`}{For `"assignments"`, column names of `data` holding
#'     known labels. The table then cross-tabulates the model's labels
#'     against them, one row per class and truth value: `assignment` (the
#'     label compared: `group_class` for a truth constant within every
#'     group of a fit with group classes, `profile` otherwise), `class`,
#'     `truth`, `value`, `n` and `proportion` (the share of units with that
#'     truth value assigned to that class). Excluded units appear under
#'     `class = NA`. Labels are arbitrary, so a truth value assigned mostly
#'     to a class with another label is relabelling, not misclassification.}
#'   \item{`by`, `adjust`}{For `"residuals"`: `by` is `"profile"` (default, each
#'     profile separately) or `"overall"` (pooled over the profiles), and
#'     `adjust` is the multiplicity correction across the pairs: `"BH"`
#'     (default), `"none"` or another `stats::p.adjust()` method.}
#'   \item{`definitions`}{For `"information_criteria"` in long format,
#'     `TRUE` adds each criterion's formula and reference; `FALSE`
#'     (default) leaves them out.}
#'   \item{`estimated`, `stable`}{For `"transitions"` of a
#'     `multilpa_transitions` fit; `NULL` (default) keeps every row.
#'     `estimated = TRUE` keeps the rows estimated from observed occupancy
#'     and `FALSE` the others; `stable = TRUE` keeps the diagonal cells
#'     (staying) and `FALSE` the off-diagonal ones (moving).}
#' }
#' `cutoff`, `using`, `criterion` and the tree arguments belong to the grid
#' tables and are described with them.
#'
#' @section Information criteria:
#' With `q` free parameters, `n` the sample size and `EN` the entropy of the
#' classification: `deviance = -2L`, `aic = -2L + 2q`, `bic = -2L + q log(n)`,
#' `sabic = -2L + q log((n + 2) / 24)`, `caic = -2L + q (log(n) + 1)`,
#' `awe = -2(L - EN) + 2q (1.5 + log(n))`, `icl = -2L + q log(n) + 2 EN`,
#' `kic = -2L + 3(q + 1)` (as tidyLPA computes it) and `clc = -2L + 2 EN` (on
#' the entropy sum, so not comparable with tidyLPA's CLC on relative entropy).
#' `EN` is the soft entropy of the complete-data classification, as in ICL, AWE
#' and CLC; in a two-level model it is the group-class entropy plus each row's
#' profile entropy within its group class.
#'
#' Every criterion with a sample size is reported under two conventions,
#' `"groups"` (`n` the number of groups) and `"individuals"` (`n` the number
#' of rows), which coincide in a single-level model; use one consistently
#' across candidates. A transition fit keeps one entropy per level instead,
#' since its profiles form a chain within a person. Rows with no observed
#' indicator, groups made only of such rows, and units of weight zero are
#' in no sample size and no entropy.
#'
#' @section Classification quality:
#' The odds of correct classification of a class is
#' `(p / (1 - p)) / (r / (1 - r))`, with `p` the average posterior in the
#' assigned class and `r` the estimated class proportion: near one, the
#' classification is no better than the class proportions; five is a
#' conventional benchmark, not a test. It is `NA` when `p` is missing or zero or
#' `r` is zero or one (so always for a single class), and `Inf` when `p` is one
#' at an interior `r`. `relative_entropy` is `1 - entropy_sum / (n log K)`, with
#' `entropy_sum` the level's own entropy, between zero and one, higher meaning
#' sharper classification; `NA` for a level with one class. A fit without group
#' classes reports individuals only.
#'
#' @section Bivariate residuals:
#' `"residuals"` checks, pair by pair, that the indicators are independent
#' within a profile, weighting each observation by its posterior. Two Gaussian
#' indicators (`kind = "gaussian"`): `observed` is the posterior-weighted
#' residual correlation, `expected` the correlation the fitted covariance
#' implies, and `statistic` the difference of their Fisher transforms times
#' `sqrt(effective_n - 3)`, referred to the normal distribution. Two categorical
#' or ordinal indicators (`"categorical"`): `statistic` is the Pearson
#' chi-square of the posterior-weighted two-way table against the model's,
#' Latent GOLD's bivariate residual times its `df = (R - 1)(C - 1)`, with
#' `observed` and `expected` the two tables' Cramer's V; the reference ignores
#' that the model's table is estimated from the same data and can be
#' conservative. A mixed pair or a pair with a count indicator (`"mixed"`,
#' `"count"`): `observed` is the posterior-weighted correlation of the values (a
#' categorical indicator entering as its category number), `expected` zero, and
#' no test. `p_adjusted` applies `adjust` across the pairs with a p-value. A
#' weighted fit reports the descriptive columns and no tests.
#'
#' @seealso [as.data.frame()] for the primary table, `summary()` for every
#'   table, [diagnostics()] for the classification tables with their checks, and
#'   the verbs [get_checks()], [get_entropy()], [get_classification()],
#'   [get_avepp()], [get_residuals()], [get_replication()], [get_candidates()],
#'   [get_best()], [get_deltas()], [get_knees()] and [get_tree()].
#' @examples
#' fit <- multilpa(
#'   course_engagement,
#'   vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
#'   id = "student", n_profiles = 2, n_group_classes = 2, n_starts = 4,
#'   seed = 1
#' )
#' get_results(fit)
#' get_results(fit, what = "profile_probabilities")
#' get_results(fit, what = "entropy")
#' names(get_results(fit, what = "all"))
#'
#' # Compare estimated assignments with the generating classes at both levels.
#' get_results(fit, what = "assignments", data = course_engagement,
#'             truth = c("engagement", "student_type"))
#' @export
get_results <- function(x, what = NULL, ...) {
  UseMethod("get_results")
}

#' @noRd
#' @export
get_results.default <-function(x, what = NULL, ...) {
  stop(errorCondition(sprintf(
    "`get_results()` has no method for an object of class %s.",
    paste(sprintf("`%s`", class(x)), collapse = ", ")),
    class = "latents_bad_argument", call = NULL))
}

#' @noRd
#' @export
get_results.multilpa <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @noRd
#' @export
get_results.multilpa_covariates <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @noRd
#' @export
get_results.multilpa_transitions <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}


#' @noRd
#' @export
get_results.multilpa_enumeration <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @noRd
#' @export
get_results.latents_comparison <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @noRd
#' @export
get_results.multilpa_bootstrap_lrt <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @noRd
#' @export
get_results.multilpa_start <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @noRd
#' @export
get_results.multilpa_diagnostics <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @noRd
#' @export
get_results.summary_multilpa <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @noRd
#' @export
get_results.summary_multilpa_covariates <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}


#' @noRd
#' @export
get_results.summary_multilpa_transitions <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @noRd
#' @export
get_results.summary_multilpa_enumeration <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' @noRd
#' @export
get_results.summary_multilpa_bootstrap_lrt <-function(x, what = NULL, ...) {
  .multilpa_dispatch_table(x, what, list(...))
}

#' Coerce an object of this package to its primary table
#'
#' The one body behind every `as.data.frame()` method here, so no class can
#' coerce to something other than the first entry of its own catalogue.
#'
#' @param x The object.
#' @param row.names The caller's `row.names`.
#' @param extra `list(...)` as the method received it, which must be empty.
#' @return A base `data.frame`.
#' @noRd
.multilpa_coerce <- function(x, row.names, extra) {
  .multilpa_reject_extra_arguments(
    extra, "as.data.frame()",
    "It coerces to the primary table; get_results(x, what = ) has the others.")
  result <- get_results(x)
  row.names(result) <- row.names
  result
}

#' Resolve one `what` against an object's catalogue and build that table
#'
#' The single place `what` is matched, arguments are checked and a table is
#' built, so every method is the same three words and no class can drift from
#' the shared contract.
#'
#' @param x The object.
#' @param what The caller's `what`, possibly `NULL`.
#' @param extra `list(...)` as the method received it.
#' @return A `data.frame`, or a named list of them for `what = "all"`.
#' @noRd
.multilpa_dispatch_table <- function(x, what, extra) {
  catalogue <- .multilpa_catalogue(x)
  if (is.null(what)) what <- .multilpa_primary_table(x, names(catalogue))
  what <- .multilpa_match_table(what, names(catalogue))
  if (identical(what, "all")) return(.multilpa_all_tables(x, catalogue, extra))
  entry <- catalogue[[what]]
  .multilpa_check_table_arguments(what, entry$args, extra)
  do.call(entry$build, c(list(x), extra[intersect(names(extra), entry$args)]))
}

#' The primary table of an object: what `get_results(x)` and
#' `as.data.frame(x)` return
#'
#' The first catalogue entry, except for a fitted model without continuous
#' indicators, whose profile-mean table would be empty: its primary table is
#' the first of its categorical, ordinal and count tables that it has.
#'
#' @param x The object.
#' @param choices The catalogue's names.
#' @return A single table name.
#' @noRd
.multilpa_primary_table <- function(x, choices) {
  stopifnot("`choices` must name at least one table" =
              is.character(choices) && length(choices) > 0L)
  first <- choices[[1L]]
  if (!identical(first, "profiles") || .multilpa_any_summary(x) ||
      !.multilpa_any_fit(x) || length(.multilpa_continuous_names(x)) > 0L) {
    return(first)
  }
  carried <- c(responses = length(x$response_probabilities) > 0L,
               ordinal = !is.null(x$extra_data$ordinal),
               count_means = !is.null(x$extra_data$count))
  available <- names(carried)[carried & names(carried) %in% choices]
  if (length(available) == 0L) first else available[[1L]]
}

#' Match a requested table name against what an object offers
#'
#' Exact matching rather than `match.arg()`'s partial matching: the catalogues
#' run to nineteen entries and share prefixes -- `starts`, `stages`,
#' `sequences`, `sequence_summary` -- so a partial match would resolve to
#' whichever happened to come first. The refusal lists every available name, so
#' nothing is lost by being exact.
#'
#' @param what The caller's `what`, possibly `NULL`.
#' @param choices The catalogue's names.
#' @return A single table name, or `"all"`.
#' @noRd
.multilpa_match_table <- function(what, choices) {
  if (is.null(what)) return(choices[1L])
  .latents_require("`what` must be a single table name, or NULL for the first one" =
                     is.character(what) && length(what) == 1L && !is.na(what))
  if (identical(what, "all") || what %in% choices) return(what)
  # An earlier name of a table still finds it, so scripts written against a
  # released version keep working.
  aliases <- c(criteria = "best")
  if (what %in% names(aliases) && aliases[[what]] %in% choices) return(aliases[[what]])
  stop(errorCondition(sprintf(
    "`what = \"%s\"` is not a table of this object. It has %s, and \"all\".",
    what, paste(sprintf("\"%s\"", choices), collapse = ", ")),
    class = "latents_bad_argument", call = NULL))
}

#' Refuse an argument the requested table does not take
#'
#' A dropped argument returns a table that answers a different question from
#' the one asked, and does so without a word. Naming it back costs the caller
#' one correction and saves them a wrong number.
#'
#' @param what The table name, for the message.
#' @param accepted The argument names that table takes.
#' @param extra `list(...)` as the method received it.
#' @return `NULL`, invisibly, when there is nothing to refuse.
#' @noRd
.multilpa_check_table_arguments <- function(what, accepted, extra) {
  if (length(extra) == 0L) return(invisible(NULL))
  supplied <- names(extra) %||% rep("", length(extra))
  unknown <- is.na(match(supplied, accepted)) | !nzchar(supplied)
  if (!any(unknown)) return(invisible(NULL))
  labels <- ifelse(nzchar(supplied[unknown]),
                   sprintf("`%s`", supplied[unknown]), "an unnamed argument")
  takes <- if (length(accepted) == 0L) "takes no further arguments" else
    sprintf("takes %s", paste(sprintf("`%s`", accepted), collapse = ", "))
  stop(errorCondition(sprintf("`what = \"%s\"` does not use %s; it %s.",
                              what, paste(labels, collapse = ", "), takes),
                      class = "latents_bad_argument", call = NULL))
}

#' Classed refusals that mean "this object has no such table"
#'
#' `what = "all"` walks the whole catalogue, and some of its entries cannot be
#' built for every object: a fit made without `time` has no sequences, a fit
#' without discrete group classes has no bivariate residuals. Those are refusals about
#' the object, not failures, and `"all"` leaves the table out. Every other
#' condition propagates, so a genuine defect is never swallowed.
#'
#' @return A character vector of condition classes.
#' @noRd
.multilpa_absent_table_conditions <- function() {
  c("latents_no_time", "latents_no_group_classes",
    "latents_no_indicator_data", "latents_incomplete_fit",
    "latents_no_continuous", "latents_no_categorical",
    # A fit whose classes cannot be separated has no three-step correction to
    # report, which is a fact about the fit and not a reason for `summary()`
    # to fail.
    "latents_inseparable_classes",
    # A noise component has no classification-error correction yet.
    "latents_unsupported_noise",
    # A grid with one solution per model has no tree.
    "latents_no_tree",
    # A fit without locally dependent pairs has no table of them.
    "latents_no_local_dependence",
    # Nor does a fit without direct effects.
    "latents_no_direct_effects")
}

#' Build every table an object can produce
#'
#' @param x The object.
#' @param catalogue Its catalogue.
#' @param extra `list(...)`, which must be empty.
#' @return A named list of `data.frame`s, without the ones this object has not.
#' @noRd
.multilpa_all_tables <- function(x, catalogue, extra) {
  .multilpa_check_table_arguments("all", character(), extra)
  absent <- .multilpa_absent_table_conditions()
  catalogue <- Filter(function(entry) !isFALSE(entry$in_all), catalogue)
  tables <- lapply(catalogue, function(entry) {
    # Errors only: a table that warns is still a table, and the warning is the
    # fit's own, so it is left to reach the caller who asked for every table.
    tryCatch(entry$build(x), error = function(condition) {
      if (inherits(condition, absent)) return(NULL)
      stop(condition)
    })
  })
  Filter(Negate(is.null), tables)
}

#' The tables an object offers, and the arguments each of them takes
#'
#' One definition per table, shared by `what = <name>` and by `what = "all"`,
#' so the two cannot return differently built versions of the same table. The
#' first entry is the object's primary table, which is what `as.data.frame()`
#' and a `NULL` `what` return.
#'
#' A table the family defines but this particular object cannot produce stays
#' in the catalogue and raises its classed condition, so the refusal names the
#' reason. A table another family defines is absent, and the refusal lists what
#' this one has.
#'
#' @param x The object.
#' @return A named list of `list(build =, args =)` entries, in reporting order:
#'   the measurement model, then the mixing structure, then the unit-level
#'   tables, then the diagnostics, then the sequence tables, then provenance.
#' @noRd
.multilpa_catalogue <- function(x) {
  if (inherits(x, "multilpa_enumeration")) {
    return(list(
      candidates = .multilpa_table(
        function(x, cutoff = 0.8) .multilpa_enumeration_grid(x, cutoff),
        "cutoff"),
      best = .multilpa_table(
        function(x) .multilpa_enumeration_minima(.multilpa_enumeration_grid(x))),
      criteria = .multilpa_table(
        function(x) .multilpa_enumeration_minima(.multilpa_enumeration_grid(x))),
      classification = .multilpa_table(.multilpa_enumeration_classification),
      deltas = .multilpa_table(
        function(x, using = "bic_groups") .latents_deltas_table(x, using), "using"),
      knees = .multilpa_table(
        function(x, using = "bic_groups") .latents_knees_table(x, using), "using"),
      tree = .multilpa_table(.latents_tree_table, .latents_tree_arguments()),
      # Every candidate's information matrix: computed when asked for, not by
      # get_results(x, "all") or summary().
      lmr = .multilpa_table(.latents_lmr_table, in_all = FALSE)))
  }
  if (inherits(x, "latents_comparison")) {
    catalogue <- list(comparison = .multilpa_table(function(x) {
      attr(x, "models") <- NULL
      x
    }))
    models <- attr(x, "models")
    if (is.null(models) || all(vapply(models, .latents_comparison_measurement_fit,
                                    logical(1)))) {
      catalogue$classification <- .multilpa_table(.latents_comparison_classification)
    }
    return(catalogue)
  }
  if (inherits(x, "multilpa_bootstrap_lrt")) {
    return(list(
      test = .multilpa_table(.multilpa_bootstrap_test_frame),
      replicates = .multilpa_table(function(x) x$replicates)))
  }
  if (inherits(x, "multilpa_start")) {
    return(list(
      profiles = .multilpa_table(.multilpa_start_profile_frame),
      responses = .multilpa_table(.multilpa_start_response_frame),
      profile_probabilities = .multilpa_table(.multilpa_start_probability_frame),
      covariances = .multilpa_table(.multilpa_start_covariance_frame)))
  }
  if (inherits(x, "multilpa_diagnostics")) {
    # The entropy table stays the default (as in 0.9.21); the checks follow.
    tables <- list(
      entropy = .multilpa_table(function(x) x$entropy),
      classification = .multilpa_table(function(x) x$classification),
      average_posteriors = .multilpa_table(function(x) x$average_posteriors),
      residuals = .multilpa_table(.multilpa_diagnostic_residuals),
      # A diagnostics object saved before 0.97.5 has no stored table.
      replication = .multilpa_table(function(x) {
        x$replication %||% .multilpa_replication_table(x$fit)
      }),
      checks = .multilpa_table(function(x) x$checks))
    # The family tables a fit did not produce (a family without the profile
    # classification tables, or an object saved before the checks) are absent.
    present <- c(checks = !is.null(x$checks), entropy = !is.null(x$entropy),
                 classification = !is.null(x$classification),
                 average_posteriors = !is.null(x$average_posteriors),
                 residuals = !is.null(x$fit) && .multilpa_any_fit(x$fit),
                 replication = TRUE)
    return(tables[present[names(tables)]])
  }
  if (.multilpa_any_summary(x)) return(.multilpa_stored_catalogue(x))
  stopifnot("`x` must be an object of this package" = .multilpa_any_fit(x))
  c(.multilpa_measurement_catalogue(),
    if (!is.null(x$gaussian_diagnostics) && nrow(x$gaussian_diagnostics) > 0L)
      list(gaussian_diagnostics = .multilpa_table(function(x) x$gaussian_diagnostics)),
    .latents_extra_catalogue(x),
    .multilpa_mixing_catalogue(x),
    .multilpa_unit_catalogue(x),
    .multilpa_diagnostic_catalogue(),
    .multilpa_sequence_catalogue(x),
    .multilpa_provenance_catalogue(x))
}

#' One catalogue entry
#' @param build A function of the object and that table's own arguments.
#' @param args The argument names it accepts.
#' @param in_all `FALSE` leaves the table out of `what = "all"` (and so of
#'   `summary()`), for a table too costly to build unasked.
#' @return A list of the three.
#' @noRd
.multilpa_table <- function(build, args = character(), in_all = TRUE) {
  list(build = build, args = args, in_all = in_all)
}

#' Is this the summary of a fitted model of this package?
#' @param x Any object.
#' @return A single logical.
#' @noRd
.multilpa_any_summary <- function(x) {
  inherits(x, c("summary_multilpa", "summary_multilpa_covariates",
                "summary_multilpa_transitions",
                "summary_multilpa_enumeration",
                "summary_multilpa_bootstrap_lrt"))
}

#' A catalogue over the tables a summary already carries
#'
#' A summary is built by asking the fit for every table, so its own catalogue
#' serves them back rather than recomputing them. That is what makes
#' `get_results(summary(fit), what)` and `get_results(fit, what)` the same table
#' rather than two computations that could diverge.
#'
#' @param x A summary object.
#' @return A named list of catalogue entries.
#' @noRd
.multilpa_stored_catalogue <- function(x) {
  tables <- x$tables
  if (!is.list(tables) || length(tables) == 0L) {
    stop(errorCondition(paste(
      "This summary carries no tables, so it was made by an older version of",
      "the package. Take the tables from the fit instead."),
      class = "latents_incomplete_fit", call = NULL))
  }
  lapply(tables, function(table) {
    force(table)
    .multilpa_table(function(x) table)
  })
}

#' The measurement tables every fitted family offers
#' @return A named list of catalogue entries.
#' @noRd
.multilpa_measurement_catalogue <- function() {
  list(
    profiles = .multilpa_table(
      function(x, data = NULL, scale = c("raw", "standardized"))
        .multilpa_profile_frame(x, data = data,
                                scale = .latents_match_arg(scale, c("raw", "standardized"),
                                                           "scale"),
                                standard_errors = TRUE),
      c("data", "scale")),
    responses = .multilpa_table(
      function(x, data = NULL) .multilpa_response_frame(x, data = data,
                                                        standard_errors = TRUE),
      "data"),
    covariances = .multilpa_table(.multilpa_covariance_frame),
    local_dependence = .multilpa_table(.latents_local_dependence_table),
    direct_effects = .multilpa_table(
      function(x, data = NULL) .latents_direct_frame(x, data = data), "data"))
}

#' The mixing-structure tables, which differ by model family
#' @param x A fitted model of this package.
#' @return A named list of catalogue entries.
#' @noRd
.multilpa_mixing_catalogue <- function(x) {
  counts <- list(counts = .multilpa_table(.multilpa_count_frame))
  if (inherits(x, "multilpa_transitions")) {
    return(c(list(
      transitions = .multilpa_table(
        function(x, estimated = NULL, stable = NULL)
          .multilpa_transitions_table(x, estimated = estimated,
                                      stable = stable),
        c("estimated", "stable")),
      initial = .multilpa_table(.multilpa_initial_frame)), counts))
  }
  if (inherits(x, "multilpa_covariates")) {
    return(c(list(coefficients = .multilpa_table(.multilpa_coefficient_frame)),
             counts))
  }
  c(list(profile_probabilities = .multilpa_table(.multilpa_probability_frame)),
    counts)
}

#' The one-row-per-unit tables
#' @param x A fitted model of this package.
#' @return A named list of catalogue entries.
#' @noRd
.multilpa_unit_catalogue <- function(x) {
  posteriors <- list(posteriors = .multilpa_table(
    function(x, format = c("long", "wide"))
      .multilpa_posterior_frame(x, .latents_match_arg(format, c("long", "wide"), "format")),
    "format"))
  assignments <- list(assignments = .multilpa_table(
    function(x, data = NULL, truth = NULL)
      .multilpa_hide_fabricated_id(
        x, .multilpa_assignments(x, data = data, truth = truth)),
    c("data", "truth")))
  # A single-level fit, or one with one group class, classifies nothing at
  # the group level: its group posteriors stay available by name but are left
  # out of "all" and summary() (API-9).
  c(posteriors,
    list(group_posteriors = .multilpa_table(
      function(x, format = c("long", "wide"))
        .multilpa_group_posterior_frame(
          x, .latents_match_arg(format, c("long", "wide"), "format")),
      "format", in_all = .multilpa_has_group_classes(x))),
    assignments)
}

#' The classification and model-fit tables every fitted family offers
#' @return A named list of catalogue entries.
#' @noRd
.multilpa_diagnostic_catalogue <- function() {
  # `level` follows the fit rather than a fixed default: a two-level model is
  # asked about at both levels, and a family with one level has only one to
  # report. This is the rule `diagnostics()` already applies.
  by_level <- function(verb) {
    .multilpa_table(function(x, level = NULL)
      verb(x, level = .multilpa_default_level(x, level)), "level")
  }
  list(
    classification = by_level(.multilpa_classification_table),
    average_posteriors = by_level(.multilpa_average_posteriors),
    classification_errors = by_level(.multilpa_classification_errors),
    bch_weights = by_level(.multilpa_bch_weights),
    entropy = .multilpa_table(.multilpa_entropy_table),
    replication = .multilpa_table(.multilpa_replication_table),
    residuals = .multilpa_table(
      function(x, data = NULL, by = c("profile", "overall"), adjust = "BH")
        .multilpa_bivariate_residuals(
          x, data = data, by = .latents_match_arg(by, c("profile", "overall"), "by"),
          adjust = .latents_match_arg(adjust, .multilpa_p_adjust_methods, "adjust")),
      c("data", "by", "adjust")),
    information_criteria = .multilpa_table(
      function(x, definitions = FALSE, format = c("wide", "long"))
        .multilpa_information_criteria(
          x, definitions = definitions,
          format = .latents_match_arg(format, c("wide", "long"), "format")),
      c("definitions", "format")),
    model = .multilpa_table(.multilpa_fit_frame))
}

#' The occasion-ordered tables
#' @param x A fitted model of this package.
#' @return A named list of catalogue entries.
#' @noRd
.multilpa_sequence_catalogue <- function(x) {
  shared <- list(
    sequences = .multilpa_table(
      function(x, format = c("long", "wide"))
        .multilpa_sequences(x, format = .latents_match_arg(format, c("long", "wide"),
                                                           "format")),
      "format"),
    sequence_summary = .multilpa_table(.multilpa_sequence_summary))
  if (!inherits(x, "multilpa_transitions")) return(shared)
  c(shared,
    list(sequence_lengths = .multilpa_table(.multilpa_sequence_length_frame)))
}

#' The tables that record how the fit was produced
#' @param x A fitted model of this package.
#' @return A named list of catalogue entries.
#' @noRd
.multilpa_provenance_catalogue <- function(x) {
  stages <- if (inherits(x, "multilpa"))
    list(stages = .multilpa_table(.multilpa_stage_frame)) else list()
  c(stages,
    list(starts = .multilpa_table(function(x) x$starts),
         data = .multilpa_table(function(x)
           .multilpa_hide_fabricated_id(x, .multilpa_model_frame(x)))))
}

#' Hide the unit column a single-level fit fabricated for itself
#'
#' `id = NULL` numbers the rows so the two-level machinery has a unit to group
#' by, and every internal rebuild of the fitting data needs that column --
#' `parameter_inference()` checks the frame against it. What the caller never
#' wrote, though, they should not be handed back, so it is dropped here, at the
#' table the caller reads, rather than at the rebuild the package relies on.
#'
#' @param x A fitted model of this package.
#' @param frame A table that may carry the fabricated column.
#' @return `frame`, without that column when the fit fabricated one.
#' @noRd
.multilpa_hide_fabricated_id <- function(x, frame) {
  if (!isTRUE(x$single_level) || !x$id %in% names(frame)) return(frame)
  frame[setdiff(names(frame), x$id)]
}

#' Which classification levels to report when the caller named none
#' @param x A fitted model of this package.
#' @param level The caller's `level`, possibly `NULL`.
#' @return A single level name.
#' @noRd
.multilpa_default_level <- function(x, level) {
  if (!is.null(level)) return(level)
  # The rule `diagnostics()` applies: a single-level fit, or one with one
  # group class, has no group-level classification to report (API-9).
  if (.multilpa_has_group_classes(x)) "both" else "individuals"
}

#' The bivariate residuals a diagnostics object holds, or the reason it has none
#' @param x A `multilpa_diagnostics` object.
#' @return The residual table.
#' @noRd
.multilpa_diagnostic_residuals <- function(x) {
  if (!is.null(x$residuals)) return(x$residuals)
  if (isTRUE(x$fit$noise)) {
    stop(errorCondition(paste(
      "This fit has no bivariate residuals: they do not account for a noise",
      "component yet."),
      class = "latents_unsupported_noise", call = NULL))
  }
  stop(errorCondition(
    "This fit has no bivariate residuals: they need a discrete group-class model.",
    class = "latents_no_group_classes", call = NULL))
}

#' How many occasions each group of a transition fit contributes
#' @param x A fitted `multilpa_transitions` model.
#' @return One row per group.
#' @noRd
.multilpa_sequence_length_frame <- function(x) {
  group_classes <- x$group_classes
  group_classes[!.multilpa_reported_groups(x)] <- NA_integer_
  data.frame(group = x$group_values, group_class = group_classes,
             occasions = unname(x$sequence_lengths),
             observations = unname(tabulate(x$group_index, nbins = x$n_groups)),
             complete = unname(x$sequence_lengths) == x$n_occasions,
             row.names = NULL)
}

#' Print every table an object carries
#'
#' A summary is the whole of what an object holds, so its `print()` shows every
#' table rather than a chosen few. A table longer than `rows` is shown to that
#' depth with its remaining row count and the call that returns it whole: the
#' alternative is either a console dump or a table the reader cannot tell was
#' truncated.
#'
#' @param tables A named list of `data.frame`s, as `get_results(x, "all")` gives.
#' @param rows How many rows of each table to show.
#' @param digits Printed significant digits.
#' @return `NULL`, invisibly. Called for the printing.
#' @noRd
.multilpa_print_tables <- function(tables, rows = 10L, digits = 4L) {
  stopifnot(
    "`rows` must be a single non-negative whole number" =
      is.numeric(rows) && length(rows) == 1L && is.finite(rows) && rows >= 0 &&
      rows == as.integer(rows),
    "the tables must be a named list of data frames" =
      is.list(tables) && !is.null(names(tables)))
  invisible(lapply(names(tables), function(name) {
    frame <- tables[[name]]
    cat(sprintf("\n-- %s %s\n", name,
                strrep("-", max(3L, 64L - nchar(name)))))
    if (nrow(frame) == 0L) {
      cat("   (no rows)\n")
      return(invisible(NULL))
    }
    # `rows = 0` asks for the headings alone; printing a zero-row head would
    # dump the column names as if the table were empty.
    if (rows > 0) print(utils::head(frame, rows), digits = digits, row.names = FALSE)
    if (nrow(frame) > rows) {
      cat(sprintf("   ... %d more rows.  get_results(x, what = \"%s\")\n",
                  nrow(frame) - rows, name))
    }
    invisible(NULL)
  }))
}

#' The closing line of a summary that has printed its tables
#' @param tables The tables that were printed.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_print_table_footer <- function(tables) {
  cat(sprintf(paste0(
    "\n%d tables above, truncated to fit. get_results(x, what = ) returns any\n",
    "of them whole, and get_results(x, what = \"all\") returns every one.\n"),
    length(tables)))
  invisible(NULL)
}

#' Validate the two printing controls every summary print method takes
#' @param digits The caller's `digits`.
#' @param rows The caller's `rows`.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_check_print_arguments <- function(digits, rows) {
  stopifnot(
    "`digits` must be a single number between 1 and 22" =
      is.numeric(digits) && length(digits) == 1L && is.finite(digits) &&
      digits >= 1 && digits <= 22,
    "`rows` must be a single non-negative whole number" =
      is.numeric(rows) && length(rows) == 1L && is.finite(rows) && rows >= 0 &&
      rows == as.integer(rows))
  invisible(NULL)
}

#' The profile means as one row per profile
#'
#' The long measurement table is one row per profile *and* indicator, which is
#' the right shape to compute on and the wrong one to read: three profiles over
#' five indicators is fifteen lines of console for a table that has three rows.
#' This is the shape the means are read in. It is built by reshaping the long
#' table rather than from `x$means`, so the two cannot report different numbers.
#'
#' @param x A fitted model of this package.
#' @return A `data.frame` with one row per profile and one column per
#'   continuous indicator.
#' @noRd
.multilpa_wide_means <- function(x) {
  long <- .multilpa_profile_frame(x)
  wide <- stats::reshape(long[c("profile", "indicator", "mean")],
                         idvar = "profile", timevar = "indicator",
                         direction = "wide")
  names(wide) <- sub("^mean\\.", "", names(wide))
  row.names(wide) <- NULL
  wide
}

#' Response probabilities one row per indicator category, one column per class
#' @param x A fitted model with categorical indicators.
#' @return A data frame: `indicator`, `category`, then `profile_1`, ...
#' @noRd
.multilpa_wide_responses <- function(x) {
  long <- .multilpa_response_frame(x)
  key <- paste(long$indicator, long$category, sep = "\r")
  keys <- unique(key)
  first <- match(keys, key)
  profiles <- sort(unique(long$profile))
  columns <- lapply(profiles, function(profile) {
    mine <- long$profile == profile
    long$probability[mine][match(keys, key[mine])]
  })
  names(columns) <- paste0("profile_", profiles)
  data.frame(indicator = long$indicator[first], category = long$category[first],
             columns, check.names = FALSE, stringsAsFactors = FALSE)
}

#' Print the estimated class sizes, and the group-class sizes of a two-level fit
#' @param x A fitted model.
#' @return `NULL`, invisibly.
#' @noRd
.multilpa_print_class_sizes <- function(x) {
  sizes <- .multilpa_count_frame(x)
  line <- function(level, label, prefix) {
    rows <- sizes[sizes$level == level, , drop = FALSE]
    if (nrow(rows) == 0L) return(invisible(NULL))
    cat(sprintf("%-14s%s\n", label, paste(sprintf("%s%s %.1f%% (%.0f)", prefix, rows$class,
                                                   100 * rows$effective_proportion,
                                                   rows$effective_count),
                                           collapse = "   ")))
  }
  cat("\n")
  line("individuals", "Class sizes", "profile_")
  if (!isTRUE(x$single_level) && !identical(as.integer(x$n_group_classes), 1L)) {
    line("groups", "Group classes", "group_class_")
  }
  invisible(NULL)
}

#' Print a fit's primary table under its header
#'
#' A fitted model auto-prints its own estimates. The alternative is a header
#' that describes a result the reader then has to fetch with a second call,
#' which is the ritual this package exists to remove: `fit` and
#' `as.data.frame(fit)` would otherwise be two calls showing one thing.
#'
#' The measurement model is the primary table, and an all-categorical fit has
#' none, so that fit prints its response probabilities instead of an empty
#' block under a heading promising means.
#'
#' @param x A fitted model of this package.
#' @param rows How many rows to print before truncating.
#' @return `NULL`, invisibly. Called for the printing.
#' @noRd
.multilpa_print_primary <- function(x, rows = 20L) {
  stopifnot("`rows` must be a single non-negative whole number" =
              is.numeric(rows) && length(rows) == 1L && is.finite(rows) &&
              rows >= 0 && rows == as.integer(rows))
  name <- names(.multilpa_catalogue(x))[1L]
  # Printing shows the estimates; the standard errors belong to the tables
  # asked for by name, and computing them on every print would be slow.
  table <- if (identical(name, "profiles")) .multilpa_profile_frame(x) else
    get_results(x)
  ## The Gaussian measurement is printed one row per profile. Everything else
  ## keeps the catalogue's own shape, which is already one row per thing.
  if (identical(name, "profiles") && nrow(table) > 0L) {
    ## How big each profile is comes first: a mean is read differently when it
    ## describes two percent of the sample. The means follow one row per
    ## indicator and one column per profile, as a class model's probabilities
    ## do, so many indicators lengthen the table instead of wrapping it.
    .multilpa_print_class_sizes(x)
    # Read by position from the wide table, never by name: an indicator may
    # itself be called "profile".
    wide <- .multilpa_wide_means(x)
    indicators <- names(wide)[-1L]
    means <- data.frame(indicator = indicators, stringsAsFactors = FALSE)
    invisible(lapply(seq_len(nrow(wide)), function(index) {
      means[[paste0("profile_", wide[[1L]][[index]])]] <<-
        round(vapply(wide[index, -1L], as.numeric, numeric(1)), 3L)
    }))
    cat("\nMeans\n")
    if (rows > 0) print(utils::head(means, rows), row.names = FALSE)
    if (nrow(means) > rows) {
      cat(sprintf("   ... %d more indicators.  get_results(x, \"profiles\")\n",
                  nrow(means) - rows))
    }
    cat("\nVariances and standard errors: get_results(x, \"profiles\").\n",
        "Every other table: get_results(x, what = ), or get_results(x, \"all\").\n", sep = "")
    return(invisible(NULL))
  }
  if (nrow(table) == 0L && length(x$categorical %||% character()) > 0L) {
    # A class model reads as a profile of response probabilities: one row per
    # indicator category, one column per class, under the class sizes.
    .multilpa_print_class_sizes(x)
    wide <- .multilpa_wide_responses(x)
    numeric_columns <- vapply(wide, is.numeric, logical(1))
    wide[numeric_columns] <- lapply(wide[numeric_columns], round, digits = 3L)
    cat("\nResponse probabilities\n")
    if (rows > 0) print(utils::head(wide, rows), row.names = FALSE)
    if (nrow(wide) > rows) {
      cat(sprintf("   ... %d more rows.  get_results(x, \"responses\")\n", nrow(wide) - rows))
    }
    cat("\nThresholds and standard errors: get_results(x, \"responses\").\n",
        "Every other table: get_results(x, what = ), or get_results(x, \"all\").\n", sep = "")
    return(invisible(NULL))
  }
  if (nrow(table) == 0L) return(invisible(NULL))
  cat("\n")
  if (rows > 0) print(utils::head(table, rows), row.names = FALSE)
  if (nrow(table) > rows) {
    cat(sprintf("   ... %d more rows.  get_results(x, \"%s\")\n",
                nrow(table) - rows, name))
  }
  cat("\nEvery other table: get_results(x, what = ), or get_results(x, \"all\").\n")
  invisible(NULL)
}
